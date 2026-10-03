"""IntraGenomeStage: 基因组内自聚类阶段。

对应 v1 Step2(docs/归档/STEP2_基因组内自聚类.md)与设计文档 2.4 节的接口契约:
    输入: Genome
    输出: GenomeIntraResult

本阶段职责:
    1. 对基因组内全部基因做全对全比对(AlignerPort.align_all_vs_all)。
    2. 用 Cutoff-I 过滤比对结果,构建图,交给 ClusteringPort(Union-Find)
       做连通分量聚类,得到局部簇。
    3. 对每个局部簇选取代表序列(双阈值决策):
       - 单成员簇(孤儿基因): 自身为代表。
       - 多成员簇: 若簇内任意两个成员的 pairwise identity >= Cutoff-II 且
         coverage >= Cutoff-C(即任意两个成员都高度相似且覆盖充分, 真正
         意义上的"整体均质"), 选 1 个"中心代表"(与簇内其他成员平均相似度
         最高者); 否则(至少存在一对成员差异较大或覆盖不足, 可能存在亚型
         分化)走 multi_representative_clustering, 用贪心覆盖(k-center
         近似)选出多个代表, 数量不超过 max_per_cluster。

       注: v1 STEP2 文档字面写的判断依据是"最大" pairwise identity, 但这在
       语义上不成立(仅凭一对高相似度成员不能代表整体均质, 且会导致贪心
       覆盖的停止条件永远无法触发), 经与用户确认后改为"最小" pairwise
       identity 作为判断依据。进一步, 第二轮漏洞3修复(2026-10-03, 与用户
       确认): 高置信判定("是否已被覆盖/是否整体均质")此前只看 identity,
       未考虑 coverage, 默认 cutoff_c=0.0(覆盖度过滤关闭)下会把
       identity 高但 coverage 很低的碎片匹配误判为高置信。修复后统一改为
       identity >= cutoff_ii **且** coverage >= cutoff_c 的双条件判定,
       与 passes_inclusion_cutoff() 的 coverage 判据口径(cutoff_c)保持
       一致; 同时 --cutoff-c 默认值由 0.0 改为 0.5, 使该过滤默认生效。

全对全比对复杂度为 O(n^2)(n 为单个基因组的基因数),这是 STEP2 文档强调的
"必须全量比对、不能增量"的固有开销;代表选取逻辑仅基于已有的相似度矩阵
做轻量级查表计算,不产生新的比对开销。
"""
import logging
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import networkx as nx

from ..core.entities import (
    COV_MODE_BIDIRECTIONAL,
    AlignmentResult,
    Genome,
    GenomeIntraResult,
    LocalCluster,
    passes_inclusion_cutoff,
)
from ..core.events import ClusterCreated, DomainEvent, EventBus
from ..ports.aligner_port import AlignerPort
from ..ports.clustering_port import ClusteringPort

logger = logging.getLogger(__name__)

_IdentityMap = Dict[Tuple[str, str], float]


class IntraGenomeStage:
    """对应 v1 Step2 的现代化实现。"""

    def __init__(
        self,
        aligner: AlignerPort,
        clustering: ClusteringPort,
        cutoff_i: float = 0.5,
        cutoff_ii: float = 0.8,
        cutoff_c: float = 0.5,
        cov_mode: int = COV_MODE_BIDIRECTIONAL,
        max_per_cluster: int = 5,
        event_bus: Optional[EventBus] = None,
    ) -> None:
        self._aligner = aligner
        self._clustering = clustering
        self._cutoff_i = cutoff_i
        self._cutoff_ii = cutoff_ii
        self._cutoff_c = cutoff_c
        self._cov_mode = cov_mode
        self._max_per_cluster = max_per_cluster
        self._event_bus = event_bus

    def run(self, genome: Genome) -> GenomeIntraResult:
        genes = genome.genes
        n_genes = len(genes)

        if n_genes == 0:
            return GenomeIntraResult(genome_id=genome.genome_id)

        # build_index 持久化建库(与用户确认,2026-10-01):建一次索引, 后续
        # GlobalJoinStage 反过来用全局库 query 该索引, 复用而不重新建库
        # (详见 docs/LINC_技术框架设计_v2.md 第2.4节"GlobalJoinStage 性能方案")。
        genome_index = self._aligner.build_index(genes, index_id=genome.genome_id)
        # self_search=True(2026-10-02, 方案3, 与用户确认): query 就是
        # genome_index 对应的同一份 genes, 让 aligner 直接复用已持久化的
        # genes.fasta/索引, 不再重复写一份内容完全相同的 query.fasta。
        hits = self._aligner.search(genes, genome_index, tag="intra_self", self_search=True)

        identity_map: _IdentityMap = {}
        coverage_map: _IdentityMap = {}
        pairwise_hits: List[AlignmentResult] = []
        graph = nx.Graph()
        graph.add_nodes_from(gene.gene_id for gene in genes)
        for hit in hits:
            if hit.query_id == hit.subject_id:
                continue   # 自身比对(all-vs-all 场景下 query/target 来自同一份序列集合)
            identity_map[(hit.query_id, hit.subject_id)] = hit.identity
            identity_map[(hit.subject_id, hit.query_id)] = hit.identity
            # 第二轮漏洞3修复(2026-10-03, 与用户确认): 同步记录 coverage,
            # 用于代表选取阶段的高置信判定不再只看 identity。
            coverage_map[(hit.query_id, hit.subject_id)] = hit.coverage
            coverage_map[(hit.subject_id, hit.query_id)] = hit.coverage
            pairwise_hits.append(hit)
            if passes_inclusion_cutoff(hit, self._cutoff_i, self._cutoff_c, self._cov_mode):
                graph.add_edge(hit.query_id, hit.subject_id)

        groups = self._clustering.cluster(graph)

        local_clusters: List[LocalCluster] = []
        local_representatives: List[str] = []
        orphan_genes: List[str] = []

        for index, member_ids in enumerate(groups):
            cluster_id = f"{genome.genome_id}|cluster_{index}"
            is_orphan = len(member_ids) == 1

            if is_orphan:
                representative_ids = list(member_ids)
                orphan_genes.extend(representative_ids)
            else:
                representative_ids = self._select_representatives(
                    member_ids, identity_map, coverage_map
                )

            local_clusters.append(
                LocalCluster(
                    cluster_id=cluster_id,
                    member_ids=set(member_ids),
                    representative_ids=representative_ids,
                    is_orphan=is_orphan,
                )
            )
            local_representatives.extend(representative_ids)

            self._publish(
                ClusterCreated(
                    timestamp=_now(),
                    payload={
                        "genome_id": genome.genome_id,
                        "cluster_id": cluster_id,
                        "member_count": len(member_ids),
                        "representative_count": len(representative_ids),
                        "is_orphan": is_orphan,
                    },
                )
            )

        result = GenomeIntraResult(
            genome_id=genome.genome_id,
            local_clusters=local_clusters,
            local_representatives=local_representatives,
            orphan_genes=orphan_genes,
            n_genes=n_genes,
            n_clusters=len(local_clusters),
            n_representatives=len(local_representatives),
            genome_index=genome_index,
            pairwise_hits=pairwise_hits,
        )
        logger.info(
            "IntraGenomeStage 完成: genome_id=%s n_genes=%d n_clusters=%d "
            "n_representatives=%d",
            genome.genome_id, n_genes, result.n_clusters, result.n_representatives,
        )
        return result

    def _select_representatives(
        self, member_ids, identity_map: _IdentityMap, coverage_map: _IdentityMap
    ) -> List[str]:
        members = sorted(member_ids)
        # 第二轮漏洞3修复(2026-10-03, 与用户确认): 是否"整体均质"需要
        # 簇内任意两个成员同时满足 identity >= cutoff_ii 且
        # coverage >= cutoff_c, 而不是只看 identity。
        all_high_confidence = True
        for i, a in enumerate(members):
            for b in members[i + 1:]:
                identity = identity_map.get((a, b), 0.0)
                coverage = coverage_map.get((a, b), 0.0)
                if identity < self._cutoff_ii or coverage < self._cutoff_c:
                    all_high_confidence = False
                    break
            if not all_high_confidence:
                break

        if all_high_confidence:
            return [self._select_central_representative(members, identity_map)]
        return self._multi_representative_clustering(members, identity_map, coverage_map)

    @staticmethod
    def _select_central_representative(
        members: List[str], identity_map: _IdentityMap
    ) -> str:
        """选取与簇内其他成员平均相似度最高的成员作为中心代表。"""
        best_member = members[0]
        best_avg = -1.0
        for member in members:
            others = [m for m in members if m != member]
            avg_identity = (
                sum(identity_map.get((member, other), 0.0) for other in others)
                / len(others)
            )
            if avg_identity > best_avg:
                best_avg = avg_identity
                best_member = member
        return best_member

    def _multi_representative_clustering(
        self, members: List[str], identity_map: _IdentityMap, coverage_map: _IdentityMap
    ) -> List[str]:
        """贪心覆盖(k-center 近似)选出多个代表。

        1. 先选 1 个"中心代表"(与簇内其他成员平均相似度最高者)。
        2. 反复选"与已选代表集合的最大相似度最低"的剩余成员作为新代表
           (这一步仍用 identity 作为贪心多样性的依据), 直到达到
           max_per_cluster 个, 或剩余成员都已被某代表"覆盖"(第二轮漏洞3
           修复, 2026-10-03, 与用户确认: 覆盖判定改为 identity >= Cutoff-II
           且 coverage >= Cutoff-C 的双条件)为止。
        """
        representatives = [
            self._select_central_representative(members, identity_map)
        ]
        remaining = [m for m in members if m not in representatives]

        while remaining and len(representatives) < self._max_per_cluster:
            def max_sim_to_reps(member: str) -> float:
                return max(
                    identity_map.get((member, rep), 0.0) for rep in representatives
                )

            def is_covered(member: str) -> bool:
                return any(
                    identity_map.get((member, rep), 0.0) >= self._cutoff_ii
                    and coverage_map.get((member, rep), 0.0) >= self._cutoff_c
                    for rep in representatives
                )

            # 修复(2026-10-03): 不能先按 identity 贪心选出候选再单独检查该
            # 候选是否被覆盖——当多个剩余成员的 identity 并列最低时,先选中
            # 的候选可能恰好"被覆盖"(双条件达标)而被提前 break,导致其他
            # identity 并列但 coverage 不达标、真正"未被覆盖"的成员被遗漏。
            # 应先筛出所有"未被覆盖"的剩余成员,只要存在未覆盖成员就必须
            # 继续选代表(按 identity 贪心挑其中相似度最低者),全部被覆盖
            # 才停止。
            uncovered = [m for m in remaining if not is_covered(m)]
            if not uncovered:
                break
            candidate = min(uncovered, key=max_sim_to_reps)
            representatives.append(candidate)
            remaining.remove(candidate)

        return representatives

    def _publish(self, event: DomainEvent) -> None:
        if self._event_bus is not None:
            self._event_bus.publish(event)


def _now() -> datetime:
    return datetime.now()
