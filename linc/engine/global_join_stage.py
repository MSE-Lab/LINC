"""GlobalJoinStage: 增量并入全局代表库(Step3,核心阶段)。

设计依据: docs/归档/STEP3_增量并入全局代表库.md(算法语义) +
docs/LINC_技术框架设计_v2.md 第2.4节"GlobalJoinStage 性能方案"(工程实现,
2026-10-01 与用户三轮讨论确认):

- Order-D 确定性排序: (基因总数 ↓, 局部代表数 ↓, genome_id 字典序 ↑)。
- query/target 方向与 IntraGenomeStage 保持一致: 全局代表库整体作为 query,
  对该基因组在 IntraGenomeStage 阶段已建好的持久化索引(target)做比对,
  复用索引、不重新建库;一个基因组一次比对调用(基因组级批处理快照)。
  正确性依据: IntraGenomeStage.multi_representative_clustering 的贪心覆盖
  保证同一基因组内所有局部代表彼此 pairwise identity 必然 < Cutoff-II,
  批处理快照不会引入"同基因组内本应共享代表却被误判为各自新建簇"的问题。
- 双阈值决策树: 命中 0/1/>=2 个簇分三路处理, 一次比对结果同时用于
  簇归属判定/合并触发/代表裁决, 不二次比对。
"""
import logging
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Optional

from ..core.entities import (
    AlignmentResult,
    Gene,
    Genome,
    GenomeIntraResult,
    GlobalRepresentativeLibrary,
    IncrementalResult,
    SSNEdge,
)
from ..core.events import (
    ClusterCreated,
    ClusterMerged,
    DomainEvent,
    EventBus,
    GenomeJoined,
    RepresentativeAdded,
)
from ..ports.aligner_port import AlignerPort

logger = logging.getLogger(__name__)


class GlobalJoinStage:
    """对应 v1 Step3: Order-D 排序 + 双阈值决策树, 增量并入全局代表库。"""

    def __init__(
        self,
        aligner: AlignerPort,
        cutoff_i: float = 0.5,
        cutoff_ii: float = 0.8,
        max_per_cluster: int = 5,
        event_bus: Optional[EventBus] = None,
    ) -> None:
        self._aligner = aligner
        self._cutoff_i = cutoff_i
        self._cutoff_ii = cutoff_ii
        self._max_per_cluster = max_per_cluster
        self._event_bus = event_bus
        self._cluster_counter = 0

    def run(
        self,
        genomes: Dict[str, Genome],
        intra_results: Dict[str, GenomeIntraResult],
    ) -> IncrementalResult:
        order = self._order_d_sort(intra_results)
        logger.info("Order-D 排序结果: %s", order)

        global_lib = GlobalRepresentativeLibrary()
        ssn_edges = set()
        total = len(order)

        for idx, genome_id in enumerate(order, start=1):
            intra = intra_results[genome_id]
            genome = genomes[genome_id]
            logger.info(
                "[%d/%d] 正在处理基因组 %s (局部代表数=%d); "
                "当前全局代表库: %d 个簇 / %d 个代表",
                idx, total, genome_id, intra.n_representatives,
                global_lib.n_active_clusters, global_lib.n_representatives,
            )
            self._process_genome(genome, intra, global_lib, ssn_edges)
            logger.info(
                "[%d/%d] 基因组 %s 处理完成; 全局代表库更新为: %d 个簇 / %d 个代表",
                idx, total, genome_id,
                global_lib.n_active_clusters, global_lib.n_representatives,
            )
            self._publish(
                GenomeJoined(
                    timestamp=_now(),
                    payload={
                        "genome_id": genome_id,
                        "n_local_representatives": intra.n_representatives,
                        "n_global_clusters_so_far": global_lib.n_active_clusters,
                        "n_global_representatives_so_far": global_lib.n_representatives,
                    },
                )
            )

        result = IncrementalResult(
            global_lib=global_lib,
            ssn_edges=ssn_edges,
            processing_order=order,
            n_clusters_final=global_lib.n_active_clusters,
            n_representatives_final=global_lib.n_representatives,
        )
        logger.info(
            "GlobalJoinStage 完成: genomes=%d n_clusters_final=%d "
            "n_representatives_final=%d",
            len(order), result.n_clusters_final, result.n_representatives_final,
        )
        return result

    @staticmethod
    def _order_d_sort(intra_results: Dict[str, GenomeIntraResult]) -> List[str]:
        """Order-D: (基因总数 ↓, 局部代表数 ↓, genome_id 字典序 ↑)。"""
        return sorted(
            intra_results.keys(),
            key=lambda gid: (
                -intra_results[gid].n_genes,
                -intra_results[gid].n_representatives,
                gid,
            ),
        )

    def _process_genome(
        self,
        genome: Genome,
        intra: GenomeIntraResult,
        global_lib: GlobalRepresentativeLibrary,
        ssn_edges: set,
    ) -> None:
        if not intra.local_representatives:
            return

        genes_by_id: Dict[str, Gene] = {g.gene_id: g for g in genome.genes}

        # 基因组级批处理快照(与用户确认,2026-10-01):全局库整体作为 query,
        # 对该基因组持久化索引做一次比对, 该基因组内所有局部代表的决策都基于
        # 同一份快照结果(处理这个基因组之前的全局库状态)。
        query_genes = global_lib.representative_genes()
        if query_genes and intra.genome_index is not None:
            hits = self._aligner.search(query_genes, intra.genome_index)
        else:
            hits = []

        hits_by_gene: Dict[str, List[AlignmentResult]] = defaultdict(list)
        for hit in hits:
            if hit.identity >= self._cutoff_i and hit.coverage >= self._cutoff_i:
                hits_by_gene[hit.subject_id].append(hit)

        for rep_id in sorted(intra.local_representatives):
            rep_gene = genes_by_id[rep_id]
            self._decide_and_apply(
                rep_gene, hits_by_gene.get(rep_id, []), global_lib, ssn_edges
            )

    def _decide_and_apply(
        self,
        rep_gene: Gene,
        rep_hits: List[AlignmentResult],
        global_lib: GlobalRepresentativeLibrary,
        ssn_edges: set,
    ) -> None:
        """双阈值决策树: 按命中簇数量(0/1/>=2)分三路处理。"""
        # 同一簇可能被多个代表命中, 取该簇的最佳(最高 identity)一条命中结果。
        hit_by_cluster: Dict[str, AlignmentResult] = {}
        for hit in rep_hits:
            cluster_id = global_lib.gene_to_cluster.get(hit.query_id)
            if cluster_id is None:
                continue
            cluster = global_lib.get_cluster(cluster_id)
            if cluster is None or not cluster.active:
                continue
            existing = hit_by_cluster.get(cluster_id)
            if existing is None or hit.identity > existing.identity:
                hit_by_cluster[cluster_id] = hit

        for hit in hit_by_cluster.values():
            ssn_edges.add(
                SSNEdge(
                    source=rep_gene.gene_id,
                    target=hit.query_id,
                    identity=hit.identity,
                    coverage=hit.coverage,
                )
            )

        if len(hit_by_cluster) == 0:
            # 路径①: 零簇命中 -> 创建新簇
            cluster_id = self._new_cluster_id()
            global_lib.add_new_cluster(cluster_id, rep_gene)
            self._publish(
                ClusterCreated(
                    timestamp=_now(),
                    payload={"cluster_id": cluster_id, "representative": rep_gene.gene_id},
                )
            )
            return

        if len(hit_by_cluster) == 1:
            # 路径②: 单簇命中
            (cluster_id, best_hit), = hit_by_cluster.items()
            if best_hit.identity >= self._cutoff_ii:
                global_lib.add_member_only(cluster_id, rep_gene)
            else:
                cluster = global_lib.get_cluster(cluster_id)
                if len(cluster.representative_ids) < self._max_per_cluster:
                    global_lib.add_representative(cluster_id, rep_gene)
                    self._publish(
                        RepresentativeAdded(
                            timestamp=_now(),
                            payload={"cluster_id": cluster_id, "representative": rep_gene.gene_id},
                        )
                    )
                else:
                    global_lib.add_member_only(cluster_id, rep_gene)
            return

        # 路径③: 跨簇命中(>=2 簇)
        cluster_ids = list(hit_by_cluster.keys())
        any_cutoff_ii = any(hit.identity >= self._cutoff_ii for hit in hit_by_cluster.values())
        merged_id = self._new_cluster_id()
        if any_cutoff_ii:
            global_lib.merge_clusters(cluster_ids, merged_id)
            global_lib.add_member_only(merged_id, rep_gene)
        else:
            global_lib.merge_clusters(cluster_ids, merged_id, new_representative=rep_gene)
        self._publish(
            ClusterMerged(
                timestamp=_now(),
                payload={
                    "merged_cluster_id": merged_id,
                    "source_cluster_ids": cluster_ids,
                    "triggered_by": rep_gene.gene_id,
                    "new_representative_added": not any_cutoff_ii,
                },
            )
        )

    def _new_cluster_id(self) -> str:
        cluster_id = f"global|cluster_{self._cluster_counter}"
        self._cluster_counter += 1
        return cluster_id

    def _publish(self, event: DomainEvent) -> None:
        if self._event_bus is not None:
            self._event_bus.publish(event)


def _now() -> datetime:
    return datetime.now()
