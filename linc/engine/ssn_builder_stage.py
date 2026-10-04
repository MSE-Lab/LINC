"""SSNBuilderStage: 序列相似性网络(SSN)构建阶段。

对应 v1 Step4(docs/归档/STEP4_稀疏SSN构建.md)+ Step5(docs/归档/STEP5_最终补边.md)
的合并现代化实现(2026-10-02 与用户确认: 不再保留 SparseSSN/CompleteSSN 两套
独立产物, 合并为一个 Stage, 输出一个 CompleteSSN)。

本阶段职责:
    1. 以 GlobalJoinStage 产出的代表-代表边(IncrementalResult.ssn_edges,
       稀疏SSN)为基础。
    2. 按 `edge_completion_mode` 补齐簇内关系, 三档可选(与用户确认):
       - none: 不补边, 直接返回稀疏SSN。
       - representative_member(默认): 补"代表-成员"边, 复用 IntraGenomeStage
         已经算出的全对全比对结果(pairwise_hits), 不产生任何新比对。
       - full: 在 representative_member 基础上, 进一步补全簇内全部基因对
         (含成员-成员)关系; 同基因组内的成员-成员关系依然复用已有结果,
         只有跨基因组的成员对才真正发起新比对。性能方案(2026-10-02 与用户
         确认, 修复"逐簇调用子进程"的性能问题; 2026-10-03 与用户确认,
         进一步修复"仍是 O(基因组对数)次子进程调用"的性能问题): 不按"簇"
         或"基因组对"循环调用比对, 而是先遍历全部活跃簇, 按"target 基因组"
         聚合出需要查询的基因集合(同一 target 基因组下, 来自不同 source
         基因组的 query 基因去重合并到同一个集合里), 再对每个 target
         基因组只调用一次 AlignerPort.search(复用 IntraGenomeStage 已
         持久化的 genome_index 作为 target, 避免重新建库; tag=
         f"full_to_{genome_b}", 与 intra 阶段的 tag="intra_self"/join 阶段的
         tag="join" 区分开, 避免三类调用写入同一个 `search/` 子目录互相
         覆盖, 2026-10-02 与用户确认), 调用次数从 O(基因组对数)进一步降到
         O(作为target出现过的基因组数)(最坏情况下从 O(G²) 降到 O(G));
         search 结果按"双方是否属于同一个簇"过滤(该过滤本身与 query 具体
         来自哪个 source 基因组无关, 因此按 target 合并查询不改变结果),
         排除同一 target 基因组下属于不同簇的无关比对结果。
    3. 巨簇处理: 仅做信息提醒(超过 `mega_cluster_warning_threshold` 时记录
       warning), 不做自动降级或拆分 —— 与用户确认(2026-10-02): 任何跳过
       真实比对的"降级"或"拆分"都会让结果偏离"全对全比对"的最初设计目标,
       产出有损/近似的 SSN, 这个决定权应交还给用户, 而不是由程序自作主张。
    4. 纳入判据一致性(2026-10-03 与用户确认, 修复第二轮漏洞2): 本阶段补边
       时无论是复用 `pairwise_hits`(代表-成员/成员-成员/同基因组跨局部簇)
       还是跨基因组新发起比对的结果, 都必须统一调用
       `passes_inclusion_cutoff(hit, cutoff_i, cutoff_c, cov_mode)` 过滤,
       与 IntraGenomeStage/GlobalJoinStage 的判据口径保持一致, 避免混入
       identity < cutoff_i 的弱边改变 MCL 细分阶段看到的权重分布。
"""
import logging
import os
from typing import Dict, List, Optional, Set

from ..core.entities import (
    COV_MODE_BIDIRECTIONAL,
    EDGE_COMPLETION_FULL,
    EDGE_COMPLETION_NONE,
    EDGE_COMPLETION_REPRESENTATIVE_MEMBER,
    AlignmentResult,
    Cluster,
    CompleteSSN,
    Gene,
    Genome,
    GenomeIntraResult,
    IncrementalResult,
    SSNEdge,
    passes_inclusion_cutoff,
)
from ..core.events import DomainEvent, EventBus
from ..ports.aligner_port import AlignerPort

logger = logging.getLogger(__name__)


class SSNBuilderStage:
    """对应 v1 Step4+5 的合并现代化实现。"""

    def __init__(
        self,
        aligner: AlignerPort,
        edge_completion_mode: str = EDGE_COMPLETION_REPRESENTATIVE_MEMBER,
        mega_cluster_warning_threshold: Optional[int] = None,
        event_bus: Optional[EventBus] = None,
        cutoff_i: float = 0.5,
        cutoff_c: float = 0.5,
        cov_mode: int = COV_MODE_BIDIRECTIONAL,
        cleanup_intermediate: bool = True,
        dump_raw_hits: Optional[str] = None,
    ) -> None:
        self._aligner = aligner
        self._edge_completion_mode = edge_completion_mode
        self._mega_cluster_warning_threshold = mega_cluster_warning_threshold
        self._event_bus = event_bus
        # 数据管理优化第3项(2026-10-03 与用户确认): full 模式下 search()
        # 落盘的中间文件(search_full_to_*/) 结果已转存进 CompleteSSN,
        # 默认用完即删。
        self._cleanup_intermediate = cleanup_intermediate
        # 与用户确认修复第二轮漏洞2(2026-10-03): 补边阶段复用/新发起比对产生
        # 的边, 必须与 IntraGenomeStage/GlobalJoinStage 一样统一走
        # passes_inclusion_cutoff() 判定, 保持 SSN 边纳入判据前后一致,
        # 避免混入 identity < cutoff_i 的弱边改变 MCL 看到的权重分布。
        self._cutoff_i = cutoff_i
        self._cutoff_c = cutoff_c
        self._cov_mode = cov_mode
        # 调试出口(2026-10-03 与用户确认): 非空时把 full 补边阶段 cutoff 过滤
        # 前的全部原始跨基因组 hit 落盘到该目录, 供离线 RBB 验证使用。
        self._dump_raw_hits = dump_raw_hits

    def run(
        self,
        genomes: Dict[str, Genome],
        intra_results: Dict[str, GenomeIntraResult],
        incremental_result: IncrementalResult,
    ) -> CompleteSSN:
        edges: Set[SSNEdge] = set(incremental_result.ssn_edges)
        logger.info(
            "SSNBuilderStage 开始: edge_completion_mode=%s, 初始(稀疏)边数=%d",
            self._edge_completion_mode, len(edges),
        )

        if self._edge_completion_mode == EDGE_COMPLETION_NONE:
            logger.info("edge_completion_mode=none, 跳过补边, 仅保留稀疏SSN")
            return CompleteSSN(edges=edges, edge_completion_mode=self._edge_completion_mode)

        self._add_known_intra_genome_edges(
            intra_results, edges,
            include_member_member=(self._edge_completion_mode == EDGE_COMPLETION_FULL),
        )

        if self._edge_completion_mode == EDGE_COMPLETION_REPRESENTATIVE_MEMBER:
            logger.info("SSNBuilderStage 完成(representative_member模式): 共 %d 条边", len(edges))
            return CompleteSSN(edges=edges, edge_completion_mode=self._edge_completion_mode)

        genes_by_id: Dict[str, Gene] = {
            gene.gene_id: gene for genome in genomes.values() for gene in genome.genes
        }
        self._complete_cross_genome_member_edges(
            intra_results, incremental_result, genes_by_id, edges
        )

        logger.info("SSNBuilderStage 完成(full模式): 共 %d 条边", len(edges))
        return CompleteSSN(edges=edges, edge_completion_mode=self._edge_completion_mode)

    def _add_known_intra_genome_edges(
        self,
        intra_results: Dict[str, GenomeIntraResult],
        edges: Set[SSNEdge],
        include_member_member: bool,
    ) -> None:
        """复用 IntraGenomeStage 已算出的 pairwise 结果, 补齐代表-成员(及
        full 模式下成员-成员)关系, 不产生任何新比对。

        依据(与用户确认, 2026-10-02): IntraGenomeStage 对单个基因组内的
        全部基因已经做过全对全比对, 同一基因组内任意两个基因之间的
        identity/coverage/bitscore 在 intra 阶段就已经算出, 补边阶段应直接
        复用, 而不是重新发起比对。

        纳入判据修复(2026-10-03, 与用户确认, 第二轮漏洞2): 复用的 hit 来自
        `pairwise_hits`, 其中包含 IntraGenomeStage 全对全比对的全部原始
        结果(不限于通过 cutoff 的), 补边前必须同样调用
        `passes_inclusion_cutoff()` 过滤, 否则会混入 identity < cutoff_i
        的弱边, 与稀疏SSN(GlobalJoinStage)的判据口径不一致。
        """
        for intra in intra_results.values():
            pair_lookup: Dict[frozenset, AlignmentResult] = {}
            for hit in intra.pairwise_hits:
                pair_lookup.setdefault(frozenset((hit.query_id, hit.subject_id)), hit)

            for local_cluster in intra.local_clusters:
                if local_cluster.is_orphan:
                    continue
                representatives = local_cluster.representative_ids
                plain_members = local_cluster.member_ids - set(representatives)

                for rep_id in representatives:
                    for member_id in plain_members:
                        hit = pair_lookup.get(frozenset((rep_id, member_id)))
                        if hit is not None and passes_inclusion_cutoff(
                            hit, self._cutoff_i, self._cutoff_c, self._cov_mode
                        ):
                            source_length, target_length = SSNBuilderStage._lengths_for(
                                hit, rep_id
                            )
                            edges.add(SSNEdge(
                                source=rep_id, target=member_id,
                                identity=hit.identity, coverage=hit.coverage,
                                bitscore=hit.bitscore,
                                source_length=source_length, target_length=target_length,
                            ))

                if include_member_member:
                    members = sorted(local_cluster.member_ids)
                    for i, gene_a in enumerate(members):
                        for gene_b in members[i + 1:]:
                            hit = pair_lookup.get(frozenset((gene_a, gene_b)))
                            if hit is not None and passes_inclusion_cutoff(
                                hit, self._cutoff_i, self._cutoff_c, self._cov_mode
                            ):
                                source_length, target_length = SSNBuilderStage._lengths_for(
                                    hit, gene_a
                                )
                                edges.add(SSNEdge(
                                    source=gene_a, target=gene_b,
                                    identity=hit.identity, coverage=hit.coverage,
                                    bitscore=hit.bitscore,
                                    source_length=source_length, target_length=target_length,
                                ))

    def _complete_cross_genome_member_edges(
        self,
        intra_results: Dict[str, GenomeIntraResult],
        incremental_result: IncrementalResult,
        genes_by_id: Dict[str, Gene],
        edges: Set[SSNEdge],
    ) -> None:
        """full 模式下真正需要新比对的部分: 同一个全局簇内, 来自不同基因组的
        成员之间从未被比对过(IntraGenomeStage 只做单基因组内全对全比对,
        GlobalJoinStage 只比对代表)。

        性能方案(2026-10-02 与用户确认, 修复"逐簇调用子进程"的性能问题;
        2026-10-03 与用户确认, 进一步修复"仍是 O(基因组对数)次子进程调用"
        的性能问题): 不按"簇"或"基因组对"循环发起比对调用, 而是先遍历全部
        活跃簇, 按"target 基因组"聚合出需要查询的基因集合(同一 target 基因组
        下, 来自不同 source 基因组的 query 基因去重合并到同一个集合里),
        再对每个 target 基因组只调用一次 `search`(复用 IntraGenomeStage 已
        持久化的 `genome_index` 作为 target, 不重新建库), 调用次数从
        O(基因组对数)进一步降到 O(作为target出现过的基因组数)(最坏情况下
        从 O(G²) 降到 O(G)); search 返回的结果可能包含"同一 target 基因组
        下但分属不同簇"的无关基因对, 需按 `gene_id -> cluster_id` 映射过滤,
        只保留双方属于同一个簇的结果(该过滤本身与 query 具体来自哪个
        source 基因组无关, 因此按 target 合并查询不会改变最终结果)。

        补边缺口修复(2026-10-03, 与用户确认, 第二轮漏洞1): 同一基因组、
        不同局部簇、但被全局合并进同一个全局簇的基因对, 之前既不在
        `_add_known_intra_genome_edges()`(只覆盖同一局部簇内)的范围内,
        也不在本方法原有逻辑(只处理 `genome_a != genome_b` 的跨基因组对,
        且 `len(genome_ids) < 2` 时整体跳过)的范围内, 导致这部分边永久
        缺失。修复方式: 同基因组内任意两个基因在 IntraGenomeStage 阶段
        已经做过全对全比对(结果都在 `pairwise_hits` 里), 这里只需按
        `gene_id -> local_cluster_id` 映射识别出"跨局部簇"的基因对,
        直接查表补边, 不发起任何新比对, 不引入额外子进程调用。

        bitscore 方向约定(第二轮漏洞7, 2026-10-03 与用户确认, 文档澄清、
        不改动逻辑): 下方"第二遍"对跨基因组基因对发起新比对时, 只按
        genome_id 字典序固定一个方向(genome_a < genome_b 时以 genome_a
        的基因作 query、genome_b 的基因作 subject), 不会反向再比对一次。
        因此产生的 SSNEdge.source/target(及其 bitscore/source_length/
        target_length)的 query/subject 身份由"谁的 genome_id 字典序更小"
        决定, 是人为约定而非语义上的"真实"方向; 详见 `SSNEdge` 类 docstring
        中的说明。
        """
        global_lib = incremental_result.global_lib
        gene_cluster_map: Dict[str, str] = {}

        # cluster_id -> Cluster 的本地索引(2026-10-04 修复: 不直接用
        # global_lib.get_cluster(), 因为其内部 _cluster_index 只在
        # add_new_cluster/merge_clusters 等"写路径"里维护; 若调用方/测试
        # 直接用 GlobalRepresentativeLibrary(clusters=[...]) 构造(未经
        # 写路径), _cluster_index 会是空的, get_cluster 恒返回 None。
        # 这里直接从 global_lib.clusters 现算一份, 保证无论 global_lib
        # 如何构造都能正确查到簇, 下面的 touching_map 扩搜索范围逻辑复用
        # 此索引。)
        cluster_by_id: Dict[str, Cluster] = {c.cluster_id: c for c in global_lib.clusters}

        # gene_id -> 所属局部簇 id(用于识别"同基因组, 跨局部簇"的基因对;
        # 同一局部簇内的基因对已经由 _add_known_intra_genome_edges() 覆盖,
        # 这里需要排除掉, 避免重复处理)。
        gene_to_local_cluster: Dict[str, str] = {}
        for intra in intra_results.values():
            for local_cluster in intra.local_clusters:
                for gene_id in local_cluster.member_ids:
                    gene_to_local_cluster[gene_id] = local_cluster.cluster_id

        # 每个基因组的 pairwise_hits 按 frozenset(gene_a, gene_b) 建索引,
        # 懒加载缓存, 供"同基因组跨局部簇"补边复用, 不重新比对。
        pair_lookup_cache: Dict[str, Dict[frozenset, AlignmentResult]] = {}

        def _pair_lookup_for(genome_id: str) -> Dict[frozenset, AlignmentResult]:
            cached = pair_lookup_cache.get(genome_id)
            if cached is not None:
                return cached
            lookup: Dict[frozenset, AlignmentResult] = {}
            intra = intra_results.get(genome_id)
            if intra is not None:
                for hit in intra.pairwise_hits:
                    lookup.setdefault(frozenset((hit.query_id, hit.subject_id)), hit)
            pair_lookup_cache[genome_id] = lookup
            return lookup

        # 第一遍: 仅做簇内基因组分布统计(不涉及任何比对调用), 用于①巨簇
        # warning 提醒 ②按"target 基因组"聚合出去重后的 query 基因集合
        # (2026-10-03 与用户确认, 从"按基因组对"改为"按 target 基因组"聚合,
        # 把调用次数从 O(基因组对数) 降到 O(作为target出现过的基因组数))
        # ③顺带补齐"同基因组跨局部簇"的边(直接查表, 不算新比对)。
        query_genes_by_target_genome: Dict[str, Set[str]] = {}
        for cluster in global_lib.clusters:
            if not cluster.active:
                continue
            members_by_genome: Dict[str, List[str]] = {}
            for gene_id in cluster.member_ids:
                gene_cluster_map[gene_id] = cluster.cluster_id
                gene = genes_by_id.get(gene_id)
                if gene is None:
                    continue
                members_by_genome.setdefault(gene.genome_id, []).append(gene_id)

            # 同基因组、跨局部簇补边: 对簇内每个基因组, 若其成员分属多个
            # 局部簇, 则把"局部簇不同"的基因对直接从该基因组的
            # pairwise_hits 里查出来补边(同局部簇内的对已由
            # _add_known_intra_genome_edges() 覆盖, 这里跳过避免重复)。
            for genome_id, gene_ids in members_by_genome.items():
                if len(gene_ids) < 2:
                    continue
                sorted_genes = sorted(gene_ids)
                distinct_local_clusters = {
                    gene_to_local_cluster.get(gid) for gid in sorted_genes
                }
                if len(distinct_local_clusters) < 2:
                    continue   # 同一局部簇内, 已被上一步覆盖
                pair_lookup = _pair_lookup_for(genome_id)
                for i, gene_a in enumerate(sorted_genes):
                    for gene_b in sorted_genes[i + 1:]:
                        if gene_to_local_cluster.get(gene_a) == gene_to_local_cluster.get(gene_b):
                            continue   # 同一局部簇内的对, 跳过
                        hit = pair_lookup.get(frozenset((gene_a, gene_b)))
                        if hit is None:
                            continue
                        if not passes_inclusion_cutoff(
                            hit, self._cutoff_i, self._cutoff_c, self._cov_mode
                        ):
                            continue
                        source_length, target_length = self._lengths_for(hit, gene_a)
                        edges.add(SSNEdge(
                            source=gene_a, target=gene_b,
                            identity=hit.identity, coverage=hit.coverage,
                            bitscore=hit.bitscore,
                            source_length=source_length, target_length=target_length,
                        ))

            genome_ids = sorted(members_by_genome.keys())
            if len(genome_ids) < 2:
                continue   # 单基因组簇, 同基因组内关系已在上两步复用完毕

            n_cross_pairs = self._count_cross_genome_pairs(members_by_genome, genome_ids)
            if (
                self._mega_cluster_warning_threshold is not None
                and n_cross_pairs > self._mega_cluster_warning_threshold
            ):
                logger.warning(
                    "SSNBuilderStage: 簇 %s 跨基因组待比对成员对数量(%d)超过阈值"
                    "(%d), full 模式补边可能耗时较长; 如需加速可考虑改用 "
                    "--edge-completion-mode representative_member, 或缩小输入规模。"
                    "(不会自动降级/拆分, 以避免产出有损/近似的SSN)",
                    cluster.cluster_id, n_cross_pairs, self._mega_cluster_warning_threshold,
                )

            for i, genome_a in enumerate(genome_ids):
                for genome_b in genome_ids[i + 1:]:
                    # 聚合键只用 target 基因组(genome_b), 不再按
                    # (genome_a, genome_b) 对分组: 多个不同的 genome_a 若
                    # 都需要查询同一个 genome_b, 它们的 query 基因会被合并
                    # 到同一个集合里, 只触发一次 search(2026-10-03 与用户
                    # 确认)。
                    query_genes_by_target_genome.setdefault(genome_b, set()).update(
                        members_by_genome[genome_a]
                    )

        # 【新增, 2026-10-04 接触判据/纳入判据解耦设计方案】扩搜索范围:
        # 上面的循环只覆盖"簇内已有成员"的基因组对, 但一个簇若只是对某个
        # 基因组有"evalue达标但纳入判据不达标"的弱命中(接触了该基因组,
        # 却从未真正合并进去), 该基因组此前永远不会被纳入 query 范围,
        # 导致"该簇其它基因组的成员其实可能比代表更强、本可成边"的情况
        # 被永久漏掉。这里按 GlobalJoinStage 记录的 touching_map(接触判据,
        # 与纳入判据/合并解耦)补充这部分 target 基因组的 query 基因,
        # query_genes_by_target_genome 是 Set, 与上面已有的成员重复添加
        # 会自动去重, 不影响已覆盖的部分。
        touching_map = incremental_result.touching_map or {}
        for genome_b, cluster_ids in touching_map.items():
            for cluster_id in cluster_ids:
                cluster = cluster_by_id.get(cluster_id)
                if cluster is None or not cluster.active:
                    continue
                for gene_id in cluster.member_ids:
                    gene_cluster_map[gene_id] = cluster.cluster_id
                    gene = genes_by_id.get(gene_id)
                    if gene is None or gene.genome_id == genome_b:
                        continue
                    query_genes_by_target_genome.setdefault(genome_b, set()).add(gene_id)

        logger.info(
            "SSNBuilderStage(full模式): 需要新发起比对的target基因组数=%d"
            "(与簇数/基因组对数无关)",
            len(query_genes_by_target_genome),
        )

        # 第二遍: 每个 target 基因组只调用一次 search(复用已持久化的
        # genome_index), 不逐簇、不逐基因组对调用子进程。
        for genome_b, query_gene_ids in query_genes_by_target_genome.items():
            query_genes = [genes_by_id[gid] for gid in sorted(query_gene_ids)]
            intra_b = intra_results.get(genome_b)
            index_b = intra_b.genome_index if intra_b is not None else None
            if index_b is not None:
                hits = self._aligner.search(
                    query_genes, index_b, tag=f"full_to_{genome_b}"
                )
                if self._cleanup_intermediate:
                    self._aligner.cleanup_search(index_b, tag=f"full_to_{genome_b}")
            else:
                # 回退: 没有持久化索引可复用(如测试/某些 adapter 实现),
                # 退化为一次性的便捷方法调用, 仍然按 target 基因组聚合,
                # 不逐簇/逐基因组对调用。
                subject_genes = [
                    gene for gid, gene in genes_by_id.items()
                    if gene.genome_id == genome_b
                ]
                hits = self._aligner.align_query_vs_subject(query_genes, subject_genes)

            self._dump_hits(hits, "full")

            for hit in hits:
                # 【改动, 2026-10-04 接触判据/纳入判据解耦设计方案】去掉
                # "同簇过滤": 此前这里要求 query/subject 双方属于同一个簇
                # 才加边, 但接触放宽后, "接触了 genome_b 但没并入 genome_b"
                # 的簇, 其成员命中 genome_b 的基因必然在别的簇(否则早被
                # 合并进同一簇了), 同簇过滤会把这部分真边全部拒绝掉,
                # 与本次改动的目的相悖。保留纳入判据(identity+coverage)
                # 作为唯一兜底, 防止混入弱边。
                if not passes_inclusion_cutoff(
                    hit, self._cutoff_i, self._cutoff_c, self._cov_mode
                ):
                    continue
                edges.add(SSNEdge(
                    source=hit.query_id, target=hit.subject_id,
                    identity=hit.identity, coverage=hit.coverage,
                    bitscore=hit.bitscore,
                    source_length=hit.query_length, target_length=hit.subject_length,
                ))

    @staticmethod
    def _lengths_for(hit: AlignmentResult, source_id: str) -> tuple:
        """把 AlignmentResult 的 query_length/subject_length 按 source/target
        的实际方向对齐返回(SSNEdge 的 source 不一定等于 hit.query_id,
        因为 pair_lookup 按 frozenset 去重, 可能存的是反向记录)。
        """
        if hit.query_id == source_id:
            return hit.query_length, hit.subject_length
        return hit.subject_length, hit.query_length

    def _dump_hits(self, hits: List[AlignmentResult], tag: str) -> None:
        """把 cutoff 过滤前的全部原始跨基因组 hit 追加落盘(调试出口)。

        列: query_id, subject_id, identity, query_cover, subject_cover,
        bitscore, query_length, subject_length。用于离线 RBB 等验证脚本。
        """
        if not self._dump_raw_hits:
            return
        os.makedirs(self._dump_raw_hits, exist_ok=True)
        path = os.path.join(self._dump_raw_hits, f"raw_{tag}_hits.tsv")
        with open(path, "a", encoding="utf-8") as fh:
            for hit in hits:
                fh.write(
                    f"{hit.query_id}\t{hit.subject_id}\t{hit.identity}\t"
                    f"{hit.query_cover}\t{hit.subject_cover}\t{hit.bitscore}\t"
                    f"{hit.query_length}\t{hit.subject_length}\n"
                )

    @staticmethod
    def _count_cross_genome_pairs(
        members_by_genome: Dict[str, List[str]], genome_ids: List[str]
    ) -> int:
        total = 0
        for i, genome_a in enumerate(genome_ids):
            for genome_b in genome_ids[i + 1:]:
                total += len(members_by_genome[genome_a]) * len(members_by_genome[genome_b])
        return total

    def _publish(self, event: DomainEvent) -> None:
        if self._event_bus is not None:
            self._event_bus.publish(event)
