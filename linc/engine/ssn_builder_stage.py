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
         只有跨基因组的成员对才真正发起新比对(调用 AlignerPort.
         align_query_vs_subject, 按"全局簇内涉及的基因组两两分组"批量调用,
         不逐 pair 调用以避免子进程开销)。
    3. 巨簇处理: 仅做信息提醒(超过 `mega_cluster_warning_threshold` 时记录
       warning), 不做自动降级或拆分 —— 与用户确认(2026-10-02): 任何跳过
       真实比对的"降级"或"拆分"都会让结果偏离"全对全比对"的最初设计目标,
       产出有损/近似的 SSN, 这个决定权应交还给用户, 而不是由程序自作主张。
"""
import logging
from typing import Dict, List, Optional, Set

from ..core.entities import (
    EDGE_COMPLETION_FULL,
    EDGE_COMPLETION_NONE,
    EDGE_COMPLETION_REPRESENTATIVE_MEMBER,
    AlignmentResult,
    CompleteSSN,
    Gene,
    Genome,
    GenomeIntraResult,
    IncrementalResult,
    SSNEdge,
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
    ) -> None:
        self._aligner = aligner
        self._edge_completion_mode = edge_completion_mode
        self._mega_cluster_warning_threshold = mega_cluster_warning_threshold
        self._event_bus = event_bus

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
        self._complete_cross_genome_member_edges(incremental_result, genes_by_id, edges)

        logger.info("SSNBuilderStage 完成(full模式): 共 %d 条边", len(edges))
        return CompleteSSN(edges=edges, edge_completion_mode=self._edge_completion_mode)

    @staticmethod
    def _add_known_intra_genome_edges(
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
                        if hit is not None:
                            edges.add(SSNEdge(
                                source=rep_id, target=member_id,
                                identity=hit.identity, coverage=hit.coverage,
                                bitscore=hit.bitscore,
                            ))

                if include_member_member:
                    members = sorted(local_cluster.member_ids)
                    for i, gene_a in enumerate(members):
                        for gene_b in members[i + 1:]:
                            hit = pair_lookup.get(frozenset((gene_a, gene_b)))
                            if hit is not None:
                                edges.add(SSNEdge(
                                    source=gene_a, target=gene_b,
                                    identity=hit.identity, coverage=hit.coverage,
                                    bitscore=hit.bitscore,
                                ))

    def _complete_cross_genome_member_edges(
        self,
        incremental_result: IncrementalResult,
        genes_by_id: Dict[str, Gene],
        edges: Set[SSNEdge],
    ) -> None:
        """full 模式下真正需要新比对的部分: 同一个全局簇内, 来自不同基因组的
        成员之间从未被比对过(IntraGenomeStage 只做单基因组内全对全比对,
        GlobalJoinStage 只比对代表), 按"全局簇内涉及的基因组两两分组"批量
        调用 align_query_vs_subject, 避免逐 pair 调用子进程。
        """
        global_lib = incremental_result.global_lib

        for cluster in global_lib.clusters:
            if not cluster.active:
                continue
            members_by_genome: Dict[str, List[str]] = {}
            for gene_id in cluster.member_ids:
                gene = genes_by_id.get(gene_id)
                if gene is None:
                    continue
                members_by_genome.setdefault(gene.genome_id, []).append(gene_id)

            genome_ids = sorted(members_by_genome.keys())
            if len(genome_ids) < 2:
                continue   # 单基因组簇, 同基因组内关系已在上一步复用完毕

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
                    query_genes = [genes_by_id[gid] for gid in members_by_genome[genome_a]]
                    subject_genes = [genes_by_id[gid] for gid in members_by_genome[genome_b]]
                    hits = self._aligner.align_query_vs_subject(query_genes, subject_genes)
                    for hit in hits:
                        edges.add(SSNEdge(
                            source=hit.query_id, target=hit.subject_id,
                            identity=hit.identity, coverage=hit.coverage,
                            bitscore=hit.bitscore,
                        ))

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
