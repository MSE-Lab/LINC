"""Pipeline: 串联 engine/ 五个 Stage 的编排层(第5步, 2026-10-02; 第8步
持久化/resume/append 扩展, 2026-10-02)。

设计依据: 与用户多轮讨论确认(见 docs/开发进度.md):
- 不单独新建 `incremental_engine.py`, "增量并入"的核心算法已完整实现在
  `engine/global_join_stage.py` 里, Stage 编排顺序这件事直接放进本文件即可。
- `Pipeline.run()` 只负责编排 5 个 Stage 的调用顺序 + 通过 `RepositoryPort`
  做每个 Stage 产出的存取, 不负责生成摘要 dict/写 json/打印等表现层逻辑
  (这些仍由 cli.py 负责, 以便未来非 CLI 调用方式 (如 Python API) 可以直接
  复用 `Pipeline`)。
- 本次不做"朴素全量比对 vs 增量结果等价性"的端到端集成测试(已与用户确认暂缓)。

第8步(2026-10-02, 与用户确认): `--resume`/`--append` 均依托 `RepositoryPort`
(实际为 `FileRepository`, 持久化到 `.linc_work/state/`)已持久化的步骤结果:
- `resume()`: 按 genomes->intra_results->join_result->complete_ssn->families 的
  顺序检查哪些 Stage 已持久化过, 从中断点接着跑完剩余 Stage;
- `append()`: 新基因组只跑 GenomeLoadStage+IntraGenomeStage, 再与历史
  `join_result` 接续运行 GlobalJoinStage(传入 `initial_global_lib` 等),
  SSNBuilderStage/FamilyRefineStage 对合并后的完整数据全量重跑
  (不做增量优化, 保证正确性优先, 与用户确认一致)。
"""
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .core.entities import (
    CompleteSSN,
    Genome,
    GenomeIntraResult,
    GenomeStatus,
    IncrementalResult,
    OrthologFamily,
)
from .engine.family_refine_stage import FamilyRefineStage
from .engine.genome_load_stage import GenomeLoadStage
from .engine.global_join_stage import GlobalJoinStage
from .engine.intra_genome_stage import IntraGenomeStage
from .engine.ssn_builder_stage import SSNBuilderStage
from .ports.storage_port import RepositoryPort


@dataclass
class PipelineConfig:
    """Pipeline 运行所需的全部算法参数, 收拢自 cli.py 的 `run` 命令选项。"""

    cutoff_i: float = 0.5
    cutoff_ii: float = 0.8
    cutoff_c: float = 0.0
    cov_mode: int = 0
    max_per_cluster: int = 5
    edge_completion_mode: str = "representative_member"
    mega_cluster_warning_threshold: Optional[int] = None
    inflation: float = 1.5
    edge_weight_mode: str = "identity_x_coverage"
    edge_weight_alpha: float = 0.5


@dataclass
class PipelineResult:
    """`Pipeline.run()` 的返回值, 汇总 5 个 Stage 的全部产出。"""

    genomes: Dict[str, Genome] = field(default_factory=dict)
    intra_results: Dict[str, GenomeIntraResult] = field(default_factory=dict)
    join_result: Optional[IncrementalResult] = None
    complete_ssn: Optional[CompleteSSN] = None
    families: List[OrthologFamily] = field(default_factory=list)


class Pipeline:
    """按顺序编排 GenomeLoadStage -> IntraGenomeStage -> GlobalJoinStage ->
    SSNBuilderStage -> FamilyRefineStage, 每个 Stage 跑完后通过
    `RepositoryPort` 保存其产出。
    """

    def __init__(
        self,
        aligner,
        clustering_for_intra,
        clustering_for_mcl,
        repository: RepositoryPort,
        config: PipelineConfig,
    ) -> None:
        self._aligner = aligner
        self._clustering_for_intra = clustering_for_intra
        self._clustering_for_mcl = clustering_for_mcl
        self._repository = repository
        self._config = config

    @staticmethod
    def _loaded_only(genomes: Dict[str, Genome]) -> Dict[str, Genome]:
        return {
            genome_id: genome
            for genome_id, genome in genomes.items()
            if genome.status == GenomeStatus.LOADED
        }

    def _run_intra(self, loaded_genomes: Dict[str, Genome]) -> Dict[str, GenomeIntraResult]:
        config = self._config
        intra_stage = IntraGenomeStage(
            aligner=self._aligner,
            clustering=self._clustering_for_intra,
            cutoff_i=config.cutoff_i,
            cutoff_ii=config.cutoff_ii,
            cutoff_c=config.cutoff_c,
            cov_mode=config.cov_mode,
            max_per_cluster=config.max_per_cluster,
        )
        return {
            genome_id: intra_stage.run(genome)
            for genome_id, genome in loaded_genomes.items()
        }

    def _run_join(
        self,
        loaded_genomes: Dict[str, Genome],
        intra_results: Dict[str, GenomeIntraResult],
        initial_global_lib=None,
        initial_ssn_edges=None,
        initial_cluster_counter: int = 0,
    ) -> IncrementalResult:
        config = self._config
        join_stage = GlobalJoinStage(
            aligner=self._aligner,
            cutoff_i=config.cutoff_i,
            cutoff_ii=config.cutoff_ii,
            cutoff_c=config.cutoff_c,
            cov_mode=config.cov_mode,
            max_per_cluster=config.max_per_cluster,
        )
        return join_stage.run(
            loaded_genomes,
            intra_results,
            initial_global_lib=initial_global_lib,
            initial_ssn_edges=initial_ssn_edges,
            initial_cluster_counter=initial_cluster_counter,
        )

    def _run_ssn(
        self,
        loaded_genomes: Dict[str, Genome],
        intra_results: Dict[str, GenomeIntraResult],
        join_result: IncrementalResult,
    ) -> CompleteSSN:
        config = self._config
        ssn_stage = SSNBuilderStage(
            aligner=self._aligner,
            edge_completion_mode=config.edge_completion_mode,
            mega_cluster_warning_threshold=config.mega_cluster_warning_threshold,
            cutoff_i=config.cutoff_i,
            cutoff_c=config.cutoff_c,
            cov_mode=config.cov_mode,
        )
        return ssn_stage.run(loaded_genomes, intra_results, join_result)

    def _run_family(
        self,
        complete_ssn: CompleteSSN,
        intra_results: Dict[str, GenomeIntraResult],
    ) -> List[OrthologFamily]:
        config = self._config
        family_stage = FamilyRefineStage(
            clustering=self._clustering_for_mcl,
            inflation=config.inflation,
            edge_weight_mode=config.edge_weight_mode,
            edge_weight_alpha=config.edge_weight_alpha,
        )
        # 传入全部基因id, 以便零边的孤立基因也能被显式补成图节点, 避免在
        # nx.connected_components() 中被静默忽略而从最终家族结果中消失
        # (2026-10-02 修复, 详见 family_refine_stage.py 模块docstring)。
        all_gene_ids = {
            gene_id
            for intra in intra_results.values()
            for cluster in intra.local_clusters
            for gene_id in cluster.member_ids
        }
        return family_stage.run(complete_ssn, all_gene_ids=all_gene_ids)

    def run(self, genome_files: List[str]) -> PipelineResult:
        """首次完整运行(无历史状态), 依次跑完 5 个 Stage 并逐一持久化。"""
        self._repository.save_input_signature(
            sorted(os.path.basename(f) for f in genome_files)
        )

        genomes = GenomeLoadStage().run(genome_files)
        self._repository.save_genomes(genomes)

        loaded_genomes = self._loaded_only(genomes)

        intra_results = self._run_intra(loaded_genomes)
        self._repository.save_intra_results(intra_results)

        join_result = self._run_join(loaded_genomes, intra_results)
        self._repository.save_join_result(join_result)

        complete_ssn = self._run_ssn(loaded_genomes, intra_results, join_result)
        self._repository.save_complete_ssn(complete_ssn)

        families = self._run_family(complete_ssn, intra_results)
        self._repository.save_families(families)

        return PipelineResult(
            genomes=genomes,
            intra_results=intra_results,
            join_result=join_result,
            complete_ssn=complete_ssn,
            families=families,
        )

    def resume(self, genome_files: List[str]) -> PipelineResult:
        """从上次中断处接着跑完剩余 Stage(2026-10-02, 第8步)。

        要求 `genome_files` 与历史 `-i` 签名一致(由 cli.py 校验后调用本方法,
        本方法内部不再重复校验, 专注于"从哪个 Stage 接着跑")。
        """
        repo = self._repository

        if repo.has_genomes():
            genomes = repo.load_genomes()
        else:
            genomes = GenomeLoadStage().run(genome_files)
            repo.save_genomes(genomes)

        loaded_genomes = self._loaded_only(genomes)

        if repo.has_intra_results():
            intra_results = repo.load_intra_results()
        else:
            intra_results = self._run_intra(loaded_genomes)
            repo.save_intra_results(intra_results)

        if repo.has_join_result():
            join_result = repo.load_join_result()
        else:
            join_result = self._run_join(loaded_genomes, intra_results)
            repo.save_join_result(join_result)

        if repo.has_complete_ssn():
            complete_ssn = repo.load_complete_ssn()
        else:
            complete_ssn = self._run_ssn(loaded_genomes, intra_results, join_result)
            repo.save_complete_ssn(complete_ssn)

        if repo.has_families():
            families = repo.load_families()
        else:
            families = self._run_family(complete_ssn, intra_results)
            repo.save_families(families)

        return PipelineResult(
            genomes=genomes,
            intra_results=intra_results,
            join_result=join_result,
            complete_ssn=complete_ssn,
            families=families,
        )

    def append(self, new_genome_files: List[str]) -> PipelineResult:
        """把 `new_genome_files` 追加到历史完整结果上(2026-10-02, 第8步)。

        要求历史状态是一次完整运行的结果(由 cli.py 在调用前校验 5 个 Stage
        均已持久化, 本方法假定这一前提已满足)。流程(与用户确认):
        - 新基因组只跑 GenomeLoadStage + IntraGenomeStage;
        - 与历史 `join_result` 接续运行 GlobalJoinStage(新基因组按 Order-D
          并入历史全局库, cluster_id 从历史计数器接续);
        - SSNBuilderStage/FamilyRefineStage 对"历史+新增"合并后的完整数据
          全量重跑(不做增量优化, 保证正确性优先)。
        """
        repo = self._repository

        historical_genomes = repo.load_genomes()
        historical_intra = repo.load_intra_results()
        historical_join = repo.load_join_result()

        existing_genome_ids = set(historical_genomes.keys())
        new_genomes = GenomeLoadStage().run(
            new_genome_files, existing_genome_ids=existing_genome_ids
        )
        new_loaded = self._loaded_only(new_genomes)

        new_intra_results = self._run_intra(new_loaded)

        merged_genomes = {**historical_genomes, **new_genomes}
        merged_intra_results = {**historical_intra, **new_intra_results}
        merged_loaded = self._loaded_only(merged_genomes)

        join_result = self._run_join(
            new_loaded,
            new_intra_results,
            initial_global_lib=historical_join.global_lib,
            initial_ssn_edges=historical_join.ssn_edges,
            initial_cluster_counter=historical_join.cluster_counter,
        )

        complete_ssn = self._run_ssn(merged_loaded, merged_intra_results, join_result)
        families = self._run_family(complete_ssn, merged_intra_results)

        # 注意: `-i` 对应的 input_signature 不随 append 变化(`-i` 始终代表
        # 历史基础数据集, 只有 `-a` 每次追加的新目录才变化, 2026-10-02 与
        # 用户确认), 因此这里不更新 input_signature, 只更新实际的 Stage 产出。
        repo.save_genomes(merged_genomes)
        repo.save_intra_results(merged_intra_results)
        repo.save_join_result(join_result)
        repo.save_complete_ssn(complete_ssn)
        repo.save_families(families)

        return PipelineResult(
            genomes=merged_genomes,
            intra_results=merged_intra_results,
            join_result=join_result,
            complete_ssn=complete_ssn,
            families=families,
        )
