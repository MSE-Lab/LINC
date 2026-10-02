"""Pipeline: 串联 engine/ 五个 Stage 的编排层(第5步, 2026-10-02)。

设计依据: 与用户多轮讨论确认(见 docs/开发进度.md):
- 不单独新建 `incremental_engine.py`, "增量并入"的核心算法已完整实现在
  `engine/global_join_stage.py` 里, Stage 编排顺序这件事直接放进本文件即可。
- `Pipeline.run()` 只负责编排 5 个 Stage 的调用顺序 + 通过 `RepositoryPort`
  做每个 Stage 产出的存取, 不负责生成摘要 dict/写 json/打印等表现层逻辑
  (这些仍由 cli.py 负责, 以便未来非 CLI 调用方式 (如 Python API) 可以直接
  复用 `Pipeline`)。
- 本次(第5步)不含 resume/--append, 只跑通"单次完整运行"主流程;
  `--append`/断点续跑留待第8步(持久化存储)一并设计实现。
- 本次不做"朴素全量比对 vs 增量结果等价性"的端到端集成测试(已与用户确认暂缓)。
"""
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

    def run(self, genome_files: List[str]) -> PipelineResult:
        config = self._config

        genomes = GenomeLoadStage().run(genome_files)
        self._repository.save_genomes(genomes)

        loaded_genomes = {
            genome_id: genome
            for genome_id, genome in genomes.items()
            if genome.status == GenomeStatus.LOADED
        }

        intra_stage = IntraGenomeStage(
            aligner=self._aligner,
            clustering=self._clustering_for_intra,
            cutoff_i=config.cutoff_i,
            cutoff_ii=config.cutoff_ii,
            cutoff_c=config.cutoff_c,
            cov_mode=config.cov_mode,
            max_per_cluster=config.max_per_cluster,
        )
        intra_results = {
            genome_id: intra_stage.run(genome)
            for genome_id, genome in loaded_genomes.items()
        }
        self._repository.save_intra_results(intra_results)

        join_stage = GlobalJoinStage(
            aligner=self._aligner,
            cutoff_i=config.cutoff_i,
            cutoff_ii=config.cutoff_ii,
            cutoff_c=config.cutoff_c,
            cov_mode=config.cov_mode,
            max_per_cluster=config.max_per_cluster,
        )
        join_result = join_stage.run(loaded_genomes, intra_results)
        self._repository.save_join_result(join_result)

        ssn_stage = SSNBuilderStage(
            aligner=self._aligner,
            edge_completion_mode=config.edge_completion_mode,
            mega_cluster_warning_threshold=config.mega_cluster_warning_threshold,
        )
        complete_ssn = ssn_stage.run(loaded_genomes, intra_results, join_result)
        self._repository.save_complete_ssn(complete_ssn)

        family_stage = FamilyRefineStage(
            clustering=self._clustering_for_mcl,
            inflation=config.inflation,
            edge_weight_mode=config.edge_weight_mode,
            edge_weight_alpha=config.edge_weight_alpha,
        )
        families = family_stage.run(complete_ssn)
        self._repository.save_families(families)

        return PipelineResult(
            genomes=genomes,
            intra_results=intra_results,
            join_result=join_result,
            complete_ssn=complete_ssn,
            families=families,
        )
