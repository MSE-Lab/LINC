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
import logging
import os
from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional

from .core.entities import (
    CompleteSSN,
    Genome,
    GenomeIntraResult,
    GenomeStatus,
    IncrementalResult,
    OrthologFamily,
)
from .engine.family_refine_stage import FamilyRefineStage
from .engine.genome_load_stage import GenomeLoadStage, build_input_signature
from .engine.global_join_stage import GlobalJoinStage
from .engine.intra_genome_stage import IntraGenomeStage
from .engine.ssn_builder_stage import SSNBuilderStage
from .ports.storage_port import RepositoryPort

logger = logging.getLogger(__name__)


class IndexMissingError(RuntimeError):
    """`--resume` 时发现历史 state 引用的比对索引已缺失(第二轮漏洞5,
    2026-10-03 与用户确认)。

    `GenomeIntraResult.genome_index` 只是 `AlignerPort.build_index` 返回的
    不透明 handle(如 mmseqs/diamond 的 DB 路径), state pkl 只持久化这个
    handle 本身, 不代表其指向的索引文件仍然存在。若用户手动删除了
    `.linc_work/mmseqs_db`(或 diamond 对应目录)而 state 还在, 直接复用该
    handle 发起比对会在子进程层面崩溃、报错不友好, 因此在 `resume()` 中
    提前校验并抛出本异常, 由调用方(cli.py)给出清晰的修复指引。
    """


class ConfigMismatchError(RuntimeError):
    """`--resume`/`--append` 时发现本次算法参数与历史运行不一致(工程审查
    问题1, 2026-10-03 与用户确认)。

    历史上 `resume()`/`append()` 会直接复用历史 state(用旧参数算出的
    `intra_results`/`join_result` 等)并用本次(可能不同的)参数继续跑后续
    Stage, 完全不报错, 产出新旧参数混搭的不一致结果。本异常在检测到参数
    不一致且调用方未显式传入 `force=True` 时抛出, 列出全部不一致的字段
    (历史值 vs 本次值), 要求用户确认后用 `--force` 显式覆盖继续, 或保持
    参数一致/重新完整运行。
    """


@dataclass
class PipelineConfig:
    """Pipeline 运行所需的全部算法参数, 收拢自 cli.py 的 `run` 命令选项。"""

    cutoff_i: float = 0.5
    cutoff_ii: float = 0.8
    cutoff_c: float = 0.5
    cov_mode: int = 0
    max_per_cluster: int = 5
    edge_completion_mode: str = "representative_member"
    mega_cluster_warning_threshold: Optional[int] = None
    inflation: float = 1.5
    edge_weight_mode: str = "identity_x_coverage"
    edge_weight_alpha: float = 0.5
    # 数据管理优化第3项(2026-10-03 与用户确认): 默认 False 即各 Stage 用完
    # search()/mcl 产生的中间文件即删, 避免 `.linc_work` 无限累积; 设为
    # True 保留全部中间产物以便调试(索引目录/state pkl 不受影响, 始终保留)。
    keep_intermediate: bool = False
    # 工程审查问题9(2026-10-03 与用户确认): 传给 mcl 的 `-te` 线程数(合并
    # 为全图一次mcl调用后, 通过 mcl 自身多线程利用多核); 0 表示不显式指定。
    # 纯性能参数, 不影响算法结果, 已在 `_check_config_consistency` 中排除,
    # 不参与 resume/append 的参数一致性校验(修改线程数不应触发强制重算)。
    threads: int = 0


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

    def _check_genome_indexes(self, intra_results: Dict[str, GenomeIntraResult]) -> None:
        """校验 `intra_results` 中每个非空 `genome_index` 对应的索引仍可用。

        仅在 `--resume` 从历史 state 加载 `intra_results` 时调用(第二轮
        漏洞5, 2026-10-03 与用户确认): 后续 Stage(GlobalJoinStage/
        SSNBuilderStage)会直接复用这些 handle 发起 `search`, 若对应索引
        已被用户手动删除, 应在此处提前失败并给出清晰指引, 而不是让比对
        工具子进程崩溃。
        """
        missing = [
            genome_id
            for genome_id, intra in intra_results.items()
            if intra.genome_index is not None
            and not self._aligner.index_exists(intra.genome_index)
        ]
        if missing:
            raise IndexMissingError(
                "--resume 失败: 以下基因组的比对索引已缺失(可能被手动删除), "
                f"而历史 state 仍引用它们: {sorted(missing)}; "
                "索引与 state 必须配套存在才能安全续跑, 请删除 -o 目录下的 "
                ".linc_work 后重新执行一次完整的 `linc run`。"
            )

    def _check_config_consistency(self, force: bool, action: str) -> bool:
        """校验本次算法参数(`self._config`)与历史持久化的参数是否一致
        (工程审查问题1, 2026-10-03 与用户确认; force 风险修复, 2026-10-03
        与用户再次确认)。

        历史 state 不存在 `config.pkl`(本次改动之前写入的旧版 state)时,
        视为"无法校验", 只打印 warning 并放行(与 `input_signature` 旧格式
        兼容处理的既有原则一致), 不阻断 resume/append, 返回 `False`
        (不视为"发生了不一致", 不触发调用方的强制重算逻辑)。

        返回值 `mismatched: bool` 表示"本次是否检测到参数不一致(无论是否
        因 --force 被放行)"。调用方(`resume`/`append`)必须在 `mismatched`
        为 `True` 时, 忽略任何已持久化的 intra_results/join_result/
        complete_ssn/families 缓存, 强制用本次(新)参数从 IntraGenomeStage
        开始对全部历史基因组重新计算——否则 --force 只是放行了检查、却仍
        然复用旧参数算出的缓存, 产出会在同一次结果里悄悄混用新旧参数
        (风险1)。force 覆盖事件会写入 `append_audit_log` 持久化留痕, 不会
        因为随后 `save_config()` 更新基线而被悄悄抹掉(风险2)。
        """
        historical_config = self._repository.load_config()
        if historical_config is None:
            logger.warning(
                "历史 state 中未找到算法参数记录(config.pkl), 可能是本次改动"
                "之前写入的旧版 state, 无法校验参数一致性, 已跳过该项校验。"
            )
            return False
        current_config = asdict(self._config)
        # 工程审查问题9(2026-10-03 与用户确认): `threads` 是纯性能参数
        # (传给 mcl 的 `-te` 线程数), 不影响算法结果, 修改它不应被视为
        # "参数不一致"而触发 ConfigMismatchError/强制重算。
        non_semantic_fields = {"threads"}
        mismatched = {
            key: (historical_config.get(key), current_config[key])
            for key in current_config
            if key not in non_semantic_fields
            and historical_config.get(key) != current_config[key]
        }
        if not mismatched:
            return False
        detail = "; ".join(
            f"{key}: 历史={hist!r}, 本次={cur!r}"
            for key, (hist, cur) in sorted(mismatched.items())
        )
        if force:
            logger.warning(
                "检测到算法参数与历史运行不一致, 但已指定 --force, 将忽略"
                "已持久化的 intra_results/join_result/complete_ssn/families "
                "缓存, 用本次新参数对全部历史基因组从 IntraGenomeStage 开始"
                "强制重算(等同于一次接近完整的重跑, 比普通 resume/append 慢)"
                ", 以避免新旧参数混搭: %s", detail,
            )
            self._repository.append_audit_log(
                {
                    "timestamp": datetime.now().isoformat(),
                    "event": "force_config_override",
                    "action": action,
                    "mismatched_fields": {
                        key: {"historical": hist, "current": cur}
                        for key, (hist, cur) in mismatched.items()
                    },
                }
            )
            return True
        raise ConfigMismatchError(
            f"--resume/--append 失败: 本次算法参数与历史运行不一致: {detail}; "
            "直接续跑会产出新旧参数混搭的不一致结果。若确认要覆盖继续, 请加 "
            "--force 显式确认(会强制从 IntraGenomeStage 开始对全部历史基因组"
            "重算, 风险自负); 否则请保持参数与历史一致, 或删除 -o 目录下的 "
            ".linc_work 后重新执行一次完整的 `linc run`。"
        )

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
            cleanup_intermediate=not config.keep_intermediate,
            edge_completion_mode=config.edge_completion_mode,
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
            cleanup_intermediate=not config.keep_intermediate,
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
            cleanup_intermediate=not config.keep_intermediate,
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
            cleanup_intermediate=not config.keep_intermediate,
            threads=config.threads,
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
            build_input_signature(genome_files)
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

        # 工程审查问题1(2026-10-03 与用户确认; 2026-10-03 再次确认挪到末尾):
        # 持久化本次运行使用的全部算法参数, 供后续 --resume/--append 校验
        # 一致性。必须等全部 Stage 产出都已成功落盘之后才保存, 否则一旦中途
        # 崩溃, config.pkl 会被"提前"更新为本次参数, 而部分 Stage 产出仍是
        # 旧参数(甚至是上一次运行遗留)的残留文件, 导致下次 resume 时一致性
        # 校验对比的是"已经被更新过的历史 config", 无法察觉这个不一致。
        self._repository.save_config(asdict(self._config))

        return PipelineResult(
            genomes=genomes,
            intra_results=intra_results,
            join_result=join_result,
            complete_ssn=complete_ssn,
            families=families,
        )

    def resume(self, genome_files: List[str], force: bool = False) -> PipelineResult:
        """从上次中断处接着跑完剩余 Stage(2026-10-02, 第8步)。

        要求 `genome_files` 与历史 `-i` 签名一致(由 cli.py 校验后调用本方法,
        本方法内部不再重复校验, 专注于"从哪个 Stage 接着跑")。

        `force`: 若本次算法参数与历史持久化的参数不一致, 默认抛出
        `ConfigMismatchError`; 传入 `force=True` 则记 warning + 持久化审计
        记录后继续执行, 并强制忽略所有已持久化的 Stage 缓存、从
        IntraGenomeStage 开始用新参数对全部基因组重算(工程审查问题1,
        2026-10-03 与用户确认; force 风险修复, 2026-10-03 与用户再次确认)。
        """
        repo = self._repository
        forced_override = self._check_config_consistency(force=force, action="resume")

        if repo.has_genomes():
            genomes = repo.load_genomes()
        else:
            genomes = GenomeLoadStage().run(genome_files)
            repo.save_genomes(genomes)

        loaded_genomes = self._loaded_only(genomes)

        # force 风险修复(2026-10-03 与用户确认, 风险1): 一旦检测到参数不
        # 一致且用户显式 --force, 不能再信任任何已持久化的 intra_results/
        # join_result/complete_ssn/families 缓存(它们可能是用旧参数算出的)
        # ——必须强制忽略缓存, 从 IntraGenomeStage 开始用本次新参数对全部
        # 历史基因组重新计算, 避免同一次产出里混用新旧参数。
        if forced_override or not repo.has_intra_results():
            intra_results = self._run_intra(loaded_genomes)
            repo.save_intra_results(intra_results)
        else:
            intra_results = repo.load_intra_results()
            self._check_genome_indexes(intra_results)

        if forced_override or not repo.has_join_result():
            join_result = self._run_join(loaded_genomes, intra_results)
            repo.save_join_result(join_result)
        else:
            join_result = repo.load_join_result()

        if forced_override or not repo.has_complete_ssn():
            complete_ssn = self._run_ssn(loaded_genomes, intra_results, join_result)
            repo.save_complete_ssn(complete_ssn)
        else:
            complete_ssn = repo.load_complete_ssn()

        if forced_override or not repo.has_families():
            families = self._run_family(complete_ssn, intra_results)
            repo.save_families(families)
        else:
            families = repo.load_families()

        # 工程审查问题1(2026-10-03 与用户确认; 2026-10-03 再次确认挪到末尾):
        # 待全部 Stage 产出都已成功落盘(或确认复用历史缓存无误)后才保存本次
        # 参数, 避免中途崩溃导致 config.pkl 与实际落盘的 Stage 结果不一致。
        repo.save_config(asdict(self._config))

        return PipelineResult(
            genomes=genomes,
            intra_results=intra_results,
            join_result=join_result,
            complete_ssn=complete_ssn,
            families=families,
        )

    def append(self, new_genome_files: List[str], force: bool = False) -> PipelineResult:
        """把 `new_genome_files` 追加到历史完整结果上(2026-10-02, 第8步)。

        要求历史状态是一次完整运行的结果(由 cli.py 在调用前校验 5 个 Stage
        均已持久化, 本方法假定这一前提已满足)。流程(与用户确认):
        - 新基因组只跑 GenomeLoadStage + IntraGenomeStage;
        - 与历史 `join_result` 接续运行 GlobalJoinStage(新基因组按 Order-D
          并入历史全局库, cluster_id 从历史计数器接续);
        - SSNBuilderStage/FamilyRefineStage 对"历史+新增"合并后的完整数据
          全量重跑(不做增量优化, 保证正确性优先)。

        `force`: 若本次算法参数与历史持久化的参数不一致, 默认抛出
        `ConfigMismatchError`; 传入 `force=True` 则记 warning + 持久化审计
        记录后继续执行, 并强制对**全部历史基因组**用新参数重算
        IntraGenomeStage、放弃增量并入优化、从头对全部基因组重跑
        GlobalJoinStage(工程审查问题1, 2026-10-03 与用户确认; force 风险
        修复, 2026-10-03 与用户再次确认)。
        """
        repo = self._repository
        forced_override = self._check_config_consistency(force=force, action="append")

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
        merged_loaded = self._loaded_only(merged_genomes)

        # force 风险修复(2026-10-03 与用户确认, 风险1): append 场景下
        # `historical_intra`/`historical_join` 是用"历史参数"算出的, 若
        # 本次参数不一致且用户 --force 放行, 不能让历史基因组继续复用旧
        # 参数产出——必须用新参数对**全部历史基因组**重新跑 IntraGenomeStage,
        # 并放弃"接续历史 join_result 增量并入"的优化路径, 改为对
        # "历史+新增"全部基因组从头跑一次 GlobalJoinStage, 避免新旧参数在
        # 同一次产出里混搭。
        if forced_override:
            historical_loaded = self._loaded_only(historical_genomes)
            historical_intra = self._run_intra(historical_loaded)
            merged_intra_results = {**historical_intra, **new_intra_results}
            join_result = self._run_join(merged_loaded, merged_intra_results)
        else:
            merged_intra_results = {**historical_intra, **new_intra_results}
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

        # 数据管理优化第1项(2026-10-03 与用户确认): append 审计轨迹——记录
        # 本次追加的原始文件/目录、实际新增的 genome_id 列表。
        # `skipped_genome_ids` 字段预留给未来可能出现的"静默跳过重复"场景,
        # 当前 `GenomeLoadStage._check_genome_id_conflicts` 遇到重复
        # genome_id 会直接抛异常中止, 不存在静默跳过, 因此该字段目前恒为空。
        repo.append_audit_log(
            {
                "timestamp": datetime.now().isoformat(),
                "new_genome_files": list(new_genome_files),
                "new_genome_ids": sorted(new_genomes.keys()),
                "skipped_genome_ids": [],
            }
        )

        # 工程审查问题1(2026-10-03 与用户确认; 2026-10-03 再次确认挪到末尾):
        # 待全部 Stage 产出都已成功落盘后才保存本次参数, 避免中途崩溃导致
        # config.pkl 与实际落盘的 Stage 结果不一致。
        repo.save_config(asdict(self._config))

        return PipelineResult(
            genomes=merged_genomes,
            intra_results=merged_intra_results,
            join_result=join_result,
            complete_ssn=complete_ssn,
            families=families,
        )
