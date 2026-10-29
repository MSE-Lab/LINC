"""只读加载 `linc run` 的 state, 供 `linc pan` 报告使用。

设计依据: docs/LINC_结果描述增强开发计划.md 第 2/3 节。本模块只负责"读",
不重跑任何 Stage,不修改任何 state 文件;缺失/损坏/版本不兼容统一抛
`ReportInputError`, 由 CLI 层转成清晰的报错信息。
"""
import logging
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from ..adapters.storage.file_repo import (
    FileRepository,
    SchemaVersionMismatchError,
    StateCorruptedError,
)
from ..core.entities import (
    CompleteSSN,
    Genome,
    GenomeIntraResult,
    IncrementalResult,
    OrthologFamily,
)

logger = logging.getLogger(__name__)


class ReportInputError(RuntimeError):
    """报告输入缺失/损坏/schema 版本不兼容。

    与 `FileRepository` 的 `StateCorruptedError` 区分: 后者表示单个 state 文件
    损坏; 本异常表示"这次报告无法从给定输出目录读到完整、兼容的 state"。
    """


@dataclass
class LoadedState:
    """一次报告所需的全部只读 state, 以及内部 ID -> 原始 ID 的映射。

    内部 `gene_id` 形如 "G00001|ORF_00001"(genome_id 是内部短编号), 对外
    输出(`Orthogroups.txt` / `gene_genome.tsv`)使用原始序列 ID(即
    `Gene.original_header`)与原始基因组名(`Genome.source_file_name`)。报告
    产出统一用原始 ID, 故在加载阶段一次性建好映射。
    """

    output_dir: str
    genomes: Dict[str, Genome]
    families: List[OrthologFamily]
    complete_ssn: Optional[CompleteSSN]
    intra_results: Dict[str, GenomeIntraResult]
    join_result: Optional[IncrementalResult]
    config: Dict[str, Any]
    # 内部 gene_id -> 原始序列 ID(如 "GCA_000005845.2|ORF_00001")。
    gene_id_to_original: Dict[str, str] = field(default_factory=dict)
    # 内部 genome_id -> 原始基因组名(如 "GCA_000005845.2")。
    genome_id_to_name: Dict[str, str] = field(default_factory=dict)


def state_dir_of(output_dir: str) -> str:
    """返回 `<output_dir>/.linc_work/state` 路径。"""
    return os.path.join(output_dir, ".linc_work", "state")


def load_state(output_dir: str) -> LoadedState:
    """只读加载 `<output_dir>/.linc_work/state` 下的全部持久化数据。

    必需文件: `genomes.pkl` 与 `families.pkl`(泛基因组与家族统计的最小依赖)。
    其余(`complete_ssn.pkl` / `intra_results.pkl` / `join_result.pkl` /
    `config.pkl`)按需加载, 缺失只告警、不失败, 相关统计由下游自行跳过。

    Raises:
        ReportInputError: state 目录不存在、必需文件缺失、或 schema 版本不兼容。
    """
    state_dir = state_dir_of(output_dir)
    if not os.path.isdir(state_dir):
        raise ReportInputError(
            f"未找到 state 目录: {state_dir}; "
            f"请确认 `{output_dir}` 是 `linc run -o` 的输出目录且已完成一次完整运行。"
        )

    repo = FileRepository(state_dir)

    # 必需文件检查(先于加载, 避免在缺文件时触发无意义的 schema/损坏报错)。
    missing = []
    if not repo.has_genomes():
        missing.append("genomes.pkl")
    if not repo.has_families():
        missing.append("families.pkl")
    if missing:
        raise ReportInputError(
            f"state 缺少必需文件: {', '.join(missing)}; "
            f"请确认 `linc run` 已完整跑完(至少产出 families.pkl 与 genomes.pkl)。"
        )

    try:
        genomes = repo.load_genomes()
        families = repo.load_families()
        complete_ssn = repo.load_complete_ssn()
        intra_results = repo.load_intra_results()
        join_result = repo.load_join_result()
        config = repo.load_config() or {}
    except SchemaVersionMismatchError as exc:
        raise ReportInputError(
            f"state schema 版本不兼容: {exc}; 请用当前代码版本重新完整运行 `linc run` 后再生成报告。"
        ) from exc
    except StateCorruptedError as exc:
        raise ReportInputError(f"state 文件损坏, 无法生成报告: {exc}") from exc

    # 可选文件缺失时只告警, 不阻断(下游据此跳过对应统计)。
    if complete_ssn is None:
        logger.warning("state 缺少 complete_ssn.pkl, SSN 相关统计将被跳过")
    if not intra_results:
        logger.warning("state 缺少 intra_results.pkl, 局部簇/代表统计将被跳过")
    if join_result is None:
        logger.warning("state 缺少 join_result.pkl, 全局簇统计将被跳过")
    if not config:
        logger.warning("state 缺少 config.pkl, 参数信息将缺失")

    gene_id_to_original: Dict[str, str] = {}
    genome_id_to_name: Dict[str, str] = {}
    for genome in genomes.values():
        name = genome.source_file_name or genome.genome_id
        genome_id_to_name[genome.genome_id] = name
        for gene in genome.genes:
            gene_id_to_original[gene.gene_id] = gene.original_header or gene.gene_id

    return LoadedState(
        output_dir=output_dir,
        genomes=genomes,
        families=families,
        complete_ssn=complete_ssn,
        intra_results=intra_results,
        join_result=join_result,
        config=config,
        gene_id_to_original=gene_id_to_original,
        genome_id_to_name=genome_id_to_name,
    )
