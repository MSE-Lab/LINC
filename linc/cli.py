"""LINC 命令行入口。

第5步(2026-10-02, 与用户确认): CLI 不再手工串联 5 个 Stage, 而是构造
`pipeline.Pipeline` 并调用其 `run()` 一次性完成全部编排(GenomeLoadStage ->
IntraGenomeStage -> GlobalJoinStage -> SSNBuilderStage -> FamilyRefineStage);
`cli.py` 只负责参数解析/校验, 以及把 `Pipeline.run()` 返回的 `PipelineResult`
转换为各阶段的摘要 json/Orthogroups.txt 并打印, 不再包含任何 Stage 编排逻辑
本身(便于未来非 CLI 调用方式直接复用 `Pipeline`)。CLI 接口形态本身
(`linc run -i <目录> -o <输出目录>`)不受影响。

IntraGenomeStage/GlobalJoinStage 的比对后端可通过 `--aligner` 选择(二者共用
同一个 aligner 实例, 以便 GlobalJoinStage 复用 IntraGenomeStage 阶段已建好的
持久化索引, 避免重复建库):
- "mmseqs"(默认): 调用系统级 mmseqs2(需要已安装并在 PATH 中, 推荐通过 conda
  `linc` 环境提供), 使用真实序列比对结果。
- "diamond": 调用系统级 diamond(需要已安装并在 PATH 中, 推荐通过 conda `linc`
  环境提供), 使用真实序列比对结果。

MockAligner(假比对实现, 不依赖外部二进制程序)仅用于单元测试内部构造
夹具, 不作为 CLI 可选项暴露(与用户确认, 2026-10-02: 避免用户误以为
默认在用真实比对, 生产环境运行必须显式使用真实比对后端)。
"""
import json
import logging
import os
import sys
from typing import Optional

import typer

from .adapters.aligners.mmseqs_adapter import MMseqsAligner, MMseqsConfig
from .adapters.aligners.diamond_adapter import DiamondAligner, DiamondConfig
from .adapters.clustering.mcl_adapter import MCLClustering, MCLConfig
from .adapters.clustering.union_find_adapter import UnionFindClustering
from .core.entities import (
    VALID_COV_MODES,
    VALID_EDGE_COMPLETION_MODES,
    VALID_EDGE_WEIGHT_MODES,
    GenomeStatus,
)
from .engine.genome_load_stage import (
    GenomeIdConflictError,
    build_input_signature,
    discover_fasta_files,
)
from .adapters.storage.file_repo import FileRepository, StateCorruptedError
from .pipeline import ConfigMismatchError, IndexMissingError, Pipeline, PipelineConfig
from .utils.lock import PipelineLockedError, acquire_pipeline_lock

app = typer.Typer(help="LINC: Linear incremental genome comparison that builds a sequence similarity network (SSN) and refines protein families via MCL.")


@app.callback()
def main() -> None:
    """LINC command-line tool."""
    pass


logger = logging.getLogger("linc")


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        stream=sys.stdout,
    )


def _build_clustering(output_dir: str):
    """构造 MCLClustering 实例, 独立成函数以便测试通过 monkeypatch 注入假的
    聚类实现(与 `_build_aligner` 同样的理由: MCLClustering 要求本机已安装
    系统级 `mcl`, 单测环境不应该强依赖该外部二进制, 2026-10-02)。
    """
    mcl_work_dir = os.path.join(output_dir, ".linc_work", "mcl")
    return MCLClustering(MCLConfig(work_dir=mcl_work_dir))


def _build_aligner(aligner_backend: str, threads: int, output_dir: str):
    """根据 --aligner 取值构造真实比对后端实例。

    独立成函数以便测试通过 monkeypatch 注入 MockAligner(2026-10-02 与用户
    确认: MockAligner 仅供单测内部构造夹具使用, 不作为 CLI 可选项暴露;
    test_cli.py 对这类集成测试改为 monkeypatch 本函数, 而不是新增隐藏的
    CLI 后门参数)。取值非法时返回 None, 由调用方负责报错退出。
    """
    if aligner_backend == "mmseqs":
        work_dir = os.path.join(output_dir, ".linc_work", "mmseqs_db")
        return MMseqsAligner(MMseqsConfig(threads=threads, work_dir=work_dir))
    if aligner_backend == "diamond":
        work_dir = os.path.join(output_dir, ".linc_work", "diamond_db")
        return DiamondAligner(DiamondConfig(threads=threads, work_dir=work_dir))
    return None


def _validate_cutoff_range(option_name: str, value: float) -> None:
    """校验阈值类参数落在 [0, 1] 区间内, 否则直接退出并提示。"""
    if not 0.0 <= value <= 1.0:
        typer.secho(
            f"{option_name} 必须在 [0, 1] 范围内, 当前为: {value}",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)


def _build_load_summary(genomes: dict) -> dict:
    loaded = [g for g in genomes.values() if g.status == GenomeStatus.LOADED]
    failed = [g for g in genomes.values() if g.status == GenomeStatus.FAILED]
    return {
        "total_genomes": len(genomes),
        "loaded_count": len(loaded),
        "failed_count": len(failed),
        "failed_genome_ids": [g.genome_id for g in failed],
        "genomes": {
            genome_id: {
                "status": genome.status.value if genome.status else None,
                "gene_count": len(genome.genes),
                "file_path": genome.file_path,
                # genome_id 为内部短编号(方案B, 2026-10-03 与用户确认), 此处
                # 附带原始文件名, 便于用户把 genome_id 对照回原始输入。
                "source_file_name": genome.source_file_name,
            }
            for genome_id, genome in genomes.items()
        },
    }


def _build_join_summary(result) -> dict:
    return {
        "processing_order": result.processing_order,
        "n_clusters_final": result.n_clusters_final,
        "n_representatives_final": result.n_representatives_final,
        "n_ssn_edges": len(result.ssn_edges),
        "clusters": [
            {
                "cluster_id": cluster.cluster_id,
                "n_representatives": len(cluster.representative_ids),
                "n_members": len(cluster.member_ids),
            }
            for cluster in result.global_lib.clusters
            if cluster.active
        ],
    }


def _build_ssn_summary(complete_ssn) -> dict:
    return {
        "edge_completion_mode": complete_ssn.edge_completion_mode,
        "n_edges": len(complete_ssn.edges),
    }


def _build_family_summary(families, inflation: float, edge_weight_mode: str) -> dict:
    return {
        "inflation": inflation,
        "edge_weight_mode": edge_weight_mode,
        "n_families": len(families),
        "families": [
            {
                "family_id": family.family_id,
                "component_id": family.component_id,
                "size": family.size,
                "edge_density": round(family.edge_density, 4),
            }
            for family in families
        ],
    }


def _build_gene_id_to_original(genomes: dict) -> dict:
    """构建 gene_id -> 原始序列 ID 的映射, 用于将 Orthogroups.txt 中的
    内部 gene_id("{genome_id}|{local_id}")回填为用户可识别的原始 ID
    (2026-10-03 与用户确认"方案B": 原始 ID 不再因内部重编码而丢失)。

    同一基因组内若存在重复的 original_header, 按出现顺序加 "#2"/"#3"...
    后缀消歧, 保证映射后的值仍然全局唯一(与 GenomeLoadStage 解析阶段的
    重复 warning 配套)。若 original_header 缺失(理论上不应出现, 除非解析
    异常), 回退使用内部 gene_id 本身, 避免丢失成员。
    """
    mapping: dict = {}
    for genome in genomes.values():
        seen_count: dict = {}
        for gene in genome.genes:
            header = gene.original_header or gene.gene_id
            seen_count[header] = seen_count.get(header, 0) + 1
            if seen_count[header] > 1:
                header = f"{header}#{seen_count[header]}"
            mapping[gene.gene_id] = header
    return mapping


def _write_orthogroups_txt(families, gene_id_to_original: dict, path: str) -> None:
    """输出 Orthogroups.txt(一行一个 ortholog group): 第一列是 OG 编号,
    第二列是该家族全部成员的原始序列 ID(不再是 LINC 内部重编码后的
    gene_id, 2026-10-03 与用户确认修复), 用空格分隔(整体格式与用户
    确认, 2026-10-02)。
    """
    with open(path, "w", encoding="utf-8") as fh:
        for family in families:
            members = [
                gene_id_to_original.get(gene_id, gene_id)
                for gene_id in sorted(family.member_genes)
            ]
            fh.write(f"{family.family_id}\t{' '.join(members)}\n")


def _write_gene_genome_tsv(genomes: dict, path: str) -> None:
    """输出 gene_genome.tsv: 每行为"原始序列ID\t基因组名"(基因组名取
    Genome.source_file_name), 便于用户把任意序列对照回其所属的原始
    基因组文件, 不依赖内部 gene_id/genome_id 编码。
    """
    with open(path, "w", encoding="utf-8") as fh:
        for genome in genomes.values():
            for gene in genome.genes:
                header = gene.original_header or gene.gene_id
                fh.write(f"{header}\t{genome.source_file_name}\n")


def _load_gene_genome_map(path: str) -> dict:
    """解析 gene_genome.tsv(每行"原始序列ID\t基因组名"), 返回
    {原始序列ID: 基因组名} 映射。"""
    mapping = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                typer.secho(
                    f"gene_genome.tsv 第 {line_no} 行格式异常(应为两列, "
                    f"实际为 {len(parts)} 列): {line!r}",
                    fg=typer.colors.RED, err=True,
                )
                raise typer.Exit(code=1)
            gene_original_id, genome_name = parts
            mapping[gene_original_id] = genome_name
    return mapping


def _parse_orthogroups_txt(path: str) -> list:
    """解析 Orthogroups.txt(每行"family_id\\tmember1 member2 ..."),
    返回 [(family_id, [原始序列ID, ...]), ...]。"""
    families = []
    with open(path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                typer.secho(
                    f"Orthogroups.txt 第 {line_no} 行格式异常(应为两列, "
                    f"实际为 {len(parts)} 列): {line!r}",
                    fg=typer.colors.RED, err=True,
                )
                raise typer.Exit(code=1)
            family_id, members_str = parts
            members = members_str.split()
            families.append((family_id, members))
    return families


def _classify_family(n_genomes_present: int, n_total_genomes: int, n_members: int) -> str:
    """按用户确认的规则对一个直系同源家族(OG)分类:
    - 家族总成员数为 1(即只出现在 1 个基因组且该基因组内只有 1 条基因,
      无任何同源序列) -> singleton, 优先于下面的比例判断。
    - 否则按"覆盖的基因组数 / 总基因组数"的比例:
      >=99% -> core | [95%, 99%) -> soft-core |
      [15%, 95%) -> shell | <15% -> cloud。
    """
    if n_members == 1:
        return "singleton"
    ratio = n_genomes_present / n_total_genomes
    if ratio >= 0.99:
        return "core"
    if ratio >= 0.95:
        return "soft-core"
    if ratio >= 0.15:
        return "shell"
    return "cloud"


_PAN_CATEGORIES = ["core", "soft-core", "shell", "cloud", "singleton"]


def _build_pan_genome_data(families: list, gene_to_genome: dict, genome_names: list) -> dict:
    """基于已解析的 Orthogroups 与 gene->genome 映射, 一次性计算:
    - presence_absence: {family_id: {genome_name: copy_number}}
    - family_category: {family_id: category}
    - genome_breakdown: {genome_name: {category: gene_count}}
    """
    n_total_genomes = len(genome_names)
    presence_absence = {}
    family_category = {}
    genome_breakdown = {
        genome_name: {category: 0 for category in _PAN_CATEGORIES}
        for genome_name in genome_names
    }

    for family_id, members in families:
        copy_counts = {}
        unknown_members = []
        for original_id in members:
            genome_name = gene_to_genome.get(original_id)
            if genome_name is None:
                unknown_members.append(original_id)
                continue
            copy_counts[genome_name] = copy_counts.get(genome_name, 0) + 1
        if unknown_members:
            typer.secho(
                f"警告: Orthogroups.txt 中家族 {family_id} 的以下序列 ID "
                f"在 gene_genome.tsv 中未找到对应基因组, 已忽略: "
                f"{unknown_members}",
                fg=typer.colors.YELLOW, err=True,
            )

        presence_absence[family_id] = copy_counts
        n_genomes_present = len(copy_counts)
        n_members = sum(copy_counts.values())
        category = _classify_family(n_genomes_present, n_total_genomes, n_members)
        family_category[family_id] = category

        for genome_name, count in copy_counts.items():
            genome_breakdown[genome_name][category] += count

    return {
        "presence_absence": presence_absence,
        "family_category": family_category,
        "genome_breakdown": genome_breakdown,
    }


def _write_presence_absence_matrix(
    presence_absence: dict, genome_names: list, path: str
) -> None:
    """输出 presence_absence_matrix.tsv: 行=family_id, 列=各基因组名,
    值=该家族在该基因组中的基因拷贝数(整数, 0 表示缺失)。"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("Orthogroup\t" + "\t".join(genome_names) + "\n")
        for family_id in sorted(presence_absence.keys()):
            copy_counts = presence_absence[family_id]
            row = [str(copy_counts.get(genome_name, 0)) for genome_name in genome_names]
            fh.write(f"{family_id}\t" + "\t".join(row) + "\n")


def _write_pan_genome_summary(family_category: dict, path: str) -> None:
    """输出 pan_genome_summary.tsv: 按分类列出 OG 数量及占比。"""
    total = len(family_category)
    counts = {category: 0 for category in _PAN_CATEGORIES}
    for category in family_category.values():
        counts[category] += 1
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("Category\tN_orthogroups\tPercentage\n")
        for category in _PAN_CATEGORIES:
            n = counts[category]
            pct = (n / total * 100) if total else 0.0
            fh.write(f"{category}\t{n}\t{pct:.2f}%\n")


def _write_per_genome_breakdown(
    genome_breakdown: dict, genome_names: list, path: str
) -> None:
    """输出 per_genome_breakdown.tsv: 每行一个基因组, 列为 core/soft-core/
    shell/cloud/singleton 五类的基因数(该基因组内属于该分类家族的基因数,
    而非家族个数)。"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("Genome\t" + "\t".join(_PAN_CATEGORIES) + "\n")
        for genome_name in genome_names:
            counts = genome_breakdown[genome_name]
            row = [str(counts[category]) for category in _PAN_CATEGORIES]
            fh.write(f"{genome_name}\t" + "\t".join(row) + "\n")


def _build_intra_summary(intra_results: dict) -> dict:
    total_genes = sum(r.n_genes for r in intra_results.values())
    total_representatives = sum(r.n_representatives for r in intra_results.values())
    compression_ratio = (
        total_representatives / total_genes if total_genes > 0 else 0.0
    )
    return {
        "total_genes": total_genes,
        "total_representatives": total_representatives,
        "compression_ratio": round(compression_ratio, 4),
        "genomes": {
            genome_id: {
                "n_genes": r.n_genes,
                "n_clusters": r.n_clusters,
                "n_representatives": r.n_representatives,
                "n_orphans": len(r.orphan_genes),
                "compression_ratio": round(
                    r.n_representatives / r.n_genes, 4
                ) if r.n_genes > 0 else 0.0,
            }
            for genome_id, r in intra_results.items()
        },
    }


@app.command()
def run(
    input_dir: str = typer.Option(
        ..., "-i", "--input-dir",
        help="Directory of input genome FASTA files (*.fasta/*.fa/*.faa). "
             "Must match the original input set (by file name) on --resume/--append.",
    ),
    output_dir: str = typer.Option(
        "output", "-o", "--output-dir", help="Output directory for result and summary files."
    ),
    append_dir: Optional[str] = typer.Option(
        None, "-a", "--append",
        help="Append mode: directory of new genome FASTA files to add to the "
             "existing results in -o. Can be used multiple times. Requires -i "
             "to still match the original input set; mutually exclusive with --resume.",
    ),
    resume: bool = typer.Option(
        False, "--resume",
        help="Resume an interrupted run from -o, continuing from the last "
             "completed stage. Mutually exclusive with --append.",
    ),
    cutoff_i: float = typer.Option(
        0.5, "--cutoff-i", help="Cutoff-I: inclusion threshold for clustering."
    ),
    cutoff_ii: float = typer.Option(
        0.8, "--cutoff-ii", help="Cutoff-II: high-confidence threshold for representative selection."
    ),
    cutoff_c: float = typer.Option(
        0.5, "--cutoff-c", help="Cutoff-C: coverage threshold used together with --cutoff-i/--cutoff-ii "
             "for clustering and representative decisions. Set to 0.0 to disable coverage filtering."
    ),
    cov_mode: int = typer.Option(
        0, "--cov-mode",
        help="Coverage mode (as in mmseqs2 --cov-mode): "
             "0=min of query/subject coverage (default), 1=subject only, 2=query only.",
    ),
    max_per_cluster: int = typer.Option(
        5, "--max-per-cluster", help="Maximum number of representatives per local cluster."
    ),
    aligner_backend: str = typer.Option(
        "mmseqs", "--aligner",
        help="Alignment backend: mmseqs (default, requires mmseqs2 installed) "
             "or diamond (requires diamond installed).",
    ),
    threads: int = typer.Option(
        0, "--threads",
        help="Number of threads for mmseqs/diamond. 0 lets the tool use its own default "
             "(all cores).",
    ),
    edge_completion_mode: str = typer.Option(
        "representative_member", "--edge-completion-mode",
        help="SSN edge completion strategy: none (representative-representative edges only), "
             "representative_member (default, adds representative-member edges at no extra "
             "alignment cost), or full (adds all member-member edges within a cluster; "
             "cross-genome pairs require new alignments and can be slow for large clusters).",
    ),
    mega_cluster_warning_threshold: Optional[int] = typer.Option(
        None, "--mega-cluster-warning-threshold",
        help="With --edge-completion-mode full, log a warning when a cluster's number of "
             "cross-genome pairs exceeds this threshold. Unset by default (no warning).",
    ),
    inflation: float = typer.Option(
        1.5, "--inflation",
        help="MCL inflation parameter controlling sub-family granularity: lower values yield "
             "coarser clusters, higher values yield finer ones. Requires mcl to be installed.",
    ),
    edge_weight_mode: str = typer.Option(
        "identity_x_coverage", "--edge-weight-mode",
        help="Edge weight used by MCL: identity | coverage | "
             "identity_x_coverage (default, identity*coverage) | "
             "identity_coverage_weighted (alpha*identity + (1-alpha)*coverage, "
             "see --edge-weight-alpha) | bitscore (raw, unnormalized).",
    ),
    edge_weight_alpha: float = typer.Option(
        0.5, "--edge-weight-alpha",
        help="Weight of identity when --edge-weight-mode identity_coverage_weighted is used "
             "(coverage weight is 1-alpha); range [0, 1].",
    ),
    keep_intermediate: bool = typer.Option(
        False, "--keep-intermediate",
        help="Keep intermediate alignment/clustering files instead of deleting them after use. "
             "Index and state files are always kept regardless of this flag.",
    ),
    force: bool = typer.Option(
        False, "--force",
        help="With --resume/--append, continue (with a warning) even if algorithm parameters "
             "differ from the historical run, instead of raising an error.",
    ),
    verbose: bool = typer.Option(False, "--verbose", help="Enable DEBUG-level logging."),
) -> None:
    """Run the LINC pipeline (GenomeLoadStage -> IntraGenomeStage -> GlobalJoinStage -> SSNBuilderStage -> FamilyRefineStage)."""
    _setup_logging(verbose)

    if not os.path.isdir(input_dir):
        typer.secho(f"输入目录不存在: {input_dir}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    _validate_cutoff_range("--cutoff-i", cutoff_i)
    _validate_cutoff_range("--cutoff-ii", cutoff_ii)
    _validate_cutoff_range("--cutoff-c", cutoff_c)
    if cutoff_i > cutoff_ii:
        typer.secho(
            f"--cutoff-i 不应大于 --cutoff-ii(当前 cutoff_i={cutoff_i}, "
            f"cutoff_ii={cutoff_ii}): cutoff_ii 是更高置信的阈值, "
            "若 cutoff_i > cutoff_ii 会产出语义错乱的结果。",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)
    if cov_mode not in VALID_COV_MODES:
        typer.secho(
            f"--cov-mode 取值非法: {cov_mode}(可选: {', '.join(str(m) for m in VALID_COV_MODES)})",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)
    if max_per_cluster < 1:
        typer.secho(
            f"--max-per-cluster 必须 >= 1, 当前为: {max_per_cluster}",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)
    if threads < 0:
        typer.secho(
            f"--threads 必须 >= 0, 当前为: {threads}",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)
    if edge_completion_mode not in VALID_EDGE_COMPLETION_MODES:
        typer.secho(
            f"--edge-completion-mode 取值非法: {edge_completion_mode}"
            f"(可选: {', '.join(VALID_EDGE_COMPLETION_MODES)})",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)
    if mega_cluster_warning_threshold is not None and mega_cluster_warning_threshold < 1:
        typer.secho(
            f"--mega-cluster-warning-threshold 必须 >= 1, 当前为: {mega_cluster_warning_threshold}",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)
    if inflation <= 0:
        typer.secho(
            f"--inflation 必须 > 0, 当前为: {inflation}",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)
    if edge_weight_mode not in VALID_EDGE_WEIGHT_MODES:
        typer.secho(
            f"--edge-weight-mode 取值非法: {edge_weight_mode}"
            f"(可选: {', '.join(VALID_EDGE_WEIGHT_MODES)})",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)
    _validate_cutoff_range("--edge-weight-alpha", edge_weight_alpha)
    if aligner_backend not in ("mmseqs", "diamond"):
        typer.secho(
            f"未知的 --aligner 取值: {aligner_backend}(可选: mmseqs | diamond)",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)
    if resume and append_dir is not None:
        typer.secho(
            "--resume 不能与 -a/--append 同时使用: --resume 是接续上次中断的"
            "运行(不引入新数据), -a 是追加新基因组(要求历史结果已完整), "
            "二者语义互斥。",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)
    if append_dir is not None and not os.path.isdir(append_dir):
        typer.secho(f"-a/--append 指向的目录不存在: {append_dir}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    genome_files = discover_fasta_files(input_dir)
    if not genome_files:
        # 工程审查问题6(2026-10-03 与用户确认): 空输入目录此前只打 warning
        # 后继续执行整个 pipeline, 最终产出空的 Orthogroups.txt, 容易让
        # 用户误以为流程正常完成; 改为直接报错退出。
        typer.secho(
            f"输入目录下未发现任何 FASTA 文件(*.fasta/*.fa/*.faa): {input_dir}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)

    # 工程审查问题7(2026-10-03 与用户确认): 并发跑同一 -o 目录会互相破坏
    # 中间文件/state(search()里的rmtree+重建、FileRepository._save()的
    # .tmp+os.replace均无并发保护)。在触及 output_dir 下任何文件(包括
    # _build_aligner/_build_clustering 创建的索引工作目录)之前, 先以
    # 非阻塞方式获取该目录的独占锁; 若已被另一进程占用, 立即报错退出,
    # 不排队等待(互斥场景而非排队场景)。整个 pipeline 执行过程持有该锁
    # 直到命令结束(正常/异常退出都在 finally 里释放)。
    try:
        lock_cm = acquire_pipeline_lock(output_dir)
        lock_cm.__enter__()
    except PipelineLockedError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    try:
        try:
            aligner = _build_aligner(aligner_backend, threads, output_dir)
        except RuntimeError as exc:
            # 工程审查问题4(2026-10-03 与用户确认): 缺失外部二进制(如 mmseqs/
            # diamond 未安装)时, 构造函数会抛出 RuntimeError, 这里统一兜底为
            # 干净的报错退出, 而不是让用户看到未捕获的 Python traceback。
            typer.secho(str(exc), fg=typer.colors.RED, err=True)
            raise typer.Exit(code=1)
        if aligner is None:
            typer.secho(
                f"未知的 --aligner 取值: {aligner_backend}(可选: mmseqs | diamond)",
                fg=typer.colors.RED,
                err=True,
            )
            raise typer.Exit(code=1)
        logger.info("比对后端: %s", aligner_backend)

        state_dir = os.path.join(output_dir, ".linc_work", "state")
        repository = FileRepository(state_dir)

        if resume or append_dir is not None:
            try:
                historical_signature = repository.load_input_signature()
            except StateCorruptedError as exc:
                typer.secho(str(exc), fg=typer.colors.RED, err=True)
                raise typer.Exit(code=1)
            if historical_signature is None:
                typer.secho(
                    f"-o {output_dir} 下未找到历史运行记录(.linc_work/state), "
                    "--resume/-a 需要先有一次 `linc run` 的历史运行。",
                    fg=typer.colors.RED, err=True,
                )
                raise typer.Exit(code=1)
            # 兼容旧版签名格式(纯 basename 字符串列表): 旧格式只能做 basename
            # 级别比对(无法感知内容是否被替换), 新格式为 (basename, size, mtime)
            # 三元组列表, 能进一步拦住"同名但内容已变"的误用场景(2026-10-03,
            # 第二轮漏洞"签名只比 basename"修复, 与用户权衡性能后采用 大小+mtime
            # 而非完整内容 MD5 的折中方案)。
            is_legacy_format = bool(historical_signature) and isinstance(historical_signature[0], str)
            if is_legacy_format:
                current_signature = sorted(os.path.basename(f) for f in genome_files)
                logger.warning(
                    "历史 input_signature 为旧版格式(仅 basename), 本次校验退化为"
                    "仅比对文件名, 无法感知文件内容是否被替换; 建议重新完整运行一次"
                    "以升级签名格式。"
                )
            else:
                current_signature = build_input_signature(genome_files)
            if current_signature != historical_signature:
                typer.secho(
                    f"-i {input_dir} 下发现的基因组文件集合与历史记录不一致"
                    "(--resume/-a 要求 -i 始终指向同一份历史基础数据集, 不允许变化): "
                    f"历史={historical_signature}, 本次={current_signature}",
                    fg=typer.colors.RED, err=True,
                )
                raise typer.Exit(code=1)

        if append_dir is not None:
            missing_stages = [
                name for name, has in (
                    ("genomes", repository.has_genomes()),
                    ("intra_results", repository.has_intra_results()),
                    ("join_result", repository.has_join_result()),
                    ("complete_ssn", repository.has_complete_ssn()),
                    ("families", repository.has_families()),
                ) if not has
            ]
            if missing_stages:
                typer.secho(
                    f"-a/--append 要求 -o {output_dir} 下是一次完整运行的结果, "
                    f"但缺少以下 Stage 的历史产出: {', '.join(missing_stages)}; "
                    "请先用 --resume 补全上次运行, 再执行 --append。",
                    fg=typer.colors.RED, err=True,
                )
                raise typer.Exit(code=1)

        try:
            clustering_for_mcl = _build_clustering(output_dir)
        except RuntimeError as exc:
            # 工程审查问题4(2026-10-03 与用户确认): 缺失外部二进制(如 mcl 未
            # 安装)时, 构造函数会抛出 RuntimeError, 这里统一兜底为干净的报错
            # 退出, 而不是让用户看到未捕获的 Python traceback。
            typer.secho(str(exc), fg=typer.colors.RED, err=True)
            raise typer.Exit(code=1)

        pipeline = Pipeline(
            aligner=aligner,
            clustering_for_intra=UnionFindClustering(),
            clustering_for_mcl=clustering_for_mcl,
            repository=repository,
            config=PipelineConfig(
                cutoff_i=cutoff_i,
                cutoff_ii=cutoff_ii,
                cutoff_c=cutoff_c,
                cov_mode=cov_mode,
                max_per_cluster=max_per_cluster,
                edge_completion_mode=edge_completion_mode,
                mega_cluster_warning_threshold=mega_cluster_warning_threshold,
                inflation=inflation,
                edge_weight_mode=edge_weight_mode,
                edge_weight_alpha=edge_weight_alpha,
                keep_intermediate=keep_intermediate,
                threads=threads,
            ),
        )

        try:
            if append_dir is not None:
                new_genome_files = discover_fasta_files(append_dir)
                if not new_genome_files:
                    # 工程审查问题6(2026-10-03 与用户确认): -a/--append 指向
                    # 空目录(或目录内无受支持的FASTA文件)时, 此前会"成功"退出
                    # 但实际未追加任何内容, 容易让用户误以为追加生效; 改为报错。
                    typer.secho(
                        f"-a/--append 指向的目录下未发现任何 FASTA 文件"
                        f"(*.fasta/*.fa/*.faa): {append_dir}",
                        fg=typer.colors.RED, err=True,
                    )
                    raise typer.Exit(code=1)
                logger.info(
                    "追加模式(--append): 发现 %d 个新基因组候选文件(目录=%s)",
                    len(new_genome_files), append_dir,
                )
                result = pipeline.append(new_genome_files, force=force)
            elif resume:
                logger.info("断点续跑(--resume): 共 %d 个候选文件", len(genome_files))
                result = pipeline.resume(genome_files, force=force)
            else:
                logger.info("开始加载基因组, 共发现 %d 个候选文件", len(genome_files))
                result = pipeline.run(genome_files)
        except GenomeIdConflictError as exc:
            typer.secho(str(exc), fg=typer.colors.RED, err=True)
            raise typer.Exit(code=1)
        except (StateCorruptedError, IndexMissingError, ConfigMismatchError) as exc:
            typer.secho(str(exc), fg=typer.colors.RED, err=True)
            raise typer.Exit(code=1)

        os.makedirs(output_dir, exist_ok=True)

        summary = _build_load_summary(result.genomes)
        summary_path = os.path.join(output_dir, "load_summary.json")
        with open(summary_path, "w", encoding="utf-8") as fh:
            json.dump(summary, fh, ensure_ascii=False, indent=2)

        typer.echo(
            f"基因组加载完成: 共 {summary['total_genomes']} 个, "
            f"成功 {summary['loaded_count']} 个, 失败 {summary['failed_count']} 个"
        )
        if summary["failed_count"]:
            typer.secho(
                f"失败的基因组: {', '.join(summary['failed_genome_ids'])}",
                fg=typer.colors.YELLOW,
            )
        typer.echo(f"摘要已写入: {summary_path}")

        intra_summary = _build_intra_summary(result.intra_results)
        intra_summary_path = os.path.join(output_dir, "intra_summary.json")
        with open(intra_summary_path, "w", encoding="utf-8") as fh:
            json.dump(intra_summary, fh, ensure_ascii=False, indent=2)

        typer.echo(
            f"基因组内自聚类完成: 共 {intra_summary['total_genes']} 个基因, "
            f"压缩为 {intra_summary['total_representatives']} 个代表 "
            f"(压缩率 {intra_summary['compression_ratio']:.2%})"
        )
        typer.echo(f"摘要已写入: {intra_summary_path}")

        join_summary = _build_join_summary(result.join_result)
        join_summary_path = os.path.join(output_dir, "join_summary.json")
        with open(join_summary_path, "w", encoding="utf-8") as fh:
            json.dump(join_summary, fh, ensure_ascii=False, indent=2)

        typer.echo(
            f"增量并入全局代表库完成: 最终 {join_summary['n_clusters_final']} 个簇, "
            f"{join_summary['n_representatives_final']} 个全局代表, "
            f"{join_summary['n_ssn_edges']} 条 SSN 边"
        )
        typer.echo(f"摘要已写入: {join_summary_path}")

        ssn_summary = _build_ssn_summary(result.complete_ssn)
        ssn_summary_path = os.path.join(output_dir, "ssn_summary.json")
        with open(ssn_summary_path, "w", encoding="utf-8") as fh:
            json.dump(ssn_summary, fh, ensure_ascii=False, indent=2)

        typer.echo(
            f"SSN 构建完成(edge_completion_mode={edge_completion_mode}): "
            f"共 {ssn_summary['n_edges']} 条边"
        )
        typer.echo(f"摘要已写入: {ssn_summary_path}")

        family_summary = _build_family_summary(result.families, inflation, edge_weight_mode)
        family_summary_path = os.path.join(output_dir, "family_summary.json")
        with open(family_summary_path, "w", encoding="utf-8") as fh:
            json.dump(family_summary, fh, ensure_ascii=False, indent=2)

        orthogroups_path = os.path.join(output_dir, "Orthogroups.txt")
        gene_id_to_original = _build_gene_id_to_original(result.genomes)
        _write_orthogroups_txt(result.families, gene_id_to_original, orthogroups_path)

        gene_genome_path = os.path.join(output_dir, "gene_genome.tsv")
        _write_gene_genome_tsv(result.genomes, gene_genome_path)

        typer.echo(f"MCL 亚家族细分完成: 共 {family_summary['n_families']} 个直系同源家族")
        typer.echo(f"摘要已写入: {family_summary_path}")
        typer.echo(f"Orthogroups 已写入: {orthogroups_path}")
        typer.echo(f"基因-基因组对照表已写入: {gene_genome_path}")
    finally:
        lock_cm.__exit__(None, None, None)


@app.command()
def pan(
    orthogroups_path: str = typer.Option(
        ..., "-g", "--orthogroups",
        help="Path to Orthogroups.txt (as produced by `linc run`).",
    ),
    gene_genome_path: str = typer.Option(
        ..., "-m", "--gene-genome",
        help="Path to gene_genome.tsv (as produced by `linc run`).",
    ),
    output_dir: str = typer.Option(
        "output", "-o", "--output-dir",
        help="Output directory for the pan-genome analysis result files.",
    ),
) -> None:
    """Pan-genome analysis based on Orthogroups.txt + gene_genome.tsv: computes
    the presence/absence matrix, classifies each orthogroup into
    core/soft-core/shell/cloud/singleton, and breaks down gene counts per
    genome by category."""
    _setup_logging(False)

    if not os.path.isfile(orthogroups_path):
        typer.secho(
            f"未找到 Orthogroups.txt 文件: {orthogroups_path}",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)
    if not os.path.isfile(gene_genome_path):
        typer.secho(
            f"未找到 gene_genome.tsv 文件: {gene_genome_path}",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)

    gene_to_genome = _load_gene_genome_map(gene_genome_path)
    genome_names = sorted(set(gene_to_genome.values()))
    if not genome_names:
        typer.secho(
            f"gene_genome.tsv 中未解析到任何基因组: {gene_genome_path}",
            fg=typer.colors.RED, err=True,
        )
        raise typer.Exit(code=1)

    families = _parse_orthogroups_txt(orthogroups_path)

    pan_data = _build_pan_genome_data(families, gene_to_genome, genome_names)

    os.makedirs(output_dir, exist_ok=True)

    matrix_path = os.path.join(output_dir, "presence_absence_matrix.tsv")
    _write_presence_absence_matrix(pan_data["presence_absence"], genome_names, matrix_path)

    summary_path = os.path.join(output_dir, "pan_genome_summary.tsv")
    _write_pan_genome_summary(pan_data["family_category"], summary_path)

    breakdown_path = os.path.join(output_dir, "per_genome_breakdown.tsv")
    _write_per_genome_breakdown(pan_data["genome_breakdown"], genome_names, breakdown_path)

    typer.echo(f"泛基因组分析完成: 共 {len(families)} 个直系同源家族, {len(genome_names)} 个基因组")
    typer.echo(f"存在/缺失矩阵已写入: {matrix_path}")
    typer.echo(f"分类汇总已写入: {summary_path}")
    typer.echo(f"各基因组分解表已写入: {breakdown_path}")


if __name__ == "__main__":
    app()
