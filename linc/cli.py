"""LINC 命令行入口。

当前阶段(第4步 engine 层开发中)已接通 GenomeLoadStage + IntraGenomeStage +
GlobalJoinStage,后续每接通一个新 Stage,直接在 `run` 命令内部继续往后扩展
流水线调用,CLI 接口形态本身(`linc run -i <目录> -o <输出目录>`)不再变动。

IntraGenomeStage/GlobalJoinStage 的比对后端可通过 `--backend` 选择(二者共用
同一个 aligner 实例, 以便 GlobalJoinStage 复用 IntraGenomeStage 阶段已建好的
持久化索引, 避免重复建库):
- "mock"(默认, 用于单测/无 mmseqs 环境下打通流程): MockAligner 的长度相似度模式。
- "mmseqs": 调用系统级 mmseqs2(需要已安装并在 PATH 中, 推荐通过 conda `linc`
  环境提供), 使用真实序列比对结果。
- "diamond": 调用系统级 diamond(需要已安装并在 PATH 中, 推荐通过 conda `linc`
  环境提供), 使用真实序列比对结果。
"""
import json
import logging
import os
import sys
import typer

from .adapters.aligners.mock_adapter import MockAligner
from .adapters.aligners.mmseqs_adapter import MMseqsAligner, MMseqsConfig
from .adapters.aligners.diamond_adapter import DiamondAligner, DiamondConfig
from .adapters.clustering.union_find_adapter import UnionFindClustering
from .core.entities import VALID_COV_MODES, GenomeStatus
from .engine.genome_load_stage import (
    GenomeIdConflictError,
    GenomeLoadStage,
    discover_fasta_files,
)
from .engine.global_join_stage import GlobalJoinStage
from .engine.intra_genome_stage import IntraGenomeStage

app = typer.Typer(help="LINC: 增量基因组比较构建 SSN 并通过 MCL 细分蛋白质家族。")


@app.callback()
def main() -> None:
    """LINC 命令行工具。"""
    pass


logger = logging.getLogger("linc")


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        stream=sys.stdout,
    )


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
        ..., "-i", "--input-dir", help="包含基因组 FASTA 文件(*.fasta/*.fa/*.faa)的目录。"
    ),
    output_dir: str = typer.Option(
        "output", "-o", "--output-dir", help="输出目录,用于写入运行摘要等结果文件。"
    ),
    cutoff_i: float = typer.Option(
        0.5, "--cutoff-i", help="Cutoff-I(纳入阈值): 判定是否成簇。"
    ),
    cutoff_ii: float = typer.Option(
        0.8, "--cutoff-ii", help="Cutoff-II(高置信阈值): 判定代表选取策略。"
    ),
    cutoff_c: float = typer.Option(
        0.0, "--cutoff-c", help="Cutoff-C(覆盖度纳入阈值): 与 --cutoff-i 共同判定是否成簇, "
             "默认 0.0 等价于不启用覆盖度过滤。"
    ),
    cov_mode: int = typer.Option(
        0, "--cov-mode",
        help="覆盖度合成方式(参考 mmseqs2 --cov-mode): "
             "0=双向(取 query/subject 覆盖度较小者, 默认) | "
             "1=仅看 subject 覆盖度 | 2=仅看 query 覆盖度。",
    ),
    max_per_cluster: int = typer.Option(
        5, "--max-per-cluster", help="单个局部簇最多选取的代表数。"
    ),
    backend: str = typer.Option(
        "mock", "--backend",
        help="比对引擎后端: mock(默认, 无需安装比对工具) | mmseqs(真实比对, 需已"
             "安装 mmseqs2) | diamond(真实比对, 需已安装 diamond)。",
    ),
    threads: int = typer.Option(
        0, "--threads",
        help="mmseqs/diamond 使用的线程数(仅 --backend mmseqs|diamond 生效)。"
             "0 = 不显式指定, 让工具使用其自身默认的全部核心数(注意: 多线程并行"
             "时结果在极少数边界情况下不保证严格可复现, 如需完全可复现结果请"
             "指定 --threads 1)。",
    ),
    verbose: bool = typer.Option(False, "--verbose", help="输出 DEBUG 级别日志。"),
) -> None:
    """运行 LINC 流水线(当前阶段已接通 GenomeLoadStage + IntraGenomeStage + GlobalJoinStage)。"""
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

    genome_files = discover_fasta_files(input_dir)
    if not genome_files:
        typer.secho(
            f"输入目录下未发现任何 FASTA 文件(*.fasta/*.fa/*.faa): {input_dir}",
            fg=typer.colors.YELLOW,
            err=True,
        )

    logger.info("开始加载基因组, 共发现 %d 个候选文件", len(genome_files))
    try:
        genomes = GenomeLoadStage().run(genome_files)
    except GenomeIdConflictError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    summary = _build_load_summary(genomes)

    os.makedirs(output_dir, exist_ok=True)
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

    loaded_genomes = {
        genome_id: genome
        for genome_id, genome in genomes.items()
        if genome.status == GenomeStatus.LOADED
    }
    if backend == "mmseqs":
        work_dir = os.path.join(output_dir, ".linc_work", "mmseqs_db")
        aligner = MMseqsAligner(MMseqsConfig(threads=threads, work_dir=work_dir))
    elif backend == "diamond":
        work_dir = os.path.join(output_dir, ".linc_work", "diamond_db")
        aligner = DiamondAligner(DiamondConfig(threads=threads, work_dir=work_dir))
    elif backend == "mock":
        aligner = MockAligner()
    else:
        typer.secho(
            f"未知的 --backend 取值: {backend}(可选: mock | mmseqs | diamond)",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)
    logger.info("比对后端: %s", backend)
    intra_stage = IntraGenomeStage(
        aligner=aligner,
        clustering=UnionFindClustering(),
        cutoff_i=cutoff_i,
        cutoff_ii=cutoff_ii,
        cutoff_c=cutoff_c,
        cov_mode=cov_mode,
        max_per_cluster=max_per_cluster,
    )
    logger.info("开始基因组内自聚类(IntraGenomeStage), 共 %d 个基因组", len(loaded_genomes))
    intra_results = {
        genome_id: intra_stage.run(genome)
        for genome_id, genome in loaded_genomes.items()
    }

    intra_summary = _build_intra_summary(intra_results)
    intra_summary_path = os.path.join(output_dir, "intra_summary.json")
    with open(intra_summary_path, "w", encoding="utf-8") as fh:
        json.dump(intra_summary, fh, ensure_ascii=False, indent=2)

    typer.echo(
        f"基因组内自聚类完成: 共 {intra_summary['total_genes']} 个基因, "
        f"压缩为 {intra_summary['total_representatives']} 个代表 "
        f"(压缩率 {intra_summary['compression_ratio']:.2%})"
    )
    typer.echo(f"摘要已写入: {intra_summary_path}")

    join_stage = GlobalJoinStage(
        aligner=aligner,
        cutoff_i=cutoff_i,
        cutoff_ii=cutoff_ii,
        cutoff_c=cutoff_c,
        cov_mode=cov_mode,
        max_per_cluster=max_per_cluster,
    )
    logger.info("开始增量并入全局代表库(GlobalJoinStage), 共 %d 个基因组", len(loaded_genomes))
    join_result = join_stage.run(loaded_genomes, intra_results)

    join_summary = _build_join_summary(join_result)
    join_summary_path = os.path.join(output_dir, "join_summary.json")
    with open(join_summary_path, "w", encoding="utf-8") as fh:
        json.dump(join_summary, fh, ensure_ascii=False, indent=2)

    typer.echo(
        f"增量并入全局代表库完成: 最终 {join_summary['n_clusters_final']} 个簇, "
        f"{join_summary['n_representatives_final']} 个全局代表, "
        f"{join_summary['n_ssn_edges']} 条 SSN 边"
    )
    typer.echo(f"摘要已写入: {join_summary_path}")


if __name__ == "__main__":
    app()
