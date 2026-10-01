"""LINC 命令行入口。

当前阶段(第4步 engine 层开发中)仅接通 GenomeLoadStage,后续每接通一个新
Stage,直接在 `run` 命令内部继续往后扩展流水线调用,CLI 接口形态本身
(`linc run -i <目录> -o <输出目录>`)不再变动。
"""
import json
import logging
import os
import sys
import typer

from .core.entities import GenomeStatus
from .engine.genome_load_stage import GenomeLoadStage, discover_fasta_files

app = typer.Typer(help="LINC: 增量基因组比较构建 SSN 并通过 MCL 细分蛋白质家族。")

logger = logging.getLogger("linc")


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        stream=sys.stdout,
    )


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


@app.command()
def run(
    input_dir: str = typer.Option(
        ..., "-i", "--input-dir", help="包含基因组 FASTA 文件(*.fasta/*.fa/*.faa)的目录。"
    ),
    output_dir: str = typer.Option(
        "output", "-o", "--output-dir", help="输出目录,用于写入运行摘要等结果文件。"
    ),
    verbose: bool = typer.Option(False, "--verbose", help="输出 DEBUG 级别日志。"),
) -> None:
    """运行 LINC 流水线(当前阶段仅接通 GenomeLoadStage)。"""
    _setup_logging(verbose)

    if not os.path.isdir(input_dir):
        typer.secho(f"输入目录不存在: {input_dir}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    genome_files = discover_fasta_files(input_dir)
    if not genome_files:
        typer.secho(
            f"输入目录下未发现任何 FASTA 文件(*.fasta/*.fa/*.faa): {input_dir}",
            fg=typer.colors.YELLOW,
            err=True,
        )

    logger.info("开始加载基因组, 共发现 %d 个候选文件", len(genome_files))
    genomes = GenomeLoadStage().run(genome_files)

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


if __name__ == "__main__":
    app()
