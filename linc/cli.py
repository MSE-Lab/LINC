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
from .engine.genome_load_stage import GenomeIdConflictError, discover_fasta_files
from .adapters.storage.file_repo import FileRepository, StateCorruptedError
from .pipeline import IndexMissingError, Pipeline, PipelineConfig

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


def _write_orthogroups_txt(families, path: str) -> None:
    """输出 Orthogroups.txt(一行一个 ortholog group): 第一列是 OG 编号,
    第二列是该家族全部成员基因ID, 用空格分隔(与用户确认, 2026-10-02)。
    """
    with open(path, "w", encoding="utf-8") as fh:
        for family in families:
            members = " ".join(sorted(family.member_genes))
            fh.write(f"{family.family_id}\t{members}\n")


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
        help="包含基因组 FASTA 文件(*.fasta/*.fa/*.faa)的目录。表示历史(首次)"
             "运行的输入集合; --resume/--append 时必须与历史记录一致"
             "(按文件名比对), 不可变化。",
    ),
    output_dir: str = typer.Option(
        "output", "-o", "--output-dir", help="输出目录,用于写入运行摘要等结果文件。"
    ),
    append_dir: Optional[str] = typer.Option(
        None, "-a", "--append",
        help="追加模式: 指定一个包含新基因组 FASTA 文件的目录, 将其追加到"
             "-o 指向的历史完整结果上(要求 -o 下已存在一次完整运行的结果)。"
             "可多次追加(每次 -a 指向不同的新目录)。-i 仍须传入且必须与历史"
             "一致, 不与 --resume 共用。",
    ),
    resume: bool = typer.Option(
        False, "--resume",
        help="断点续跑: 从 -o 历史中断的 Stage 接着跑完剩余步骤, -i 必须与"
             "历史一致, 不得同时传入 -a。",
    ),
    cutoff_i: float = typer.Option(
        0.5, "--cutoff-i", help="Cutoff-I(纳入阈值): 判定是否成簇。"
    ),
    cutoff_ii: float = typer.Option(
        0.8, "--cutoff-ii", help="Cutoff-II(高置信阈值): 判定代表选取策略。"
    ),
    cutoff_c: float = typer.Option(
        0.5, "--cutoff-c", help="Cutoff-C(覆盖度纳入阈值): 与 --cutoff-i 共同判定是否成簇, "
             "同时也与 --cutoff-ii 共同判定高置信代表选取/并入决策(第二轮漏洞3修复, "
             "2026-10-03); 默认 0.5。设为 0.0 等价于不启用覆盖度过滤(回退为旧行为)。"
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
    aligner_backend: str = typer.Option(
        "mmseqs", "--aligner",
        help="比对引擎后端: mmseqs(默认, 真实比对, 需已安装 mmseqs2) | "
             "diamond(真实比对, 需已安装 diamond)。",
    ),
    threads: int = typer.Option(
        0, "--threads",
        help="mmseqs/diamond 使用的线程数(仅 --aligner mmseqs|diamond 生效)。"
             "0 = 不显式指定, 让工具使用其自身默认的全部核心数(注意: 多线程并行"
             "时结果在极少数边界情况下不保证严格可复现, 如需完全可复现结果请"
             "指定 --threads 1)。",
    ),
    edge_completion_mode: str = typer.Option(
        "representative_member", "--edge-completion-mode",
        help="SSN 补边模式: none(不补边, 仅保留代表-代表边) | "
             "representative_member(默认, 补代表-成员边, 复用已有比对结果,"
             "零新增比对开销) | full(簇内全量补边, 同基因组内部分依旧免费复用,"
             "仅跨基因组成员对会真正发起新比对, 大规模簇可能耗时较长)。",
    ),
    mega_cluster_warning_threshold: Optional[int] = typer.Option(
        None, "--mega-cluster-warning-threshold",
        help="仅 --edge-completion-mode full 生效: 当某个簇的跨基因组待比对"
             "成员对数量超过该阈值时, 记录 warning 提醒(不会自动降级或拆分,"
             "流程照常完整执行)。默认不设阈值, 不做任何提醒。",
    ),
    inflation: float = typer.Option(
        1.5, "--inflation",
        help="MCL(Markov Clustering)膨胀参数, 控制亚家族细分的精细程度: "
             "值越低聚类越粗(更多基因归入同一家族), 值越高聚类越细。需要已"
             "安装系统级 mcl(推荐: conda install -c bioconda mcl)。",
    ),
    edge_weight_mode: str = typer.Option(
        "identity_x_coverage", "--edge-weight-mode",
        help="喂给 MCL 的边权重计算方式: identity | coverage | "
             "identity_x_coverage(默认, identity*coverage) | "
             "identity_coverage_weighted(alpha*identity+(1-alpha)*coverage,"
             "配合 --edge-weight-alpha) | bitscore(原始值, 不做归一化)。",
    ),
    edge_weight_alpha: float = typer.Option(
        0.5, "--edge-weight-alpha",
        help="仅 --edge-weight-mode identity_coverage_weighted 生效: "
             "identity 的权重系数, coverage 权重为 (1-alpha), 范围 [0, 1]。",
    ),
    keep_intermediate: bool = typer.Option(
        False, "--keep-intermediate",
        help="默认(False)下各 Stage 用完 search()/mcl 产生的中间文件即删"
             "(如 search_join/search_full_to_*/component_*), 避免 "
             ".linc_work 无限累积(数据管理优化第3项, 2026-10-03); "
             "索引目录/state pkl 不受影响, 始终保留。设为 True 保留全部中间"
             "产物以便调试。",
    ),
    verbose: bool = typer.Option(False, "--verbose", help="输出 DEBUG 级别日志。"),
) -> None:
    """运行 LINC 流水线(当前阶段已接通 GenomeLoadStage + IntraGenomeStage +
    GlobalJoinStage + SSNBuilderStage + FamilyRefineStage)。"""
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
        typer.secho(
            f"输入目录下未发现任何 FASTA 文件(*.fasta/*.fa/*.faa): {input_dir}",
            fg=typer.colors.YELLOW,
            err=True,
        )

    aligner = _build_aligner(aligner_backend, threads, output_dir)
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
        current_signature = sorted(os.path.basename(f) for f in genome_files)
        if historical_signature is None:
            typer.secho(
                f"-o {output_dir} 下未找到历史运行记录(.linc_work/state), "
                "--resume/-a 需要先有一次 `linc run` 的历史运行。",
                fg=typer.colors.RED, err=True,
            )
            raise typer.Exit(code=1)
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

    pipeline = Pipeline(
        aligner=aligner,
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=_build_clustering(output_dir),
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
        ),
    )

    try:
        if append_dir is not None:
            new_genome_files = discover_fasta_files(append_dir)
            logger.info(
                "追加模式(--append): 发现 %d 个新基因组候选文件(目录=%s)",
                len(new_genome_files), append_dir,
            )
            result = pipeline.append(new_genome_files)
        elif resume:
            logger.info("断点续跑(--resume): 共 %d 个候选文件", len(genome_files))
            result = pipeline.resume(genome_files)
        else:
            logger.info("开始加载基因组, 共发现 %d 个候选文件", len(genome_files))
            result = pipeline.run(genome_files)
    except GenomeIdConflictError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except (StateCorruptedError, IndexMissingError) as exc:
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
    _write_orthogroups_txt(result.families, orthogroups_path)

    typer.echo(f"MCL 亚家族细分完成: 共 {family_summary['n_families']} 个直系同源家族")
    typer.echo(f"摘要已写入: {family_summary_path}")
    typer.echo(f"Orthogroups 已写入: {orthogroups_path}")


if __name__ == "__main__":
    app()
