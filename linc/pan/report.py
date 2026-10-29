"""报告编排: 写 run_summary.json、report.html, 并串联 TSV/SVG/SSN 输出。

设计依据: docs/LINC_结果描述增强开发计划.md 第 3/4 节。本模块是 `linc pan`
的编排入口: 加载 -> 分析 -> 写 TSV -> 写 SVG -> 写 SSN -> 写 run_summary.json
-> 写 report.html, 全程只读 state、只写 report 目录。
"""
import json
import os
from datetime import datetime
from typing import Any, Dict, List, Optional

from .analysis import ReportData, analyze, write_tsvs
from .charts import render_charts
from .loader import LoadedState, load_state
from .ssn_html import render_ssn_html

SUMMARY_SCHEMA = "linc-run-summary/1"


def _fmt(value: Any) -> str:
    """数值展示: 浮点保留 3 位, 整数/字符串原样。"""
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def _warnings(state: LoadedState) -> List[str]:
    warnings: List[str] = []
    if state.complete_ssn is None:
        warnings.append("缺少 complete_ssn.pkl, SSN 相关统计被跳过")
    if not state.intra_results:
        warnings.append("缺少 intra_results.pkl, 局部簇/代表统计被跳过")
    if state.join_result is None:
        warnings.append("缺少 join_result.pkl, 全局簇统计被跳过")
    if not state.config:
        warnings.append("缺少 config.pkl, 参数信息缺失")
    return warnings


def build_run_summary(state: LoadedState, data: ReportData,
                      visualization: Dict[str, Any]) -> Dict[str, Any]:
    """组装 run_summary.json 的顶层结构。"""
    return {
        "schema": SUMMARY_SCHEMA,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "output_dir": state.output_dir,
        "parameters": state.config,
        "counts": data.counts,
        "derived": data.derived,
        "ssn": data.ssn,
        "families": data.families,
        "pan_genome": data.pan_genome,
        "mcl": data.mcl,
        "cluster_size_distribution": data.cluster_size_distribution,
        "genomes": data.genomes,
        "top_families": data.top_families,
        "top_components": data.top_components,
        "visualization": visualization,
        "warnings": _warnings(state),
    }


# ---------------------------------------------------------------------------
# report.html
# ---------------------------------------------------------------------------
def _cards_html(data: ReportData) -> str:
    c = data.counts
    d = data.derived
    s = data.ssn
    f = data.families
    m = data.mcl
    p = data.pan_genome["categories"]
    cards = [
        ("Genomes", c["n_genomes"]),
        ("Genes", c["n_genes"]),
        ("Avg genes per genome", d["avg_genes_per_genome"]),
        ("Local representatives", c["n_local_representatives"]),
        ("Avg local reps per genome", d["avg_local_reps_per_genome"]),
        ("Global clusters", c["n_global_clusters"]),
        ("Global representatives", c["n_global_representatives"]),
        ("Avg global reps per cluster", d["avg_global_reps_per_cluster"]),
        ("SSN nodes", s["n_nodes"]),
        ("SSN edges", s["n_edges"]),
        ("SSN components", s["n_components"]),
        ("Largest component", d["largest_component_size"]),
        ("Isolated genes", s["n_isolated_nodes"]),
        ("MCL families", c["n_families"]),
        ("Singleton families", f["n_singleton"]),
        ("Median family size", f["median_size"]),
        ("Largest family", f["max_size"]),
        ("Components split by MCL", m["n_components_split"]),
        ("Core genes", p["core"]["n_genes"]),
        ("Soft-core genes", p["soft-core"]["n_genes"]),
        ("Shell genes", p["shell"]["n_genes"]),
        ("Cloud genes", p["cloud"]["n_genes"]),
    ]
    items = "\n".join(
        f'<div class="card"><div class="card-value">{_fmt(v)}</div>'
        f'<div class="card-label">{label}</div></div>'
        for label, v in cards
    )
    return f'<div class="cards">{items}</div>'


def _table(headers: List[str], rows: List[List[Any]]) -> str:
    head = "".join(f"<th>{h}</th>" for h in headers)
    body = "\n".join(
        "<tr>" + "".join(f"<td>{_fmt(c)}</td>" for c in row) + "</tr>"
        for row in rows
    )
    return f'<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def _top_families_table(data: ReportData) -> str:
    headers = ["family", "size", "genomes", "density", "category", "representative"]
    rows = [[f["family_id"], f["size"], f["n_genomes"], f["edge_density"],
             f["category"], f["representative"]] for f in data.top_families]
    return _table(headers, rows)


def _top_components_table(data: ReportData) -> str:
    headers = ["component", "genes", "families"]
    rows = [[c["component_id"], c["n_genes"], c["n_families"]]
            for c in data.top_components]
    return _table(headers, rows)


def _genome_table(state: LoadedState, data: ReportData) -> str:
    headers = ["genome", "genes", "families", "singleton genes"]
    rows = []
    for gid in sorted(state.genomes.keys()):
        name = state.genome_id_to_name.get(gid, gid)
        g = data.genomes.get(gid, {})
        rows.append([name, g.get("n_genes", 0), g.get("n_families", 0),
                     g.get("n_singleton_genes", 0)])
    return _table(headers, rows[:100])


def _params_table(state: LoadedState) -> str:
    headers = ["parameter", "value"]
    rows = [[k, v] for k, v in sorted(state.config.items())]
    return _table(headers, rows)


def build_report_html(state: LoadedState, data: ReportData) -> str:
    """生成自包含 report.html(除 ssn.html 的 Cytoscape CDN 外无外部依赖)。"""
    cards = _cards_html(data)
    top_families = _top_families_table(data)
    top_components = _top_components_table(data)
    genome_table = _genome_table(state, data)
    params_table = _params_table(state)

    links = "\n".join(
        f'<li><a href="{name}">{name}</a></li>'
        for name in (
            "run_summary.json", "presence_absence_matrix.tsv",
            "per_genome_breakdown.tsv", "pan_genome_summary.tsv",
            "ssn_nodes.tsv", "ssn_edges.tsv", "family_members.tsv",
            "family_stats.tsv", "mcl_components.tsv",
            "cluster_size_distribution.tsv",
        )
    )

    warnings = "".join(f"<li>{w}</li>" for w in _warnings(state)) or "<li>无</li>"

    template = """<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8">
<title>LINC 运行报告</title>
<style>
  body { font-family: sans-serif; margin: 20px; color: #222; }
  h1 { font-size: 22px; } h2 { font-size: 16px; border-bottom: 1px solid #ccc; padding-bottom: 4px; }
  .cards { display: flex; flex-wrap: wrap; gap: 10px; }
  .card { border: 1px solid #ddd; border-radius: 6px; padding: 10px 14px; min-width: 130px; }
  .card-value { font-size: 20px; font-weight: bold; }
  .card-label { font-size: 12px; color: #666; }
  table { border-collapse: collapse; margin: 8px 0; }
  th, td { border: 1px solid #ddd; padding: 4px 8px; font-size: 12px; text-align: left; }
  th { background: #f5f5f5; }
  img { max-width: 100%; border: 1px solid #eee; margin: 6px 0; }
  iframe { width: 100%; height: 600px; border: 1px solid #ccc; }
  .warn { color: #b8860b; }
</style>
</head>
<body>
<h1>LINC 运行报告</h1>

<h2>关键指标</h2>
__CARDS__

<h2>泛基因组分析</h2>
<img src="pan_genome_families.svg" alt="pan-genome families">
<img src="pan_genome_genes.svg" alt="pan-genome genes">

<h2>家族交互子图</h2>
<iframe src="ssn.html"></iframe>

<h2>簇大小分布</h2>
<img src="cluster_size_distribution.svg" alt="cluster size distribution">

<h2>MCL 前后对比</h2>
<img src="mcl_before_after.svg" alt="MCL before/after">
<img src="mcl_top_components.svg" alt="MCL top components">

<h2>Top families</h2>
__TOP_FAMILIES__

<h2>被 MCL 拆分最多的连通分量</h2>
__TOP_COMPONENTS__

<h2>各基因组概况(前 100)</h2>
__GENOME_TABLE__

<h2>运行参数</h2>
__PARAMS_TABLE__

<h2>结果文件</h2>
<ul>__LINKS__</ul>

<h2>警告</h2>
<ul class="warn">__WARNINGS__</ul>
</body>
</html>
"""
    return (
        template
        .replace("__CARDS__", cards)
        .replace("__TOP_FAMILIES__", top_families)
        .replace("__TOP_COMPONENTS__", top_components)
        .replace("__GENOME_TABLE__", genome_table)
        .replace("__PARAMS_TABLE__", params_table)
        .replace("__LINKS__", links)
        .replace("__WARNINGS__", warnings)
    )


# ---------------------------------------------------------------------------
# 编排入口
# ---------------------------------------------------------------------------
def generate_report(output_dir: str, report_dir: Optional[str] = None,
                    family_ids: Optional[List[str]] = None,
                    graphml: bool = False) -> str:
    """完整生成报告, 返回报告目录路径。

    流程: load_state -> analyze -> write_tsvs -> render_charts ->
    render_ssn_html -> run_summary.json(两处) -> report.html。
    """
    state = load_state(output_dir)
    data = analyze(state)

    target = report_dir or os.path.join(output_dir, "report")
    os.makedirs(target, exist_ok=True)

    tsvs = write_tsvs(state, data, target)
    charts = render_charts(state, data, target)
    ssn = render_ssn_html(state, data, target, family_ids=family_ids, graphml=graphml)

    visualization = {
        "report_dir": target,
        "tsv": tsvs,
        "charts": charts,
        "ssn": ssn,
    }
    summary = build_run_summary(state, data, visualization)

    # run_summary.json 同时写到 report 目录与输出根目录。
    summary_payload = json.dumps(summary, ensure_ascii=False, indent=2)
    with open(os.path.join(target, "run_summary.json"), "w", encoding="utf-8") as fh:
        fh.write(summary_payload)
    with open(os.path.join(output_dir, "run_summary.json"), "w", encoding="utf-8") as fh:
        fh.write(summary_payload)

    html = build_report_html(state, data)
    with open(os.path.join(target, "report.html"), "w", encoding="utf-8") as fh:
        fh.write(html)

    return target
