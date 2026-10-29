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

COLLAPSED_MODE = "collapsed"
COLLAPSED_NOTE = (
    "family_mode=collapsed: SSN statistics are on the representative graph "
    "(edges of folded genes were removed), so SSN node/edge counts are smaller "
    "than the number of genes in families."
)


def _fmt(value: Any) -> str:
    """数值展示: 整数加千分位, 浮点保留 3 位, 其余原样。"""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return f"{value:,.3f}"
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def _warnings(state: LoadedState) -> List[str]:
    warnings: List[str] = []
    if state.complete_ssn is None:
        warnings.append("complete_ssn.pkl is missing; SSN statistics are skipped.")
    if not state.intra_results:
        warnings.append("intra_results.pkl is missing; local-cluster/representative statistics are skipped.")
    if state.join_result is None:
        warnings.append("join_result.pkl is missing; global-cluster statistics are skipped.")
    if not state.config:
        warnings.append("config.pkl is missing; parameter information is unavailable.")
    if state.config.get("family_mode") == COLLAPSED_MODE:
        warnings.append(COLLAPSED_NOTE)
    return warnings


def build_run_summary(state: LoadedState, data: ReportData,
                      visualization: Dict[str, Any]) -> Dict[str, Any]:
    """组装 run_summary.json 的顶层结构。"""
    ssn = dict(data.ssn)
    ssn["scope"] = (
        "representative_graph"
        if state.config.get("family_mode") == COLLAPSED_MODE
        else "full_graph"
    )
    return {
        "schema": SUMMARY_SCHEMA,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "output_dir": state.output_dir,
        "parameters": state.config,
        "counts": data.counts,
        "derived": data.derived,
        "ssn": ssn,
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
def _card(label: str, value: Any) -> str:
    return (
        f'<div class="card"><div class="card-value">{_fmt(value)}</div>'
        f'<div class="card-label">{label}</div></div>'
    )


def _card_group(title: str, items: List[tuple]) -> str:
    cards = "\n".join(_card(label, value) for label, value in items)
    return (
        f'<div class="card-group"><h3>{title}</h3>'
        f'<div class="cards">{cards}</div></div>'
    )


def _cards_html(data: ReportData) -> str:
    c = data.counts
    d = data.derived
    s = data.ssn
    f = data.families
    m = data.mcl
    p = data.pan_genome["categories"]

    groups = [
        ("Dataset", [
            ("Genomes", c["n_genomes"]),
            ("Genes", c["n_genes"]),
            ("Avg genes per genome", d["avg_genes_per_genome"]),
        ]),
        ("Representatives", [
            ("Local representatives", c["n_local_representatives"]),
            ("Avg local reps per genome", d["avg_local_reps_per_genome"]),
            ("Global clusters", c["n_global_clusters"]),
            ("Global representatives", c["n_global_representatives"]),
            ("Avg global reps per cluster", d["avg_global_reps_per_cluster"]),
        ]),
        ("Sequence similarity network", [
            ("SSN nodes", s["n_nodes"]),
            ("SSN edges", s["n_edges"]),
            ("SSN components", s["n_components"]),
            ("Largest component", d["largest_component_size"]),
            ("Isolated genes", s["n_isolated_nodes"]),
        ]),
        ("Families (MCL)", [
            ("MCL families", c["n_families"]),
            ("Singleton families", f["n_singleton"]),
            ("Median family size", f["median_size"]),
            ("Largest family", f["max_size"]),
            ("Components split by MCL", m["n_components_split"]),
        ]),
        ("Pan-genome", [
            ("Core genes", p["core"]["n_genes"]),
            ("Soft-core genes", p["soft-core"]["n_genes"]),
            ("Shell genes", p["shell"]["n_genes"]),
            ("Cloud genes", p["cloud"]["n_genes"]),
        ]),
    ]
    return "\n".join(_card_group(title, items) for title, items in groups)


def _table(headers: List[str], rows: List[List[Any]],
           num_cols: Optional[List[int]] = None) -> str:
    """渲染表格; num_cols 指定需要右对齐的列索引(从 0 开始)。"""
    num_cols = num_cols or []
    head = "".join(
        f'<th class="num">{h}</th>' if i in num_cols else f"<th>{h}</th>"
        for i, h in enumerate(headers)
    )
    body = []
    for row in rows:
        cells = []
        for i, cell in enumerate(row):
            cls = ' class="num"' if i in num_cols else ""
            cells.append(f"<td{cls}>{_fmt(cell)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return (
        '<table><thead><tr>' + head + '</tr></thead><tbody>'
        + "\n".join(body) + '</tbody></table>'
    )


def _pan_genome_table(data: ReportData) -> str:
    headers = ["Category", "Families", "Families %", "Genes", "Genes %"]
    # 派生行预计算: single-copy-core 是 core 的子集, cloud+singleton 是合并视图。
    n_core = data.pan_genome["categories"]["core"]["n_families"]
    n_sc = data.pan_genome["n_single_copy_core"]
    n_genomes = data.counts["n_genomes"]
    n_genes = data.counts["n_genes"]
    sc_families_pct = f"{n_sc / n_core * 100:.2f}% of core" if n_core else "—"
    sc_genes = n_sc * n_genomes
    sc_genes_pct = f"{sc_genes / n_genes * 100:.2f}%" if n_genes else "—"
    sc_row = ["single-copy-core*", n_sc, sc_families_pct, sc_genes, sc_genes_pct]

    ci = data.pan_genome["cloud_including_singleton"]
    ci_row = [
        "cloud+singleton*", ci["n_families"],
        f"{ci['family_pct'] * 100:.2f}%",
        ci["n_genes"], f"{ci['gene_pct'] * 100:.2f}%",
    ]

    rows = []
    for cat in ("core", "soft-core", "shell", "cloud", "singleton"):
        c = data.pan_genome["categories"][cat]
        rows.append([
            cat, c["n_families"],
            f"{c['family_pct'] * 100:.2f}%",
            c["n_genes"],
            f"{c['gene_pct'] * 100:.2f}%",
        ])
        # 派生行紧跟其"父"分类: single-copy-core 在 core 后, cloud+singleton 在 singleton 后。
        if cat == "core":
            rows.append(list(sc_row))
        if cat == "singleton":
            rows.append(list(ci_row))

    table = _table(headers, rows, num_cols=[1, 2, 3, 4])
    note = ('<p class="note">* single-copy-core is a subset of core (core families '
            'with exactly one copy per genome); its Families % is relative to core. '
            'cloud+singleton is a derived merged view of cloud and singleton.</p>')
    return table + note


def _top_families_table(data: ReportData) -> str:
    headers = ["Family", "Size", "Genomes", "Density", "Category", "Representative"]
    rows = [[f["family_id"], f["size"], f["n_genomes"], f["edge_density"],
             f["category"], f["representative"]] for f in data.top_families]
    return _table(headers, rows, num_cols=[1, 2, 3])


def _top_components_table(data: ReportData) -> str:
    headers = ["Component", "Genes", "Families"]
    rows = [[c["component_id"], c["n_genes"], c["n_families"]]
            for c in data.top_components]
    return _table(headers, rows, num_cols=[1, 2])


def _genome_table(state: LoadedState, data: ReportData) -> str:
    headers = ["Genome", "Genes", "Families", "Singleton genes"]
    rows = []
    for gid in sorted(state.genomes.keys()):
        name = state.genome_id_to_name.get(gid, gid)
        g = data.genomes.get(gid, {})
        rows.append([name, g.get("n_genes", 0), g.get("n_families", 0),
                     g.get("n_singleton_genes", 0)])
    return _table(headers, rows[:100], num_cols=[1, 2, 3])


def _params_table(state: LoadedState) -> str:
    headers = ["Parameter", "Value"]
    rows = [[k, v] for k, v in sorted(state.config.items())]
    return _table(headers, rows)


def _figure(src: str, caption: str) -> str:
    return (
        f'<figure><img src="{src}" alt="{caption}">'
        f'<figcaption>{caption}</figcaption></figure>'
    )


def build_report_html(state: LoadedState, data: ReportData) -> str:
    """生成自包含 report.html(除 ssn.html 的 Cytoscape CDN 外无外部依赖)。"""
    cards = _cards_html(data)
    pan_table = _pan_genome_table(data)
    top_families = _top_families_table(data)
    top_components = _top_components_table(data)
    genome_table = _genome_table(state, data)
    params_table = _params_table(state)

    generated_at = datetime.now().isoformat(timespec="seconds")

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

    warnings = "".join(f"<li>{w}</li>" for w in _warnings(state))

    template = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>LINC Run Report</title>
<style>
  :root {
    --accent: #1b4f72;
    --accent-light: #2980b9;
    --bg: #f5f7fa;
    --card-bg: #ffffff;
    --border: #dfe3e8;
    --text: #1a1d21;
    --muted: #5f6b7a;
  }
  * { box-sizing: border-box; }
  body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI",
         Roboto, "Helvetica Neue", Arial, sans-serif; background: var(--bg);
         color: var(--text); line-height: 1.55; }
  header { background: linear-gradient(120deg, #14334f, #1b4f72 55%, #2980b9);
           color: #fff; padding: 30px 40px; }
  header h1 { margin: 0 0 6px; font-size: 26px; font-weight: 600; letter-spacing: -0.01em; }
  header .subtitle { font-size: 14px; opacity: 0.88; }
  header .subtitle code { font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
                          background: rgba(255,255,255,.15); padding: 1px 6px; border-radius: 4px; }
  main { max-width: 1080px; margin: 0 auto; padding: 28px 40px 60px; }
  section { margin: 0 0 40px; }
  h2 { font-size: 17px; font-weight: 600; color: var(--accent);
       border-bottom: 2px solid var(--accent-light); padding-bottom: 6px; margin: 0 0 16px; }
  h3 { font-size: 12px; font-weight: 600; color: var(--muted);
       text-transform: uppercase; letter-spacing: 0.06em; margin: 0 0 10px; }
  .card-group { margin: 0 0 18px; }
  .cards { display: grid; grid-template-columns: repeat(auto-fill, minmax(160px, 1fr)); gap: 12px; }
  .card { background: var(--card-bg); border: 1px solid var(--border);
          border-radius: 8px; padding: 14px 16px;
          box-shadow: 0 1px 2px rgba(16,24,40,.04); }
  .card-value { font-size: 24px; font-weight: 650; font-variant-numeric: tabular-nums; }
  .card-label { font-size: 12px; color: var(--muted); margin-top: 3px; }
  .figure-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-bottom: 20px; }
  .figure-grid.wide { grid-template-columns: 1fr; }
  figure { margin: 0; background: var(--card-bg); border: 1px solid var(--border);
           border-radius: 8px; padding: 12px; }
  figure img { width: 100%; height: auto; display: block; }
  figure figcaption { font-size: 12px; color: var(--muted); margin-top: 8px; text-align: center; }
  .iframe-wrap { background: var(--card-bg); border: 1px solid var(--border);
                 border-radius: 8px; padding: 8px; }
  iframe { width: 100%; height: 620px; border: none; display: block; }
  .table-wrap { background: var(--card-bg); border: 1px solid var(--border);
                border-radius: 8px; overflow: hidden; }
  table { width: 100%; border-collapse: collapse; font-size: 13px; }
  thead th { background: #eef2f6; color: var(--accent); text-align: left;
             padding: 9px 12px; font-weight: 600; border-bottom: 1px solid var(--border); }
  tbody td { padding: 7px 12px; border-bottom: 1px solid #f0f2f5; }
  tbody tr:nth-child(even) { background: #fafbfc; }
  tbody tr:hover { background: #f0f5fa; }
  th.num, td.num { text-align: right; font-variant-numeric: tabular-nums; }
  p.note { font-size: 12px; color: var(--muted); margin-top: 8px; }
  ul.links { columns: 2; list-style: none; padding: 0; margin: 0; }
  ul.links li { margin: 4px 0; }
  ul.links a { color: var(--accent-light); text-decoration: none;
               font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 13px; }
  ul.links a:hover { text-decoration: underline; }
  .warnings { background: #fff8e6; border: 1px solid #f0d98c; border-radius: 8px;
              padding: 12px 16px; color: #7a5c00; font-size: 13px; }
  .warnings ul { margin: 0; padding-left: 18px; }
  footer { text-align: center; color: var(--muted); font-size: 12px; padding: 20px 0; }
  @media (max-width: 720px) { .figure-grid { grid-template-columns: 1fr; } }
</style>
</head>
<body>
<header>
  <h1>LINC Run Report</h1>
  <div class="subtitle">Generated at __GENERATED_AT__ &middot; output <code>__OUTPUT_DIR__</code></div>
</header>

<main>
  <section>
    <h2>Summary</h2>
    __CARDS__
  </section>

  <section>
    <h2>Pan-genome</h2>
    <div class="figure-grid">
      __PAN_FAMILIES__
      __PAN_GENES__
    </div>
    <div class="table-wrap">__PAN_TABLE__</div>
  </section>

  <section>
    <h2>Sequence similarity network</h2>
    <div class="iframe-wrap"><iframe src="ssn.html"></iframe></div>
  </section>

  <section>
    <h2>Distributions</h2>
    <div class="figure-grid wide">__CLUSTER_SIZE__</div>
    <div class="figure-grid">
      __MCL_BEFORE_AFTER__
      __MCL_TOP_COMPONENTS__
    </div>
  </section>

  <section>
    <h2>Top families</h2>
    <div class="table-wrap">__TOP_FAMILIES__</div>
  </section>

  <section>
    <h2>Top connected components</h2>
    <div class="table-wrap">__TOP_COMPONENTS__</div>
  </section>

  <section>
    <h2>Genomes</h2>
    <div class="table-wrap">__GENOME_TABLE__</div>
  </section>

  <section>
    <h2>Parameters</h2>
    <div class="table-wrap">__PARAMS_TABLE__</div>
  </section>

  <section>
    <h2>Output files</h2>
    <ul class="links">__LINKS__</ul>
  </section>

  <section>
    <h2>Warnings</h2>
    <div class="warnings"><ul>__WARNINGS__</ul></div>
  </section>
</main>

<footer>Generated by LINC pan</footer>
</body>
</html>
"""
    return (
        template
        .replace("__GENERATED_AT__", generated_at)
        .replace("__OUTPUT_DIR__", state.output_dir)
        .replace("__CARDS__", cards)
        .replace("__PAN_FAMILIES__", _figure("pan_genome_families.svg", "Families per pan-genome category"))
        .replace("__PAN_GENES__", _figure("pan_genome_genes.svg", "Genes per pan-genome category"))
        .replace("__PAN_TABLE__", pan_table)
        .replace("__CLUSTER_SIZE__", _figure("cluster_size_distribution.svg", "MCL family size distribution"))
        .replace("__MCL_BEFORE_AFTER__", _figure("mcl_before_after.svg", "Components before vs families after MCL"))
        .replace("__MCL_TOP_COMPONENTS__", _figure("mcl_top_components.svg", "Top components split by MCL"))
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
                    graphml: bool = False, top_n: int = 30) -> str:
    """完整生成报告, 返回报告目录路径。

    流程: load_state -> analyze -> write_tsvs -> render_charts ->
    render_ssn_html -> run_summary.json(两处) -> report.html。
    """
    state = load_state(output_dir)
    data = analyze(state, top_n=top_n)

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
