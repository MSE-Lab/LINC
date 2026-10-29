"""纯标准库 SVG 图表生成(不依赖 matplotlib)。

设计依据: docs/LINC_结果描述增强开发计划.md 第 4.4 节。产出 5 张 SVG:
cluster_size_distribution / mcl_before_after / mcl_top_components /
pan_genome_families / pan_genome_genes。所有图都是自包含 SVG, 可直接嵌入
HTML、可缩放。
"""
import os
from typing import Any, Dict, List, Tuple

from .analysis import PAN_CATEGORIES, ReportData, bin_sizes
from .loader import LoadedState

# 五类泛基因组分类的固定配色(确定性, 不随数据变化)。
_CATEGORY_COLORS = {
    "core": "#4c78a8",
    "soft-core": "#59a14f",
    "shell": "#edc948",
    "cloud": "#ff9da7",
    "singleton": "#b07aa1",
}
# 家族子图堆叠用的一串可区分颜色。
_STACK_COLORS = [
    "#4c78a8", "#f58518", "#e45756", "#72b7b2", "#54a24b",
    "#eeca3b", "#b279a2", "#ff9da6", "#9d755d", "#bab0ac",
]


def _esc(text: Any) -> str:
    """转义 XML 特殊字符(标签多为普通文本, 防御性处理)。"""
    return (
        str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    )


def _svg_doc(width: int, height: int, title: str, body: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        f'  <rect width="{width}" height="{height}" fill="white"/>\n'
        f'  <text x="{width / 2}" y="30" text-anchor="middle" font-size="18" '
        f'font-family="sans-serif">{_esc(title)}</text>\n'
        f"{body}"
        f"</svg>\n"
    )


def _bar_chart(title: str, labels: List[str], values: List[float],
               color: str = "#4c78a8") -> str:
    """单系列竖直柱状图。"""
    width, height = 800, 480
    ml, mr, mt, mb = 70, 20, 60, 70
    plot_w = width - ml - mr
    plot_h = height - mt - mb
    n = len(labels)
    if n == 0:
        return _svg_doc(width, height, title, "")
    max_v = max(values) if values else 1
    slot = plot_w / n
    bar_w = slot * 0.7

    body = []
    for i, (label, value) in enumerate(zip(labels, values)):
        h = (value / max_v * plot_h) if max_v else 0.0
        x = ml + i * slot + (slot - bar_w) / 2
        y = mt + plot_h - h
        cx = x + bar_w / 2
        body.append(f'  <rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" '
                    f'height="{h:.1f}" fill="{color}"/>')
        body.append(f'  <text x="{cx:.1f}" y="{y - 5:.1f}" text-anchor="middle" '
                    f'font-size="12" font-family="sans-serif">{_esc(value)}</text>')
        body.append(f'  <text x="{cx:.1f}" y="{height - mb + 20:.1f}" text-anchor="middle" '
                    f'font-size="11" font-family="sans-serif">{_esc(label)}</text>')
    return _svg_doc(width, height, title, "\n".join(body) + "\n")


def _grouped_bar_chart(title: str, labels: List[str],
                       series_a: List[float], series_b: List[float],
                       label_a: str, label_b: str) -> str:
    """双系列分组柱状图(MCL 前后对比)。"""
    width, height = 800, 480
    ml, mr, mt, mb = 70, 20, 60, 90
    plot_w = width - ml - mr
    plot_h = height - mt - mb
    n = len(labels)
    if n == 0:
        return _svg_doc(width, height, title, "")
    max_v = max(max(series_a, default=0), max(series_b, default=0)) or 1
    slot = plot_w / n
    group_w = slot * 0.8
    bar_w = group_w / 2

    body = []
    for i, label in enumerate(labels):
        gx = ml + i * slot + (slot - group_w) / 2
        for j, (series, color) in enumerate(((series_a, "#4c78a8"), (series_b, "#f58518"))):
            value = series[i]
            h = (value / max_v * plot_h) if max_v else 0.0
            x = gx + j * bar_w
            y = mt + plot_h - h
            body.append(f'  <rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" '
                        f'height="{h:.1f}" fill="{color}"/>')
        body.append(f'  <text x="{gx + group_w / 2:.1f}" y="{height - mb + 20:.1f}" '
                    f'text-anchor="middle" font-size="11" font-family="sans-serif">{_esc(label)}</text>')

    # 图例
    ly = height - mb + 40
    body.append(f'  <rect x="{ml:.0f}" y="{ly - 12:.0f}" width="14" height="14" fill="#4c78a8"/>')
    body.append(f'  <text x="{ml + 20:.0f}" y="{ly:.0f}" font-size="12" font-family="sans-serif">{_esc(label_a)}</text>')
    body.append(f'  <rect x="{ml + 140:.0f}" y="{ly - 12:.0f}" width="14" height="14" fill="#f58518"/>')
    body.append(f'  <text x="{ml + 160:.0f}" y="{ly:.0f}" font-size="12" font-family="sans-serif">{_esc(label_b)}</text>')
    return _svg_doc(width, height, title, "\n".join(body) + "\n")


def _stacked_hbar(title: str, rows: List[Tuple[str, List[Tuple[str, float]]]]) -> str:
    """水平堆叠条: 每行 = 一个分量, 段 = 其内部亚家族。"""
    width, height = 800, max(200, 80 + 40 * len(rows))
    ml, mr, mt, mb = 60, 220, 60, 20
    plot_w = width - ml - mr
    row_h = 30

    body = []
    for i, (label, segments) in enumerate(rows):
        total = sum(v for _, v in segments) or 1
        y = mt + i * (row_h + 10)
        x = ml
        for j, (seg_label, seg_value) in enumerate(segments):
            w = seg_value / total * plot_w
            color = _STACK_COLORS[j % len(_STACK_COLORS)]
            body.append(f'  <rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" '
                        f'height="{row_h:.1f}" fill="{color}"/>')
            x += w
        body.append(f'  <text x="{ml - 8:.1f}" y="{y + row_h / 2 + 4:.1f}" text-anchor="end" '
                    f'font-size="11" font-family="sans-serif">{_esc(label)}</text>')
        body.append(f'  <text x="{ml + plot_w + 8:.1f}" y="{y + row_h / 2 + 4:.1f}" '
                    f'font-size="11" font-family="sans-serif">{total} genes</text>')
    return _svg_doc(width, height, title, "\n".join(body) + "\n")


# ---------------------------------------------------------------------------
# 5 张图
# ---------------------------------------------------------------------------
def _cluster_size_distribution(data: ReportData) -> str:
    bins = data.cluster_size_distribution["bins"]
    return _bar_chart(
        "MCL family size distribution",
        [b["label"] for b in bins],
        [b["n"] for b in bins],
    )


def _mcl_before_after(state: LoadedState, data: ReportData) -> str:
    component_sizes = sorted((len(c) for c in data.components.values()))
    family_sizes = sorted((f.size for f in state.families))
    labels = [label for label, _ in bin_sizes(component_sizes)]
    comp_counts = [n for _, n in bin_sizes(component_sizes)]
    fam_counts = [n for _, n in bin_sizes(family_sizes)]
    return _grouped_bar_chart(
        "MCL before (components) vs after (families) by size",
        labels, comp_counts, fam_counts,
        "components", "families",
    )


def _mcl_top_components(state: LoadedState, data: ReportData) -> str:
    # 按分量大小取前 20 个, 每个分量内部按家族拆分(最多 20 个亚家族)。
    comp_to_fams: Dict[int, List[Tuple[str, int]]] = {}
    for fam in state.families:
        comp_to_fams.setdefault(fam.component_id, []).append((fam.family_id, fam.size))

    top = sorted(data.components.items(), key=lambda kv: -len(kv[1]))[:20]
    rows = []
    for cid, genes in top:
        segs = sorted(comp_to_fams.get(cid, []), key=lambda s: -s[1])[:20]
        rows.append((f"comp {cid}", [(fid, size) for fid, size in segs]))
    return _stacked_hbar("Top components split into subfamilies by MCL", rows)


def _pan_genome_multi_color(data: ReportData, key: str, title: str) -> str:
    """泛基因组五类柱状图, 每类用固定配色。"""
    cats = data.pan_genome["categories"]
    labels = PAN_CATEGORIES
    values = [cats[cat][key] for cat in labels]

    width, height = 800, 480
    ml, mr, mt, mb = 70, 20, 60, 70
    plot_w = width - ml - mr
    plot_h = height - mt - mb
    n = len(labels)
    max_v = max(values) if values else 1
    slot = plot_w / n
    bar_w = slot * 0.7

    body = []
    for i, (label, value) in enumerate(zip(labels, values)):
        h = (value / max_v * plot_h) if max_v else 0.0
        x = ml + i * slot + (slot - bar_w) / 2
        y = mt + plot_h - h
        cx = x + bar_w / 2
        color = _CATEGORY_COLORS[label]
        body.append(f'  <rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" '
                    f'height="{h:.1f}" fill="{color}"/>')
        body.append(f'  <text x="{cx:.1f}" y="{y - 5:.1f}" text-anchor="middle" '
                    f'font-size="12" font-family="sans-serif">{_esc(value)}</text>')
        body.append(f'  <text x="{cx:.1f}" y="{height - mb + 20:.1f}" text-anchor="middle" '
                    f'font-size="11" font-family="sans-serif">{_esc(label)}</text>')
    return _svg_doc(width, height, title, "\n".join(body) + "\n")


def _pan_genome_families(data: ReportData) -> str:
    return _pan_genome_multi_color(data, "n_families", "Pan-genome families by category")


def _pan_genome_genes(data: ReportData) -> str:
    return _pan_genome_multi_color(data, "n_genes", "Pan-genome genes by category")


# ---------------------------------------------------------------------------
# 写盘
# ---------------------------------------------------------------------------
def render_charts(state: LoadedState, data: ReportData, report_dir: str) -> List[str]:
    """把 5 张 SVG 写入 report_dir, 返回文件名列表。"""
    os.makedirs(report_dir, exist_ok=True)
    charts = [
        ("cluster_size_distribution.svg", _cluster_size_distribution(data)),
        ("mcl_before_after.svg", _mcl_before_after(state, data)),
        ("mcl_top_components.svg", _mcl_top_components(state, data)),
        ("pan_genome_families.svg", _pan_genome_families(data)),
        ("pan_genome_genes.svg", _pan_genome_genes(data)),
    ]
    for name, svg in charts:
        with open(os.path.join(report_dir, name), "w", encoding="utf-8") as fh:
            fh.write(svg)
    return [name for name, _ in charts]
