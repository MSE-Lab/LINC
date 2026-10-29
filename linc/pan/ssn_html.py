"""生成 Cytoscape.js 家族交互子图(ssn.html)与 GraphML。

设计依据: docs/LINC_结果描述增强开发计划.md 第 4.5 节。只做"按 `--family`
指定的家族诱导子图", 不做全图渲染(全图太大且呈毛线球, 价值低)。全量 SSN
数据由 analysis.py 的 `ssn_nodes.tsv` / `ssn_edges.tsv` 承载。
"""
import json
import os
from typing import Any, Dict, List, Optional, Set, Tuple

from ..core.entities import SSNEdge
from .analysis import ReportData
from .loader import LoadedState

# 基因组配色(确定性, 循环使用)。
_GENOME_PALETTE = [
    "#4c78a8", "#f58518", "#e45756", "#72b7b2", "#54a24b", "#eeca3b",
    "#b279a2", "#ff9da6", "#9d755d", "#bab0ac", "#9b6ff0", "#17becf",
]

_DEFAULT_TOP_FAMILIES = 10

_CYTOSCAPE_CDN = "https://unpkg.com/cytoscape@3.26.0/dist/cytoscape.min.js"


def _genome_colors(genome_names: List[str]) -> Dict[str, str]:
    return {
        name: _GENOME_PALETTE[i % len(_GENOME_PALETTE)]
        for i, name in enumerate(sorted(genome_names))
    }


def select_families(state: LoadedState, family_ids: Optional[List[str]] = None,
                    top_n: int = _DEFAULT_TOP_FAMILIES) -> List[Any]:
    """选出要展示的家族: 显式指定则按给定顺序, 否则取 size 最大的 top_n。"""
    if family_ids:
        by_id = {f.family_id: f for f in state.families}
        return [by_id[fid] for fid in family_ids if fid in by_id]
    return sorted(state.families, key=lambda f: (-f.size, f.family_id))[:top_n]


def _build_elements(state: LoadedState, data: ReportData,
                    families: List[Any], genome_colors: Dict[str, str]) -> Tuple[List[dict], List[dict]]:
    """构建 Cytoscape 元素(nodes/edges)。node id 用原始序列 ID, 全局唯一。"""
    node_ids: Set[str] = set()
    for fam in families:
        node_ids |= fam.member_genes

    nodes = []
    for gene_id in sorted(node_ids):
        orig = state.gene_id_to_original.get(gene_id, gene_id)
        gid = data.gene_to_genome.get(gene_id, "")
        genome_name = state.genome_id_to_name.get(gid, gid)
        nodes.append({
            "data": {
                "id": orig,
                "label": orig,
                "genome": genome_name,
                "family": data.gene_to_family.get(gene_id, ""),
                "degree": data.degree.get(gene_id, 0),
                "is_rep": 3 if gene_id in data.representative_genes else 0,
                "color": genome_colors.get(genome_name, "#999"),
            },
        })

    edges = []
    if state.complete_ssn is not None:
        for e in sorted(state.complete_ssn.edges, key=lambda e: (e.source, e.target)):
            if e.source in node_ids and e.target in node_ids:
                edges.append({
                    "data": {
                        "id": f"{e.source}__{e.target}",
                        "source": state.gene_id_to_original.get(e.source, e.source),
                        "target": state.gene_id_to_original.get(e.target, e.target),
                        "identity": round(e.identity, 4),
                        "coverage": round(e.coverage, 4),
                    },
                })
    return nodes, edges


def _cytoscape_html(title: str, nodes: List[dict], edges: List[dict],
                    genome_names: List[str]) -> str:
    """生成自包含 Cytoscape.js HTML。"""
    elements = json.dumps([{"data": n["data"]} for n in nodes]
                          + [{"data": e["data"]} for e in edges],
                          ensure_ascii=False)
    genome_options = "\n".join(
        f'<option value="{g}">{g}</option>' for g in sorted(genome_names)
    )
    legend_items = "\n".join(
        f'<span style="display:inline-block;width:12px;height:12px;'
        f'background:{_GENOME_PALETTE[i % len(_GENOME_PALETTE)]};margin-right:4px"></span>{g}'
        for i, g in enumerate(sorted(genome_names))
    )

    return """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>%s</title>
<script src="%s"></script>
<style>
  body { font-family: sans-serif; margin: 0; }
  #toolbar { padding: 8px; border-bottom: 1px solid #ccc; background: #fafafa; }
  #toolbar button, #toolbar select { margin-right: 10px; }
  #cy { position: absolute; top: 52px; bottom: 0; left: 0; right: 0; }
  #info { position: absolute; top: 60px; right: 10px; background: white;
          border: 1px solid #ccc; padding: 8px; max-width: 360px; font-size: 12px; }
  #legend { position: absolute; bottom: 10px; left: 10px; background: white;
            border: 1px solid #ccc; padding: 8px; font-size: 12px; line-height: 1.6; }
</style>
</head>
<body>
<div id="toolbar">
  <button onclick="cy.fit()">Fit</button>
  <button onclick="cy.layout({name:'cose'}).run()">Re-layout</button>
  <label>Genome:
    <select id="genome-filter" onchange="filterGenome(this.value)">
      <option value="">all</option>
      %s
    </select>
  </label>
  <label><input type="checkbox" id="label-toggle" checked onchange="toggleLabels(this.checked)"> labels</label>
</div>
<div id="cy"></div>
<div id="info">Click a node for details</div>
<div id="legend">%s</div>
<script>
var elements = %s;
var cy = cytoscape({
  container: document.getElementById('cy'),
  elements: elements,
  style: [
    { selector: 'node', style: {
        'background-color': 'data(color)',
        'width': 'mapData(degree, 0, 20, 8, 30)',
        'border-width': 'data(is_rep)',
        'border-color': '#222',
        'label': 'data(label)',
        'font-size': 8,
        'text-valign': 'bottom',
        'text-margin-y': 2
    }},
    { selector: 'edge', style: {
        'width': 1,
        'line-color': '#ccc',
        'target-arrow-color': '#ccc',
        'target-arrow-shape': 'triangle',
        'curve-style': 'bezier'
    }}
  ],
  layout: { name: 'cose' }
});

function filterGenome(g) {
  cy.batch(function() {
    cy.nodes().forEach(function(n) {
      var show = (g === '' || n.data('genome') === g);
      n.style('display', show ? 'element' : 'none');
    });
  });
}

function toggleLabels(on) {
  cy.nodes().style('label', on ? 'data(label)' : '');
}

cy.on('tap', 'node', function(evt) {
  var n = evt.target;
  document.getElementById('info').innerHTML =
    '<b>' + n.data('id') + '</b><br>' +
    'genome: ' + n.data('genome') + '<br>' +
    'family: ' + n.data('family') + '<br>' +
    'degree: ' + n.data('degree') + '<br>' +
    'representative: ' + (n.data('is_rep') ? 'yes' : 'no');
});
</script>
</body>
</html>
""" % (title, _CYTOSCAPE_CDN, genome_options, legend_items, elements)


def _graphml_xml(state: LoadedState, families: List[Any],
                 nodes: List[dict], edges: List[dict]) -> str:
    """生成所选家族诱导子图的 GraphML。"""
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">',
        '  <key id="genome" for="node" attr.name="genome" attr.type="string"/>',
        '  <key id="family" for="node" attr.name="family" attr.type="string"/>',
        '  <key id="identity" for="edge" attr.name="identity" attr.type="double"/>',
        '  <key id="coverage" for="edge" attr.name="coverage" attr.type="double"/>',
        '  <graph id="G" edgedefault="undirected">',
    ]
    for n in nodes:
        d = n["data"]
        lines.append(
            f'    <node id="{d["id"]}">'
            f'<data key="genome">{d["genome"]}</data>'
            f'<data key="family">{d["family"]}</data></node>'
        )
    for e in edges:
        d = e["data"]
        lines.append(
            f'    <edge source="{d["source"]}" target="{d["target"]}">'
            f'<data key="identity">{d["identity"]}</data>'
            f'<data key="coverage">{d["coverage"]}</data></edge>'
        )
    lines.append("  </graph>")
    lines.append("</graphml>")
    return "\n".join(lines) + "\n"


def render_ssn_html(state: LoadedState, data: ReportData, report_dir: str,
                    family_ids: Optional[List[str]] = None,
                    graphml: bool = False) -> List[str]:
    """生成 ssn.html(及可选 GraphML), 返回写入的文件名列表。"""
    os.makedirs(report_dir, exist_ok=True)
    written: List[str] = []

    families = select_families(state, family_ids)
    genome_names = sorted(state.genome_id_to_name.values())
    genome_colors = _genome_colors(genome_names)
    nodes, edges = _build_elements(state, data, families, genome_colors)

    title = "LINC SSN - " + (
        ", ".join(f.family_id for f in families) if family_ids else "top families"
    )
    html = _cytoscape_html(title, nodes, edges, genome_names)
    with open(os.path.join(report_dir, "ssn.html"), "w", encoding="utf-8") as fh:
        fh.write(html)
    written.append("ssn.html")

    if graphml:
        xml = _graphml_xml(state, families, nodes, edges)
        with open(os.path.join(report_dir, "ssn_view.graphml"), "w", encoding="utf-8") as fh:
            fh.write(xml)
        written.append("ssn_view.graphml")

    return written
