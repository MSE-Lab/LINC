# linc pan — 结果描述/可视化

`linc pan` 从一次已完成的 `linc run` 输出中只读生成报告与可视化,不重跑任何
Stage,不修改任何 `.linc_work/state`。

## 用法

```bash
linc pan -o <output_dir>
```

`<output_dir>` 是 `linc run -o` 的输出目录(其中含 `.linc_work/state/`)。

可选参数:

```text
--report-dir <dir>   自定义报告子目录(默认 <output_dir>/report)
--family <family_id> 只渲染指定家族的交互子图(可重复,默认取最大的家族)
--graphml            额外导出所选家族的 ssn_view.graphml
```

## 输出

报告写入 `<output_dir>/report/`(另在 `<output_dir>/` 根写一份
`run_summary.json`):

| 文件 | 内容 |
|---|---|
| `run_summary.json` | 统一统计(参数、counts、derived、ssn、families、pan_genome、mcl、top 列表等) |
| `report.html` | 综合 HTML 报告(指标卡片 + 图 + 表格) |
| `ssn.html` | Cytoscape.js 家族交互子图(需联网加载 CDN) |
| `pan_genome_families.svg` / `pan_genome_genes.svg` | 泛基因组五类占比 |
| `cluster_size_distribution.svg` / `mcl_before_after.svg` / `mcl_top_components.svg` | 家族大小 / MCL 前后 / 分量拆分 |
| `presence_absence_matrix.tsv` / `per_genome_breakdown.tsv` / `pan_genome_summary.tsv` | 泛基因组 TSV |
| `ssn_nodes.tsv` / `ssn_edges.tsv` | 全量 SSN 节点/边 |
| `family_members.tsv` / `family_stats.tsv` | 家族成员/统计 |
| `mcl_components.tsv` / `cluster_size_distribution.tsv` | 分量/分箱 |

## 模块

- `loader.py` — 只读加载 state + ID 映射;
- `analysis.py` — 泛基因组分类(迁移自 cli.py)+ SSN/家族/MCL 统计 + TSV;
- `charts.py` — 纯标准库 SVG;
- `ssn_html.py` — Cytoscape.js 家族子图 + GraphML;
- `report.py` — 编排,写 `run_summary.json` / `report.html`。

## 已知限制

- `ssn.html` 的 Cytoscape.js 来自 CDN,离线时交互图不可用(TSV/GraphML 仍可用);
- 报告会加载较大的 pickle,数百基因组时内存与耗时都不低;
- 交互视图只覆盖 `--family` 指定的家族,不做全图;完整 SSN 以 TSV 为准;
- `run_summary.json` 由 `linc pan` 生成,不是 `linc run` 自动生成。
