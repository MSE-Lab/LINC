# LINC 结果描述增强开发计划

> 状态: **待实现**(2026-10-29)。本文是开发计划,尚未编写任何代码。
> 目标: 扩展现有 `linc pan` 子程序,使其成为统一的结果描述/可视化入口,只读 `.linc_work/state`,不重跑任何 Stage。
> 实现载体: 复用现有 `linc pan` 子命令(而非新建 `linc/reporting/` 模块),代码落到 `linc/pan/` 包。

---

## 1. 目标与边界

### 1.1 现有产出的不足

当前 `linc run` 的主要结果是:

- `Orthogroups.txt`:一行一个 MCL 家族,扁平列表;
- `family_summary.json`:只有 `family_id / component_id / size / edge_density`;
- `ssn_summary.json`:只有 `edge_completion_mode + n_edges`;
- `join_summary.json`:有全局簇和代表数,但没有节点/边明细和可视化。

由此带来的问题:

1. 没有 SSN 图,无法直观地看簇之间的关系;
2. 没有 family 成员、基因组分布、identity/coverage 等统计;
3. MCL 前后没有对比,无法判断"亚家族细分"是否发生;
4. 没有统一的泛基因组 core/soft-core/shell/cloud 分类(现有 `linc pan` 有,但不在统一报告里);
5. 没有统一的 `run_summary.json` 和 HTML 报告。

### 1.2 目标

扩展现有 `linc pan` 子程序,使其成为统一的结果描述/可视化入口,只读取 `.linc_work/state`,不重跑任何 Stage,不修改任何 state 文件。

### 1.3 边界

- 不修改 `linc run` 主流程(`linc/pipeline.py`、`linc/engine/*`、以及 `linc/cli.py` 里的 run/append/resume);
- `linc pan` 子程序是本次的实现载体,会被扩展:现有泛基因组分类逻辑迁入 `linc/pan/` 包,`pan` 命令本身改为薄封装;
- `linc pan` **移除 `-g/--orthogroups`、`-m/--gene-genome` 参数**,直接读取 `<output_dir>/.linc_work/state/`;不保留"对任意 Orthogroups.txt 独立做泛基因组"的兼容模式;
- 报告生成失败不影响 `linc run`(报告是独立命令,可失败、可重跑)。

---

## 2. 总体架构

新增包(承载实现):

```text
linc/pan/
  __init__.py
  loader.py
  analysis.py
  charts.py
  ssn_html.py
  report.py
  README.md
```

其中 `analysis.py` 承接从 `linc/cli.py` 迁入的泛基因组分类逻辑(`_classify_family` / `_PAN_CATEGORIES` / `_build_pan_genome_data` 及三个 TSV writer)。

数据流:

```text
linc run 的 state(<output_dir>/.linc_work/state)
  → loader.py(FileRepository 只读加载)
  → analysis.py(泛基因组分类 + 统计 + 连通分量)
  → charts.py(SVG) / ssn_html.py(家族交互子图)
  → report.py(run_summary.json / TSV / report.html)
```

入口(复用现有子命令,`linc pan` 命令改为薄封装):

```bash
linc pan -o <output_dir>
```

---

## 3. 模块职责

| 文件 | 职责 |
|---|---|
| `loader.py` | 通过 `FileRepository` 只读加载 state;缺文件抛 `ReportInputError`;处理 `SchemaVersionMismatchError`(见 8.1) |
| `analysis.py` | 承接自 `cli.py` 迁入的泛基因组分类;构建 `GeneMeta` / `ReportData`;union-find 连通分量;统计;TSV;交互子图选择;GraphML |
| `charts.py` | 纯标准库 SVG:家族大小分布、MCL 前后、最大分量拆分、泛基因组 family/gene 组成 |
| `ssn_html.py` | 生成 Cytoscape.js 家族交互子图(`--family` 指定;缩放、点击、过滤、图例) |
| `report.py` | 编排全部输出,写 `run_summary.json`、`report.html`、`ssn.html` |
| `README.md` | 使用说明和输出清单 |

> 说明:`linc/cli.py` 里的 `pan` 命令改为薄封装(解析参数 + 调 `linc/pan/` 的 `report.py`),不再内联泛基因组/统计逻辑。

---

## 4. 输出规格

### 4.1 `run_summary.json`

同时写到 `<output>/run_summary.json` 和 `<output>/report/run_summary.json`。

顶层字段:

```text
schema
generated_at
output_dir
parameters        # 含 family_mode / edge_completion_mode / inflation / cutoffs 等
counts
derived
ssn
families
pan_genome
mcl
cluster_size_distribution
genomes
top_families
top_components
visualization
warnings
```

#### 4.1.1 `counts`

- genome / gene / local cluster / local representative;
- global cluster / global representative;
- SSN edge;
- MCL component / family;
- 有 family / 无 family 基因数。

#### 4.1.2 `derived`

- 平均 genes/genome、平均 local reps/genome、平均 global reps/cluster;
- 平均 family size;
- 最大连通分量大小及占比;
- 最大 MCL family;
- median degree;
- isolated gene fraction;
- cross-genome edge fraction;
- 被 MCL 拆分的分量数;
- 单个分量最大亚家族数。

#### 4.1.3 `ssn`

- nodes / edges / components / density;
- degree 分布;
- isolated nodes;
- cross-genome edge fraction;
- same-family edge fraction;
- identity / coverage / weight 的 mean / median / min / max;
- edge completion mode;
- **family_mode(见 8.2, collapsed 模式下注明"代表图口径")**。

#### 4.1.4 `mcl`

- inflation、edge weight mode;
- MCL 前连通分量统计;
- MCL 后 family 统计;
- 被 MCL 拆分的分量数;
- 单分量最大亚家族数。

#### 4.1.5 `top_families` / `top_components`

- 最大 30 个 family:大小、基因组数、边密度、identity、coverage、代表基因;
- 最大 30 个连通分量:gene 数、family 数、最大 family、family 拆分。

#### 4.1.6 `genomes`

每个基因组:gene 数、family 数、singleton gene 数、global cluster 数。

### 4.2 报告指标卡片(目标 22 个)

- Genomes / Genes / Avg genes per genome;
- Local representatives / Avg local reps per genome;
- Global clusters / Global representatives / Avg global reps per cluster;
- SSN nodes / SSN edges / SSN components / Largest component / Isolated genes;
- MCL families / Singleton families / Median family size / Largest family / Components split by MCL;
- Core / Soft-core / Shell / Cloud genes。

> 派生指标(如 largest component %、median degree、mean family size)保留在 `run_summary.json` 的 `derived` 中,卡片是否展示由实现阶段决定,不预先写死。

### 4.3 泛基因组分类

**复用现有 `linc pan` 的分类逻辑**(迁到 `linc/pan/analysis.py`),采用其 5 类互斥约定:

```text
core      : 覆盖基因组比例 >= 99%
soft-core : 95% – 99%
shell     : 15% – 95%
cloud     : < 15%
singleton : 家族总成员数为 1(优先于上面的比例判断)
```

> 即沿用 `linc pan` 现有的 `_PAN_CATEGORIES = ["core","soft-core","shell","cloud","singleton"]`,不另立"singleton 作为 cloud 子集"的约定。

在此基础上的额外派生统计:

- `single_copy_core`:core 中每个出现基因组拷贝数都为 1 的 family,占 core 的比例;
- `cloud_including_singleton`:cloud + singleton 合并视图(派生列,不改变基础 5 类分类);
- 每类的 `n_families / family_pct / n_genes / gene_pct / avg_genes_per_genome / avg_families_per_genome`。

输出(前三个与 `linc pan` 现有输出保持一致,后两个为新增图):

- `presence_absence_matrix.tsv`、`per_genome_breakdown.tsv`、`pan_genome_summary.tsv`;
- `pan_genome_families.svg`、`pan_genome_genes.svg`。

### 4.4 SVG 图表(纯标准库,不依赖 matplotlib)

1. `cluster_size_distribution.svg`:MCL family 大小分布;singleton 比例、median、max;
2. `mcl_before_after.svg`:MCL 前连通分量 vs MCL 后 family 的分箱对比;
3. `mcl_top_components.svg`:最大连通分量的内部亚家族堆叠;单分量最多画 20 个亚家族,其余聚合为 other;
4. `pan_genome_families.svg`:Core/Soft-core/Shell/Cloud/Singleton 的 family 占比,柱上标注 family 数和百分比;
5. `pan_genome_genes.svg`:同样五类的 gene 占比。

### 4.5 家族交互子图(`ssn.html`)

> **设计取舍**:SSN 全图通常有几十万~几百万条边,直接可视化既渲染不动、又是"毛线球",肉眼提取不出生物学信息。因此**不做全图交互渲染**,只做**按家族的子图**——一个家族通常几十~几百个基因,可读、有生物学意义。全量 SSN 通过 `ssn_nodes.tsv` / `ssn_edges.tsv` 导出,交给 Cytoscape 桌面版 / Gephi 等外部工具。

- 按 `--family <family_id>` 指定一个或多个家族,生成这些家族的诱导子图(家族内基因 + 家族内边);
- 节点 = 基因;颜色 = 基因组;节点大小 = degree;代表节点有边框;
- 边 = 家族内的 `CompleteSSN.edges`,标注 identity / coverage;
- 支持:滚轮缩放、拖拽平移、点击节点/边查看详情、genome 过滤、edge weight 过滤、标签开关、Fit/Reset、图例;
- 分量级总览由静态图 `mcl_top_components.svg` 承担(见 4.4),不做基因级全图;
- `--graphml` 额外导出所选家族的 `ssn_view.graphml`;
- 全量 `ssn_nodes.tsv` / `ssn_edges.tsv` 照常导出(供外部工具使用)。

### 4.6 `report.html`

1. 指标卡片(22 个);
2. 泛基因组分析(family/gene 两张 SVG + `pan_genome_summary.tsv` 表);
3. 家族交互子图 iframe;
4. 簇大小分布图;
5. MCL 前后对比图;
6. Top families 表;
7. 被 MCL 拆分最多的连通分量表;
8. 每个基因组概况(前 100 个);
9. 参数表;
10. 所有结果文件链接;
11. 警告提示。

CSS 自包含;除 `ssn.html` 的 Cytoscape.js CDN 外,不依赖外部资源。

### 4.7 TSV 导出

| 文件 | 内容 |
|---|---|
| `presence_absence_matrix.tsv` | 行=family,列=基因组,值=拷贝数(沿用现有 `linc pan`) |
| `per_genome_breakdown.tsv` | 每个基因组的五类基因数(沿用现有 `linc pan`) |
| `pan_genome_summary.tsv` | 泛基因组五类 + single-copy core(沿用现有 `linc pan`) |
| `ssn_nodes.tsv` | 全部基因节点:family、global cluster、local cluster、degree、代表标记、长度 |
| `ssn_edges.tsv` | 全部 SSN 边:identity、coverage、bitscore、weight、same_family |
| `family_members.tsv` | 每个 MCL family 的成员明细 |
| `family_stats.tsv` | family 大小、基因组数、边密度、平均 identity/coverage、代表 |
| `mcl_components.tsv` | MCL 前连通分量及亚家族拆分 |
| `cluster_size_distribution.tsv` | 分量/family 分箱计数 |

---

## 5. 命令行入口

```bash
linc pan -o <output_dir>
```

其中 `<output_dir>` 为 `linc run -o` 的输出目录,state 位于 `<output_dir>/.linc_work/state/`。

> **输入方式(已定)**:`linc pan` 直接读取 `linc run` 的输出(state),**移除现有 `-g/--orthogroups`、`-m/--gene-genome` 两个参数**。SSN/家族/MCL 等数据只在 state 里,文本文件无法支撑,故不保留"对任意 Orthogroups.txt 独立做泛基因组"的兼容模式。

可选参数:

```text
--report-dir <dir>     # 自定义报告子目录
--top-families <int>   # top family 数量(默认 30)
--family <family_id>   # 可重复,生成指定 family 的交互子图
--graphml              # 额外导出所选家族的 GraphML
--quiet
```

---

## 6. 关键设计决策

1. **独立 pan 命令**:报告/可视化逻辑放在 `linc/pan/`,不把 matplotlib/Cytoscape 依赖塞进 `linc run` 主流程;报告失败不影响 `linc run`。
2. **纯 SVG**:不依赖 matplotlib;SVG 可直接嵌入 HTML、可缩放。
3. **Cytoscape.js CDN**:交互能力成熟;代价是离线时交互图不可用,TSV/GraphML 仍可用。
4. **不做全图交互渲染,只做家族子图**:SSN 全图太大且呈"毛线球",全图渲染价值低;交互视图只做 `--family` 指定的家族子图,全量 SSN 通过 TSV/GraphML 交给外部工具。
5. **泛基因组分类复用现有 `linc pan` 逻辑**:把 `cli.py` 里的 `_classify_family` 等迁到 `linc/pan/analysis.py`,报告与 `linc pan` 共用同一套 5 类互斥分类,不另立约定。
6. **报告只读 state**:不重新计算 intra/join/SSN/MCL,不修改任何 `.linc_work/state`。

---

## 7. 实现步骤(按依赖排序)

| Step | 内容 | 产出 | 依赖 |
|---|---|---|---|
| 1 | `loader.py`:FileRepository 只读加载;缺文件抛 `ReportInputError`;处理 `SchemaVersionMismatchError`(见 8.1) | loader + 错误处理 | — |
| 2 | `analysis.py`:迁移自 `cli.py` 的泛基因组分类 + `GeneMeta`/`ReportData` + union-find 连通分量 + 统计 + TSV | analysis + 全部 TSV | Step 1 |
| 3 | `charts.py`:5 张纯标准库 SVG | 5 张 SVG | Step 2 |
| 4 | `ssn_html.py`:Cytoscape.js 家族交互子图(`--family`)+ GraphML | ssn.html + ssn_view.graphml | Step 2 |
| 5 | `report.py`:编排,写 `run_summary.json` / `report.html` | report.html + run_summary.json | Step 2–4 |
| 6 | `cli.py` 的 `pan` 命令改为薄封装 + `README.md` | `linc pan` 可命令行调用 | Step 5 |
| 7 | 单元测试(见 9) | tests/ | Step 1–6 |

建议实现顺序:先 1→2→3(拿到统计和静态图,可尽早验证口径),再做 4→5(交互图 + HTML 编排)、6(接线到 `linc pan`),最后补 7。

---

## 8. 实现前必须想清楚的坑

### 8.1 SCHEMA_VERSION 兼容

- `loader.py` 走 `FileRepository` 读 state,会触发 schema 版本检查(当前 `SCHEMA_VERSION = 3`)。
- **决策**:遇到 `SchemaVersionMismatchError` 时,报告层**报错退出**,并提示"请用当前代码版本重跑 `linc run` 再生成报告"。
- 第一版**不引入** `--ignore-schema` 之类的绕过开关(保持安全、简单)。

### 8.2 与 collapsed 模式的交互(最容易被忽略)

collapsed 模式下:

- `complete_ssn.pkl` 是**代表图**(折叠基因的边已被 `_finalize` 删掉);
- `families.pkl` 是**展开后的全部基因**。

因此报告会看到 "SSN 节点数 << 家族基因数" 的错位。**必须显式处理**,否则误导用户:

- `run_summary.json` 增加 `family_mode` 字段;
- `ssn` / `top_components` 等 SSN 口径的统计,在 collapsed 模式下明确标注"代表图口径(折叠模式)";
- 第一版**不尝试识别单个折叠基因**(需要重算 `collapse_map`,而它不落盘),只标注口径;后续如需,可在报告里加"折叠成员"视图,从 `intra_results + join_result` 重算 `collapse_map`。

### 8.3 gene ID 格式

- state 内部的 `gene_id` 与对外输出(`Orthogroups.txt` / `gene_genome.tsv`)的 ID 必须保持一致。
- 实现前先确认当前 Scheme B 重编码的映射关系(`genome_id` 内部是 `G00001`,对外是原始 `GCA_...`),避免报告里的基因 ID 和 `Orthogroups.txt` 对不上。

### 8.4 内存与惰性加载

- 几百基因组时,一次加载 `genomes.pkl` / `intra_results.pkl` / `join_result.pkl` / `complete_ssn.pkl` / `families.pkl` 内存和时间都不低。
- **建议**:loader 分层加载——大部分统计只需 `families` + `complete_ssn`;只有需要节点元数据(基因长度、基因组归属)时才读 `genomes`;`intra_results` 仅在需要 local cluster / local representative 信息时读。

### 8.5 边权口径

- `SSNEdge` 存 identity / coverage / bitscore 三个原值,`weight` 需按 `edge_weight_mode`(默认 identity×coverage)**现算**,与 `FamilyRefineStage` 的算法保持一致,不能直接假设 weight 已存。

### 8.6 `run_summary.json` 的生成时机

- `run_summary.json` 由报告命令生成,不是 `linc run` 自动生成。
- **建议(可选,不阻塞)**:在 `linc run` 结束时打印一行提示"可用 `linc pan -o ...` 生成可视化报告",避免用户以为跑完就有完整报告。

---

## 9. 测试计划

- `loader`:缺文件抛 `ReportInputError`;schema 不匹配抛错;正常加载;
- `analysis`:连通分量;泛基因组分类边界(99% / 95% / 15%);single-copy core 判定;TSV 字段;
- `charts`:5 张 SVG 生成;空数据、单家族、无 edge 的边界;
- `ssn_html`:`--family` 子图生成(单个/多个家族);GraphML 导出;未指定 `--family` 时的行为;
- `report`:输出文件完整性;与 collapsed 模式同跑时口径标注正确;
- 回归:确认 `linc run` 主流程不受影响(报告层独立,不 import 主流程模块之外的副作用)。

---

## 10. 验收标准

1. `linc pan -o <output_dir>` 能在现有输出上生成完整报告(JSON + TSV + SVG + HTML);
2. 报告过程**不修改**任何 `.linc_work/state` 文件;
3. 在 Streptomyces(跨物种)和 E. coli 32(同物种)两套数据上都能跑通;
4. collapsed 模式下 SSN 口径正确标注"代表图",不产生误导;
5. 全部单测通过,`linc run` 主流程无回归。

---

## 11. 已知限制(实现后需在 README 中说明)

1. `ssn.html` 的 Cytoscape.js 来自 CDN,离线时交互图不可用;TSV/GraphML 仍可用;
2. 报告生成会加载较大的 pickle,几百基因组时内存和时间都不低;
3. 交互视图只覆盖 `--family` 指定的家族,不做全图;完整 SSN 数据以 `ssn_nodes.tsv` / `ssn_edges.tsv` 为准;
4. `run_summary.json` 由报告命令生成,不是 `linc run` 自动生成;
5. 泛基因组分类采用 core/soft-core/shell/cloud/singleton 5 类互斥;single-copy core 是 core 的子集;`cloud+singleton` 仅是派生合并视图,不改变基础分类。
