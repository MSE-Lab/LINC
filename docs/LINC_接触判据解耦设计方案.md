# LINC 接触判据 / 纳入判据解耦设计方案

> 状态:设计稿(待评审)
> 日期:2026-10-03
> 关联:`docs/LINC_vs_OrthoFinder_gap根因分析报告.md`(gap 与顺序偏差同源:SSN 不完整度)
> 目标:在现有增量框架内,把"接触"(要不要搜成员)与"纳入"(够不够强成边)解耦,
> 用独立的宽松接触判据提高跨基因组比对完备性,进而降低顺序偏差、改善低相似 recall。

---

## 1. 背景与目标

### 1.1 现状问题

LINC 的跨基因组比对完备性由"接触"决定,而"接触"当前是 join 阶段"纳入判据放行合并"的
**副产品**——它们被同一个严格阈值 `passes_inclusion_cutoff`(identity ≥ cutoff_i AND
coverage ≥ cutoff_c)绑死:

- 一个簇的代表若对某基因组只有"低 identity 但 evalue 显著"的弱命中,则:
  1. 弱命中过不了纳入判据 → **不合并**;
  2. 不合并 → 该簇**不"接触"该基因组**;
  3. 不接触 → full 补边阶段**不去搜该簇的成员**;
  4. 于是"成员其实比代表强、本可成边"的情况被永久漏掉。

这直接导致(见根因报告):Streptomyces 跨物种 F1 低、RBB 覆盖率仅 3.3%、随机分批 append
的顺序偏差——**三者的共同病根都是"增量搜索不完整"**。

### 1.2 目标

在不放弃增量、不做全对全的前提下,引入一个**独立的、宽松的"接触判据"**(evalue),让
"哪些簇的成员被搜"不再被"是否够强合并"卡死;同时保持"纳入判据"(identity+coverage)严格,
确保进 SSN 的边仍是高置信的。

---

## 2. 两个判据的职责划分

| 判据 | 判据函数 | 松/严 | 回答的问题 | 用在哪 |
|---|---|---|---|---|
| **接触判据** | `evalue ≤ touch_evalue` | 宽松 | 这个簇要不要被搜成员? | join(打标记) |
| **纳入判据** | `passes_inclusion_cutoff`(identity+coverage) | 严格 | 这条 hit 够不够格合并/成边? | join(合并)、full(加边) |

**关键**:接触判据只"打标记",不参与合并;cutoff_ii 的"成员 vs 新代表"逻辑完全不变(它只
吃强 hit)。

---

## 3. 改动总览

| # | 文件 | 改动 | 性质 |
|---|---|---|---|
| 1 | `linc/core/entities.py` | `IncrementalResult` 加 `touching_map` | 数据模型 |
| 2 | `linc/pipeline.py` | `PipelineConfig` 加 `touch_evalue`;透传两个 stage | 配置 |
| 3 | `linc/cli.py` | 加 `--touch-evalue` 选项 | 配置 |
| 4 | `linc/engine/global_join_stage.py` | join 加"接触标记"(与合并解耦) | 核心 |
| 5 | `linc/engine/ssn_builder_stage.py` | full 扩搜索范围 + 去掉"同簇过滤" | 核心 |
| 6 | `linc/pipeline.py` / `file_repo.py` | 一致性 + SCHEMA_VERSION | 持久化 |

---

## 4. 详细改动

### 4.1 数据模型:`IncrementalResult` 加 `touching_map`

文件:`linc/core/entities.py`(约 373 行)

```python
@dataclass
class IncrementalResult:
    global_lib: GlobalRepresentativeLibrary = field(default_factory=GlobalRepresentativeLibrary)
    ssn_edges: Set["SSNEdge"] = field(default_factory=set)
    processing_order: List[str] = field(default_factory=list)
    n_clusters_final: int = 0
    n_representatives_final: int = 0
    cluster_counter: int = 0
    # 新增(2026-10-03): 接触关系映射。key = target 基因组 id,
    # value = "接触了该基因组"的全局簇 id 集合。GlobalJoinStage 按宽松接触判据
    # 记录, SSNBuilderStage full 补边据此决定"哪些簇的成员需要被搜",
    # 与"合并(纳入判据)"解耦。
    touching_map: Dict[str, Set[str]] = field(default_factory=dict)
```

> 确认 `Dict` 已 import(entities.py 通常已有 `Set`;缺 `Dict` 则补 `from typing import Dict`)。

---

### 4.2 配置:`PipelineConfig` 加 `touch_evalue`

文件:`linc/pipeline.py`

```python
    dump_raw_hits: Optional[str] = None
    # 新增: 接触判据的 evalue 阈值。None = 最宽松(所有 search 返回的 hit 都算
    # 接触, 即 evalue <= 搜索阈值); 设为更小值(如 1e-10)可收紧"接触"、
    # 减少 full 补边的搜索量(性能旋钮)。语义参数, 影响结果, 纳入一致性校验。
    touch_evalue: Optional[float] = None
```

---

### 4.3 CLI:加 `--touch-evalue`

文件:`linc/cli.py`(在 `dump_raw_hits` 附近)

```python
    touch_evalue: Optional[float] = typer.Option(
        None, "--touch-evalue",
        help="E-value threshold for the 'touching' criterion used by --edge-completion-mode full: "
             "a global cluster whose representative hits a genome below this e-value is marked as "
             "'touching' that genome, so its members get cross-compared even if the hit is too weak "
             "to merge. None (default) treats every search hit as touching (maximally relaxed); "
             "set a smaller value (e.g. 1e-10) to tighten it and reduce full-stage cost.",
    ),
```

透传进 `PipelineConfig`:

```python
            config=PipelineConfig(
                ...
                dump_raw_hits=dump_raw_hits,
                touch_evalue=touch_evalue,
            ),
```

`Pipeline._run_join` / `_run_ssn` 把 `config.touch_evalue` 传给 `GlobalJoinStage` /
`SSNBuilderStage`。

---

### 4.4 核心改动 1:`GlobalJoinStage`

#### `__init__`

```python
    def __init__(self, ..., dump_raw_hits=None, touch_evalue=None):
        ...
        self._dump_raw_hits = dump_raw_hits
        self._touch_evalue = touch_evalue
        self._touching_map: Dict[str, Set[str]] = defaultdict(set)
```

#### `run()` 初始化并写回结果

```python
    def run(self, ...):
        ...
        self._touching_map = defaultdict(set)

        for idx, genome_id in enumerate(order, start=1):
            ...
            self._process_genome(genome, intra, global_lib, ssn_edges)
            ...

        result = IncrementalResult(
            ...
            cluster_counter=self._cluster_counter,
            touching_map=dict(self._touching_map),   # 新增
        )
```

#### `_process_genome` 加"接触标记"(与合并解耦)

现在(约 196-199 行):

```python
        hits_by_gene: Dict[str, List[AlignmentResult]] = defaultdict(list)
        for hit in hits:
            if passes_inclusion_cutoff(hit, self._cutoff_i, self._cutoff_c, self._cov_mode):
                hits_by_gene[hit.subject_id].append(hit)
```

改成(在纳入判据**之前**加独立接触标记循环):

```python
        # 【新增】接触判据: 只打标记, 不参与 merge。默认(None)所有返回 hit 都算接触。
        for hit in hits:
            if self._touch_evalue is not None and hit.evalue > self._touch_evalue:
                continue
            cluster_id = global_lib.gene_to_cluster.get(hit.query_id)
            if cluster_id is not None:
                self._touching_map[genome.genome_id].add(cluster_id)

        # 【现有, 不变】纳入判据: 只有强 hit 参与合并(进而决定 cutoff_ii 成员/代表)
        hits_by_gene: Dict[str, List[AlignmentResult]] = defaultdict(list)
        for hit in hits:
            if passes_inclusion_cutoff(hit, self._cutoff_i, self._cutoff_c, self._cov_mode):
                hits_by_gene[hit.subject_id].append(hit)
```

> **cutoff_ii 的成员/代表逻辑一行不改**——`_decide_and_apply` 只吃 `hits_by_gene`(强 hit),
> 接触标记完全在其之外、与之正交。

> **方向说明**:join 的 search 是 `query=全局代表, subject=本基因组基因`,故
> `hit.query_id` 是代表 → `gene_to_cluster[hit.query_id]` 是该簇;语义是"该簇接触了本基因组"。

---

### 4.5 核心改动 2:`SSNBuilderStage`

#### 第一遍:扩搜索范围

现有第一遍只处理"横跨 ≥2 基因组的簇"。在它**之后**追加:

```python
        # 【新增】接触判据放宽后的搜索范围: 对每个 target 基因组 B, 把"接触了 B
        # (但未必在 B 里有成员)"的簇的成员(来自其它基因组)也加进 query。
        touching_map = getattr(incremental_result, "touching_map", None) or {}
        for genome_b, cluster_ids in touching_map.items():
            for cluster_id in cluster_ids:
                cluster = global_lib.get_cluster(cluster_id)
                if cluster is None or not cluster.active:
                    continue
                for gene_id in cluster.member_ids:
                    gene = genes_by_id.get(gene_id)
                    if gene is None or gene.genome_id == genome_b:
                        continue
                    query_genes_by_target_genome.setdefault(genome_b, set()).add(gene_id)
```

> `query_genes_by_target_genome` 是 `Set`,与既有 merge 逻辑重复添加会自动去重。

> 需确认 `_complete_cross_genome_member_edges` 能拿到 `incremental_result`;若没有,
> 从 `run()` 传入 `incremental_result.touching_map`。

#### 第二遍:去掉"同簇过滤"

现在(约 379-390 行):

```python
            for hit in hits:
                if gene_cluster_map.get(hit.query_id) != gene_cluster_map.get(hit.subject_id):
                    continue                       # 同簇过滤
                if not passes_inclusion_cutoff(...):
                    continue
                edges.add(SSNEdge(...))
```

改成:

```python
            for hit in hits:
                # 【改动】去掉"同簇过滤": "接触了 B 但没并入 B"的簇, 其成员命中的
                # B 基因必然在别的簇, 同簇过滤会把真边全拒掉。只保留纳入判据兜底。
                if not passes_inclusion_cutoff(
                    hit, self._cutoff_i, self._cutoff_c, self._cov_mode
                ):
                    continue
                edges.add(SSNEdge(
                    source=hit.query_id, target=hit.subject_id,
                    identity=hit.identity, coverage=hit.coverage,
                    bitscore=hit.bitscore,
                    source_length=hit.query_length, target_length=hit.subject_length,
                ))
```

---

## 5. 持久化与一致性

1. **`touching_map` 随 `join_result.pkl` 自动持久化**:`IncrementalResult` 是 dataclass,
   走既有 gzip pickle,无需额外代码。
2. **`touch_evalue` 是语义参数**:影响结果,**不加**进 `non_semantic_fields`;默认走
   `_check_config_consistency` 校验(改它触发重算/报错)。
3. **`SCHEMA_VERSION` 递增**:`IncrementalResult` 新增字段,旧 state 反序列化会缺字段,
   `linc/adapters/storage/file_repo.py` 的 `SCHEMA_VERSION` 从 2 递增到 3。

---

## 6. 边界情况与风险

| 情况 | 处理 |
|---|---|
| 接触标记与 merge 重叠 | `query_genes_by_target_genome` 是 `Set`,自动去重 |
| 单基因组簇(只接触、不跨基因组) | 第一遍 skip;新逻辑补上其成员 → target |
| 去掉"同簇过滤"后边变多 | 预期内(补回跨簇真边);纳入判据兜底防弱边 |
| append 模式 | `touching_map` 需随历史 `join_result` 复用,`run()` 加 `initial_touching_map` 参数(本方案先做全量,append 留后续) |
| 性能 | 接触越宽松 → full 搜越多 → 越慢;`--touch-evalue` 收紧可降成本 |

---

## 7. 验证计划

1. **单测**:
   - join:弱 hit(evalue 达标、identity 不达标)→ 只打接触标记、不合并、不影响 cutoff_ii;
   - full:接触但不同簇的 hit → 能加边(同簇过滤已去掉);
   - 边界:`touch_evalue=None` 时所有 hit 都标记接触。

2. **离线收益验证(复用 `--dump-raw-hits` 已落盘的原始 hit)**:
   - 模拟"接触放宽后多搜出哪些边、RBB 覆盖率从 3.3% 升到多少、家族 F1 抬多少";
   - 与"全对全 RBB = F1 1.0"上界对比,量化旋钮的真实收益。

3. **回归**:
   - `touch_evalue=None` + full 模式,在 E. coli 上仍复现等价性 F1 0.995;
   - 随机分批 append 的顺序偏差随接触放宽而下降(量化)。

---

## 8. 落地顺序

1. 先做 4.1~4.5(数据模型 + 配置 + join 标记 + full 扩范围/去过滤)→ 最小可跑闭环;
2. 跑通后立即做"离线收益验证"(7.2),确认收益再继续;
3. 再补第 5 节(持久化/一致性)、append 接续、单测、回归。

---

## 附:与既有结论的关系

- 本方案**不引入 RBB、不做全对全、不改 cutoff_ii**,保持增量范式;
- 它是把根因报告里"SSN 不完整度"这个病根,用一个**可调旋钮**在"速度 ↔ 完备性 ↔ 顺序偏差"
  三者之间做权衡;
- `touch_evalue=None`(默认)是"最大松弛"端,等价于"所有 search 命中都触发接触";
  `--touch-evalue 1e-10` 是"收紧"端,用于降成本。
