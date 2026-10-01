# ============================================================
# LINC 算法设计文档 (工程开发参考)
# ============================================================

> 本文档基于 `STEP1`~`STEP6`、`MMseqs2_Aligner.md` 等现有技术文档汇总提炼而成,
> 面向**工程实现**,聚焦模块划分、数据结构、接口契约、算法规格与边界情况,
> 作为后续编码开发的设计依据。生物学原理与概念讲解请参考原始文档。

---

## 目录

1. [总体架构](#1-总体架构)
2. [核心数据结构定义](#2-核心数据结构定义)
3. [关键算法逻辑规格化](#3-关键算法逻辑规格化)
4. [配置参数汇总表](#4-配置参数汇总表)
5. [模块间接口契约](#5-模块间接口契约)
6. [边界情况与工程注意事项](#6-边界情况与工程注意事项)
7. [待实现/待确认事项清单](#7-待实现待确认事项清单)

---

## 1. 总体架构

### 1.1 六阶段流水线总览

```
FASTA 文件
   │
   ▼
Step 1  多基因组输入        genome_loader.py        → Dict[genome_id, GenomeInfo]
   │
   ▼
Step 2  基因组内自聚类(并行)  intra_genomic_clusterer  → Dict[genome_id, GenomeIntraResult]
   │
   ▼
Step 3  增量并入全局代表库(核心) incremental_joiner     → IncrementalResult
   │
   ▼
Step 4  稀疏 SSN 构建        sparse_ssn_builder      → SparseSSN
   │
   ▼
Step 5  最终补边             edge_completer          → CompleteSSN
   │
   ▼
Step 6  MCL 亚家族细分       mcl_refiner             → MCLResult
   │
   ▼
输出: OrthologFamily 列表 / SSNEdge 表 / GlobalRepresentativeLibrary / logs / run_summary.json
```

### 1.2 建议模块划分 (包结构)

```
src/
├── pipeline.py                     # LINC 主入口,串联六个 Step
├── config/
│   └── config_loader.py            # YAML 配置加载与校验
├── step1_input/
│   └── genome_loader.py            # GenomeLoader
├── step2_intra/
│   ├── intra_genomic_clusterer.py  # IntraGenomicClusterer (对外主类)
│   ├── union_find.py               # Union-Find 数据结构
│   └── representative_selector.py  # 代表选取 / 多代表决策树
├── step3_incremental/
│   ├── order_d_sorter.py           # Order-D 排序
│   ├── incremental_joiner.py       # IncrementalJoiner (核心)
│   └── decision_tree.py            # 双阈值决策树 (路径①②③)
├── step4_sparse_ssn/
│   └── sparse_ssn_builder.py       # SparseSSNBuilder
├── step5_edge_completion/
│   ├── edge_completer.py           # EdgeCompleter
│   └── pruning.py                  # 长度预过滤 + compared_pairs 去重
├── step6_mcl/
│   └── mcl_refiner.py              # MCLRefiner
├── aligner/
│   ├── base_aligner.py             # Aligner 抽象基类 (统一接口)
│   ├── mmseqs_aligner.py           # MMseqs2 实现
│   ├── diamond_aligner.py          # DIAMOND 实现 (可选)
│   └── mock_aligner.py             # Mock 实现 (测试用)
├── models/
│   └── data_structures.py          # 所有 @dataclass 定义 (见第 2 节)
└── utils/
    ├── fasta_parser.py             # FASTA 文件解析
    └── logger.py                   # 日志工具
```

> 说明: 目前项目仅有 `main.py`,尚未建立 `src/` 包结构。上述划分基于各 STEP 文档中出现的
> `pipeline.py (LINC.run)`、`step1_input/genome_loader.py` 等既有调用关系路径名推断整理,
> 实际落地时文件名/包名可按团队习惯微调,但建议保持"每个 Step 一个子包"的组织原则,
> 便于单元测试和职责隔离。

---

## 2. 核心数据结构定义

> 以下定义统一汇总自各 Step 文档,建议集中放在 `models/data_structures.py`,
> 避免跨模块重复定义导致类型不一致。

```python
from dataclasses import dataclass, field
from typing import List, Dict, Set, Optional
import networkx as nx

# ---------- Step 1 ----------
@dataclass
class GenomeInfo:
    genome_id: str           # 基因组唯一标识符
    file_path: str           # FASTA 文件绝对路径
    n_sequences: int         # 蛋白编码基因数量
    total_length: int        # 序列总长度
    avg_length: float        # 平均序列长度
    min_length: int          # 最短序列长度
    max_length: int          # 最长序列长度
    status: str              # pending / loaded / failed


# ---------- Step 2 ----------
@dataclass
class LocalCluster:
    cluster_id: str
    member_ids: List[str]
    representative_ids: List[str]
    is_orphan: bool = False


@dataclass
class GenomeIntraResult:
    genome_id: str
    local_clusters: List[LocalCluster]
    local_representatives: List[str]
    orphan_genes: List[str]
    n_genes: int
    n_clusters: int
    n_representatives: int
    processing_time: float


# ---------- 比对结果 (MMseqs2/DIAMOND/Mock 统一输出) ----------
@dataclass
class AlignmentResult:
    query_id: str
    subject_id: str
    identity: float          # 序列一致性 (0-1)
    coverage: float          # min(query_cover, subject_cover)
    evalue: float
    bitscore: float
    query_cover: float
    subject_cover: float
    alignment_length: int


# ---------- Step 3 ----------
@dataclass
class SSNEdge:
    source: str
    target: str
    identity: float
    coverage: float
    weight: float = 0.0      # 默认 identity * coverage,计算方式见 6.3


@dataclass
class GlobalCluster:
    cluster_id: str
    representative_ids: Set[str]
    member_ids: Set[str]


@dataclass
class GlobalRepresentativeLibrary:
    clusters: Dict[str, GlobalCluster]
    representative_to_cluster: Dict[str, str]   # 反查表: rep_id -> cluster_id


class ComparedPairsTracker:
    """已比对基因对集合,去重用。内部用 Set[Tuple[str, str]] (排序后的二元组) 实现。"""
    def __init__(self):
        self._pairs: Set[tuple] = set()

    def mark(self, a: str, b: str) -> None: ...
    def has_compared(self, a: str, b: str) -> bool: ...


@dataclass
class IncrementalResult:
    global_lib: GlobalRepresentativeLibrary
    ssn_edges: Set[SSNEdge]
    compared_pairs: ComparedPairsTracker
    processing_order: List[str]
    n_clusters_final: int
    n_representatives_final: int


# ---------- Step 4 ----------
@dataclass
class SparseSSN:
    graph: nx.Graph
    edges: Set[SSNEdge]
    nodes: Set[str]
    representative_nodes: Set[str]
    n_edges: int
    n_nodes: int
    connected_components: List[Set[str]]


# ---------- Step 5 ----------
@dataclass
class CompleteSSN:
    edges: Set[SSNEdge]
    n_edges: int
    n_nodes: int
    completion_ratio: float
    skipped_pairs: int
    processed_pairs: int


# ---------- Step 6 ----------
@dataclass
class OrthologFamily:
    family_id: str
    component_id: int
    member_genes: Set[str]
    size: int
    edge_density: float


@dataclass
class MCLResult:
    ortholog_families: List[OrthologFamily]
    n_families: int
    inflation_used: float
    processing_time: float
```

---

## 3. 关键算法逻辑规格化

### 3.1 Order-D 确定性排序 (Step 3 前置)

**排序键 (降序/升序见标注),用于对基因组列表排序,决定并入全局代表库的处理顺序:**

```python
def order_d_sort_key(genome: GenomeInfo, intra: GenomeIntraResult):
    return (
        -genome.n_sequences,              # 优先级1: 基因总数 ↓ (大基因组优先)
        -intra.n_representatives,         # 优先级2: 局部代表数 ↓ (高多样性优先)
        genome.genome_id,                 # 优先级3: 文件名/ID 字典序 ↑ (兜底唯一)
    )

def order_d_sort(genomes, intra_results):
    return sorted(genomes, key=lambda g: order_d_sort_key(g, intra_results[g.genome_id]))
```

> 要求: 结果必须**完全确定**,相同输入任意次运行排序结果一致,是全流程可复现性的基础。

### 3.2 Union-Find 聚类 (Step 2.3)

```python
class UnionFind:
    def __init__(self, elements: List[str]):
        self.parent = {e: e for e in elements}

    def find(self, x: str) -> str:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]  # 路径压缩
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb

def cluster_by_union_find(hits: List[AlignmentResult], cutoff_i: float) -> Dict[str, List[str]]:
    """
    hits: 全对全比对结果 (已含 identity/coverage)
    cutoff_i: 纳入阈值, 判定条件: identity >= cutoff_i AND min(qcov, scov) >= cutoff_i
    返回: {cluster_root: [member_ids]}
    """
    ...
```

### 3.3 代表选取 / 多代表决策 (Step 2.4)

```python
def select_representatives(member_ids: List[str], hit_matrix, cutoff_ii: float, max_per_cluster: int) -> List[str]:
    if len(member_ids) == 1:
        return member_ids  # 孤儿基因,自身为代表, is_orphan=True

    max_identity = max_pairwise_identity(member_ids, hit_matrix)

    if max_identity >= cutoff_ii:
        # 簇内高度均质 -> 单代表: 选取与簇内其他成员平均相似性最高者 (中心代表)
        return [select_central_representative(member_ids, hit_matrix)]
    else:
        # 存在亚型分化 -> 多代表聚类 (上限 max_per_cluster)
        return multi_representative_clustering(member_ids, hit_matrix, max_per_cluster)
```

> `multi_representative_clustering` 的具体子算法当前文档未展开实现细节,
> 需要开发时明确 (见第 7 节待确认事项)。建议思路: 对簇内再次做层次聚类或贪心最大覆盖,
> 直到所有成员都被某代表以 `cutoff_ii` 覆盖,或达到 `max_per_cluster` 上限。

### 3.4 双阈值决策树 (Step 3.3,核心)

```python
def decide(new_rep: str, hits: List[AlignmentResult], cutoff_i: float, cutoff_ii: float,
           global_lib: GlobalRepresentativeLibrary, max_per_cluster: int):
    # hits: new_rep 与全局代表库所有代表的比对结果 (一次比对,后续三处复用)
    passing_hits = [h for h in hits if h.identity >= cutoff_i and h.coverage >= cutoff_i]
    hit_clusters = {global_lib.representative_to_cluster[h.subject_id] for h in passing_hits}

    # ---- 路径①: 零簇命中 ----
    if len(hit_clusters) == 0:
        create_new_cluster(new_rep)
        return

    # ---- 路径②: 单簇命中 ----
    if len(hit_clusters) == 1:
        cluster_id = next(iter(hit_clusters))
        reaches_ii = any(h.identity >= cutoff_ii and h.coverage >= cutoff_ii for h in passing_hits)
        if reaches_ii:
            attach_to_cluster_no_new_rep(new_rep, cluster_id)      # 沿用已有代表
        else:
            add_as_representative(new_rep, cluster_id, max_per_cluster)  # 增补为代表
        return

    # ---- 路径③: 跨簇命中 (>=2 簇) ----
    reaches_ii_any = any(h.identity >= cutoff_ii and h.coverage >= cutoff_ii for h in passing_hits)
    merged_cluster_id = merge_clusters(hit_clusters)
    if reaches_ii_any:
        attach_to_cluster_no_new_rep(new_rep, merged_cluster_id)   # 只合并,不新增代表
    else:
        add_as_representative(new_rep, merged_cluster_id, max_per_cluster)  # 合并且新增代表
```

**关键原则 — 一次比对,全局复用:** `hits` 只计算一次,同时用于①簇归属判定、②合并触发、③代表裁决,
不允许为这三个决策分别发起二次比对 (避免额外开销与不确定性)。

### 3.5 补边阶段剪枝策略 (Step 5.2)

```python
def length_filter(len_a: int, len_b: int, threshold: float) -> bool:
    """返回 True 表示应跳过该比对对"""
    ratio = max(len_a, len_b) / min(len_a, len_b)
    return ratio > threshold

def should_compare(a: str, b: str, compared_pairs: ComparedPairsTracker) -> bool:
    return not compared_pairs.has_compared(a, b)

def generate_pairs_to_complete(cluster: GlobalCluster, compared_pairs, length_lookup, threshold):
    pairs = []
    members = list(cluster.member_ids)
    for i in range(len(members)):
        for j in range(i + 1, len(members)):
            a, b = members[i], members[j]
            if not should_compare(a, b, compared_pairs):
                continue
            if length_filter(length_lookup[a], length_lookup[b], threshold):
                continue
            pairs.append((a, b))
    return pairs
```

### 3.6 MCL Inflation 自适应选择 (Step 6.2)

```python
def determine_inflation(graph: nx.Graph, default=2.0, min_v=1.5, max_v=4.0, step=0.5) -> float:
    weights = [d["weight"] for _, _, d in graph.edges(data=True)]
    if not weights:
        return default
    variance = statistics_variance(weights)   # 如 numpy.var
    if variance > 0.1:
        return min(max_v, default + step)    # 权重分散 -> 细分
    else:
        return max(min_v, default - step)    # 权重集中 -> 合并
```

### 3.7 边权重计算 (Step 6.1)

```python
def calculate_edge_weight(edge: AlignmentResult, mode: str = "identity_x_coverage") -> float:
    if mode == "identity_only":
        return edge.identity
    if mode == "coverage_only":
        return edge.coverage
    return edge.identity * edge.coverage  # 默认 & identity_x_coverage
```

---

## 4. 配置参数汇总表

| 参数 | 默认值 | 所属阶段 | 说明 |
|------|--------|---------|------|
| `cutoff_i` | 0.50 | Step 2/3/5 | 纳入阈值: identity & coverage 双向判定是否成簇/建边 |
| `cutoff_ii` | 0.80 | Step 2/3 | 高置信阈值: 判定是否新增代表 |
| `max_per_cluster` | 5 | Step 2/3 | 单簇最大代表数上限 |
| `length_ratio_threshold` | 10.0 (STEP2文档示例) | Step 5 | 长度预过滤比值上限,超过则跳过比对 |
| `dedup_enabled` | true | Step 5 | 是否启用 `compared_pairs` 去重 |
| `aligner.backend` | `"mmseqs"` | 全局 | 比对引擎: mmseqs / diamond / mock |
| `aligner.mmseqs.sensitivity` | 4.5 | MMseqs2 | 1.0(快)~7.0(灵敏),平衡模式推荐 4.5 |
| `aligner.mmseqs.evalue` | 1e-5 | MMseqs2 | E-value 阈值 |
| `aligner.mmseqs.max_accept` | 500 | MMseqs2 | 最大保留比对数 |
| `aligner.mmseqs.min_alignment_score` | 20 | MMseqs2 | 最小比对分数 |
| `parallel.n_threads` | 0 (全部核心) | 全局 | 并行线程数 |
| `parallel.inter_genome_parallel` | true | Step 2 | 基因组间并行开关 |
| `parallel.finalize_parallel` | true | Step 5 | 补边并行开关 |
| `mcl.inflation` | 2.0 | Step 6 | MCL 膨胀参数,默认值;自适应范围建议 [1.5, 4.0] |
| `mcl.auto_inflation` | true/false (可配) | Step 6 | 是否根据边权重方差自适应选择 inflation |
| `mcl.edge_weight_mode` | `"identity_x_coverage"` | Step 6 | 边权重计算模式 |

> 建议统一放入 `config.yaml`,并在 `config_loader.py` 中做 schema 校验
> (类型、取值范围、必填项),避免非法配置在流程中段才暴露问题。

---

## 5. 模块间接口契约

| 阶段 | 输入类型 | 输出类型 | 主要入口类/函数 (推断命名) |
|------|---------|---------|---------------------------|
| Step 1 | `config.yaml` | `Dict[str, GenomeInfo]` | `GenomeLoader.load()` |
| Step 2 | `Dict[str, GenomeInfo]` | `Dict[str, GenomeIntraResult]` | `IntraGenomicClusterer.run()` |
| Step 3 | `GenomeInfo` + `GenomeIntraResult` (全部基因组) | `IncrementalResult` | `IncrementalJoiner.run()` |
| Step 4 | `IncrementalResult` | `SparseSSN` | `SparseSSNBuilder.build()` |
| Step 5 | `SparseSSN` + `IncrementalResult` | `CompleteSSN` | `EdgeCompleter.complete()` |
| Step 6 | `CompleteSSN` | `MCLResult` | `MCLRefiner.refine()` |

**调用顺序约束:** 严格串行 Step1→Step2→Step3,Step2 内部基因组间可并行;
Step3 内部必须按 Order-D 排序后**串行**处理基因组 (存在全局状态依赖,不可并行);
Step5 内部各簇可并行 (无共享可变状态冲突,`compared_pairs` 建议按簇分区或加锁)。

---

## 6. 边界情况与工程注意事项

### 6.1 数据/输入层面
- **空基因组文件** (0 条序列): Step 1 应标记 `status="failed"` 并在摘要中告警,而非直接崩溃。
- **单基因组运行** (只有 1 个基因组): Step 3 的"全局代表库"退化为该基因组自身的局部代表,
  Order-D 排序结果恒定,应有专门测试用例覆盖。
- **重复 genome_id**: 需要在 Step 1 做唯一性校验,冲突时报错而非静默覆盖。

### 6.2 算法逻辑层面
- **孤儿基因 (ORFans) 不可丢弃**: Step 2 必须保留单成员簇,`is_orphan=True`,
  否则会在全局层面永久丧失跨基因组同源发现机会。
- **空簇 / 合并后空簇**: Step 3 路径③合并操作后,原簇对象应被标记失效或从
  `global_lib.clusters` 中移除,避免 Step 4/5 遍历到悬空引用。
- **`compared_pairs` 一致性**: Step 3 增量阶段产生的比对记录必须被 Step 5 完整复用,
  任何"重复比对同一对"都违反数学等价性承诺,需要单元测试显式验证
  (对照"全量两两比对"结果做一致性断言)。
- **MCL 单节点连通分量**: `size=1` 的连通分量不应调用 MCL (无边可聚),
  直接作为单基因家族输出。

### 6.3 并行化层面
- Step 2 各基因组间**无数据依赖**,适合 `ThreadPoolExecutor`,但注意 I/O (FASTA 读取)
  和 CPU 密集型比对混合时的线程数配置。
- Step 3 **不可并行** (Order-D 顺序依赖 + 全局代表库共享可变状态),误用并行会破坏
  "一次比对全局复用"和决策树的正确性。
- Step 5 并行时,`compared_pairs` 若为全局共享集合,需要加锁或改为"先并行计算结果,
  最后统一合并写入"的模式,避免竞态。

### 6.4 比对引擎层面
- `base_aligner.py` 应定义统一抽象接口 (如 `run_all_vs_all()`、`run_query()`),
  使 MMseqs2 / DIAMOND / Mock 三种实现可互换,便于测试时用 Mock 替代真实比对。
- MMseqs2 未安装时应给出清晰报错并提示安装方式 (`conda install -c bioconda mmseqs2`),
  而不是抛出难以理解的子进程异常。

---

## 7. 待实现/待确认事项清单

以下内容在现有 `docs/` 文档中提及但未展开具体算法细节,建议开发前与相关方确认:

1. **`multi_representative_clustering` 的具体算法**: Step 2.4 提到"存在亚型分化,需要多代表聚类",
   但未说明具体子聚类算法 (层次聚类?贪心覆盖?),需要明确后再实现。
2. **`select_central_representative` 的打分公式**: 文档描述为"与簇内其他成员平均相似性最高者",
   需要明确是用 identity 均值、coverage 均值,还是两者加权。
3. **`edge_density` (OrthologFamily) 的计算方式**: Step 6 输出数据结构中出现该字段,
   但未给出计算公式,推测为 `实际边数 / C(n,2)`,需要确认。
4. **`合并簇` (merge_clusters) 的代表集合处理**: 路径③合并两个簇后,原两组代表是否全部保留,
   还是需要去重/重新评估,文档未明确。
5. **多线程下 `compared_pairs` 的并发安全实现**: Step 5 文档提到"GIL 保护"下共享集合写入,
   但 Python 的 `set.add()` 是否严格原子、是否需要显式锁,建议实现时用锁或线程安全结构兜底。
6. **MCL 库选择**: `markov-clustering` (Python 库) vs 系统级 `mcl` 命令行工具,
   两者是否都要支持,或者以哪个为主,需要在 `step6_mcl/mcl_refiner.py` 设计时确定
   (当前 conda 环境两者都已安装)。
