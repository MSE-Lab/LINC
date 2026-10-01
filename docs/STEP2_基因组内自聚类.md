# ============================================================
# Step 2 技术文档: 基因组内自聚类
# ============================================================

## 概述

Step 2 对每个基因组独立执行全对全蛋白比对和局部同源簇聚类,提取局部代表序列。

---

## 输入

- **数据**: Step 1 输出的 `GenomeInfo` 对象列表
- **参数**:
  - Cutoff-I (默认 0.50)
  - Cutoff-II (默认 0.80)
  - max_per_cluster (默认 5)

---

## 流程图

```
┌─────────────────────────────────────────────────────────┐
│            Step 2: 基因组内自聚类 (并行)                  │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  2.1 读取基因组序列 (FASTA Parser)                       │
│     - 逐个读取蛋白序列                                    │
│     - 构建序列字典 {gene_id: FastaSequence}              │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  2.2 全对全比对 (All-vs-All)                            │
│     ┌─────────────────────────────────────────────────┐ │
│     │  genome01: gene_0001 ──比对── gene_0002        │ │
│     │        │                    │                   │ │
│     │        ├──比对──┐          │                   │ │
│     │        │       │          │                   │ │
│     │        ▼       ▼          ▼                   │ │
│     │    [hit1]  [hit2]    [hit3]                  │ │
│     │                                               │ │
│     │  gene_0003 ───────────────────────────► gene_0004 │
│     └─────────────────────────────────────────────────┘ │
│     - 并行处理多个基因组 (ThreadPoolExecutor)            │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  2.3 Union-Find 连通分量聚类                            │
│     ┌─────────────────────────────────────────────────┐ │
│     │  hits (Cutoff-I 过滤后)                         │ │
│     │  hit1: gene_0001 <-> gene_0002 (identity=0.95) │ │
│     │  hit2: gene_0001 <-> gene_0003 (identity=0.60) │ │
│     │  hit3: gene_0002 <-> gene_0003 (identity=0.58) │ │
│     │                                                   │ │
│     │  Union-Find 合并:                               │ │
│     │  gene_0001, gene_0002, gene_0003 → cluster_A   │ │
│     │  gene_0004 → cluster_B (孤儿基因)               │ │
│     └─────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  2.4 提取局部代表序列                                   │
│     ┌─────────────────────────────────────────────────┐ │
│     │  cluster_A (多成员):                            │ │
│     │    - 选取 1 个中心代表 (与簇内其他成员             │ │
│     │      平均相似性最高)                              │ │
│     │                                                   │ │
│     │  cluster_B (单成员/孤儿):                       │ │
│     │    - 自身为代表                                   │ │
│     │    - is_orphan = True                           │ │
│     └─────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  输出: Dict[genome_id, GenomeIntraResult]               │
│        ├── local_clusters: 9 个局部簇                  │
│        ├── local_representatives: 9 个代表              │
│        └── orphan_genes: 孤儿基因列表                   │
└─────────────────────────────────────────────────────────┘
```

---

## 关键设计决策

### 1. 必须使用全对全比对而非增量聚类

**为什么基因组内部不能增量?**

```
问题场景: 基因组内基因 a, b, c, d

增量方式 (错误):
  1. 添加 a -> 创建簇 {a}, 代表 {a}
  2. 添加 b -> 与 {a} 比对, 合并为 {a,b}
  3. 添加 c -> 与 {a,b} 比对, 发现 c 与 b 相似度更高
     但此时 a 已经是代表, b 不是...

全对全方式 (正确):
  - 一次性比较所有 (a,b), (a,c), (a,d), (b,c), (b,d), (c,d)
  - 获得完整的相似性矩阵
  - 根据全局最优做聚类决策
```

### 2. 孤儿基因不丢弃

**生物学依据: ORFans (孤儿基因)**

```
孤儿基因: 在单个基因组内没有同源伙伴,但可能在其他物种中存在同源基因

丢弃风险: 如果在自聚类阶段丢弃,这些基因将永久丧失在全局层面
          发现跨基因组同源关系的机会

LINC 处理: 保留孤儿基因作为单成员簇,其自身为代表
```

### 3. 代表选取策略

```python
def _select_representatives(member_ids, hit_dict):
    if len(members) == 1:
        return [member]  # 孤儿基因,自身为代表

    # 检查是否需要多代表
    max_identity = max(all_pairwise_identities)

    if max_identity >= Cutoff-II:
        # 簇内高度均质,单代表即可
        return [select_central_representative(members)]
    else:
        # 存在亚型分化,需要多代表
        return multi_representative_clustering(members)
```

---

## 局部簇类型

| 类型 | 描述 | 代表数 | 示例 |
|------|------|--------|------|
| 多成员簇 | 有同源伙伴 | 1-5 | 同源基因家族 |
| 单成员簇 | 孤儿基因 | 1 | ORFans |

---

## 输出数据结构

### GenomeIntraResult

```python
@dataclass
class GenomeIntraResult:
    genome_id: str
    local_clusters: List[LocalCluster]    # 局部簇列表
    local_representatives: List[str]      # 代表基因 ID 列表
    orphan_genes: List[str]              # 孤儿基因列表
    n_genes: int                        # 基因总数
    n_clusters: int                     # 簇数量
    n_representatives: int              # 代表数量
    processing_time: float               # 处理耗时
```

---

## 并行化策略

```python
# 所有基因组同时处理,无数据依赖
with ThreadPoolExecutor(max_workers=n_cpus) as executor:
    futures = {
        executor.submit(_process_single_genome, g): g
        for g in genomes
    }
```

- 基因组间**无数据依赖**,天然可并行
- 预期加速比:**接近线性**

---

## 压缩效果

```
输入: 20 个基因
输出: 9 个局部代表
压缩率: 9/20 = 45%

这意味着后续全局比对只需处理 45% 的序列!
```

---

## 日志输出示例

```
============================================================
Step 2: Intra-Genome Clustering
============================================================
Running 3 genomes in parallel with 3 workers
Completed genome01: 3 clusters, 3 representatives (0.01s)
Completed genome02: 3 clusters, 3 representatives (0.01s)
Completed genome03: 3 clusters, 3 representatives (0.01s)
============================================================
Intra-Genome Clustering Summary:
  Total genomes processed: 3
  Total genes: 20
  Total local clusters: 9
  Total local representatives: 9
  Overall compression ratio: 45.00%
============================================================
```

---

## 下一步

```
Step 2 Output                              Step 3 Input
┌──────────────────────┐              ┌──────────────────────┐
│ GenomeIntraResult    │ ──────────► │ IncrementalJoiner    │
│  - local_clusters    │              │  - genome_infos      │
│  - local_reps       │              │  - intra_results     │
│  - n_representatives│              │  - Order-D 排序      │
└──────────────────────┘              └──────────────────────┘
```
