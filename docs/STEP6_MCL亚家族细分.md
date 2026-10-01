# ============================================================
# Step 6 技术文档: MCL 亚家族细分
# ============================================================

## 概述

Step 6 在完整 SSN 连通分量上运行 MCL (Markov Cluster Algorithm),将聚类粒度从"超家族级"推进到"直系同源家族级"。

---

## 输入

- **数据**:
  - `CompleteSSN`: 完整 SSN
- **参数**:
  - MCL Inflation 参数 (默认 2.0)
  - edge_weight_mode (边权重计算模式)

---

## 流程图

```
┌─────────────────────────────────────────────────────────┐
│              Step 6: MCL 亚家族细分                      │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  6.1 构建 NetworkX 图                                   │
│     ┌─────────────────────────────────────────────┐   │
│     │  将 CompleteSSN 转换为 NetworkX 图             │   │
│     │                                               │   │
│     │  边权重计算:                                  │   │
│     │  weight = identity × coverage               │   │
│     │  (或仅 identity / 仅 coverage)              │   │
│     └─────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  6.2 确定 Inflation 参数                               │
│     ┌─────────────────────────────────────────────┐   │
│     │  auto_inflation = True:                    │   │
│     │  ┌─────────────────────────────────────┐  │   │
│     │  │ 根据边权重分布自适应选择:             │  │   │
│     │  │ - 权重分散 → 较高 inflation          │  │   │
│     │  │ - 权重集中 → 较低 inflation          │  │   │
│     │  └─────────────────────────────────────┘  │   │
│     │                                               │   │
│     │  auto_inflation = False:                    │   │
│     │  使用配置值 (默认 2.0)                    │   │
│     └─────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  6.3 对每个连通分量运行 MCL                           │
│     ┌─────────────────────────────────────────────┐   │
│     │  连通分量:                                   │   │
│     │  {A, B, C, D, E, F, G, H}                 │   │
│     │                                               │   │
│     │  MCL 聚类结果:                              │   │
│     │  - family_1: {A, B, C}                    │   │
│     │  - family_2: {D, E}                      │   │
│     │  - family_3: {F, G}                      │   │
│     │  - family_4: {H}                        │   │
│     └─────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  输出: MCLResult                                      │
│        ├── ortholog_families: 直系同源家族列表         │
│        ├── n_families: 家族总数                       │
│        └── inflation_used: 使用的 Inflation 参数        │
└─────────────────────────────────────────────────────────┘
```

---

## MCL 算法简介

### Markov Cluster Algorithm

MCL 通过模拟网络中的**随机游走过程**来发现密集连接的子区域。

```
核心思想:
  - 随机游走在密集区域停留时间长
  - 在稀疏区域容易跳出
  - 多次迭代后,随机游走趋向收敛到稠密子图

数学过程:
  1. 转移矩阵: 表示从节点到邻居的概率
  2. 膨胀 (Inflation): 强化强连接,弱化弱连接
  3. 扩张 (Expansion): 增加新连接
  4. 重复直到收敛
```

### Inflation 参数的影响

```
Inflation 参数控制聚类的精细程度:

低 Inflation (1.5):
  ├─ 聚类较粗
  ├─ 更多基因被归入同一家族
  └─ 适合保守的家族划分

高 Inflation (4.0):
  ├─ 聚类较细
  ├─ 家族被细分更多
  └─ 适合精细的功能分型

LINC 默认: 2.0
```

---

## 聚类粒度演进

```
┌─────────────────────────────────────────────────────────┐
│  粗聚类: SSN 连通分量 (超家族级)                        │
│  - 基于序列相似性网络连通性                            │
│  - 可能包含功能分化的旁系同源                          │
│  - 例如: 所有血红蛋白相关基因被归在一起                │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼ MCL 细分
┌─────────────────────────────────────────────────────────┐
│  细聚类: MCL 亚簇 (直系同源家族级)                     │
│  - 基于边权重的密集子区域                              │
│  - 反映功能特异性分化                                  │
│  - 例如: α-球蛋白 vs β-球蛋白 分离                    │
└─────────────────────────────────────────────────────────┘
```

---

## 边权重计算

```python
def calculate_edge_weight(edge, mode="identity_x_coverage"):
    if mode == "identity_x_coverage":
        return edge.identity * edge.coverage
    elif mode == "identity_only":
        return edge.identity
    elif mode == "coverage_only":
        return edge.coverage
    else:
        return edge.identity * edge.coverage
```

---

## 自适应 Inflation

```python
def determine_inflation(graph):
    """
    根据边权重分布自动选择 Inflation 参数
    """
    weights = [d["weight"] for _, _, d in graph.edges(data=True)]

    if not weights:
        return DEFAULT_INFLATION

    variance = np.var(weights)

    if variance > 0.1:
        # 权重分散 → 较高 inflation 细分
        return min(MAX_INFLATION, INFLATION + 0.5)
    else:
        # 权重集中 → 较低 inflation 合并
        return max(MIN_INFLATION, INFLATION - 0.5)
```

---

## 输出数据结构

### OrthologFamily

```python
@dataclass
class OrthologFamily:
    family_id: str           # 家族 ID
    component_id: int        # 所属连通分量
    member_genes: Set[str]   # 成员基因集合
    size: int                # 家族大小
    edge_density: float      # 内部边密度
```

### MCLResult

```python
@dataclass
class MCLResult:
    ortholog_families: List[OrthologFamily]  # 家族列表
    n_families: int                         # 家族总数
    inflation_used: float                    # Inflation 参数
    processing_time: float                    # 处理耗时
```

---

## 生物学意义

### 从超家族到直系同源家族

```
超家族 (连通分量):
  "所有与血红蛋白相关的基因"
  - α-球蛋白基因
  - β-球蛋白基因
  - 肌红蛋白基因
  - 细胞球蛋白基因
  ...

直系同源家族 (MCL 细分):
  "物种分化后保留的直系同源基因"
  - α-球蛋白家族 (人类、小鼠、鸡...)
  - β-球蛋白家族 (人类、小鼠、鸡...)
  - 肌红蛋白家族 (脊椎动物共同祖先)
```

### 为什么需要这一步?

```
问题: 传统的直系同源推断工具只输出粗粒度结果

LINC 的改进:
  1. Step 1-3: 构建全局同源关系网络
  2. Step 4-5: 补全边集合
  3. Step 6: MCL 细分 → 直系同源家族

优势:
  - 既能处理大规模数据 (增量架构)
  - 又能精细区分功能分化 (MCL 细分)
  - 不需要额外的外部注释或物种树
```

---

## 日志输出示例

```
============================================================
Step 6: MCL Subfamily Refinement
============================================================
MCL Subfamily Refinement Summary:
  Inflation parameter: 2.00
  Total ortholog families: 45
  Processing time: 0.12s
============================================================
```

---

## 最终输出汇总

```
LINC 完整流程输出:
  1. 直系同源基因家族 (OrthologFamily 列表)
  2. SSN 边表 (SSNEdge 列表)
  3. 全局代表库 (GlobalRepresentativeLibrary)
  4. 运行日志 (logs/linc_*.log)
  5. 运行摘要 (output/run_summary.json)
```
