# ============================================================
# Step 4 技术文档: 稀疏 SSN 构建
# ============================================================

## 概述

Step 4 基于增量阶段的结果构建稀疏序列相似性网络 (Sparse SSN),为后续补边和 MCL 聚类做准备。

---

## 输入

- **数据**:
  - `IncrementalResult`: 增量并入结果
    - `global_lib`: 全局代表库
    - `ssn_edges`: 代表级边集合
- **参数**: 无新增参数

---

## 流程图

```
┌─────────────────────────────────────────────────────────┐
│              Step 4: 构建稀疏 SSN                       │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  4.1 添加代表-代表边                                    │
│     ┌─────────────────────────────────────────────┐   │
│     │  全局簇结构:                                  │   │
│     │  cluster_A: reps={A1, A2}, members={...}   │   │
│     │  cluster_B: reps={B1}, members={...}        │   │
│     │  cluster_C: reps={C1, C2, C3}, members={...}│   │
│     │                                               │   │
│     │  边添加:                                      │   │
│     │  - cluster_A: A1-A2 (代表边)                │   │
│     │  - cluster_B: 无 (单代表)                   │   │
│     │  - cluster_C: C1-C2, C1-C3, C2-C3          │   │
│     └─────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  4.2 添加所有节点                                      │
│     ┌─────────────────────────────────────────────┐   │
│     │  节点类型:                                    │   │
│     │  - 代表节点 (大圆): 全局代表库的成员         │   │
│     │  - 普通节点 (小圆): 簇内其他成员             │   │
│     │                                               │   │
│     │  边类型:                                      │   │
│     │  - 代表-代表边 (实线): 已添加                 │   │
│     │  - 代表-普通边 (虚线): 尚未添加              │   │
│     │  - 普通-普通边 (虚线): 尚未添加              │   │
│     └─────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  4.3 计算连通分量                                      │
│     ┌─────────────────────────────────────────────┐   │
│     │  NetworkX 连通分量分析                        │   │
│     │  component_1: {A1, A2, B1, C1}             │   │
│     │  component_2: {C2, C3}                     │   │
│     │  component_3: {D1}                          │   │
│     └─────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────┐
│  输出: SparseSSN                                      │
│        ├── graph: NetworkX 图对象                      │
│        ├── nodes: 所有节点集合                          │
│        ├── representative_nodes: 代表节点集合            │
│        ├── edges: SSN 边集合                           │
│        └── connected_components: 连通分量列表           │
└─────────────────────────────────────────────────────────┘
```

---

## 稀疏 SSN 的特征

### 与完整 SSN 的区别

| 特性 | 稀疏 SSN | 完整 SSN |
|------|---------|---------|
| 边数量 | 少 | 多 |
| 边类型 | 仅代表-代表边 | 所有边 |
| 构建时机 | 增量阶段 | 补边阶段 |
| 用途 | 快速初步分析 | 最终结果 |

### 为什么先构建稀疏 SSN?

```
效率考量:
  - 增量阶段的核心任务是决定簇归属,只需代表级信息
  - 普通成员间的边不影响归属判定
  - 延迟到补边阶段再添加,避免 O(N²) 的边维护开销

完备性保证:
  - 稀疏 SSN 是中间产物
  - 通过补边阶段可恢复完整 SSN
  - 数学上等价于全量计算
```

---

## 图结构

```python
# 使用 NetworkX 构建无向加权图
import networkx as nx

graph = nx.Graph()

# 添加节点
graph.add_node("genome1|gene_0001")
graph.add_node("genome2|gene_0002")
# ...

# 添加边 (带权重)
graph.add_edge("A", "B", weight=0.85)
graph.add_edge("B", "C", weight=0.72)
# ...
```

---

## 连通分量分析

```python
# 计算所有连通分量
components = list(nx.connected_components(graph))

# 每个连通分量代表一个超家族
for i, comp in enumerate(components):
    print(f"Component {i}: {len(comp)} nodes")
```

---

## 输出数据结构

### SparseSSN

```python
@dataclass
class SparseSSN:
    graph: nx.Graph                       # NetworkX 图对象
    edges: Set[SSNEdge]                  # SSN 边集合
    nodes: Set[str]                      # 节点集合
    representative_nodes: Set[str]        # 代表节点集合
    n_edges: int                         # 边数量
    n_nodes: int                         # 节点数量
    connected_components: List[Set[str]] # 连通分量列表
```

---

## 统计信息

```python
def analyze_components(ssn):
    sizes = [len(c) for c in ssn.connected_components]
    return {
        "n_components": len(sizes),
        "largest_size": max(sizes),
        "singleton_count": sum(1 for s in sizes if s == 1),
        "mean_size": sum(sizes) / len(sizes)
    }
```

---

## 日志输出示例

```
============================================================
Sparse SSN Summary:
  Total nodes: 9
  Representative nodes: 9
  Total edges: 0
  Connected components: 9
  Largest component: 1 nodes
  Singleton components: 9
============================================================
```

---

## 下一步

```
Step 4 Output                              Step 5 Input
┌──────────────────────┐              ┌──────────────────────┐
│ SparseSSN           │ ──────────► │ EdgeCompleter       │
│  - graph            │              │  - sparse_ssn       │
│  - nodes            │              │  - incremental_result│
│  - components       │              │                      │
└──────────────────────┘              └──────────────────────┘
```
