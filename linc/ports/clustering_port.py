"""ClusteringPort: 聚类算法的抽象接口 (Protocol)。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.2 节。
统一 Union-Find (Step2 基因组内聚类) 与 MCL (Step6 亚家族细分) 的调用接口,
具体实现见 adapters/clustering/union_find_adapter.py 与 adapters/clustering/mcl_adapter.py。
"""
from typing import Protocol, runtime_checkable

import networkx as nx


@runtime_checkable
class ClusteringPort(Protocol):
    """聚类算法抽象接口。"""

    def cluster(self, graph: nx.Graph, **params) -> list[set[str]]:
        """输入图, 输出聚类结果 (每个 set 为一个家族的节点集合)。"""
        ...
