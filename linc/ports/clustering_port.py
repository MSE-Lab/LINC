"""ClusteringPort: 聚类算法的抽象接口 (Protocol)。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.2 节。
统一 Union-Find (Step2 基因组内聚类) 与 MCL (Step6 亚家族细分) 的调用接口,
具体实现见 adapters/clustering/union_find_adapter.py 与 adapters/clustering/mcl_adapter.py。
"""
from typing import Optional, Protocol, runtime_checkable

import networkx as nx


@runtime_checkable
class ClusteringPort(Protocol):
    """聚类算法抽象接口。"""

    def cluster(self, graph: nx.Graph, **params) -> list[set[str]]:
        """输入图, 输出聚类结果 (每个 set 为一个家族的节点集合)。"""
        ...

    def cleanup(self, index_id: Optional[str] = None) -> None:
        """删除 `cluster(index_id=...)` 产生的中间文件(数据管理优化第3项,
        2026-10-03 与用户确认): 如 `MCLClustering` 落盘的 `.abc` 输入/mcl
        原生输出, 其结果已被调用方解析并转存进 `OrthologFamily` 列表, 中间
        文件本身不再被后续步骤读取, 调用方应在消费完 `cluster()` 的返回值
        后立即调用本方法清理。不持久化中间文件的实现(如
        `UnionFindClustering`)可实现为空操作。
        """
        ...
