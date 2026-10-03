"""UnionFindClustering: ClusteringPort 的并查集(Union-Find)实现。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.3 节, 用于 Step2 基因组内自聚类:
输入一张已按相似度阈值过滤好的图(节点=基因, 边=满足阈值的相似关系),
输出若干个互不相交的节点集合(每个集合即一个簇)。

算法本身是标准的带路径压缩 + 按秩合并的并查集, 复杂度接近 O(n * α(n))。
注意: 阈值过滤(Cutoff-I 判定)不是本 adapter 的职责, 应在构造 graph 时
(或调用方 Stage 中)已经完成, 本 adapter 只做"连通性"计算。
"""
from typing import Dict, List, Optional, Set

import networkx as nx


class _DisjointSet:
    """内部使用的并查集数据结构(路径压缩 + 按秩合并)。"""

    def __init__(self, nodes) -> None:
        self._parent: Dict[str, str] = {n: n for n in nodes}
        self._rank: Dict[str, int] = {n: 0 for n in nodes}

    def find(self, x: str) -> str:
        root = x
        while self._parent[root] != root:
            root = self._parent[root]
        # 路径压缩
        while self._parent[x] != root:
            self._parent[x], x = root, self._parent[x]
        return root

    def union(self, x: str, y: str) -> None:
        root_x, root_y = self.find(x), self.find(y)
        if root_x == root_y:
            return
        # 按秩合并
        if self._rank[root_x] < self._rank[root_y]:
            root_x, root_y = root_y, root_x
        self._parent[root_y] = root_x
        if self._rank[root_x] == self._rank[root_y]:
            self._rank[root_x] += 1

    def groups(self) -> List[Set[str]]:
        buckets: Dict[str, Set[str]] = {}
        for node in self._parent:
            root = self.find(node)
            buckets.setdefault(root, set()).add(node)
        return list(buckets.values())


class UnionFindClustering:
    """ClusteringPort 的并查集实现。"""

    def cluster(self, graph: nx.Graph, **params) -> list[set[str]]:
        """以图的边为合并依据, 返回每个连通分量的节点集合。

        孤立节点(无任何边)也会作为单独的簇返回, 对应"孤儿基因"场景。
        """
        dsu = _DisjointSet(graph.nodes)
        for u, v in graph.edges:
            dsu.union(u, v)
        return dsu.groups()

    def cleanup(self, index_id: Optional[str] = None) -> None:
        """UnionFind 实现不落盘任何中间文件, 本方法为满足 ClusteringPort
        接口而提供, 无需做任何事(2026-10-03 与用户确认)。
        """
        return None
