"""ClusteringPort 结构化子类型测试。"""
import networkx as nx

from linc.ports.clustering_port import ClusteringPort


class FakeConnectedComponentsClustering:
    """最小化假实现: 直接用连通分量作为聚类结果。"""

    def cluster(self, graph: nx.Graph, **params) -> list[set[str]]:
        return [set(component) for component in nx.connected_components(graph)]

    def cleanup(self, index_id=None) -> None:
        return None


def test_fake_clustering_satisfies_protocol_via_isinstance():
    assert isinstance(FakeConnectedComponentsClustering(), ClusteringPort)


def test_cluster_returns_connected_components():
    graph = nx.Graph()
    graph.add_edge("a", "b")
    graph.add_node("c")
    clustering: ClusteringPort = FakeConnectedComponentsClustering()
    result = clustering.cluster(graph)
    assert {"a", "b"} in result
    assert {"c"} in result


def test_cluster_accepts_kwargs():
    graph = nx.Graph()
    graph.add_edge("x", "y")
    clustering: ClusteringPort = FakeConnectedComponentsClustering()
    result = clustering.cluster(graph, inflation=2.0)
    assert result == [{"x", "y"}]


def test_incomplete_implementation_fails_isinstance_check():
    class IncompleteClustering:
        def not_cluster(self, graph):
            return []

    assert not isinstance(IncompleteClustering(), ClusteringPort)
