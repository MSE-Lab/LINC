"""UnionFindClustering 单元测试。"""
import networkx as nx

from linc.adapters.clustering.union_find_adapter import UnionFindClustering
from linc.ports.clustering_port import ClusteringPort


def test_union_find_satisfies_clustering_port():
    assert isinstance(UnionFindClustering(), ClusteringPort)


def test_single_connected_component():
    graph = nx.Graph()
    graph.add_edge("a", "b")
    graph.add_edge("b", "c")
    result = UnionFindClustering().cluster(graph)
    assert result == [{"a", "b", "c"}]


def test_multiple_disjoint_components():
    graph = nx.Graph()
    graph.add_edge("a", "b")
    graph.add_edge("c", "d")
    result = UnionFindClustering().cluster(graph)
    assert {"a", "b"} in result
    assert {"c", "d"} in result
    assert len(result) == 2


def test_isolated_node_forms_its_own_singleton_cluster():
    graph = nx.Graph()
    graph.add_edge("a", "b")
    graph.add_node("orphan")
    result = UnionFindClustering().cluster(graph)
    assert {"orphan"} in result
    assert {"a", "b"} in result


def test_empty_graph_returns_empty_list():
    graph = nx.Graph()
    result = UnionFindClustering().cluster(graph)
    assert result == []


def test_accepts_arbitrary_kwargs_without_error():
    graph = nx.Graph()
    graph.add_edge("x", "y")
    result = UnionFindClustering().cluster(graph, inflation=2.0, threshold=0.5)
    assert result == [{"x", "y"}]


def test_chain_merges_into_one_component_regardless_of_edge_order():
    graph = nx.Graph()
    graph.add_edges_from([("d", "c"), ("a", "b"), ("b", "c")])
    result = UnionFindClustering().cluster(graph)
    assert result == [{"a", "b", "c", "d"}]
