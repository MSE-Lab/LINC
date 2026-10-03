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


def test_cleanup_is_noop_since_no_intermediate_files_are_written():
    """UnionFindClustering 不落盘任何中间文件, `cleanup()` 应为空操作
    (数据管理优化第3项, 2026-10-03 与用户确认)。"""
    # 不应抛出任何异常。
    UnionFindClustering().cleanup(index_id="anything")
    UnionFindClustering().cleanup()
