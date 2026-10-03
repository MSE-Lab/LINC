"""`linc/engine/family_refine_stage.py` 的单元测试。

测试策略: 直接手工构造 CompleteSSN 夹具(而非跑完整上游流水线), 用
UnionFindClustering 作为假的 ClusteringPort 实现(不依赖本机安装的外部
mcl 二进制, 对每个连通分量/子图直接返回其自身作为唯一一簇), 聚焦验证
FamilyRefineStage 本身的连通分量拆分/孤立节点处理/family_id 排序/
edge_density 计算是否符合预期。
"""
from typing import List, Set

import networkx as nx
import pytest

from linc.adapters.clustering.union_find_adapter import UnionFindClustering
from linc.core.entities import CompleteSSN, SSNEdge
from linc.engine.family_refine_stage import FamilyRefineStage


class _OneClusterPerComponent:
    """假 ClusteringPort: 把传入图里每个连通分量各自当作一个簇返回(不做
    真正的 MCL 细分), 用于验证"连通分量拆分 -> 调用 cluster() -> 汇总"这
    条主干逻辑, 不依赖真实 MCL 效果。

    工程审查问题9(2026-10-03 与用户确认)把"每个连通分量一次cluster()
    调用"合并为"全图一次cluster()调用"后, 本假实现必须按连通分量分别
    返回簇(而不是盲目把传入的全部节点合并成一个簇), 才能在"传入的图
    可能包含多个不连通分量"的新调用方式下正确模拟真实 mcl 的行为
    (mcl 对不连通的子图天然不会产生交叉簇)。
    """

    def cluster(self, graph: nx.Graph, **params) -> List[Set[str]]:
        return [set(component) for component in nx.connected_components(graph)]

    def cleanup(self, index_id=None) -> None:
        return None


class _SplitEveryEdgeIntoOwnCluster:
    """假 ClusteringPort: 把每条边拆成独立的簇(用于验证一个连通分量可以
    产出多个 OrthologFamily 的场景)。"""

    def cluster(self, graph: nx.Graph, **params) -> List[Set[str]]:
        return [set(edge) for edge in graph.edges]

    def cleanup(self, index_id=None) -> None:
        return None


def _edge(a: str, b: str, identity: float = 0.9, coverage: float = 0.9, bitscore: float = 100.0) -> SSNEdge:
    return SSNEdge(source=a, target=b, identity=identity, coverage=coverage, bitscore=bitscore)


def test_single_connected_component_produces_one_family():
    complete_ssn = CompleteSSN(edges={_edge("A", "B"), _edge("B", "C")})
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())

    families = stage.run(complete_ssn)

    assert len(families) == 1
    assert families[0].member_genes == {"A", "B", "C"}
    assert families[0].component_id == 1
    assert families[0].family_id == "OG0000000"


def test_multiple_connected_components_each_get_own_component_id():
    complete_ssn = CompleteSSN(edges={_edge("A", "B"), _edge("C", "D")})
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())

    families = stage.run(complete_ssn)

    assert len(families) == 2
    component_ids = {f.component_id for f in families}
    assert component_ids == {1, 2}
    member_sets = {frozenset(f.member_genes) for f in families}
    assert member_sets == {frozenset({"A", "B"}), frozenset({"C", "D"})}


def test_isolated_node_in_ssn_is_not_sent_to_mcl():
    """CompleteSSN 本身只含边, 孤立基因不会出现在图中(与 SSNBuilderStage
    的既有边界一致), 这里验证的是"连通分量恰好只有2个节点"这种最小连通
    场景能正确产出一个家族, 不依赖 mcl 调用次数(由 mcl_adapter 自己的
    测试覆盖"无边图不调用 mcl"这一行为)。
    """
    complete_ssn = CompleteSSN(edges={_edge("A", "B")})
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())

    families = stage.run(complete_ssn)

    assert len(families) == 1
    assert families[0].member_genes == {"A", "B"}


def test_family_id_assignment_is_deterministic_and_sorted():
    complete_ssn = CompleteSSN(edges={
        _edge("Z", "Y"), _edge("B", "A"),
    })
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())

    families = stage.run(complete_ssn)

    # component_id 按连通分量最小节点排序分配(A/B分量先于Y/Z分量)。
    ids_by_component = {f.component_id: f.family_id for f in families}
    assert len(families) == 2
    family_ids = sorted(f.family_id for f in families)
    assert family_ids == ["OG0000000", "OG0000001"]
    # 包含 A 的家族应排在包含 Y 的家族之前。
    family_with_a = next(f for f in families if "A" in f.member_genes)
    family_with_y = next(f for f in families if "Y" in f.member_genes)
    assert family_with_a.family_id < family_with_y.family_id


def test_edge_density_computed_for_induced_subgraph():
    # A-B-C-A 三角形(3条边, 3个节点): 密度 = 2*3/(3*2) = 1.0(完全图)。
    complete_ssn = CompleteSSN(edges={_edge("A", "B"), _edge("B", "C"), _edge("A", "C")})
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())

    families = stage.run(complete_ssn)

    assert len(families) == 1
    assert families[0].edge_density == pytest.approx(1.0)


def test_edge_density_less_than_one_for_sparse_subgraph():
    # A-B-C-D 链(3条边, 4个节点): 密度 = 2*3/(4*3) = 0.5。
    complete_ssn = CompleteSSN(edges={_edge("A", "B"), _edge("B", "C"), _edge("C", "D")})
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())

    families = stage.run(complete_ssn)

    assert len(families) == 1
    assert families[0].edge_density == pytest.approx(0.5)


def test_connected_component_can_be_split_into_multiple_families_by_mcl():
    # A-B-C 链, 假聚类器把每条边拆成独立簇 -> 应产出2个家族, 但同属一个
    # component_id(因为它们在同一个 SSN 连通分量内)。
    complete_ssn = CompleteSSN(edges={_edge("A", "B"), _edge("B", "C")})
    stage = FamilyRefineStage(clustering=_SplitEveryEdgeIntoOwnCluster())

    families = stage.run(complete_ssn)

    assert len(families) == 2
    assert {f.component_id for f in families} == {1}
    member_sets = {frozenset(f.member_genes) for f in families}
    assert member_sets == {frozenset({"A", "B"}), frozenset({"B", "C"})}


def test_edge_weight_mode_is_passed_through_to_graph_weight():
    """edge_weight_mode 应影响图上的 weight 属性(供真实 MCL 消费),
    这里用一个会检查 weight 的假聚类器验证传递链路正确。
    """
    captured_weights = {}

    class _CaptureWeights:
        def cluster(self, graph: nx.Graph, **params) -> List[Set[str]]:
            for u, v, data in graph.edges(data=True):
                captured_weights[frozenset((u, v))] = data["weight"]
            return [set(graph.nodes)]

        def cleanup(self, index_id=None) -> None:
            pass

    complete_ssn = CompleteSSN(edges={_edge("A", "B", identity=0.5, coverage=0.5, bitscore=999.0)})
    stage = FamilyRefineStage(clustering=_CaptureWeights(), edge_weight_mode="bitscore")

    stage.run(complete_ssn)

    assert captured_weights[frozenset(("A", "B"))] == pytest.approx(999.0)


def test_works_with_real_union_find_clustering_as_fallback_implementation():
    """UnionFindClustering 同样满足 ClusteringPort 接口, 用真实实现(而非
    假聚类器)跑一遍完整流程, 确认没有对假聚类器产生隐性依赖。"""
    complete_ssn = CompleteSSN(edges={_edge("A", "B"), _edge("B", "C")})
    stage = FamilyRefineStage(clustering=UnionFindClustering())

    families = stage.run(complete_ssn)

    assert len(families) == 1
    assert families[0].member_genes == {"A", "B", "C"}


def test_empty_complete_ssn_produces_no_families():
    complete_ssn = CompleteSSN(edges=set())
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())

    families = stage.run(complete_ssn)

    assert families == []


def test_gene_with_zero_edges_anywhere_is_not_silently_dropped():
    """2026-10-02 修复: 若某基因在 CompleteSSN.edges 中完全不出现(局部
    孤儿+全局也未找到任何可连代表), 传入 all_gene_ids 后应被显式补成
    size=1 的独立家族, 而不是在 nx.connected_components() 中彻底消失。
    """
    complete_ssn = CompleteSSN(edges={_edge("A", "B")})
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())

    # "C" 从未出现在任何边里, 但仍属于全部基因集合。
    families = stage.run(complete_ssn, all_gene_ids={"A", "B", "C"})

    assert len(families) == 2
    member_sets = {frozenset(f.member_genes) for f in families}
    assert member_sets == {frozenset({"A", "B"}), frozenset({"C"})}
    lone_family = next(f for f in families if f.member_genes == {"C"})
    assert lone_family.edge_density == pytest.approx(0.0)


def test_all_gene_ids_omitted_keeps_backward_compatible_behavior():
    """未传入 all_gene_ids 时(默认None), 行为与修复前完全一致:
    零边基因不会被补充, 保持向后兼容。"""
    complete_ssn = CompleteSSN(edges={_edge("A", "B")})
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())

    families = stage.run(complete_ssn)

    assert len(families) == 1
    assert families[0].member_genes == {"A", "B"}


def test_multiple_isolated_genes_each_become_own_size_one_family():
    """多个零边孤立基因应分别成为各自独立的 size=1 家族, 而不是被合并
    成一个家族或互相干扰。"""
    complete_ssn = CompleteSSN(edges=set())
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())

    families = stage.run(complete_ssn, all_gene_ids={"X", "Y", "Z"})

    assert len(families) == 3
    member_sets = {frozenset(f.member_genes) for f in families}
    assert member_sets == {frozenset({"X"}), frozenset({"Y"}), frozenset({"Z"})}
    assert all(f.edge_density == pytest.approx(0.0) for f in families)


class _IncompletePartitionClustering:
    """假 ClusteringPort: 模拟 mcl 异常返回不完整分区(丢失其中一个节点),
    用于验证工程审查问题5的兜底补齐逻辑(2026-10-03 与用户确认)。"""

    def cluster(self, graph: nx.Graph, **params) -> List[Set[str]]:
        nodes = sorted(graph.nodes)
        # 故意丢弃最后一个节点, 模拟 mcl 返回不完整分区。
        return [set(nodes[:-1])]

    def cleanup(self, index_id=None) -> None:
        return None


def test_incomplete_mcl_partition_missing_genes_are_backfilled_as_singletons(caplog):
    """工程审查问题5(2026-10-03 与用户确认): 若 mcl 对某连通分量返回的
    分区未覆盖全部节点(工具异常/退化输入等), 缺失的基因应被各自补回为
    size=1 的独立家族, 而不是静默从最终结果中消失; 同时应打印 warning
    日志便于事后排查。
    """
    complete_ssn = CompleteSSN(edges={_edge("A", "B"), _edge("B", "C")})
    stage = FamilyRefineStage(clustering=_IncompletePartitionClustering())

    with caplog.at_level("WARNING"):
        families = stage.run(complete_ssn)

    # A/B 被 mcl 正常分到一簇, C 被 mcl"异常丢失", 应被补回为独立家族。
    assert len(families) == 2
    member_sets = {frozenset(f.member_genes) for f in families}
    assert member_sets == {frozenset({"A", "B"}), frozenset({"C"})}
    assert any("缺失" in record.message for record in caplog.records)


class _CountingClustering:
    """假 ClusteringPort: 记录 cluster() 被调用的次数及每次传入的参数,
    用于验证工程审查问题9(2026-10-03 与用户确认)的"全图一次调用"优化:
    无论图里有多少个连通分量, cluster() 都应只被调用1次(而非N次)。"""

    def __init__(self) -> None:
        self.call_count = 0
        self.received_threads = []

    def cluster(self, graph: nx.Graph, **params) -> List[Set[str]]:
        self.call_count += 1
        self.received_threads.append(params.get("threads"))
        return [set(component) for component in nx.connected_components(graph)]

    def cleanup(self, index_id=None) -> None:
        return None


def test_multiple_components_trigger_only_one_cluster_call():
    """工程审查问题9(2026-10-03 与用户确认): 把"每个连通分量一次cluster()
    调用"合并为"全图一次调用", 调用次数应从 O(分量数) 降到 O(1), 即使有
    3个独立的非孤立连通分量, 也应只触发1次 cluster() 调用; 且每个分量
    仍各自保有独立的 component_id(与改造前语义一致)。"""
    complete_ssn = CompleteSSN(edges={
        _edge("A", "B"), _edge("C", "D"), _edge("E", "F"),
    })
    fake = _CountingClustering()
    stage = FamilyRefineStage(clustering=fake)

    families = stage.run(complete_ssn)

    assert fake.call_count == 1
    assert len(families) == 3
    component_ids = {f.component_id for f in families}
    assert component_ids == {1, 2, 3}
    member_sets = {frozenset(f.member_genes) for f in families}
    assert member_sets == {
        frozenset({"A", "B"}), frozenset({"C", "D"}), frozenset({"E", "F"}),
    }


def test_threads_param_is_passed_through_to_cluster_call():
    """工程审查问题9(2026-10-03 与用户确认): FamilyRefineStage 的
    `threads` 构造参数应原样传给 cluster(threads=...), 供
    MCLClustering 转换为 mcl 的 `-te` 参数。"""
    complete_ssn = CompleteSSN(edges={_edge("A", "B")})
    fake = _CountingClustering()
    stage = FamilyRefineStage(clustering=fake, threads=8)

    stage.run(complete_ssn)

    assert fake.received_threads == [8]


def test_isolated_genes_do_not_trigger_any_cluster_call():
    """若全部基因都是零边孤立节点(无任何非孤立连通分量), 不应触发任何
    cluster() 调用(与改造前"孤立节点不调用mcl"的既有行为一致)。"""
    complete_ssn = CompleteSSN(edges=set())
    fake = _CountingClustering()
    stage = FamilyRefineStage(clustering=fake)

    families = stage.run(complete_ssn, all_gene_ids={"X", "Y"})

    assert fake.call_count == 0
    assert len(families) == 2


def test_incomplete_mcl_partition_all_missing_backfills_every_node(caplog):
    """mcl 返回完全空分区(最极端的异常情况)时, 连通分量内的全部节点都
    应被补回为各自独立的 size=1 家族。"""

    class _EmptyPartitionClustering:
        def cluster(self, graph: nx.Graph, **params) -> List[Set[str]]:
            return []

        def cleanup(self, index_id=None) -> None:
            return None

    complete_ssn = CompleteSSN(edges={_edge("A", "B")})
    stage = FamilyRefineStage(clustering=_EmptyPartitionClustering())

    with caplog.at_level("WARNING"):
        families = stage.run(complete_ssn)

    assert len(families) == 2
    member_sets = {frozenset(f.member_genes) for f in families}
    assert member_sets == {frozenset({"A"}), frozenset({"B"})}
