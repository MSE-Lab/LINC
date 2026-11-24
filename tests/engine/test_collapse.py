"""`linc/engine/collapse.py` 与 FamilyRefineStage 折叠展开逻辑的单元测试。

使用手工构造的夹具, 不依赖外部 mcl/比对工具。
"""
from typing import List, Set

import networkx as nx

from linc.core.entities import (
    Cluster,
    CompleteSSN,
    GenomeIntraResult,
    GlobalRepresentativeLibrary,
    IncrementalResult,
    LocalCluster,
    SSNEdge,
)
from linc.engine.collapse import compute_collapse_map
from linc.engine.family_refine_stage import FamilyRefineStage


class _OneClusterPerComponent:
    def cluster(self, graph: nx.Graph, **params) -> List[Set[str]]:
        return [set(c) for c in nx.connected_components(graph)]

    def cleanup(self, index_id=None) -> None:
        return None


def _edge(a, b, identity=0.9, coverage=0.9):
    return SSNEdge(source=a, target=b, identity=identity, coverage=coverage, bitscore=100.0)


def _join(clusters, edges):
    lib = GlobalRepresentativeLibrary(clusters=clusters)
    return IncrementalResult(global_lib=lib, ssn_edges=set(edges))


def test_global_member_collapsed_to_highest_identity_rep():
    # 簇 c1 有两个代表 R1/R2, 成员 M 与两者都高相似, 依附 identity 更高的 R2
    cluster = Cluster(
        cluster_id="c1",
        representative_ids={"R1", "R2"},
        member_ids={"R1", "R2", "M"},
    )
    join = _join([cluster], [_edge("M", "R1", 0.85), _edge("M", "R2", 0.95)])
    result = compute_collapse_map(join, {}, cutoff_ii=0.8, cutoff_c=0.5)
    assert result == {"M": "R2"}


def test_low_identity_or_low_coverage_member_not_collapsed():
    cluster = Cluster(
        cluster_id="c1",
        representative_ids={"R1"},
        member_ids={"R1", "M1", "M2"},
    )
    join = _join(
        [cluster],
        [_edge("M1", "R1", identity=0.7), _edge("M2", "R1", identity=0.95, coverage=0.3)],
    )
    assert compute_collapse_map(join, {}, 0.8, 0.5) == {}


def test_inactive_cluster_ignored():
    cluster = Cluster(
        cluster_id="c1", representative_ids={"R1"},
        member_ids={"R1", "M"}, active=False,
    )
    join = _join([cluster], [_edge("M", "R1")])
    assert compute_collapse_map(join, {}, 0.8, 0.5) == {}


def test_local_members_follow_collapsed_local_representative():
    cluster = Cluster(
        cluster_id="c1", representative_ids={"R"}, member_ids={"R", "L"},
    )
    join = _join([cluster], [_edge("L", "R", 0.9)])
    # L 是局部代表(已折叠), X 是 L 的局部成员
    intra = GenomeIntraResult(
        genome_id="g",
        local_clusters=[
            LocalCluster(cluster_id="lc", member_ids={"L", "X"}, representative_ids=["L"]),
        ],
    )
    result = compute_collapse_map(join, {"g": intra}, 0.8, 0.5)
    assert result == {"L": "R", "X": "R"}


def test_local_members_not_collapsed_if_any_local_rep_not_collapsed():
    cluster = Cluster(
        cluster_id="c1", representative_ids={"R", "L2"}, member_ids={"R", "L1", "L2"},
    )
    join = _join([cluster], [_edge("L1", "R", 0.9)])
    intra = GenomeIntraResult(
        genome_id="g",
        local_clusters=[
            LocalCluster(
                cluster_id="lc", member_ids={"L1", "L2", "X"},
                representative_ids=["L1", "L2"],
            ),
        ],
    )
    result = compute_collapse_map(join, {"g": intra}, 0.8, 0.5)
    assert result == {"L1": "R"}


def test_local_member_tie_break_is_lexicographic_with_prefix_ids():
    # identity 并列(无 pairwise 记录)时取字典序最小的局部代表: "g|1" 是 "g|10" 的前缀, 应选 "g|1"
    cluster = Cluster(
        cluster_id="c1", representative_ids={"R1", "R2"},
        member_ids={"R1", "R2", "g|1", "g|10"},
    )
    join = _join([cluster], [_edge("g|1", "R1", 0.9), _edge("g|10", "R2", 0.9)])
    intra = GenomeIntraResult(
        genome_id="g",
        local_clusters=[
            LocalCluster(
                cluster_id="lc", member_ids={"g|1", "g|10", "X"},
                representative_ids=["g|1", "g|10"],
            ),
        ],
    )
    result = compute_collapse_map(join, {"g": intra}, 0.8, 0.5)
    assert result["X"] == "R1"


def test_family_refine_expands_collapsed_genes_into_anchor_family():
    ssn = CompleteSSN(edges={_edge("A", "B"), _edge("M", "A")})
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())
    families = stage.run(
        ssn,
        all_gene_ids=["A", "B", "M", "Z"],
        collapse_map={"M": "A"},
    )
    by_members = sorted(sorted(f.member_genes) for f in families)
    assert by_members == [["A", "B", "M"], ["Z"]]


def test_family_refine_without_collapse_map_unchanged():
    ssn = CompleteSSN(edges={_edge("A", "B")})
    stage = FamilyRefineStage(clustering=_OneClusterPerComponent())
    families = stage.run(ssn, all_gene_ids=["A", "B", "C"])
    assert sorted(sorted(f.member_genes) for f in families) == [["A", "B"], ["C"]]
