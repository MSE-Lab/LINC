"""`linc/engine/ssn_builder_stage.py` 的单元测试。

测试策略: 直接手工构造 GenomeIntraResult/IncrementalResult 夹具(而非跑完整
GenomeLoadStage/IntraGenomeStage/GlobalJoinStage 流水线), 聚焦验证
SSNBuilderStage 本身在三种 edge_completion_mode 下的补边行为是否符合预期。
"""
import logging

from linc.adapters.aligners.mock_adapter import MockAligner
from linc.core.entities import (
    AlignmentResult,
    Gene,
    Genome,
    GenomeIntraResult,
    GlobalRepresentativeLibrary,
    IncrementalResult,
    LocalCluster,
    SSNEdge,
)
from linc.core.entities import Cluster as GlobalCluster
from linc.engine.ssn_builder_stage import SSNBuilderStage


def _gene(gene_id: str, genome_id: str, length: int = 100) -> Gene:
    return Gene(gene_id=gene_id, genome_id=genome_id, length=length)


def _hit(a: str, b: str, identity: float = 1.0, coverage: float = 1.0, bitscore: float = 150.0) -> AlignmentResult:
    return AlignmentResult(
        query_id=a, subject_id=b, identity=identity,
        query_cover=coverage, subject_cover=coverage,
        evalue=1e-50, bitscore=bitscore, alignment_length=100,
    )


def _build_two_genome_fixture():
    """两个基因组, 各自一个局部簇(1代表+2成员), 全局合并为一个簇。"""
    a1, a2, a3 = _gene("A|1", "A"), _gene("A|2", "A"), _gene("A|3", "A")
    b1, b2, b3 = _gene("B|1", "B"), _gene("B|2", "B"), _gene("B|3", "B")

    genomes = {
        "A": Genome(genome_id="A", file_path="a.faa", genes=[a1, a2, a3]),
        "B": Genome(genome_id="B", file_path="b.faa", genes=[b1, b2, b3]),
    }

    intra_results = {
        "A": GenomeIntraResult(
            genome_id="A",
            local_clusters=[
                LocalCluster(
                    cluster_id="A|cluster_0",
                    member_ids={"A|1", "A|2", "A|3"},
                    representative_ids=["A|1"],
                    is_orphan=False,
                )
            ],
            local_representatives=["A|1"],
            n_genes=3, n_clusters=1, n_representatives=1,
            pairwise_hits=[
                _hit("A|1", "A|2"), _hit("A|1", "A|3"), _hit("A|2", "A|3"),
            ],
        ),
        "B": GenomeIntraResult(
            genome_id="B",
            local_clusters=[
                LocalCluster(
                    cluster_id="B|cluster_0",
                    member_ids={"B|1", "B|2", "B|3"},
                    representative_ids=["B|1"],
                    is_orphan=False,
                )
            ],
            local_representatives=["B|1"],
            n_genes=3, n_clusters=1, n_representatives=1,
            pairwise_hits=[
                _hit("B|1", "B|2"), _hit("B|1", "B|3"), _hit("B|2", "B|3"),
            ],
        ),
    }

    global_lib = GlobalRepresentativeLibrary(
        clusters=[
            GlobalCluster(
                cluster_id="global|cluster_0",
                representative_ids={"A|1", "B|1"},
                member_ids={"A|1", "A|2", "A|3", "B|1", "B|2", "B|3"},
                active=True,
            )
        ]
    )
    incremental_result = IncrementalResult(
        global_lib=global_lib,
        ssn_edges={SSNEdge(source="A|1", target="B|1", identity=1.0, coverage=1.0, bitscore=150.0)},
        processing_order=["A", "B"],
        n_clusters_final=1,
        n_representatives_final=2,
    )
    return genomes, intra_results, incremental_result


def test_none_mode_returns_sparse_ssn_unchanged():
    genomes, intra_results, incremental_result = _build_two_genome_fixture()
    stage = SSNBuilderStage(aligner=MockAligner(), edge_completion_mode="none")

    result = stage.run(genomes, intra_results, incremental_result)

    assert result.edge_completion_mode == "none"
    assert result.edges == incremental_result.ssn_edges


def test_representative_member_mode_adds_rep_member_edges_only():
    genomes, intra_results, incremental_result = _build_two_genome_fixture()
    stage = SSNBuilderStage(aligner=MockAligner(), edge_completion_mode="representative_member")

    result = stage.run(genomes, intra_results, incremental_result)

    expected_pairs = {
        frozenset(("A|1", "B|1")),
        frozenset(("A|1", "A|2")),
        frozenset(("A|1", "A|3")),
        frozenset(("B|1", "B|2")),
        frozenset(("B|1", "B|3")),
    }
    actual_pairs = {frozenset((e.source, e.target)) for e in result.edges}
    assert actual_pairs == expected_pairs
    # 不应包含同局部簇内的成员-成员边(representative_member 模式不补这部分)
    assert frozenset(("A|2", "A|3")) not in actual_pairs
    # 不应包含跨基因组的成员-成员边
    assert frozenset(("A|2", "B|2")) not in actual_pairs


def test_full_mode_adds_member_member_and_cross_genome_edges():
    genomes, intra_results, incremental_result = _build_two_genome_fixture()
    stage = SSNBuilderStage(aligner=MockAligner(), edge_completion_mode="full")

    result = stage.run(genomes, intra_results, incremental_result)

    actual_pairs = {frozenset((e.source, e.target)) for e in result.edges}
    # 同基因组内成员-成员边(免费复用)应该存在
    assert frozenset(("A|2", "A|3")) in actual_pairs
    assert frozenset(("B|2", "B|3")) in actual_pairs
    # 跨基因组的成员-成员边(需要真正发起新比对)也应该存在
    assert frozenset(("A|2", "B|2")) in actual_pairs
    assert frozenset(("A|3", "B|3")) in actual_pairs
    # 代表-成员边依旧保留
    assert frozenset(("A|1", "A|2")) in actual_pairs
    assert frozenset(("A|1", "B|1")) in actual_pairs


def test_full_mode_reuses_intra_hits_without_recomputing(monkeypatch):
    """同基因组内的边应直接复用 pairwise_hits, 不应触发新的比对调用。"""
    genomes, intra_results, incremental_result = _build_two_genome_fixture()

    call_count = {"n": 0}
    aligner = MockAligner()
    original = aligner.align_query_vs_subject

    def _counting_align(query, subject):
        call_count["n"] += 1
        return original(query, subject)

    monkeypatch.setattr(aligner, "align_query_vs_subject", _counting_align)

    stage = SSNBuilderStage(aligner=aligner, edge_completion_mode="full")
    stage.run(genomes, intra_results, incremental_result)

    # 只有跨基因组(A vs B)这一组基因组对会触发一次批量比对调用
    assert call_count["n"] == 1


def test_single_genome_cluster_skips_cross_genome_completion():
    """全局簇只涉及单一基因组时, full 模式不应触发任何新比对。"""
    a1, a2 = _gene("A|1", "A"), _gene("A|2", "A")
    genomes = {"A": Genome(genome_id="A", file_path="a.faa", genes=[a1, a2])}
    intra_results = {
        "A": GenomeIntraResult(
            genome_id="A",
            local_clusters=[
                LocalCluster(
                    cluster_id="A|cluster_0", member_ids={"A|1", "A|2"},
                    representative_ids=["A|1"], is_orphan=False,
                )
            ],
            local_representatives=["A|1"],
            n_genes=2, n_clusters=1, n_representatives=1,
            pairwise_hits=[_hit("A|1", "A|2")],
        ),
    }
    global_lib = GlobalRepresentativeLibrary(
        clusters=[GlobalCluster(cluster_id="global|cluster_0",
                                 representative_ids={"A|1"}, member_ids={"A|1", "A|2"}, active=True)]
    )
    incremental_result = IncrementalResult(global_lib=global_lib, ssn_edges=set())

    class _FailingAligner(MockAligner):
        def align_query_vs_subject(self, query, subject):
            raise AssertionError("单基因组簇不应触发跨基因组比对")

    stage = SSNBuilderStage(aligner=_FailingAligner(), edge_completion_mode="full")
    result = stage.run(genomes, intra_results, incremental_result)

    assert frozenset(("A|1", "A|2")) in {frozenset((e.source, e.target)) for e in result.edges}


def test_mega_cluster_warning_is_logged_but_does_not_skip_completion(caplog):
    genomes, intra_results, incremental_result = _build_two_genome_fixture()
    stage = SSNBuilderStage(
        aligner=MockAligner(),
        edge_completion_mode="full",
        mega_cluster_warning_threshold=1,
    )

    with caplog.at_level(logging.WARNING):
        result = stage.run(genomes, intra_results, incremental_result)

    assert any("超过阈值" in record.message for record in caplog.records)
    # 即使触发了 warning, 补边依然完整执行(不自动降级/拆分)
    actual_pairs = {frozenset((e.source, e.target)) for e in result.edges}
    assert frozenset(("A|2", "B|2")) in actual_pairs


def test_mega_cluster_warning_not_triggered_when_threshold_not_exceeded():
    genomes, intra_results, incremental_result = _build_two_genome_fixture()
    stage = SSNBuilderStage(
        aligner=MockAligner(),
        edge_completion_mode="full",
        mega_cluster_warning_threshold=1000,
    )
    import logging as _logging
    logger = _logging.getLogger("linc.engine.ssn_builder_stage")
    records = []
    handler = _logging.Handler()
    handler.emit = lambda record: records.append(record)
    logger.addHandler(handler)
    try:
        stage.run(genomes, intra_results, incremental_result)
    finally:
        logger.removeHandler(handler)
    assert not any("超过阈值" in r.getMessage() for r in records)


def test_inactive_global_clusters_are_skipped():
    genomes, intra_results, incremental_result = _build_two_genome_fixture()
    incremental_result.global_lib.clusters[0].active = False

    stage = SSNBuilderStage(aligner=MockAligner(), edge_completion_mode="full")
    result = stage.run(genomes, intra_results, incremental_result)

    actual_pairs = {frozenset((e.source, e.target)) for e in result.edges}
    # 非活跃簇不应触发跨基因组补边, 但同基因组内的免费复用边依旧保留
    assert frozenset(("A|2", "B|2")) not in actual_pairs
    assert frozenset(("A|2", "A|3")) in actual_pairs
