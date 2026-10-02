"""`linc/core/entities.py` 的单元测试。"""
from linc.core.entities import (
    AlignmentResult,
    Cluster,
    CompleteSSN,
    Gene,
    Genome,
    GenomeIntraResult,
    GenomeStatus,
    LocalCluster,
    OrthologFamily,
    SSNEdge,
)


def test_gene_is_frozen_and_immutable():
    gene = Gene(gene_id="g1|1", genome_id="g1", length=300)
    assert gene.gene_id == "g1|1"
    assert gene.sequence is None
    try:
        gene.length = 400  # type: ignore[misc]
    except Exception:
        pass
    else:
        raise AssertionError("Gene 应为不可变(frozen=True),赋值应抛出异常")


def test_genome_default_genes_list_is_independent():
    g1 = Genome(genome_id="g1", file_path="g1.fasta")
    g2 = Genome(genome_id="g2", file_path="g2.fasta")
    g1.genes.append(Gene(gene_id="g1|1", genome_id="g1", length=100))
    assert g1.genes != g2.genes
    assert g1.status is None
    assert GenomeStatus.LOADED.value == "loaded"


def test_alignment_result_coverage_takes_min():
    hit = AlignmentResult(
        query_id="g1|1",
        subject_id="g2|1",
        identity=0.9,
        query_cover=0.95,
        subject_cover=0.80,
        evalue=1e-50,
        bitscore=200.0,
        alignment_length=150,
    )
    assert hit.coverage == 0.80


def test_cluster_default_sets_are_independent_and_active_by_default():
    c1 = Cluster(cluster_id="c1")
    c2 = Cluster(cluster_id="c2")
    c1.member_ids.add("g1|1")
    assert c1.member_ids != c2.member_ids
    assert c1.active is True


def test_ssn_edge_weight_modes():
    edge = SSNEdge(source="g1|1", target="g2|1", identity=0.8, coverage=0.5, bitscore=200.0)
    # 旧别名(向后兼容)
    assert edge.weight("identity_only") == 0.8
    assert edge.weight("coverage_only") == 0.5
    # 新命名
    assert edge.weight("identity") == 0.8
    assert edge.weight("coverage") == 0.5
    assert edge.weight("identity_x_coverage") == 0.4
    assert edge.weight() == 0.4  # 默认模式与显式指定一致
    # 新增: identity/coverage 加权平均
    assert edge.weight("identity_coverage_weighted") == 0.65  # 0.5*0.8+0.5*0.5
    assert edge.weight("identity_coverage_weighted", alpha=0.0) == 0.5
    assert edge.weight("identity_coverage_weighted", alpha=1.0) == 0.8
    # 新增: bitscore 原始值, 不归一化
    assert edge.weight("bitscore") == 200.0


def test_ssn_edge_bitscore_defaults_to_zero_for_backward_compat():
    edge = SSNEdge(source="g1|1", target="g2|1", identity=0.8, coverage=0.5)
    assert edge.bitscore == 0.0


def test_alignment_result_length_fields_default_to_zero_for_backward_compat():
    hit = AlignmentResult(
        query_id="g1|1", subject_id="g2|1", identity=0.8,
        query_cover=0.9, subject_cover=0.9, evalue=1e-30,
        bitscore=150.0, alignment_length=100,
    )
    assert hit.query_length == 0
    assert hit.subject_length == 0


def test_ssn_edge_length_fields_default_to_zero_and_can_be_set():
    # 默认 0(向后兼容既有构造调用)
    edge_default = SSNEdge(source="g1|1", target="g2|1", identity=0.8, coverage=0.5)
    assert edge_default.source_length == 0
    assert edge_default.target_length == 0
    # 显式设置, 为未来 MCL 权重归一化预留(2026-10-02 与用户确认)
    edge = SSNEdge(
        source="g1|1", target="g2|1", identity=0.8, coverage=0.5,
        source_length=120, target_length=150,
    )
    assert edge.source_length == 120
    assert edge.target_length == 150


def test_complete_ssn_defaults_are_independent():
    ssn1 = CompleteSSN()
    ssn2 = CompleteSSN()
    ssn1.edges.add(SSNEdge(source="g1|1", target="g2|1", identity=0.8, coverage=0.5))
    assert ssn1.edges != ssn2.edges
    assert ssn1.edge_completion_mode == "representative_member"


def test_ssn_edge_is_hashable_and_usable_in_set():
    edge = SSNEdge(source="g1|1", target="g2|1", identity=0.8, coverage=0.5)
    # 可哈希, 可放入 set(验证 RepositoryPort.save_ssn_edges(edges: set) 场景可用)
    assert {edge} == {edge}


def test_ssn_edge_equality_ignores_direction():
    forward = SSNEdge(source="g1|1", target="g2|1", identity=0.8, coverage=0.5)
    backward = SSNEdge(source="g2|1", target="g1|1", identity=0.8, coverage=0.5)
    assert forward == backward
    assert hash(forward) == hash(backward)


def test_ssn_edge_equality_ignores_identity_coverage_float_noise():
    edge_a = SSNEdge(source="g1|1", target="g2|1", identity=0.8, coverage=0.5)
    edge_b = SSNEdge(source="g1|1", target="g2|1", identity=0.8000001, coverage=0.5)
    # 同一对基因, 即使 identity/coverage 因重新计算产生微小浮点误差, 仍视为同一条边
    assert edge_a == edge_b
    # 放入 set 后应去重为一条边(验证补边阶段"按基因对去重"的业务语义)
    assert len({edge_a, edge_b}) == 1


def test_ssn_edge_not_equal_to_different_gene_pair():
    edge_1 = SSNEdge(source="g1|1", target="g2|1", identity=0.8, coverage=0.5)
    edge_2 = SSNEdge(source="g1|1", target="g3|1", identity=0.8, coverage=0.5)
    assert edge_1 != edge_2


def test_ortholog_family_size_property():
    family = OrthologFamily(
        family_id="fam1",
        component_id=0,
        member_genes={"g1|1", "g1|2", "g2|1"},
        edge_density=0.75,
    )
    assert family.size == 3


def test_local_cluster_defaults_are_independent():
    c1 = LocalCluster(cluster_id="g1|cluster_0")
    c2 = LocalCluster(cluster_id="g1|cluster_1")
    c1.member_ids.add("g1|1")
    c1.representative_ids.append("g1|1")
    assert c1.member_ids != c2.member_ids
    assert c1.representative_ids != c2.representative_ids
    assert c1.is_orphan is False


def test_genome_intra_result_defaults_are_independent():
    r1 = GenomeIntraResult(genome_id="g1")
    r2 = GenomeIntraResult(genome_id="g2")
    r1.local_clusters.append(LocalCluster(cluster_id="g1|cluster_0"))
    assert r1.local_clusters != r2.local_clusters
    assert r1.n_genes == 0
    assert r1.n_clusters == 0
    assert r1.pairwise_hits == []


def test_genome_intra_result_pairwise_hits_defaults_are_independent():
    r1 = GenomeIntraResult(genome_id="g1")
    r2 = GenomeIntraResult(genome_id="g2")
    r1.pairwise_hits.append(AlignmentResult(
        query_id="g1|1", subject_id="g1|2", identity=0.9,
        query_cover=0.9, subject_cover=0.9, evalue=1e-30,
        bitscore=150.0, alignment_length=100,
    ))
    assert r1.pairwise_hits != r2.pairwise_hits
