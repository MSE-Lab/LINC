"""`linc/core/entities.py` 的单元测试。"""
from linc.core.entities import (
    AlignmentResult,
    Cluster,
    Gene,
    Genome,
    GenomeStatus,
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
    edge = SSNEdge(source="g1|1", target="g2|1", identity=0.8, coverage=0.5)
    assert edge.weight("identity_only") == 0.8
    assert edge.weight("coverage_only") == 0.5
    assert edge.weight("identity_x_coverage") == 0.4
    assert edge.weight() == 0.4  # 默认模式与显式指定一致


def test_ortholog_family_size_property():
    family = OrthologFamily(
        family_id="fam1",
        component_id=0,
        member_genes={"g1|1", "g1|2", "g2|1"},
        edge_density=0.75,
    )
    assert family.size == 3
