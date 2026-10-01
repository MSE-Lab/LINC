"""IntraGenomeStage 单元测试。"""
from linc.adapters.aligners.mock_adapter import MockAligner
from linc.adapters.clustering.union_find_adapter import UnionFindClustering
from linc.core.entities import AlignmentResult, Gene, Genome
from linc.engine.intra_genome_stage import IntraGenomeStage


def _make_gene(genome_id: str, local_id: int, length: int = 100) -> Gene:
    return Gene(
        gene_id=f"{genome_id}|{local_id}",
        genome_id=genome_id,
        length=length,
        sequence="M" * length,
    )


def _result(a: Gene, b: Gene, identity: float) -> AlignmentResult:
    return AlignmentResult(
        query_id=a.gene_id,
        subject_id=b.gene_id,
        identity=identity,
        query_cover=identity,
        subject_cover=identity,
        evalue=0.0,
        bitscore=identity * 100,
        alignment_length=min(a.length, b.length),
    )


def _stage(aligner, cutoff_i=0.5, cutoff_ii=0.8, max_per_cluster=5) -> IntraGenomeStage:
    return IntraGenomeStage(
        aligner=aligner,
        clustering=UnionFindClustering(),
        cutoff_i=cutoff_i,
        cutoff_ii=cutoff_ii,
        max_per_cluster=max_per_cluster,
    )


def test_empty_genome_returns_empty_result():
    genome = Genome(genome_id="g0", file_path="g0.faa", genes=[])
    result = _stage(MockAligner()).run(genome)

    assert result.n_genes == 0
    assert result.n_clusters == 0
    assert result.local_clusters == []
    assert result.orphan_genes == []


def test_single_gene_is_orphan_and_self_representative():
    g1 = _make_gene("g1", 1)
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=[g1])

    result = _stage(MockAligner()).run(genome)

    assert result.n_clusters == 1
    assert result.orphan_genes == [g1.gene_id]
    assert result.local_clusters[0].is_orphan is True
    assert result.local_clusters[0].representative_ids == [g1.gene_id]


def test_two_dissimilar_genes_become_two_orphan_clusters():
    g1, g2 = _make_gene("g1", 1), _make_gene("g1", 2)
    aligner = MockAligner()
    aligner.set_result(g1.gene_id, g2.gene_id, _result(g1, g2, 0.1))
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=[g1, g2])

    result = _stage(aligner, cutoff_i=0.5).run(genome)

    assert result.n_clusters == 2
    assert set(result.orphan_genes) == {g1.gene_id, g2.gene_id}


def test_homogeneous_cluster_selects_single_central_representative():
    genes = [_make_gene("g1", i) for i in range(1, 4)]
    aligner = MockAligner()
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            aligner.set_result(genes[i].gene_id, genes[j].gene_id, _result(genes[i], genes[j], 0.9))
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=genes)

    result = _stage(aligner, cutoff_i=0.5, cutoff_ii=0.8).run(genome)

    assert result.n_clusters == 1
    cluster = result.local_clusters[0]
    assert cluster.is_orphan is False
    assert len(cluster.representative_ids) == 1
    assert cluster.representative_ids[0] in {g.gene_id for g in genes}


def test_divergent_cluster_uses_multi_representative_clustering():
    # 构造一个簇: a-b-c-d 两两 identity >= cutoff_i(保证成簇),
    # 但 a 与 b 高度相似(0.95), c 与 d 高度相似(0.95),
    # 跨子群(a/b vs c/d)相似度仅 0.55(< cutoff_ii)。
    # 判断依据为"最小" pairwise identity(=0.55) < cutoff_ii=0.8,
    # 应触发多代表逻辑;覆盖阈值同为 0.8,a/c 应分别覆盖各自子群,
    # 贪心应恰好选出 2 个代表(a 覆盖 {a,b},c 覆盖 {c,d}),不应退化为全选。
    a, b, c, d = (_make_gene("g1", i) for i in range(1, 5))
    aligner = MockAligner()
    aligner.set_result(a.gene_id, b.gene_id, _result(a, b, 0.95))
    aligner.set_result(c.gene_id, d.gene_id, _result(c, d, 0.95))
    aligner.set_result(a.gene_id, c.gene_id, _result(a, c, 0.55))
    aligner.set_result(a.gene_id, d.gene_id, _result(a, d, 0.55))
    aligner.set_result(b.gene_id, c.gene_id, _result(b, c, 0.55))
    aligner.set_result(b.gene_id, d.gene_id, _result(b, d, 0.55))
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=[a, b, c, d])

    result = _stage(aligner, cutoff_i=0.5, cutoff_ii=0.8).run(genome)

    assert result.n_clusters == 1
    cluster = result.local_clusters[0]
    assert cluster.is_orphan is False
    assert len(cluster.member_ids) == 4
    # 贪心覆盖应恰好选出 2 个代表,而不是退化为全体成员都是代表
    assert len(cluster.representative_ids) == 2
    reps = set(cluster.representative_ids)
    # 代表应能覆盖两个子群(每个子群恰有一个代表)
    assert reps & {a.gene_id, b.gene_id}
    assert reps & {c.gene_id, d.gene_id}


def test_representative_count_never_exceeds_max_per_cluster():
    # 构造 6 个两两低相似度(刚好达到 cutoff_i 成簇)的基因, 验证代表数不超过上限
    genes = [_make_gene("g1", i) for i in range(1, 7)]
    aligner = MockAligner()
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            identity = 0.5 if (i, j) != (0, 1) else 0.6
            aligner.set_result(genes[i].gene_id, genes[j].gene_id, _result(genes[i], genes[j], identity))
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=genes)

    result = _stage(aligner, cutoff_i=0.4, cutoff_ii=0.99, max_per_cluster=3).run(genome)

    assert result.n_clusters == 1
    cluster = result.local_clusters[0]
    assert len(cluster.representative_ids) <= 3


def test_n_representatives_matches_sum_across_clusters():
    g1, g2 = _make_gene("g1", 1), _make_gene("g1", 2)
    aligner = MockAligner()
    aligner.set_result(g1.gene_id, g2.gene_id, _result(g1, g2, 0.1))
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=[g1, g2])

    result = _stage(aligner).run(genome)

    assert result.n_representatives == sum(
        len(c.representative_ids) for c in result.local_clusters
    )
    assert result.n_representatives == len(result.local_representatives)
