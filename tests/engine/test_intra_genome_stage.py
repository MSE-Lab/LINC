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


def _result(
    a: Gene, b: Gene, identity: float, query_cover=None, subject_cover=None
) -> AlignmentResult:
    return AlignmentResult(
        query_id=a.gene_id,
        subject_id=b.gene_id,
        identity=identity,
        query_cover=identity if query_cover is None else query_cover,
        subject_cover=identity if subject_cover is None else subject_cover,
        evalue=0.0,
        bitscore=identity * 100,
        alignment_length=min(a.length, b.length),
    )


def _stage(
    aligner, cutoff_i=0.5, cutoff_ii=0.8, cutoff_c=0.0, cov_mode=0, max_per_cluster=5
) -> IntraGenomeStage:
    return IntraGenomeStage(
        aligner=aligner,
        clustering=UnionFindClustering(),
        cutoff_i=cutoff_i,
        cutoff_ii=cutoff_ii,
        cutoff_c=cutoff_c,
        cov_mode=cov_mode,
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


def test_high_identity_but_low_coverage_pair_triggers_multi_representative():
    # 第二轮漏洞3回归测试: a-b-c 两两 identity=0.9(均>=cutoff_ii=0.8), 但
    # a-c 这一对 coverage 很低(0.1, < cutoff_c=0.5), a-b/b-c 的 coverage
    # 都达标(0.9)。a-c 因 coverage 不达标不会直接连边(仍可经由 b 间接
    # 连通成一个簇), 但代表选取阶段若仍只看"最小 identity"会误判为"整体
    # 均质"而只选 1 个中心代表; 修复后 _select_representatives 要求簇内
    # 任意两个成员同时满足 identity>=cutoff_ii 且 coverage>=cutoff_c, a-c
    # 这一对不达标应触发 multi_representative_clustering, 选出 >= 2 个代表。
    a, b, c = (_make_gene("g1", i) for i in range(1, 4))
    aligner = MockAligner()
    aligner.set_result(a.gene_id, b.gene_id, _result(a, b, identity=0.9, query_cover=0.9, subject_cover=0.9))
    aligner.set_result(b.gene_id, c.gene_id, _result(b, c, identity=0.9, query_cover=0.9, subject_cover=0.9))
    aligner.set_result(a.gene_id, c.gene_id, _result(a, c, identity=0.9, query_cover=0.1, subject_cover=0.1))
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=[a, b, c])

    result = _stage(aligner, cutoff_i=0.5, cutoff_ii=0.8, cutoff_c=0.5).run(genome)

    assert result.n_clusters == 1
    cluster = result.local_clusters[0]
    assert len(cluster.member_ids) == 3
    # 修复前会错误地只选出 1 个中心代表; 修复后应为多代表(>= 2)。
    assert len(cluster.representative_ids) >= 2


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


def test_cutoff_c_filters_out_low_coverage_hit_even_if_identity_passes():
    # identity 达标但覆盖度很低, 设置 cutoff_c=0.5 应让两个基因各自独立成孤儿簇,
    # 而不是因 identity 达标而被连边。
    g1, g2 = _make_gene("g1", 1), _make_gene("g1", 2)
    aligner = MockAligner()
    aligner.set_result(
        g1.gene_id, g2.gene_id,
        _result(g1, g2, identity=0.9, query_cover=0.1, subject_cover=0.1),
    )
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=[g1, g2])

    result = _stage(aligner, cutoff_i=0.5, cutoff_c=0.5).run(genome)

    assert result.n_clusters == 2
    assert set(result.orphan_genes) == {g1.gene_id, g2.gene_id}


def test_cov_mode_subject_only_uses_subject_cover():
    # cov_mode=1(仅看 subject_cover): query_cover 很低但 subject_cover 达标,
    # 在该模式下应被连边成簇(与 cov_mode=0 双向模式的结果不同)。
    g1, g2 = _make_gene("g1", 1), _make_gene("g1", 2)
    aligner = MockAligner()
    aligner.set_result(
        g1.gene_id, g2.gene_id,
        _result(g1, g2, identity=0.9, query_cover=0.1, subject_cover=0.9),
    )
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=[g1, g2])

    result = _stage(aligner, cutoff_i=0.5, cutoff_c=0.5, cov_mode=1).run(genome)

    assert result.n_clusters == 1
    assert set(result.orphan_genes) == set()


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


def test_pairwise_hits_includes_all_non_self_comparisons_for_ssn_builder_reuse():
    # pairwise_hits 应包含全部非自身比对结果(不论是否通过 cutoff 成簇),
    # 供 SSNBuilderStage 补边阶段直接复用, 避免重复发起比对
    # (2026-10-02 与用户确认)。
    g1, g2, g3 = (_make_gene("g1", i) for i in range(1, 4))
    aligner = MockAligner()
    aligner.set_result(g1.gene_id, g2.gene_id, _result(g1, g2, 0.9))
    # g1-g3 低于 cutoff_i, 不会成簇, 但仍应出现在 pairwise_hits 里
    aligner.set_result(g1.gene_id, g3.gene_id, _result(g1, g3, 0.1))
    aligner.set_result(g2.gene_id, g3.gene_id, _result(g2, g3, 0.1))
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=[g1, g2, g3])

    result = _stage(aligner, cutoff_i=0.5).run(genome)

    pairs = {frozenset((h.query_id, h.subject_id)) for h in result.pairwise_hits}
    assert pairs == {
        frozenset((g1.gene_id, g2.gene_id)),
        frozenset((g1.gene_id, g3.gene_id)),
        frozenset((g2.gene_id, g3.gene_id)),
    }


def test_cleans_up_search_intermediate_by_default(monkeypatch):
    """数据管理优化第3项(2026-10-03): search() 的中间文件(search_intra_self/)
    结果已转存进 pairwise_hits, 默认(cleanup_intermediate=True)应在消费完
    结果后立即调用 `cleanup_search` 清理, index/tag 需与 search() 一致。
    """
    g1, g2 = _make_gene("g1", 1), _make_gene("g1", 2)
    aligner = MockAligner()
    aligner.set_result(g1.gene_id, g2.gene_id, _result(g1, g2, 0.9))
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=[g1, g2])

    cleanup_calls = []
    original_cleanup = aligner.cleanup_search

    def _recording_cleanup(index, tag=None):
        cleanup_calls.append((index, tag))
        return original_cleanup(index, tag=tag)

    monkeypatch.setattr(aligner, "cleanup_search", _recording_cleanup)

    stage = IntraGenomeStage(aligner=aligner, clustering=UnionFindClustering())
    result = stage.run(genome)

    assert cleanup_calls == [(result.genome_index, "intra_self")]


def test_keep_intermediate_skips_cleanup(monkeypatch):
    """`cleanup_intermediate=False` 时不应调用 `cleanup_search`(保留中间
    文件以便调试, 2026-10-03 与用户确认)。
    """
    g1, g2 = _make_gene("g1", 1), _make_gene("g1", 2)
    aligner = MockAligner()
    aligner.set_result(g1.gene_id, g2.gene_id, _result(g1, g2, 0.9))
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=[g1, g2])

    cleanup_calls = {"n": 0}
    original_cleanup = aligner.cleanup_search

    def _counting_cleanup(index, tag=None):
        cleanup_calls["n"] += 1
        return original_cleanup(index, tag=tag)

    monkeypatch.setattr(aligner, "cleanup_search", _counting_cleanup)

    stage = IntraGenomeStage(
        aligner=aligner, clustering=UnionFindClustering(), cleanup_intermediate=False
    )
    stage.run(genome)

    assert cleanup_calls["n"] == 0
