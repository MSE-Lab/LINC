"""GlobalJoinStage 单元测试。"""
from linc.adapters.aligners.mock_adapter import MockAligner
from linc.core.entities import AlignmentResult, Gene, Genome, GenomeIntraResult
from linc.engine.global_join_stage import GlobalJoinStage


def _gene(genome_id: str, local_id: int, length: int = 100) -> Gene:
    return Gene(gene_id=f"{genome_id}|{local_id}", genome_id=genome_id, length=length)


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


def _make_genome_and_intra(genome_id: str, genes, representatives):
    """构造一个只有一个局部簇、全部基因都是局部代表的简单基因组。"""
    genome = Genome(genome_id=genome_id, file_path=f"{genome_id}.faa", genes=genes)
    intra = GenomeIntraResult(
        genome_id=genome_id,
        local_representatives=representatives,
        n_genes=len(genes),
        n_representatives=len(representatives),
        genome_index=genes,  # MockAligner.build_index 返回 genes 列表本身, 这里直接模拟
    )
    return genome, intra


def _stage(
    aligner, cutoff_i=0.5, cutoff_ii=0.8, cutoff_c=0.0, cov_mode=0, max_per_cluster=5,
    touch_evalue=None,
) -> GlobalJoinStage:
    return GlobalJoinStage(
        aligner=aligner,
        cutoff_i=cutoff_i,
        cutoff_ii=cutoff_ii,
        cutoff_c=cutoff_c,
        cov_mode=cov_mode,
        max_per_cluster=max_per_cluster,
        touch_evalue=touch_evalue,
    )


def test_initial_cluster_counter_is_respected_for_new_cluster_ids():
    # 2026-10-02 第8步(--append): initial_cluster_counter 传入非0值时,
    # 本次运行产生的新 cluster_id 应从该值接续, 不从0重新开始。
    g1 = _gene("g1", 1)
    genome, intra = _make_genome_and_intra("g1", [g1], [g1.gene_id])

    result = _stage(MockAligner()).run(
        {"g1": genome}, {"g1": intra}, initial_cluster_counter=7,
    )

    assert result.n_clusters_final == 1
    assert result.cluster_counter >= 7
    [cluster] = [c for c in result.global_lib.clusters if c.active]
    assert cluster.cluster_id == "global|cluster_7"


def test_initial_global_lib_is_extended_not_replaced():
    # 2026-10-02 第8步(--append): 传入历史 initial_global_lib 时, 新基因组
    # 应在该历史库基础上继续并入, 而不是从空库重新开始(历史簇应保留)。
    from linc.core.entities import GlobalRepresentativeLibrary

    historical_lib = GlobalRepresentativeLibrary()
    historical_rep = _gene("old_genome", 1, length=999)
    historical_lib.add_new_cluster("global|cluster_0", historical_rep)

    new_gene = _gene("g1", 1, length=10)
    genome, intra = _make_genome_and_intra("g1", [new_gene], [new_gene.gene_id])

    result = _stage(MockAligner()).run(
        {"g1": genome}, {"g1": intra},
        initial_global_lib=historical_lib,
        initial_cluster_counter=1,
    )

    # 历史簇应仍然存在(未被从空库重新开始覆盖), 外加新基因组产生的簇。
    active_cluster_ids = {c.cluster_id for c in result.global_lib.clusters if c.active}
    assert "global|cluster_0" in active_cluster_ids
    assert result.n_clusters_final == 2


def test_cleans_up_search_intermediate_by_default(monkeypatch):
    """数据管理优化第3项(2026-10-03): search() 的中间文件(search_join/)
    结果已转存进 IncrementalResult, 默认(cleanup_intermediate=True)应在
    消费完结果后立即调用 `cleanup_search` 清理, index/tag 需与 search()
    一致。第一个基因组处理时全局库为空、不会发起 search, 第二个基因组
    才会真正调用 search(复用 global_lib 中第一个基因组的代表作为 query)。
    """
    a = _gene("g1", 1, length=100)
    b = _gene("g2", 1, length=10)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])

    aligner = MockAligner()
    cleanup_calls = []
    original_cleanup = aligner.cleanup_search

    def _recording_cleanup(index, tag=None):
        cleanup_calls.append((index, tag))
        return original_cleanup(index, tag=tag)

    monkeypatch.setattr(aligner, "cleanup_search", _recording_cleanup)

    _stage(aligner, cutoff_i=0.5).run(
        {"g1": genome_a, "g2": genome_b}, {"g1": intra_a, "g2": intra_b},
    )

    assert cleanup_calls == [(intra_b.genome_index, "join")]


def test_keep_intermediate_skips_cleanup(monkeypatch):
    """`cleanup_intermediate=False` 时不应调用 `cleanup_search`(保留中间
    文件以便调试, 2026-10-03 与用户确认)。
    """
    a = _gene("g1", 1, length=100)
    b = _gene("g2", 1, length=10)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])

    aligner = MockAligner()
    cleanup_calls = {"n": 0}
    original_cleanup = aligner.cleanup_search

    def _counting_cleanup(index, tag=None):
        cleanup_calls["n"] += 1
        return original_cleanup(index, tag=tag)

    monkeypatch.setattr(aligner, "cleanup_search", _counting_cleanup)

    stage = GlobalJoinStage(
        aligner=aligner, cutoff_i=0.5, cleanup_intermediate=False
    )
    stage.run({"g1": genome_a, "g2": genome_b}, {"g1": intra_a, "g2": intra_b})

    assert cleanup_calls["n"] == 0


def test_single_genome_single_representative_creates_one_cluster():
    g1 = _gene("g1", 1)
    genome, intra = _make_genome_and_intra("g1", [g1], [g1.gene_id])

    result = _stage(MockAligner()).run({"g1": genome}, {"g1": intra})

    assert result.n_clusters_final == 1
    assert result.n_representatives_final == 1
    assert result.processing_order == ["g1"]


def test_second_genome_zero_hit_creates_new_cluster():
    # 两个基因组各有1个局部代表, 彼此长度差异很大(低相似度), 应各自独立成簇。
    a = _gene("g1", 1, length=100)
    b = _gene("g2", 1, length=10)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])

    aligner = MockAligner()
    result = _stage(aligner, cutoff_i=0.5).run(
        {"g1": genome_a, "g2": genome_b}, {"g1": intra_a, "g2": intra_b}
    )

    assert result.n_clusters_final == 2
    assert result.n_representatives_final == 2


def test_second_genome_cutoff_ii_hit_reuses_existing_representative():
    # g1 的代表与 g2 的代表长度相同(MockAligner 默认规则 identity=1.0,
    # 必然 >= cutoff_ii), 应沿用已有代表, 不新增。
    a = _gene("g1", 1, length=100)
    b = _gene("g2", 1, length=100)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])

    aligner = MockAligner()
    result = _stage(aligner, cutoff_i=0.5, cutoff_ii=0.8).run(
        {"g1": genome_a, "g2": genome_b}, {"g1": intra_a, "g2": intra_b}
    )

    assert result.n_clusters_final == 1
    assert result.n_representatives_final == 1   # g2 的代表未新增, 沿用 g1 的代表
    # SSN 边应记录这次高置信命中
    assert len(result.ssn_edges) == 1
    # 边应携带 source_length/target_length(为未来 MCL 权重归一化预留,
    # 2026-10-02 与用户确认), MockAligner 按 Gene.length 填充, 两基因长度均为100。
    [edge] = result.ssn_edges
    assert edge.source_length == 100
    assert edge.target_length == 100


def test_second_genome_cutoff_i_only_hit_adds_new_representative():
    # 构造一对长度相似但不到 1.0 的基因(identity=0.7, 介于 cutoff_i=0.5 与
    # cutoff_ii=0.9 之间), 应增补为新代表而不是创建新簇或直接复用。
    a = _gene("g1", 1, length=100)
    b = _gene("g2", 1, length=100)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])

    aligner = MockAligner()
    aligner.set_result(a.gene_id, b.gene_id, _result(a, b, 0.7))

    result = _stage(aligner, cutoff_i=0.5, cutoff_ii=0.9).run(
        {"g1": genome_a, "g2": genome_b}, {"g1": intra_a, "g2": intra_b}
    )

    assert result.n_clusters_final == 1
    assert result.n_representatives_final == 2   # b 增补为新代表, 共享同一个簇
    cluster = result.global_lib.clusters[0]
    assert cluster.representative_ids == {a.gene_id, b.gene_id}


def test_third_genome_hits_two_clusters_merges_without_new_representative():
    # 先用两个彼此不相似的基因组各自建立独立簇 A、B;
    # 第三个基因组的代表与 A、B 的代表都达到 cutoff_ii, 触发合并且不新增代表。
    a = _gene("g1", 1, length=10)
    b = _gene("g2", 1, length=1000)
    c = _gene("g3", 1, length=50)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])
    genome_c, intra_c = _make_genome_and_intra("g3", [c], [c.gene_id])

    aligner = MockAligner()
    # a/b 长度差异极大, 默认规则下彼此 identity 很低, 不会先合并
    aligner.set_result(a.gene_id, c.gene_id, _result(a, c, 0.95))
    aligner.set_result(b.gene_id, c.gene_id, _result(b, c, 0.95))

    result = _stage(aligner, cutoff_i=0.5, cutoff_ii=0.8).run(
        {"g1": genome_a, "g2": genome_b, "g3": genome_c},
        {"g1": intra_a, "g2": intra_b, "g3": intra_c},
    )

    # a, b 各自独立成簇(n_clusters_final 中间态为2), 第三个基因组触发合并,
    # 最终只剩 1 个存活簇, c 本身不新增代表(因为命中达到 cutoff_ii)。
    assert result.n_clusters_final == 1
    assert result.n_representatives_final == 2   # 仍是 a、b 两个代表, c 未新增
    active_clusters = [cl for cl in result.global_lib.clusters if cl.active]
    assert len(active_clusters) == 1
    assert active_clusters[0].representative_ids == {a.gene_id, b.gene_id}


def test_third_genome_hits_two_clusters_merges_with_new_representative():
    # 与上一个用例类似, 但第三个基因组的代表与 A、B 都仅达到 cutoff_i(较低相似度),
    # 应触发合并且该代表本身也新增为合并后簇的代表("连接远缘家族的新变体")。
    a = _gene("g1", 1, length=10)
    b = _gene("g2", 1, length=1000)
    c = _gene("g3", 1, length=50)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])
    genome_c, intra_c = _make_genome_and_intra("g3", [c], [c.gene_id])

    aligner = MockAligner()
    aligner.set_result(a.gene_id, c.gene_id, _result(a, c, 0.6))
    aligner.set_result(b.gene_id, c.gene_id, _result(b, c, 0.6))

    result = _stage(aligner, cutoff_i=0.5, cutoff_ii=0.8).run(
        {"g1": genome_a, "g2": genome_b, "g3": genome_c},
        {"g1": intra_a, "g2": intra_b, "g3": intra_c},
    )

    assert result.n_clusters_final == 1
    assert result.n_representatives_final == 3   # a、b、c 均为代表
    active_clusters = [cl for cl in result.global_lib.clusters if cl.active]
    assert len(active_clusters) == 1
    assert active_clusters[0].representative_ids == {a.gene_id, b.gene_id, c.gene_id}


def test_order_d_sorts_by_gene_count_desc_then_rep_count_desc_then_id_asc():
    g_small, intra_small = _make_genome_and_intra(
        "genome_b", [_gene("genome_b", 1)], ["genome_b|1"]
    )
    g_big, intra_big = _make_genome_and_intra(
        "genome_a",
        [_gene("genome_a", i) for i in range(3)],
        [f"genome_a|{i}" for i in range(3)],
    )
    intra_small.n_genes = 1
    intra_big.n_genes = 3

    result = _stage(MockAligner()).run(
        {"genome_a": g_big, "genome_b": g_small},
        {"genome_a": intra_big, "genome_b": intra_small},
    )

    assert result.processing_order == ["genome_a", "genome_b"]


def test_cutoff_c_filters_out_low_coverage_hit_even_if_identity_passes():
    # identity 达到 cutoff_i, 但覆盖度很低(0.1), 设置 cutoff_c=0.5 应过滤掉该命中,
    # 两个基因组各自独立成簇。
    a = _gene("g1", 1, length=100)
    b = _gene("g2", 1, length=100)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])

    aligner = MockAligner()
    aligner.set_result(
        a.gene_id, b.gene_id,
        _result(a, b, identity=0.95, query_cover=0.1, subject_cover=0.1),
    )

    result = _stage(aligner, cutoff_i=0.5, cutoff_ii=0.8, cutoff_c=0.5).run(
        {"g1": genome_a, "g2": genome_b}, {"g1": intra_a, "g2": intra_b}
    )

    assert result.n_clusters_final == 2
    assert result.n_representatives_final == 2


def test_cov_mode_query_only_ignores_low_subject_coverage():
    # cov_mode=2(仅看 query_cover): subject_cover 很低但 query_cover 达标,
    # 在该模式下"是否纳入SSN"的 inclusion 过滤应通过(与 cov_mode=0 双向
    # 模式的结果不同), 故两个基因仍会落入同一个簇(n_clusters_final==1)。
    #
    # 但第二轮漏洞3修复(2026-10-03, 与用户确认)后, "高置信(cutoff_ii)"
    # 判定统一改用双向 coverage(= min(query_cover, subject_cover))与
    # cutoff_c 比较, 不受 cov_mode 影响——这里双向 coverage=0.1<cutoff_c=0.5,
    # 不达标, 不会被误判为高置信, 所以应新增代表而非只做 member
    # (n_representatives_final==2)。
    a = _gene("g1", 1, length=100)
    b = _gene("g2", 1, length=100)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])

    aligner = MockAligner()
    aligner.set_result(
        a.gene_id, b.gene_id,
        _result(a, b, identity=0.95, query_cover=0.9, subject_cover=0.1),
    )

    result = _stage(
        aligner, cutoff_i=0.5, cutoff_ii=0.8, cutoff_c=0.5, cov_mode=2
    ).run({"g1": genome_a, "g2": genome_b}, {"g1": intra_a, "g2": intra_b})

    assert result.n_clusters_final == 1
    assert result.n_representatives_final == 2


def test_cutoff_ii_high_identity_but_low_coverage_adds_new_representative():
    # 第二轮漏洞3回归测试: 路径②(单簇命中) identity >= cutoff_ii, 但双向
    # coverage(query/subject 中较小者)低于 cutoff_c(碎片匹配)。
    #
    # 用 cov_mode=2(仅看 query_cover)让该命中先通过"是否纳入SSN"的
    # inclusion 过滤(query_cover=0.9 达标), 但 subject_cover 很低(0.1),
    # 双向 coverage(= min(0.9, 0.1) = 0.1)不达标。修复前"高置信"判定只看
    # identity, 会被误判为高置信直接 add_member_only; 修复后高置信判定统一
    # 用双向 coverage 与 cutoff_c 比较(与 inclusion 过滤的 cov_mode 相对
    # 独立, 更保守), 应改为新增代表, 不丢失家族内亚型分辨率。
    a = _gene("g1", 1, length=100)
    b = _gene("g2", 1, length=100)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])

    aligner = MockAligner()
    aligner.set_result(
        a.gene_id, b.gene_id,
        _result(a, b, identity=0.95, query_cover=0.9, subject_cover=0.1),
    )

    result = _stage(
        aligner, cutoff_i=0.5, cutoff_ii=0.8, cutoff_c=0.5, cov_mode=2
    ).run({"g1": genome_a, "g2": genome_b}, {"g1": intra_a, "g2": intra_b})

    assert result.n_clusters_final == 1
    assert result.n_representatives_final == 2   # b 应增补为新代表, 而非只做 member
    cluster = result.global_lib.clusters[0]
    assert cluster.representative_ids == {a.gene_id, b.gene_id}


def test_cutoff_ii_merge_path_high_identity_but_low_coverage_adds_representative():
    # 第二轮漏洞3回归测试: 路径③(跨簇命中, >=2簇) 同样的高置信判定修复,
    # 用法同上(cov_mode=2 让命中先通过 inclusion 过滤, 但双向 coverage 不
    # 达标), 所有命中都不应被判定为 any_cutoff_ii, 新代表应被新增
    # (merge_clusters 带 new_representative)。
    a = _gene("g1", 1, length=10)
    b = _gene("g2", 1, length=1000)
    c = _gene("g3", 1, length=50)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])
    genome_c, intra_c = _make_genome_and_intra("g3", [c], [c.gene_id])

    aligner = MockAligner()
    aligner.set_result(
        a.gene_id, c.gene_id,
        _result(a, c, identity=0.95, query_cover=0.9, subject_cover=0.1),
    )
    aligner.set_result(
        b.gene_id, c.gene_id,
        _result(b, c, identity=0.95, query_cover=0.9, subject_cover=0.1),
    )

    result = _stage(
        aligner, cutoff_i=0.5, cutoff_ii=0.8, cutoff_c=0.5, cov_mode=2
    ).run(
        {"g1": genome_a, "g2": genome_b, "g3": genome_c},
        {"g1": intra_a, "g2": intra_b, "g3": intra_c},
    )

    assert result.n_clusters_final == 1
    assert result.n_representatives_final == 3   # a、b、c 均为代表(c 未被跳过)
    active_clusters = [cl for cl in result.global_lib.clusters if cl.active]
    assert len(active_clusters) == 1
    assert active_clusters[0].representative_ids == {a.gene_id, b.gene_id, c.gene_id}


def test_max_per_cluster_caps_single_cluster_path_representatives():
    # 验证路径②(单簇命中、仅达 cutoff_i)增补代表时不超过 max_per_cluster。
    base = _gene("g1", 1, length=100)
    genome_base, intra_base = _make_genome_and_intra("g1", [base], [base.gene_id])

    aligner = MockAligner()
    genomes = {"g1": genome_base}
    intras = {"g1": intra_base}
    # 依次加入3个与 base 恰好 cutoff_i<=identity<cutoff_ii 的代表, max_per_cluster=2
    for i in range(2, 5):
        gi = _gene(f"g{i}", 1, length=100)
        aligner.set_result(base.gene_id, gi.gene_id, _result(base, gi, 0.6))
        genome_i, intra_i = _make_genome_and_intra(f"g{i}", [gi], [gi.gene_id])
        genomes[f"g{i}"] = genome_i
        intras[f"g{i}"] = intra_i

    result = _stage(aligner, cutoff_i=0.5, cutoff_ii=0.8, max_per_cluster=2).run(genomes, intras)

    assert result.n_clusters_final == 1
    cluster = result.global_lib.clusters[0]
    assert len(cluster.representative_ids) <= 2


# ---------------------------------------------------------------------------
# 接触判据/纳入判据解耦设计方案(2026-10-04): touching_map 相关单测。
# ---------------------------------------------------------------------------

def test_touching_map_records_weak_hit_even_when_inclusion_cutoff_fails():
    """核心场景: 弱 hit(evalue 达标、identity 不达标纳入判据)不应合并,
    但应被记录进 touching_map(接触判据与纳入判据解耦的关键验证)。
    """
    a = _gene("g1", 1, length=100)
    b = _gene("g2", 1, length=100)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])

    aligner = MockAligner()
    weak_hit = AlignmentResult(
        query_id=a.gene_id, subject_id=b.gene_id,
        identity=0.1, query_cover=0.1, subject_cover=0.1,
        evalue=1e-5, bitscore=10.0, alignment_length=10,
    )
    aligner.set_result(a.gene_id, b.gene_id, weak_hit)

    result = _stage(aligner, cutoff_i=0.5, cutoff_c=0.0).run(
        {"g1": genome_a, "g2": genome_b}, {"g1": intra_a, "g2": intra_b}
    )

    # 纳入判据(identity>=0.5)不达标, 不应合并 -> 各自独立成簇。
    assert result.n_clusters_final == 2
    [a_cluster_id] = [
        c.cluster_id for c in result.global_lib.clusters
        if c.active and a.gene_id in c.representative_ids
    ]
    # 但接触判据(touch_evalue=None, 最宽松)应记录 g2 接触了 a 所在的簇。
    assert a_cluster_id in result.touching_map.get("g2", set())


def test_touch_evalue_filters_out_weak_hits_from_touching_map():
    """收紧 touch_evalue 后, evalue 超过阈值的 hit 不应被记录为接触。"""
    a = _gene("g1", 1, length=100)
    b = _gene("g2", 1, length=100)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])

    aligner = MockAligner()
    weak_hit = AlignmentResult(
        query_id=a.gene_id, subject_id=b.gene_id,
        identity=0.1, query_cover=0.1, subject_cover=0.1,
        evalue=1e-5, bitscore=10.0, alignment_length=10,
    )
    aligner.set_result(a.gene_id, b.gene_id, weak_hit)

    result = _stage(aligner, cutoff_i=0.5, cutoff_c=0.0, touch_evalue=1e-10).run(
        {"g1": genome_a, "g2": genome_b}, {"g1": intra_a, "g2": intra_b}
    )

    assert result.touching_map.get("g2", set()) == set()


def test_touch_evalue_none_also_marks_strong_hits_as_touching():
    """接触标记与合并逻辑正交: 一条足够强、正常合并的 hit, 同样应被记录进
    touching_map(接触判据只是"打标记", 不影响、也不被合并逻辑排斥)。
    """
    a = _gene("g1", 1, length=100)
    b = _gene("g2", 1, length=100)
    genome_a, intra_a = _make_genome_and_intra("g1", [a], [a.gene_id])
    genome_b, intra_b = _make_genome_and_intra("g2", [b], [b.gene_id])

    aligner = MockAligner()
    aligner.set_result(a.gene_id, b.gene_id, _result(a, b, identity=0.9))

    result = _stage(aligner, cutoff_i=0.5, cutoff_ii=0.8, cutoff_c=0.0).run(
        {"g1": genome_a, "g2": genome_b}, {"g1": intra_a, "g2": intra_b}
    )

    assert result.n_clusters_final == 1   # 正常合并, cutoff_ii 逻辑不受影响
    [cluster_id] = [c.cluster_id for c in result.global_lib.clusters if c.active]
    assert cluster_id in result.touching_map.get("g2", set())
