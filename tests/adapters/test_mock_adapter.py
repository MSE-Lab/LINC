"""MockAligner 单元测试。"""
from linc.adapters.aligners.mock_adapter import MockAligner
from linc.core.entities import AlignmentResult, Gene
from linc.ports.aligner_port import AlignerPort


def _gene(gene_id: str, length: int) -> Gene:
    return Gene(gene_id=gene_id, genome_id="g1", length=length)


def test_mock_aligner_satisfies_aligner_port():
    assert isinstance(MockAligner(), AlignerPort)


def test_default_result_identical_length_gives_identity_one():
    aligner = MockAligner()
    a, b = _gene("g1|1", 100), _gene("g1|2", 100)
    [result] = aligner.align_all_vs_all([a, b])
    assert result.identity == 1.0
    assert result.coverage == 1.0
    assert result.evalue == 0.0


def test_default_result_different_length_gives_lower_identity():
    aligner = MockAligner()
    a, b = _gene("g1|1", 100), _gene("g1|2", 50)
    [result] = aligner.align_all_vs_all([a, b])
    assert result.identity == 0.5
    assert result.alignment_length == 50


def test_align_all_vs_all_no_self_and_no_duplicate_pairs():
    aligner = MockAligner()
    genes = [_gene(f"g1|{i}", 100) for i in range(4)]
    results = aligner.align_all_vs_all(genes)
    # C(4,2) = 6 对, 不含自身, 不重复
    assert len(results) == 6


def test_preset_result_takes_priority_over_default():
    preset = AlignmentResult(
        query_id="g1|1",
        subject_id="g1|2",
        identity=0.99,
        query_cover=0.99,
        subject_cover=0.99,
        evalue=1e-50,
        bitscore=999.0,
        alignment_length=100,
    )
    aligner = MockAligner(preset_results={("g1|1", "g1|2"): preset})
    a, b = _gene("g1|1", 100), _gene("g1|2", 10)  # 长度差异很大, 若走默认规则 identity 应很低
    [result] = aligner.align_all_vs_all([a, b])
    assert result.identity == 0.99


def test_set_result_is_bidirectional():
    preset = AlignmentResult(
        query_id="qx",
        subject_id="sx",
        identity=0.77,
        query_cover=0.77,
        subject_cover=0.77,
        evalue=0.0,
        bitscore=10.0,
        alignment_length=5,
    )
    aligner = MockAligner()
    aligner.set_result("qx", "sx", preset)
    query = [_gene("sx", 5)]
    subject = [_gene("qx", 999)]  # 调换 query/subject 顺序, 应仍命中同一条预设
    [result] = aligner.align_query_vs_subject(query, subject)
    assert result.identity == 0.77


def test_align_query_vs_subject_full_cartesian_product():
    aligner = MockAligner()
    query = [_gene("q1", 10), _gene("q2", 20)]
    subject = [_gene("s1", 10), _gene("s2", 20), _gene("s3", 30)]
    results = aligner.align_query_vs_subject(query, subject)
    assert len(results) == 6
