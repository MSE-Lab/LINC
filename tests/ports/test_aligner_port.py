"""AlignerPort 结构化子类型测试。

验证: 任何实现了对应方法签名的类, 均可被当作 AlignerPort 使用(鸭子类型),
不需要显式继承, 符合 typing.Protocol 的设计初衷。
"""
from typing import List, Sequence

from linc.core.entities import AlignmentResult, Gene
from linc.ports.aligner_port import AlignerPort


class FakeAligner:
    """最小化假实现, 未显式继承 AlignerPort。"""

    def align_all_vs_all(self, genes: Sequence[Gene]) -> List[AlignmentResult]:
        return [
            AlignmentResult(
                query_id=genes[i].gene_id,
                subject_id=genes[j].gene_id,
                identity=1.0,
                query_cover=1.0,
                subject_cover=1.0,
                evalue=0.0,
                bitscore=100.0,
                alignment_length=genes[i].length,
            )
            for i in range(len(genes))
            for j in range(i + 1, len(genes))
        ]

    def align_query_vs_subject(
        self, query: Sequence[Gene], subject: Sequence[Gene]
    ) -> List[AlignmentResult]:
        return [
            AlignmentResult(
                query_id=q.gene_id,
                subject_id=s.gene_id,
                identity=0.5,
                query_cover=0.5,
                subject_cover=0.5,
                evalue=1e-5,
                bitscore=50.0,
                alignment_length=min(q.length, s.length),
            )
            for q in query
            for s in subject
        ]

    def build_index(self, genes: Sequence[Gene], index_id=None):
        return list(genes)

    def search(self, query: Sequence[Gene], index) -> List[AlignmentResult]:
        return self.align_query_vs_subject(query, index)

    def index_exists(self, index) -> bool:
        return bool(index)


def _make_genes() -> List[Gene]:
    return [
        Gene(gene_id="g1|1", genome_id="g1", length=100),
        Gene(gene_id="g1|2", genome_id="g1", length=200),
    ]


def test_fake_aligner_satisfies_protocol_via_isinstance():
    aligner = FakeAligner()
    assert isinstance(aligner, AlignerPort)


def test_fake_aligner_align_all_vs_all():
    aligner: AlignerPort = FakeAligner()
    genes = _make_genes()
    results = aligner.align_all_vs_all(genes)
    assert len(results) == 1
    assert results[0].query_id == "g1|1"
    assert results[0].subject_id == "g1|2"


def test_fake_aligner_align_query_vs_subject():
    aligner: AlignerPort = FakeAligner()
    query = [Gene(gene_id="q1", genome_id="q", length=50)]
    subject = _make_genes()
    results = aligner.align_query_vs_subject(query, subject)
    assert len(results) == 2
    assert all(r.query_id == "q1" for r in results)


def test_incomplete_implementation_fails_isinstance_check():
    class IncompleteAligner:
        def align_all_vs_all(self, genes):
            return []

    assert not isinstance(IncompleteAligner(), AlignerPort)
