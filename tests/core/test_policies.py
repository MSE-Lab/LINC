"""`linc/core/policies.py` 的单元测试。

Protocol 是结构化子类型,这里验证:
1) 符合签名的类即被视为满足接口(isinstance 对 runtime_checkable Protocol 的行为除外,
   本测试只做"调用层面"的行为验证,不依赖 @runtime_checkable)。
2) `ThresholdDecision` 作为普通 dataclass 可被正常构造和使用。
"""
from typing import List, Sequence

from linc.core.entities import AlignmentResult, Gene, Genome
from linc.core.policies import ThresholdDecision, ThresholdPolicy


class _FakeOrderingPolicy:
    """未显式继承 OrderingPolicy,仅靠结构化子类型满足接口。"""

    def sort(self, genomes: Sequence[Genome]) -> List[Genome]:
        return sorted(genomes, key=lambda g: g.genome_id)


class _FakeThresholdPolicy:
    cutoff_i: float = 0.5
    cutoff_ii: float = 0.8

    def decide(self, new_rep: Gene, hits: List[AlignmentResult], lib: object) -> ThresholdDecision:
        if not hits:
            return ThresholdDecision(action="create_new", target_clusters=[])
        return ThresholdDecision(action="attach_no_rep", target_clusters=["c1"])


def test_ordering_policy_structural_typing():
    policy = _FakeOrderingPolicy()
    genomes = [Genome(genome_id="b", file_path="b.fasta"), Genome(genome_id="a", file_path="a.fasta")]
    ordered = policy.sort(genomes)
    assert [g.genome_id for g in ordered] == ["a", "b"]


def test_threshold_policy_decide_no_hits_creates_new_cluster():
    policy: ThresholdPolicy = _FakeThresholdPolicy()
    gene = Gene(gene_id="g3|1", genome_id="g3", length=200)
    decision = policy.decide(gene, [], lib=None)
    assert decision.action == "create_new"
    assert decision.target_clusters == []


def test_threshold_policy_decide_with_hits_attaches_to_cluster():
    policy: ThresholdPolicy = _FakeThresholdPolicy()
    gene = Gene(gene_id="g3|1", genome_id="g3", length=200)
    hit = AlignmentResult(
        query_id="g3|1",
        subject_id="g1|1",
        identity=0.6,
        query_cover=0.9,
        subject_cover=0.9,
        evalue=1e-30,
        bitscore=100.0,
        alignment_length=180,
    )
    decision = policy.decide(gene, [hit], lib=None)
    assert decision.action == "attach_no_rep"
    assert decision.target_clusters == ["c1"]


def test_threshold_decision_default_target_clusters_is_independent():
    d1 = ThresholdDecision(action="create_new")
    d2 = ThresholdDecision(action="create_new")
    d1.target_clusters.append("c1")
    assert d2.target_clusters == []
