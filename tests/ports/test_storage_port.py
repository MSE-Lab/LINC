"""RepositoryPort 结构化子类型测试。

2026-10-02 第5步(pipeline打通)与用户确认重新设计: 接口按 engine/ 五个
Stage 的实际产出类型逐一对应提供 save/load 方法(取代早期版本里粗粒度的
checkpoint/global_library/ssn_edges)。
"""
from typing import Dict, List, Optional

from linc.core.entities import CompleteSSN, Genome, GenomeIntraResult, IncrementalResult, OrthologFamily
from linc.ports.storage_port import RepositoryPort


class FakeRepository:
    """最小化假实现, 未显式继承 RepositoryPort。"""

    def __init__(self):
        self._genomes: Dict[str, Genome] = {}
        self._intra_results: Dict[str, GenomeIntraResult] = {}
        self._join_result: Optional[IncrementalResult] = None
        self._complete_ssn: Optional[CompleteSSN] = None
        self._families: List[OrthologFamily] = []

    def save_genomes(self, genomes):
        self._genomes = genomes

    def load_genomes(self):
        return self._genomes

    def save_intra_results(self, results):
        self._intra_results = results

    def load_intra_results(self):
        return self._intra_results

    def save_join_result(self, result):
        self._join_result = result

    def load_join_result(self):
        return self._join_result

    def save_complete_ssn(self, ssn):
        self._complete_ssn = ssn

    def load_complete_ssn(self):
        return self._complete_ssn

    def save_families(self, families):
        self._families = families

    def load_families(self):
        return self._families


def test_fake_repository_satisfies_protocol_via_isinstance():
    assert isinstance(FakeRepository(), RepositoryPort)


def test_genomes_roundtrip():
    repo: RepositoryPort = FakeRepository()
    assert repo.load_genomes() == {}
    genome = Genome(genome_id="g1", file_path="g1.faa")
    repo.save_genomes({"g1": genome})
    assert repo.load_genomes()["g1"] is genome


def test_intra_results_roundtrip():
    repo: RepositoryPort = FakeRepository()
    assert repo.load_intra_results() == {}
    result = GenomeIntraResult(genome_id="g1")
    repo.save_intra_results({"g1": result})
    assert repo.load_intra_results()["g1"] is result


def test_join_result_roundtrip():
    repo: RepositoryPort = FakeRepository()
    assert repo.load_join_result() is None
    result = IncrementalResult()
    repo.save_join_result(result)
    assert repo.load_join_result() is result


def test_complete_ssn_roundtrip():
    repo: RepositoryPort = FakeRepository()
    assert repo.load_complete_ssn() is None
    ssn = CompleteSSN()
    repo.save_complete_ssn(ssn)
    assert repo.load_complete_ssn() is ssn


def test_families_roundtrip():
    repo: RepositoryPort = FakeRepository()
    assert repo.load_families() == []
    family = OrthologFamily(family_id="OG0000000", component_id=0, member_genes={"g1|1"})
    repo.save_families([family])
    assert repo.load_families() == [family]


def test_incomplete_implementation_fails_isinstance_check():
    class IncompleteRepository:
        def save_genomes(self, genomes):
            pass

    assert not isinstance(IncompleteRepository(), RepositoryPort)
