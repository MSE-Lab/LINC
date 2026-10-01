"""RepositoryPort 结构化子类型测试。"""
from typing import Optional

from linc.ports.storage_port import RepositoryPort


class FakeRepository:
    """最小化假实现, 未显式继承 RepositoryPort。"""

    def __init__(self):
        self._checkpoints: dict = {}
        self._library: dict = {}
        self._edges: set = set()

    def save_checkpoint(self, stage: str, state: dict) -> None:
        self._checkpoints[stage] = state

    def load_checkpoint(self, stage: str) -> Optional[dict]:
        return self._checkpoints.get(stage)

    def save_global_library(self, clusters: dict) -> None:
        self._library = clusters

    def load_global_library(self) -> dict:
        return self._library

    def save_ssn_edges(self, edges: set) -> None:
        self._edges = edges


def test_fake_repository_satisfies_protocol_via_isinstance():
    assert isinstance(FakeRepository(), RepositoryPort)


def test_checkpoint_roundtrip():
    repo: RepositoryPort = FakeRepository()
    assert repo.load_checkpoint("step2") is None
    repo.save_checkpoint("step2", {"done": True})
    assert repo.load_checkpoint("step2") == {"done": True}


def test_global_library_roundtrip():
    repo: RepositoryPort = FakeRepository()
    assert repo.load_global_library() == {}
    repo.save_global_library({"c1": object()})
    assert "c1" in repo.load_global_library()


def test_save_ssn_edges():
    repo: RepositoryPort = FakeRepository()
    repo.save_ssn_edges({("a", "b")})
    assert ("a", "b") in repo._edges


def test_incomplete_implementation_fails_isinstance_check():
    class IncompleteRepository:
        def save_checkpoint(self, stage, state):
            pass

    assert not isinstance(IncompleteRepository(), RepositoryPort)
