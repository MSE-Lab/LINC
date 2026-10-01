"""MemoryRepository 单元测试。"""
from linc.adapters.storage.memory_repo import MemoryRepository
from linc.core.entities import Cluster
from linc.ports.storage_port import RepositoryPort


def test_memory_repo_satisfies_repository_port():
    assert isinstance(MemoryRepository(), RepositoryPort)


def test_checkpoint_roundtrip_and_missing_returns_none():
    repo = MemoryRepository()
    assert repo.load_checkpoint("stage1") is None
    repo.save_checkpoint("stage1", {"progress": 3})
    assert repo.load_checkpoint("stage1") == {"progress": 3}


def test_checkpoint_save_is_deep_copied_not_aliased():
    repo = MemoryRepository()
    state = {"progress": 1, "items": [1, 2, 3]}
    repo.save_checkpoint("stage1", state)
    state["items"].append(4)  # 修改原对象, 不应影响已保存的快照
    saved = repo.load_checkpoint("stage1")
    assert saved["items"] == [1, 2, 3]


def test_checkpoint_load_returns_independent_copy():
    repo = MemoryRepository()
    repo.save_checkpoint("stage1", {"items": [1, 2]})
    loaded_once = repo.load_checkpoint("stage1")
    loaded_once["items"].append(99)
    loaded_twice = repo.load_checkpoint("stage1")
    assert loaded_twice["items"] == [1, 2]


def test_multiple_stage_checkpoints_are_independent():
    repo = MemoryRepository()
    repo.save_checkpoint("stage1", {"a": 1})
    repo.save_checkpoint("stage2", {"b": 2})
    assert repo.load_checkpoint("stage1") == {"a": 1}
    assert repo.load_checkpoint("stage2") == {"b": 2}


def test_global_library_roundtrip():
    repo = MemoryRepository()
    assert repo.load_global_library() == {}
    cluster = Cluster(cluster_id="c1", representative_ids={"g1"}, member_ids={"g1", "g2"})
    repo.save_global_library({"c1": cluster})
    loaded = repo.load_global_library()
    assert loaded["c1"] is cluster


def test_global_library_save_overwrites_previous():
    repo = MemoryRepository()
    repo.save_global_library({"c1": Cluster(cluster_id="c1")})
    repo.save_global_library({"c2": Cluster(cluster_id="c2")})
    assert list(repo.load_global_library().keys()) == ["c2"]


def test_save_ssn_edges_stores_as_set():
    repo = MemoryRepository()
    repo.save_ssn_edges({("a", "b"), ("b", "c")})
    assert repo._ssn_edges == {("a", "b"), ("b", "c")}
