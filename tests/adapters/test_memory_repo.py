"""MemoryRepository 单元测试。

2026-10-02 第5步(pipeline打通)与用户确认重新设计接口: 按 engine/ 五个
Stage 的实际产出类型逐一对应提供 save/load 方法(取代早期版本里粗粒度的
checkpoint/global_library/ssn_edges)。
"""
from linc.adapters.storage.memory_repo import MemoryRepository
from linc.core.entities import CompleteSSN, Genome, GenomeIntraResult, IncrementalResult, OrthologFamily, SSNEdge
from linc.ports.storage_port import RepositoryPort


def test_memory_repo_satisfies_repository_port():
    assert isinstance(MemoryRepository(), RepositoryPort)


def test_genomes_roundtrip_and_missing_returns_empty_dict():
    repo = MemoryRepository()
    assert repo.load_genomes() == {}
    genome = Genome(genome_id="g1", file_path="g1.faa")
    repo.save_genomes({"g1": genome})
    loaded = repo.load_genomes()
    assert loaded["g1"] == genome


def test_genomes_save_is_deep_copied_not_aliased():
    repo = MemoryRepository()
    genome = Genome(genome_id="g1", file_path="g1.faa")
    genomes = {"g1": genome}
    repo.save_genomes(genomes)
    genomes["g1"].file_path = "mutated.faa"  # 修改原对象, 不应影响已保存的快照
    assert repo.load_genomes()["g1"].file_path == "g1.faa"


def test_genomes_load_returns_independent_copy():
    repo = MemoryRepository()
    repo.save_genomes({"g1": Genome(genome_id="g1", file_path="g1.faa")})
    loaded_once = repo.load_genomes()
    loaded_once["g1"].file_path = "mutated.faa"
    loaded_twice = repo.load_genomes()
    assert loaded_twice["g1"].file_path == "g1.faa"


def test_intra_results_roundtrip_and_missing_returns_empty_dict():
    repo = MemoryRepository()
    assert repo.load_intra_results() == {}
    result = GenomeIntraResult(genome_id="g1", n_genes=3)
    repo.save_intra_results({"g1": result})
    assert repo.load_intra_results()["g1"].n_genes == 3


def test_join_result_roundtrip_and_missing_returns_none():
    repo = MemoryRepository()
    assert repo.load_join_result() is None
    result = IncrementalResult(processing_order=["g1", "g2"], n_clusters_final=2)
    repo.save_join_result(result)
    loaded = repo.load_join_result()
    assert loaded.processing_order == ["g1", "g2"]
    assert loaded.n_clusters_final == 2


def test_join_result_save_overwrites_previous():
    repo = MemoryRepository()
    repo.save_join_result(IncrementalResult(n_clusters_final=1))
    repo.save_join_result(IncrementalResult(n_clusters_final=2))
    assert repo.load_join_result().n_clusters_final == 2


def test_complete_ssn_roundtrip_and_missing_returns_none():
    repo = MemoryRepository()
    assert repo.load_complete_ssn() is None
    edge = SSNEdge(source="g1|1", target="g2|1", identity=0.8, coverage=0.5)
    ssn = CompleteSSN(edges={edge}, edge_completion_mode="full")
    repo.save_complete_ssn(ssn)
    loaded = repo.load_complete_ssn()
    assert loaded.edge_completion_mode == "full"
    assert len(loaded.edges) == 1


def test_families_roundtrip_and_missing_returns_empty_list():
    repo = MemoryRepository()
    assert repo.load_families() == []
    family = OrthologFamily(family_id="OG0000000", component_id=0, member_genes={"g1|1", "g2|1"})
    repo.save_families([family])
    loaded = repo.load_families()
    assert len(loaded) == 1
    assert loaded[0].family_id == "OG0000000"
