"""FileRepository 单元测试(2026-10-02, 第8步: 持久化存储 + resume/append)。"""
import pytest

from linc.adapters.storage.file_repo import FileRepository, StateCorruptedError
from linc.core.entities import CompleteSSN, Genome, GenomeIntraResult, IncrementalResult, OrthologFamily
from linc.ports.storage_port import RepositoryPort


def test_file_repository_satisfies_protocol(tmp_path):
    repo = FileRepository(str(tmp_path / "state"))
    assert isinstance(repo, RepositoryPort)


def test_has_methods_false_before_save_true_after(tmp_path):
    repo = FileRepository(str(tmp_path / "state"))
    assert repo.has_genomes() is False
    repo.save_genomes({"g1": Genome(genome_id="g1", file_path="g1.faa")})
    assert repo.has_genomes() is True


def test_genomes_roundtrip_across_instances(tmp_path):
    state_dir = str(tmp_path / "state")
    repo1 = FileRepository(state_dir)
    genome = Genome(genome_id="g1", file_path="g1.faa", genes=[])
    repo1.save_genomes({"g1": genome})

    # 新建一个实例模拟跨进程重新加载
    repo2 = FileRepository(state_dir)
    loaded = repo2.load_genomes()
    assert loaded["g1"].genome_id == "g1"


def test_intra_results_roundtrip(tmp_path):
    repo = FileRepository(str(tmp_path / "state"))
    assert repo.load_intra_results() == {}
    result = GenomeIntraResult(genome_id="g1")
    repo.save_intra_results({"g1": result})
    assert repo.has_intra_results()
    assert repo.load_intra_results()["g1"].genome_id == "g1"


def test_join_result_roundtrip(tmp_path):
    repo = FileRepository(str(tmp_path / "state"))
    assert repo.load_join_result() is None
    result = IncrementalResult(cluster_counter=5)
    repo.save_join_result(result)
    assert repo.has_join_result()
    assert repo.load_join_result().cluster_counter == 5


def test_complete_ssn_roundtrip(tmp_path):
    repo = FileRepository(str(tmp_path / "state"))
    assert repo.load_complete_ssn() is None
    ssn = CompleteSSN()
    repo.save_complete_ssn(ssn)
    assert repo.has_complete_ssn()
    assert repo.load_complete_ssn() is not None


def test_families_roundtrip(tmp_path):
    repo = FileRepository(str(tmp_path / "state"))
    assert repo.load_families() == []
    family = OrthologFamily(family_id="OG0000000", component_id=0, member_genes={"g1|1"})
    repo.save_families([family])
    assert repo.has_families()
    assert repo.load_families()[0].family_id == "OG0000000"


def test_input_signature_roundtrip(tmp_path):
    repo = FileRepository(str(tmp_path / "state"))
    assert repo.load_input_signature() is None
    repo.save_input_signature(["a.faa", "b.faa"])
    assert repo.load_input_signature() == ["a.faa", "b.faa"]


def test_corrupted_state_file_raises_not_treated_as_missing(tmp_path):
    state_dir = tmp_path / "state"
    state_dir.mkdir()
    (state_dir / "genomes.pkl").write_bytes(b"not a valid pickle stream !!!")

    repo = FileRepository(str(state_dir))
    with pytest.raises(StateCorruptedError):
        repo.load_genomes()


def test_save_uses_atomic_replace_leaves_no_tmp_file(tmp_path):
    state_dir = tmp_path / "state"
    repo = FileRepository(str(state_dir))
    repo.save_genomes({"g1": Genome(genome_id="g1", file_path="g1.faa")})
    assert not (state_dir / "genomes.pkl.tmp").exists()
    assert (state_dir / "genomes.pkl").exists()
