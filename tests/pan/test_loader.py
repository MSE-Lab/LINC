"""`linc/pan/loader.py` 的单元测试: 只读加载 + 错误处理 + ID 映射。"""
import gzip
import pickle

import pytest

from linc.adapters.storage.file_repo import FileRepository
from linc.core.entities import CompleteSSN, Gene, Genome, OrthologFamily
from linc.pan.loader import LoadedState, ReportInputError, load_state


def _make_genome(genome_id: str, source_name: str, n_genes: int = 2) -> Genome:
    genes = [
        Gene(
            gene_id=f"{genome_id}|g{i}",
            genome_id=genome_id,
            length=100,
            original_header=f"{source_name}|g{i}",
        )
        for i in range(n_genes)
    ]
    return Genome(
        genome_id=genome_id,
        file_path=f"/fake/{source_name}.faa",
        genes=genes,
        source_file_name=source_name,
    )


def _save_minimal_state(out_dir, genomes=None, families=None, ssn=None, config=None):
    repo = FileRepository(str(out_dir / ".linc_work" / "state"))
    repo.save_genomes(genomes if genomes is not None else {})
    repo.save_families(families if families is not None else [])
    if ssn is not None:
        repo.save_complete_ssn(ssn)
    if config is not None:
        repo.save_config(config)


def _write_wrong_schema_genomes(state_dir):
    """手动写入一个 schema_version 不匹配的 genomes.pkl, 用于触发
    SchemaVersionMismatchError。"""
    state_dir.mkdir(parents=True, exist_ok=True)
    path = state_dir / "genomes.pkl"
    with gzip.open(path, "wb") as fh:
        pickle.dump({"__schema_version__": 999, "payload": {}}, fh)


def test_missing_state_dir_raises(tmp_path):
    with pytest.raises(ReportInputError, match="未找到 state 目录"):
        load_state(str(tmp_path / "no_such_run"))


def test_missing_required_files_raises(tmp_path):
    out = tmp_path / "out"
    # 只有空的 state 目录, 没有任何必需文件。
    (out / ".linc_work" / "state").mkdir(parents=True)
    with pytest.raises(ReportInputError, match="缺少必需文件"):
        load_state(str(out))


def test_load_state_builds_id_mappings(tmp_path):
    out = tmp_path / "out"
    genome = _make_genome("G00001", "GCA_0001", n_genes=2)
    fam = OrthologFamily(
        family_id="OG0000000",
        component_id=0,
        member_genes={"G00001|g0", "G00001|g1"},
    )
    _save_minimal_state(
        out,
        genomes={"G00001": genome},
        families=[fam],
        ssn=CompleteSSN(edges=set()),
        config={"cutoff_i": 0.5, "family_mode": "full_graph"},
    )

    state = load_state(str(out))

    assert isinstance(state, LoadedState)
    assert state.genome_id_to_name == {"G00001": "GCA_0001"}
    assert state.gene_id_to_original == {
        "G00001|g0": "GCA_0001|g0",
        "G00001|g1": "GCA_0001|g1",
    }
    assert len(state.families) == 1
    assert state.families[0].member_genes == {"G00001|g0", "G00001|g1"}
    assert state.complete_ssn is not None
    assert state.config == {"cutoff_i": 0.5, "family_mode": "full_graph"}


def test_load_state_without_optional_files_still_succeeds(tmp_path):
    out = tmp_path / "out"
    _save_minimal_state(out, genomes={"G00001": _make_genome("G00001", "GCA_0001")})

    # 只有 genomes + families(空), 无 complete_ssn/intra/join/config。
    state = load_state(str(out))
    assert state.complete_ssn is None
    assert state.intra_results == {}
    assert state.join_result is None
    assert state.config == {}


def test_schema_version_mismatch_raises(tmp_path):
    out = tmp_path / "out"
    state_dir = out / ".linc_work" / "state"
    _write_wrong_schema_genomes(state_dir)
    # families.pkl 也要存在, 否则会先触发"缺少必需文件"而非 schema 错误。
    repo = FileRepository(str(state_dir))
    repo.save_families([])

    with pytest.raises(ReportInputError, match="schema 版本不兼容"):
        load_state(str(out))
