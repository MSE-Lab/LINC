"""FileRepository 单元测试(2026-10-02, 第8步: 持久化存储 + resume/append)。"""
import gzip
import pickle

import pytest

from linc.adapters.storage.file_repo import (
    FileRepository,
    SchemaVersionMismatchError,
    StateCorruptedError,
)
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


def test_saved_state_file_is_gzip_compressed(tmp_path):
    """2026-10-02 方案A: 状态文件应以 gzip 压缩写入, 而非原始 pickle。"""
    state_dir = tmp_path / "state"
    repo = FileRepository(str(state_dir))
    # 构造一段高重复度的长序列, 验证压缩后文件体积确实小于未压缩pickle。
    long_sequence = "MK" * 5000
    genome = Genome(
        genome_id="g1",
        file_path="g1.faa",
        genes=[],
    )
    repo.save_genomes({"g1": genome, "_padding": long_sequence})

    path = state_dir / "genomes.pkl"
    assert path.exists()
    # 文件头应是 gzip magic number, 而不是 pickle 的 opcode。
    with open(path, "rb") as fh:
        magic = fh.read(2)
    assert magic == b"\x1f\x8b"

    raw_pickle_size = len(pickle.dumps({"g1": genome, "_padding": long_sequence}))
    compressed_size = path.stat().st_size
    assert compressed_size < raw_pickle_size


def test_load_falls_back_to_plain_pickle_for_legacy_uncompressed_files(tmp_path):
    """向后兼容: 本次改动之前写入的旧版非压缩pickle文件仍应能正常加载,
    不应被误判为损坏。"""
    state_dir = tmp_path / "state"
    state_dir.mkdir()
    genome = Genome(genome_id="g1", file_path="g1.faa")
    legacy_path = state_dir / "genomes.pkl"
    with open(legacy_path, "wb") as fh:
        pickle.dump({"g1": genome}, fh)

    repo = FileRepository(str(state_dir))
    loaded = repo.load_genomes()
    assert loaded["g1"].genome_id == "g1"


def test_gzip_roundtrip_matches_manual_gzip_write(tmp_path):
    """确认 FileRepository 写入的文件能被标准 gzip 工具正常解压, 不是自定义格式。

    注: 文件内容是 `{"__schema_version__": N, "payload": obj}` 包装格式
    (工程审查问题2, 2026-10-03 与用户确认引入的 schema 版本标记), 这里需要
    先解出 "payload" 才是真正的 `IncrementalResult` 对象。
    """
    state_dir = tmp_path / "state"
    repo = FileRepository(str(state_dir))
    result = IncrementalResult(cluster_counter=7)
    repo.save_join_result(result)

    path = state_dir / "join_result.pkl"
    with gzip.open(path, "rb") as fh:
        raw = pickle.load(fh)
    assert raw["__schema_version__"] == FileRepository.SCHEMA_VERSION
    assert raw["payload"].cluster_counter == 7


def test_append_audit_log_accumulates_across_calls(tmp_path):
    """数据管理优化第1项(2026-10-03 与用户确认): append_audit_log 应只追加,
    不覆盖历史记录。"""
    repo = FileRepository(str(tmp_path / "state"))
    assert repo.load_append_history() == []

    repo.append_audit_log({"timestamp": "t1", "new_genome_ids": ["g2"]})
    repo.append_audit_log({"timestamp": "t2", "new_genome_ids": ["g3", "g4"]})

    history = repo.load_append_history()
    assert history == [
        {"timestamp": "t1", "new_genome_ids": ["g2"]},
        {"timestamp": "t2", "new_genome_ids": ["g3", "g4"]},
    ]


def test_append_audit_log_persists_across_instances(tmp_path):
    """跨进程持久化验证: 重新构造一个 FileRepository 实例仍能读到历史记录。"""
    state_dir = str(tmp_path / "state")
    repo1 = FileRepository(state_dir)
    repo1.append_audit_log({"timestamp": "t1", "new_genome_ids": ["g2"]})

    repo2 = FileRepository(state_dir)
    assert repo2.load_append_history() == [{"timestamp": "t1", "new_genome_ids": ["g2"]}]


def test_append_audit_log_is_human_readable_json_lines(tmp_path):
    """存储格式应为 JSON Lines(每行一条独立 JSON 对象), 而非 pickle。"""
    state_dir = tmp_path / "state"
    repo = FileRepository(str(state_dir))
    repo.append_audit_log({"timestamp": "t1", "new_genome_ids": ["g2"]})
    repo.append_audit_log({"timestamp": "t2", "new_genome_ids": ["g3"]})

    path = state_dir / "append_history.jsonl"
    assert path.exists()
    lines = path.read_text(encoding="utf-8").strip().split("\n")
    assert len(lines) == 2
    import json
    assert json.loads(lines[0]) == {"timestamp": "t1", "new_genome_ids": ["g2"]}


def test_schema_version_mismatch_raises_clear_error(tmp_path):
    """工程审查问题2(2026-10-03 与用户确认): 持久化文件标记的 schema_version
    与当前代码不一致时, 应抛出 `SchemaVersionMismatchError`(而不是让程序
    继续用缺字段/错类型的旧数据跑下去, 在下游某个随机位置抛出语义不清的
    AttributeError)。"""
    state_dir = tmp_path / "state"
    state_dir.mkdir()
    genome = Genome(genome_id="g1", file_path="g1.faa")
    wrapped = {"__schema_version__": 999, "payload": {"g1": genome}}
    with gzip.open(state_dir / "genomes.pkl", "wb") as fh:
        pickle.dump(wrapped, fh)

    repo = FileRepository(str(state_dir))
    with pytest.raises(SchemaVersionMismatchError):
        repo.load_genomes()


def test_schema_version_mismatch_is_subclass_of_state_corrupted_error(tmp_path):
    """`SchemaVersionMismatchError` 应是 `StateCorruptedError` 的子类,
    使 cli.py 中既有的 `except StateCorruptedError` 分支能一并捕获,
    无需额外修改异常处理代码。"""
    assert issubclass(SchemaVersionMismatchError, StateCorruptedError)


def test_legacy_file_without_schema_version_marker_loads_with_warning(tmp_path, caplog):
    """向后兼容(2026-10-03 与用户确认): 本次引入 schema_version 标记之前
    写入的旧版文件(无 `__schema_version__` 包装)不应被判定为版本不兼容,
    仅打印 warning 并原样按当前版本继续使用。"""
    state_dir = tmp_path / "state"
    state_dir.mkdir()
    genome = Genome(genome_id="g1", file_path="g1.faa")
    # 模拟本次改动之前的写入方式: 不包 __schema_version__ 包装, 直接 pickle
    # 原始数据。
    with gzip.open(state_dir / "genomes.pkl", "wb") as fh:
        pickle.dump({"g1": genome}, fh)

    repo = FileRepository(str(state_dir))
    loaded = repo.load_genomes()
    assert loaded["g1"].genome_id == "g1"


def test_run_id_is_persisted_and_reused_across_saves(tmp_path):
    """工程审查问题3(2026-10-03 与用户确认): 同一个 FileRepository 实例内
    多次 save_*() 应共用同一个 run_id, 且跨进程(新实例)复用 resume 场景下
    会读到同一个已落盘的 run_id, 而不是每次都重新生成。"""
    state_dir = str(tmp_path / "state")
    repo1 = FileRepository(state_dir)
    repo1.save_genomes({})
    run_id_1 = repo1._get_or_create_run_id()

    repo1.save_input_signature([("a.faa", 1, 1.0)])
    assert repo1._get_or_create_run_id() == run_id_1

    repo2 = FileRepository(state_dir)
    assert repo2._get_or_create_run_id() == run_id_1


def test_run_id_mismatch_raises_state_corrupted_error(tmp_path):
    """工程审查问题3(2026-10-03 与用户确认): 若 state 目录中某个文件的
    run_id 与目录当前的 run_id 不一致(典型场景: 文件被从别处拷贝拼接过),
    加载时应抛出 `StateCorruptedError`, 而不是静默当作同一次运行使用。"""
    state_dir = tmp_path / "state"
    state_dir.mkdir()
    genome = Genome(genome_id="g1", file_path="g1.faa")
    wrapped = {
        "__schema_version__": FileRepository.SCHEMA_VERSION,
        "__run_id__": "run-from-another-execution",
        "payload": {"g1": genome},
    }
    with gzip.open(state_dir / "genomes.pkl", "wb") as fh:
        pickle.dump(wrapped, fh)

    repo = FileRepository(str(state_dir))
    with pytest.raises(StateCorruptedError):
        repo.load_genomes()


def test_run_id_consistent_across_multiple_files_in_same_run(tmp_path):
    """同一次运行内保存的多个文件应携带相同的 run_id, 正常加载不应报错。"""
    state_dir = str(tmp_path / "state")
    repo = FileRepository(state_dir)
    repo.save_genomes({})
    repo.save_config({"cutoff_i": 0.5})

    # 正常加载不应抛出任何异常。
    assert repo.load_genomes() == {}
    assert repo.load_config() == {"cutoff_i": 0.5}


def test_config_roundtrip_across_instances(tmp_path):
    """工程审查问题1(2026-10-03 与用户确认): 算法参数应可持久化并跨进程
    加载, 供 --resume/--append 校验一致性。"""
    state_dir = str(tmp_path / "state")
    repo1 = FileRepository(state_dir)
    assert repo1.load_config() is None
    repo1.save_config({"cutoff_i": 0.5, "inflation": 1.5})

    repo2 = FileRepository(state_dir)
    assert repo2.load_config() == {"cutoff_i": 0.5, "inflation": 1.5}
