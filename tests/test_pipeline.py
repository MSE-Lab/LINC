"""Pipeline 的 resume()/append() 单元测试(2026-10-02, 第8步)。"""
import os

import pytest

from linc.adapters.aligners.mock_adapter import MockAligner
from linc.adapters.clustering.union_find_adapter import UnionFindClustering
from linc.adapters.storage.file_repo import FileRepository
from linc.pipeline import IndexMissingError, Pipeline, PipelineConfig


class _IndexAlwaysMissingAligner(MockAligner):
    """包装 MockAligner, 模拟"索引已被手动删除"的场景(第二轮漏洞5回归测试,
    2026-10-03): `build_index`/`search` 等行为不变, 只有 `index_exists`
    恒返回 False。
    """

    def index_exists(self, index) -> bool:
        return False


def _write_fasta(path: str, records: list) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for header, sequence in records:
            fh.write(f">{header}\n{sequence}\n")


def _make_pipeline(repository) -> Pipeline:
    return Pipeline(
        aligner=MockAligner(),
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repository,
        config=PipelineConfig(),
    )


def test_resume_when_fully_saved_reuses_persisted_results_without_rerunning(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_path = input_dir / "g1.faa"
    _write_fasta(str(file_path), [("x", "MKV")])

    repo = FileRepository(str(tmp_path / "state"))
    pipeline = _make_pipeline(repo)
    first = pipeline.run([str(file_path)])

    second_pipeline = _make_pipeline(repo)
    second = second_pipeline.resume([str(file_path)])

    assert second.join_result.n_clusters_final == first.join_result.n_clusters_final
    assert len(second.families) == len(first.families)


def test_resume_from_partial_state_continues_remaining_stages(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_path = input_dir / "g1.faa"
    _write_fasta(str(file_path), [("x", "MKV")])

    repo = FileRepository(str(tmp_path / "state"))
    pipeline = _make_pipeline(repo)
    full = pipeline.run([str(file_path)])

    # 模拟"只跑到 intra_results 就中断": 清掉之后的3个状态文件
    state_dir = tmp_path / "state"
    for name in ("join_result.pkl", "complete_ssn.pkl", "families.pkl"):
        p = state_dir / name
        if p.exists():
            os.remove(p)

    resumed_pipeline = _make_pipeline(repo)
    resumed = resumed_pipeline.resume([str(file_path)])

    assert resumed.join_result.n_clusters_final == full.join_result.n_clusters_final
    assert len(resumed.families) == len(full.families)


def test_append_merges_new_genome_into_existing_result(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_a = input_dir / "gA.faa"
    _write_fasta(str(file_a), [("x", "MKV")])

    repo = FileRepository(str(tmp_path / "state"))
    pipeline = _make_pipeline(repo)
    first = pipeline.run([str(file_a)])
    assert set(first.genomes.keys()) == {"gA"}

    append_dir = tmp_path / "append"
    append_dir.mkdir()
    file_b = append_dir / "gB.faa"
    _write_fasta(str(file_b), [("y", "MKL")])

    append_pipeline = _make_pipeline(repo)
    appended = append_pipeline.append([str(file_b)])

    assert set(appended.genomes.keys()) == {"gA", "gB"}
    # 历史簇的 cluster_counter 应被延续, 不从0重新开始(新簇id应>=历史counter)。
    assert appended.join_result.cluster_counter >= first.join_result.cluster_counter


def test_append_cluster_ids_do_not_collide_with_historical_ones(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_a = input_dir / "gA.faa"
    _write_fasta(str(file_a), [("x", "MKV"), ("x2", "AAAAAAAAAAAAAAAAAAAA")])

    repo = FileRepository(str(tmp_path / "state"))
    pipeline = _make_pipeline(repo)
    first = pipeline.run([str(file_a)])
    historical_cluster_ids = {c.cluster_id for c in first.join_result.global_lib.clusters if c.active}

    append_dir = tmp_path / "append"
    append_dir.mkdir()
    file_b = append_dir / "gB.faa"
    _write_fasta(str(file_b), [("y", "CCCCCCCCCCCCCCCCCCCC")])

    append_pipeline = _make_pipeline(repo)
    appended = append_pipeline.append([str(file_b)])

    new_cluster_ids = {c.cluster_id for c in appended.join_result.global_lib.clusters if c.active}
    # 历史簇id应全部仍然存在于合并后的库中(未被覆盖)。
    assert historical_cluster_ids.issubset(new_cluster_ids)


def test_resume_raises_clear_error_when_index_missing(tmp_path):
    # 第二轮漏洞5回归测试: --resume 只校验 state 是否存在, 但未校验
    # intra_results.genome_index 指向的索引是否仍可用。若索引被手动删除
    # (这里用 _IndexAlwaysMissingAligner 模拟), resume 应在复用索引前
    # 主动抛出 IndexMissingError(而非让比对工具子进程崩溃), 错误信息中
    # 应点名缺失的 genome_id。
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_path = input_dir / "g1.faa"
    _write_fasta(str(file_path), [("x", "MKV")])

    repo = FileRepository(str(tmp_path / "state"))
    first_pipeline = Pipeline(
        aligner=MockAligner(),
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(),
    )
    first_pipeline.run([str(file_path)])

    # 模拟"只跑到 intra_results 就中断"(resume 需要真正走到
    # `repo.has_intra_results()` 分支才会触发索引校验), 同时换上一个
    # "索引恒缺失"的 aligner 模拟用户手动删除了索引目录。
    state_dir = tmp_path / "state"
    for name in ("join_result.pkl", "complete_ssn.pkl", "families.pkl"):
        p = state_dir / name
        if p.exists():
            os.remove(p)

    resume_pipeline = Pipeline(
        aligner=_IndexAlwaysMissingAligner(),
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(),
    )
    with pytest.raises(IndexMissingError) as exc_info:
        resume_pipeline.resume([str(file_path)])
    assert "g1" in str(exc_info.value)


def test_append_writes_audit_log_entry(tmp_path):
    """数据管理优化第1项(append 审计轨迹, 2026-10-03 与用户确认):
    每次 append() 成功后应写入一条审计记录, 包含本次追加的原始文件列表和
    实际新增的 genome_id 列表; 多次 append 应累加而非覆盖历史记录。"""
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_a = input_dir / "gA.faa"
    _write_fasta(str(file_a), [("x", "MKV")])

    repo = FileRepository(str(tmp_path / "state"))
    pipeline = _make_pipeline(repo)
    pipeline.run([str(file_a)])
    assert repo.load_append_history() == []  # run() 本身不写审计记录

    append_dir = tmp_path / "append"
    append_dir.mkdir()
    file_b = append_dir / "gB.faa"
    _write_fasta(str(file_b), [("y", "MKL")])

    append_pipeline = _make_pipeline(repo)
    append_pipeline.append([str(file_b)])

    history = repo.load_append_history()
    assert len(history) == 1
    assert history[0]["new_genome_ids"] == ["gB"]
    assert history[0]["new_genome_files"] == [str(file_b)]
    assert "timestamp" in history[0]

    # 再追加一次, 验证历史记录是累加而非覆盖。
    file_c = append_dir / "gC.faa"
    _write_fasta(str(file_c), [("z", "AAAA")])
    second_append_pipeline = _make_pipeline(repo)
    second_append_pipeline.append([str(file_c)])

    history_after_second = repo.load_append_history()
    assert len(history_after_second) == 2
    assert history_after_second[1]["new_genome_ids"] == ["gC"]
