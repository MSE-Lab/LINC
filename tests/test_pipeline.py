"""Pipeline 的 resume()/append() 单元测试(2026-10-02, 第8步)。"""
import os

from linc.adapters.aligners.mock_adapter import MockAligner
from linc.adapters.clustering.union_find_adapter import UnionFindClustering
from linc.adapters.storage.file_repo import FileRepository
from linc.pipeline import Pipeline, PipelineConfig


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
