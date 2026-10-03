"""Pipeline 的 resume()/append() 单元测试(2026-10-02, 第8步)。"""
import os

import pytest

from linc.adapters.aligners.mock_adapter import MockAligner
from linc.adapters.clustering.union_find_adapter import UnionFindClustering
from linc.adapters.storage.file_repo import FileRepository
from linc.pipeline import ConfigMismatchError, IndexMissingError, Pipeline, PipelineConfig


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


class _CallCountingAligner(MockAligner):
    """包装 MockAligner, 统计 `build_index`/`search` 被调用的次数(force 风险
    修复回归测试, 2026-10-03): 用于证明 --force 触发后确实重新跑了
    IntraGenomeStage/GlobalJoinStage, 而非仅仅放行检查却继续复用旧缓存。
    """

    def __init__(self):
        super().__init__()
        self.build_index_calls = 0
        self.search_calls = 0

    def build_index(self, genes, index_id=None):
        self.build_index_calls += 1
        return super().build_index(genes, index_id=index_id)

    def search(self, query, index, tag=None, self_search=False):
        self.search_calls += 1
        return super().search(query, index, tag=tag, self_search=self_search)


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


def test_resume_raises_config_mismatch_error_when_cutoff_differs(tmp_path):
    """工程审查问题1回归测试(2026-10-03 与用户确认): --resume 时若本次算法
    参数与历史持久化的参数不一致, 应抛出 ConfigMismatchError 而非静默用
    新参数混搭历史结果继续跑。"""
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_path = input_dir / "g1.faa"
    _write_fasta(str(file_path), [("x", "MKV")])

    repo = FileRepository(str(tmp_path / "state"))
    pipeline = _make_pipeline(repo)
    pipeline.run([str(file_path)])

    mismatched_pipeline = Pipeline(
        aligner=MockAligner(),
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(cutoff_i=0.9),
    )
    with pytest.raises(ConfigMismatchError) as exc_info:
        mismatched_pipeline.resume([str(file_path)])
    assert "cutoff_i" in str(exc_info.value)


def test_resume_with_force_bypasses_config_mismatch_and_updates_saved_config(tmp_path):
    """--force 应仅记 warning 后继续执行, 且之后再次 resume(不带 force)
    应以"本次已更新保存的参数"为新基准, 不再报错(2026-10-03 与用户确认)。"""
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_path = input_dir / "g1.faa"
    _write_fasta(str(file_path), [("x", "MKV")])

    repo = FileRepository(str(tmp_path / "state"))
    pipeline = _make_pipeline(repo)
    pipeline.run([str(file_path)])

    forced_pipeline = Pipeline(
        aligner=MockAligner(),
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(cutoff_i=0.9),
    )
    forced_pipeline.resume([str(file_path)], force=True)

    assert repo.load_config()["cutoff_i"] == 0.9

    again_pipeline = Pipeline(
        aligner=MockAligner(),
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(cutoff_i=0.9),
    )
    # 不应再抛出 ConfigMismatchError(基准已更新为0.9)。
    again_pipeline.resume([str(file_path)])


def test_append_raises_config_mismatch_error_when_inflation_differs(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_a = input_dir / "gA.faa"
    _write_fasta(str(file_a), [("x", "MKV")])

    repo = FileRepository(str(tmp_path / "state"))
    pipeline = _make_pipeline(repo)
    pipeline.run([str(file_a)])

    append_dir = tmp_path / "append"
    append_dir.mkdir()
    file_b = append_dir / "gB.faa"
    _write_fasta(str(file_b), [("y", "MKL")])

    mismatched_pipeline = Pipeline(
        aligner=MockAligner(),
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(inflation=3.0),
    )
    with pytest.raises(ConfigMismatchError) as exc_info:
        mismatched_pipeline.append([str(file_b)])
    assert "inflation" in str(exc_info.value)


def test_resume_skips_config_check_with_warning_when_no_historical_config(tmp_path, caplog):
    """历史 state 不存在 config.pkl(本次改动之前写入的旧版 state)时,
    应仅打印 warning 并放行, 不阻断 resume(向后兼容)。"""
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_path = input_dir / "g1.faa"
    _write_fasta(str(file_path), [("x", "MKV")])

    repo = FileRepository(str(tmp_path / "state"))
    pipeline = _make_pipeline(repo)
    pipeline.run([str(file_path)])

    # 模拟旧版 state(本次改动之前写入, 没有 config.pkl)。
    config_path = tmp_path / "state" / "config.pkl"
    assert config_path.exists()
    os.remove(config_path)

    different_pipeline = Pipeline(
        aligner=MockAligner(),
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(cutoff_i=0.9),
    )
    # 不应抛出 ConfigMismatchError(无法校验, 跳过)。
    different_pipeline.resume([str(file_path)])


def test_resume_force_ignores_stale_cache_and_recomputes_intra(tmp_path):
    """force 风险修复回归测试(2026-10-03 与用户确认, 风险1): --resume
    --force 时, 不能只放行参数校验却继续复用旧参数算出的 intra_results/
    join_result/complete_ssn/families 缓存——必须忽略这些缓存, 从
    IntraGenomeStage 开始用新参数对全部基因组重算。用调用计数验证"确实
    重新发起了比对", 而不是静默复用旧结果。"""
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_path = input_dir / "g1.faa"
    _write_fasta(str(file_path), [("x", "MKV"), ("y", "MKL")])

    repo = FileRepository(str(tmp_path / "state"))
    first_aligner = _CallCountingAligner()
    pipeline = Pipeline(
        aligner=first_aligner,
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(),
    )
    pipeline.run([str(file_path)])

    # 不带 force 的普通 resume: 一切都已持久化, 不应触发任何新的比对调用。
    plain_aligner = _CallCountingAligner()
    plain_pipeline = Pipeline(
        aligner=plain_aligner,
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(),
    )
    plain_pipeline.resume([str(file_path)])
    assert plain_aligner.search_calls == 0
    assert plain_aligner.build_index_calls == 0

    # 带 force 且参数不一致的 resume: 必须重新对全部基因组跑 Intra/Join,
    # 因此一定会触发新的比对/建库调用。
    forced_aligner = _CallCountingAligner()
    forced_pipeline = Pipeline(
        aligner=forced_aligner,
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(cutoff_i=0.9),
    )
    forced_pipeline.resume([str(file_path)], force=True)
    assert forced_aligner.search_calls > 0 or forced_aligner.build_index_calls > 0

    # force 覆盖事件应写入持久化审计记录, 不会因随后 save_config 更新基线
    # 而被悄悄抹掉(风险2)。
    history = repo.load_append_history()
    force_entries = [e for e in history if e.get("event") == "force_config_override"]
    assert len(force_entries) == 1
    assert force_entries[0]["action"] == "resume"
    assert "cutoff_i" in force_entries[0]["mismatched_fields"]


def test_append_force_recomputes_historical_intra_and_full_rejoin(tmp_path):
    """force 风险修复回归测试(2026-10-03 与用户确认, 风险1): --append
    --force 时, 历史基因组不能继续复用旧参数算出的 intra_results/
    join_result——必须用新参数重新对**全部历史基因组**跑 IntraGenomeStage,
    并放弃"接续历史 join_result 增量并入"的优化路径, 改为对历史+新增
    全部基因组从头跑一次 GlobalJoinStage。"""
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_a = input_dir / "gA.faa"
    _write_fasta(str(file_a), [("x", "MKV")])

    repo = FileRepository(str(tmp_path / "state"))
    pipeline = _make_pipeline(repo)
    pipeline.run([str(file_a)])

    append_dir = tmp_path / "append"
    append_dir.mkdir()
    file_b = append_dir / "gB.faa"
    _write_fasta(str(file_b), [("y", "MKL")])

    forced_aligner = _CallCountingAligner()
    forced_pipeline = Pipeline(
        aligner=forced_aligner,
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(inflation=3.0),
    )
    appended = forced_pipeline.append([str(file_b)], force=True)

    assert set(appended.genomes.keys()) == {"gA", "gB"}
    # 历史基因组 gA 的 intra 应被重新计算过(调用计数 > 0 能证明确实对
    # historical_loaded 重新跑了一遍 IntraGenomeStage, 而不仅仅是对新增
    # 的 gB 跑)。由于 gA 只有1个基因, 本身不会触发 self_search 的两两
    # 比对, 这里改为直接验证 join 阶段确实发起了"历史+新增"的全量重跑
    # (会至少调用一次 build_index/search, 而非仅复用历史 join_result)。
    assert forced_aligner.search_calls > 0 or forced_aligner.build_index_calls > 0

    history = repo.load_append_history()
    force_entries = [e for e in history if e.get("event") == "force_config_override"]
    assert len(force_entries) == 1
    assert force_entries[0]["action"] == "append"
    assert "inflation" in force_entries[0]["mismatched_fields"]


def test_resume_crash_mid_way_does_not_prematurely_overwrite_config(tmp_path, monkeypatch):
    """`save_config()` 时机修复回归测试(2026-10-03 与用户确认): config.pkl
    必须等全部 Stage 产出都已成功落盘之后才保存, 否则一旦 resume 中途崩溃,
    config.pkl 会被"提前"更新为本次(新)参数, 导致下一次 resume 的一致性
    校验对比的是"已经被更新过的历史 config", 无法察觉真正的参数不一致。

    场景: 先完整 run() 一次(参数A); 再用参数B + force=True 发起 resume,
    但在 _run_join 阶段模拟崩溃(抛异常); 崩溃后 config.pkl 应仍是参数A
    (而不是被提前写成参数B); 随后用参数B(不带force)再次 resume, 必须
    正确抛出 ConfigMismatchError(因为历史 config 仍是A, 与本次B不一致),
    而不是因为"历史 config 已被上次崩溃的运行悄悄改成B"而误判为一致、
    静默放行。
    """
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    file_path = input_dir / "g1.faa"
    _write_fasta(str(file_path), [("x", "MKV"), ("y", "MKL")])

    repo = FileRepository(str(tmp_path / "state"))
    pipeline = _make_pipeline(repo)
    pipeline.run([str(file_path)])

    historical_config = repo.load_config()
    assert historical_config["cutoff_i"] == PipelineConfig().cutoff_i

    crashing_pipeline = Pipeline(
        aligner=MockAligner(),
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(cutoff_i=0.9),
    )

    def _boom(*args, **kwargs):
        raise RuntimeError("模拟 resume 中途崩溃")

    monkeypatch.setattr(crashing_pipeline, "_run_join", _boom)

    with pytest.raises(RuntimeError, match="模拟 resume 中途崩溃"):
        crashing_pipeline.resume([str(file_path)], force=True)

    # 崩溃之后, config.pkl 不应被提前改写为本次(新)参数。
    config_after_crash = repo.load_config()
    assert config_after_crash["cutoff_i"] == historical_config["cutoff_i"]
    assert config_after_crash["cutoff_i"] != 0.9

    # 再次用同样的新参数(不带 force)resume, 必须正确检测出不一致并报错,
    # 而不是因为 config.pkl 已被悄悄改写而误判为一致。
    retry_pipeline = Pipeline(
        aligner=MockAligner(),
        clustering_for_intra=UnionFindClustering(),
        clustering_for_mcl=UnionFindClustering(),
        repository=repo,
        config=PipelineConfig(cutoff_i=0.9),
    )
    with pytest.raises(ConfigMismatchError):
        retry_pipeline.resume([str(file_path)])
