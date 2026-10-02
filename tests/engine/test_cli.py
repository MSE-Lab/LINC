"""linc.cli 的 CLI 集成测试(使用 typer.testing.CliRunner)。

本文件全部用例默认通过 autouse fixture 将 `linc.cli._build_aligner` monkeypatch 为
返回 MockAligner(2026-10-02 与用户确认: MockAligner 不再作为 CLI 可选项暴露,
不引入隐藏的 --aligner 后门参数, 这类集成测试本意是验证 CLI 参数校验/日志/
摘要文件格式, 不需要依赖本机安装真实 mmseqs/diamond, 也不需要调用子进程)。"""
import json
import os

import pytest
from typer.testing import CliRunner

import linc.cli as cli_module
from linc.adapters.aligners.mock_adapter import MockAligner
from linc.adapters.clustering.union_find_adapter import UnionFindClustering
from linc.cli import app

runner = CliRunner()


@pytest.fixture(autouse=True)
def _use_mock_aligner(monkeypatch):
    monkeypatch.setattr(cli_module, "_build_aligner", lambda *a, **kw: MockAligner())


@pytest.fixture(autouse=True)
def _use_fake_clustering(monkeypatch):
    # FamilyRefineStage 的 MCL 细分步骤要求本机已安装系统级 mcl, 单测环境
    # 不应该强依赖该外部二进制; UnionFindClustering 同样满足 ClusteringPort
    # 接口(cluster(graph, **params) -> List[Set[str]]), 对每个连通分量
    # (本身已经是连通的子图)返回其自身作为唯一簇, 足以验证 CLI 层面的
    # 参数校验/日志/摘要文件格式, 不需要真实 MCL 聚类效果(2026-10-02)。
    monkeypatch.setattr(cli_module, "_build_clustering", lambda *a, **kw: UnionFindClustering())


def _write_fasta(path: str, records: list) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for header, sequence in records:
            fh.write(f">{header}\n{sequence}\n")


def test_run_loads_genomes_and_writes_summary(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV"), ("g2", "MKVLA")])

    output_dir = tmp_path / "out"

    result = runner.invoke(
        app, ["run", "-i", str(input_dir), "-o", str(output_dir)]
    )

    assert result.exit_code == 0
    assert "成功 1 个" in result.stdout

    summary_path = output_dir / "load_summary.json"
    assert summary_path.exists()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["total_genomes"] == 1
    assert summary["loaded_count"] == 1
    assert summary["failed_count"] == 0
    assert summary["genomes"]["genome01"]["gene_count"] == 2


def test_run_reports_failed_genomes_without_crashing(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "good.faa"), [("g1", "MK")])
    (input_dir / "empty.faa").write_text("", encoding="utf-8")

    output_dir = tmp_path / "out"

    result = runner.invoke(
        app, ["run", "-i", str(input_dir), "-o", str(output_dir)]
    )

    assert result.exit_code == 0
    assert "失败 1 个" in result.stdout
    assert "empty" in result.stdout

    summary = json.loads((output_dir / "load_summary.json").read_text(encoding="utf-8"))
    assert summary["failed_genome_ids"] == ["empty"]


def test_run_cutoff_i_greater_than_cutoff_ii_exits_with_error(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV")])
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        [
            "run", "-i", str(input_dir), "-o", str(output_dir),
            "--cutoff-i", "0.9", "--cutoff-ii", "0.5",
        ],
    )

    assert result.exit_code != 0
    assert not output_dir.exists()


def test_run_cutoff_i_equal_to_cutoff_ii_is_allowed(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV")])
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        [
            "run", "-i", str(input_dir), "-o", str(output_dir),
            "--cutoff-i", "0.5", "--cutoff-ii", "0.5",
        ],
    )

    assert result.exit_code == 0


def test_run_nonexistent_input_dir_exits_with_error(tmp_path):
    missing_dir = str(tmp_path / "does_not_exist")
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app, ["run", "-i", missing_dir, "-o", str(output_dir)]
    )

    assert result.exit_code != 0
    assert not output_dir.exists()


def test_run_empty_input_dir_still_succeeds_with_zero_genomes(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app, ["run", "-i", str(input_dir), "-o", str(output_dir)]
    )

    assert result.exit_code == 0
    summary = json.loads((output_dir / "load_summary.json").read_text(encoding="utf-8"))
    assert summary["total_genomes"] == 0


def test_run_writes_intra_summary_for_loaded_genomes(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    # g1/g2 序列长度相同 -> MockAligner 默认模式下 identity=1.0, 会被聚为同一簇;
    # g3 长度差异很大 -> 与 g1/g2 不成簇, 为孤儿基因。
    _write_fasta(
        str(input_dir / "genome01.faa"),
        [("g1", "MKV"), ("g2", "MKL"), ("g3", "MKVLAQWERTYUIOPASDFGH")],
    )

    output_dir = tmp_path / "out"

    result = runner.invoke(
        app, ["run", "-i", str(input_dir), "-o", str(output_dir)]
    )

    assert result.exit_code == 0
    assert "基因组内自聚类完成" in result.stdout

    intra_summary_path = output_dir / "intra_summary.json"
    assert intra_summary_path.exists()
    intra_summary = json.loads(intra_summary_path.read_text(encoding="utf-8"))
    assert intra_summary["total_genes"] == 3
    genome_summary = intra_summary["genomes"]["genome01"]
    assert genome_summary["n_genes"] == 3
    assert genome_summary["n_orphans"] >= 1


def test_run_intra_summary_skipped_for_failed_genomes(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    (input_dir / "empty.faa").write_text("", encoding="utf-8")

    output_dir = tmp_path / "out"

    result = runner.invoke(
        app, ["run", "-i", str(input_dir), "-o", str(output_dir)]
    )

    assert result.exit_code == 0
    intra_summary = json.loads(
        (output_dir / "intra_summary.json").read_text(encoding="utf-8")
    )
    assert intra_summary["total_genes"] == 0
    assert "empty" not in intra_summary["genomes"]


def test_run_writes_join_summary_across_genomes(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV"), ("g2", "MKL")])
    _write_fasta(str(input_dir / "genome02.faa"), [("g3", "MKV")])

    output_dir = tmp_path / "out"

    result = runner.invoke(
        app, ["run", "-i", str(input_dir), "-o", str(output_dir)]
    )

    assert result.exit_code == 0
    assert "增量并入全局代表库完成" in result.stdout

    join_summary_path = output_dir / "join_summary.json"
    assert join_summary_path.exists()
    join_summary = json.loads(join_summary_path.read_text(encoding="utf-8"))
    assert set(join_summary["processing_order"]) == {"genome01", "genome02"}
    assert join_summary["n_clusters_final"] >= 1
    assert join_summary["n_representatives_final"] >= 1
    assert isinstance(join_summary["clusters"], list)


def test_run_join_summary_empty_when_no_genomes_loaded(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    (input_dir / "empty.faa").write_text("", encoding="utf-8")

    output_dir = tmp_path / "out"

    result = runner.invoke(
        app, ["run", "-i", str(input_dir), "-o", str(output_dir)]
    )

    assert result.exit_code == 0
    join_summary = json.loads(
        (output_dir / "join_summary.json").read_text(encoding="utf-8")
    )
    assert join_summary["processing_order"] == []
    assert join_summary["n_clusters_final"] == 0
    assert join_summary["n_representatives_final"] == 0


def test_run_writes_ssn_summary_with_default_mode(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV"), ("g2", "MKL")])
    _write_fasta(str(input_dir / "genome02.faa"), [("g3", "MKV")])

    output_dir = tmp_path / "out"

    result = runner.invoke(
        app, ["run", "-i", str(input_dir), "-o", str(output_dir)]
    )

    assert result.exit_code == 0
    assert "SSN 构建完成" in result.stdout

    ssn_summary_path = output_dir / "ssn_summary.json"
    assert ssn_summary_path.exists()
    ssn_summary = json.loads(ssn_summary_path.read_text(encoding="utf-8"))
    assert ssn_summary["edge_completion_mode"] == "representative_member"
    assert isinstance(ssn_summary["n_edges"], int)


def test_run_edge_completion_mode_none_is_accepted(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV")])
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        ["run", "-i", str(input_dir), "-o", str(output_dir), "--edge-completion-mode", "none"],
    )

    assert result.exit_code == 0
    ssn_summary = json.loads((output_dir / "ssn_summary.json").read_text(encoding="utf-8"))
    assert ssn_summary["edge_completion_mode"] == "none"


def test_run_invalid_edge_completion_mode_exits_with_error(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV")])
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        ["run", "-i", str(input_dir), "-o", str(output_dir), "--edge-completion-mode", "bogus"],
    )

    assert result.exit_code != 0
    assert not output_dir.exists()


def test_run_invalid_mega_cluster_warning_threshold_exits_with_error(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV")])
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        [
            "run", "-i", str(input_dir), "-o", str(output_dir),
            "--mega-cluster-warning-threshold", "0",
        ],
    )

    assert result.exit_code != 0
    assert not output_dir.exists()


def test_run_invalid_aligner_exits_with_error(tmp_path):
    """--aligner 取值非法(含已移除的 "mock")应报错退出, 不产出任何摘要文件
    (2026-10-02 与用户确认: mock 不再作为 CLI 可选项暴露; 由于本用例要
    验证"真实 _build_aligner 对非法取值返回 None"的分支, 这里不使用
    autouse fixture 注入的 MockAligner, 而是单独 monkeypatch 还原真实
    _build_aligner 行为)。
    """
    import linc.cli as cli_module
    from linc.cli import _build_aligner

    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV")])
    output_dir = tmp_path / "out"

    import pytest as _pytest
    monkeypatch = _pytest.MonkeyPatch()
    monkeypatch.setattr(cli_module, "_build_aligner", _build_aligner)
    try:
        result = runner.invoke(
            app,
            ["run", "-i", str(input_dir), "-o", str(output_dir), "--aligner", "mock"],
        )
    finally:
        monkeypatch.undo()

    assert result.exit_code != 0
    assert not output_dir.exists()


def test_run_full_edge_completion_mode_end_to_end(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV"), ("g2", "MKL")])
    _write_fasta(str(input_dir / "genome02.faa"), [("g3", "MKV")])
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        ["run", "-i", str(input_dir), "-o", str(output_dir), "--edge-completion-mode", "full"],
    )

    assert result.exit_code == 0
    ssn_summary = json.loads((output_dir / "ssn_summary.json").read_text(encoding="utf-8"))
    assert ssn_summary["edge_completion_mode"] == "full"


def test_run_writes_family_summary_and_orthogroups(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV"), ("g2", "MKL")])
    _write_fasta(str(input_dir / "genome02.faa"), [("g3", "MKV")])
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app, ["run", "-i", str(input_dir), "-o", str(output_dir)]
    )

    assert result.exit_code == 0
    assert "MCL 亚家族细分完成" in result.stdout

    family_summary_path = output_dir / "family_summary.json"
    assert family_summary_path.exists()
    family_summary = json.loads(family_summary_path.read_text(encoding="utf-8"))
    assert family_summary["inflation"] == pytest.approx(1.5)
    assert family_summary["edge_weight_mode"] == "identity_x_coverage"
    assert isinstance(family_summary["n_families"], int)
    assert family_summary["n_families"] >= 1
    for family in family_summary["families"]:
        assert family["family_id"].startswith("OG")
        assert family["size"] >= 1

    orthogroups_path = output_dir / "Orthogroups.txt"
    assert orthogroups_path.exists()
    lines = orthogroups_path.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == family_summary["n_families"]
    for line in lines:
        og_id, members = line.split("\t")
        assert og_id.startswith("OG")
        assert len(members.split(" ")) >= 1


def test_run_custom_inflation_and_edge_weight_mode(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV")])
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        [
            "run", "-i", str(input_dir), "-o", str(output_dir),
            "--inflation", "2.0", "--edge-weight-mode", "bitscore",
        ],
    )

    assert result.exit_code == 0
    family_summary = json.loads(
        (output_dir / "family_summary.json").read_text(encoding="utf-8")
    )
    assert family_summary["inflation"] == pytest.approx(2.0)
    assert family_summary["edge_weight_mode"] == "bitscore"


def test_run_invalid_inflation_exits_with_error(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV")])
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        ["run", "-i", str(input_dir), "-o", str(output_dir), "--inflation", "0"],
    )

    assert result.exit_code != 0
    assert not output_dir.exists()


def test_run_invalid_edge_weight_mode_exits_with_error(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV")])
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        [
            "run", "-i", str(input_dir), "-o", str(output_dir),
            "--edge-weight-mode", "bogus",
        ],
    )

    assert result.exit_code != 0
    assert not output_dir.exists()


def test_run_invalid_edge_weight_alpha_exits_with_error(tmp_path):
    input_dir = tmp_path / "genomes"
    input_dir.mkdir()
    _write_fasta(str(input_dir / "genome01.faa"), [("g1", "MKV")])
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        [
            "run", "-i", str(input_dir), "-o", str(output_dir),
            "--edge-weight-alpha", "1.5",
        ],
    )

    assert result.exit_code != 0
    assert not output_dir.exists()
