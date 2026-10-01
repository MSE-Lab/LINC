"""linc.cli 的 CLI 集成测试(使用 typer.testing.CliRunner)。"""
import json
import os

from typer.testing import CliRunner

from linc.cli import app

runner = CliRunner()


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
