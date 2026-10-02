"""DiamondAligner 单元测试。

通过 monkeypatch 替换 shutil.which/subprocess.run, 不依赖真实安装的 diamond
二进制程序, 可在任意环境下运行(包括未安装 diamond 的 CI)。风格对齐
`test_mmseqs_adapter.py`。
"""
import os
import subprocess

import pytest

from linc.adapters.aligners.diamond_adapter import DiamondAligner, DiamondConfig
from linc.core.entities import Gene
from linc.ports.aligner_port import AlignerPort


def _tmp_work_dir(tmp_path):
    return str(tmp_path / "linc_work")


def _gene(gene_id: str, length: int) -> Gene:
    return Gene(gene_id=gene_id, genome_id="g1", length=length, sequence="M" * length)


def test_missing_binary_raises_clear_error(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(RuntimeError, match="未找到 diamond"):
        DiamondAligner()


def test_satisfies_aligner_port(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")
    assert isinstance(DiamondAligner(), AlignerPort)


def _fake_run_factory(tsv_lines):
    """构造一个假的 subprocess.run, 把预设的 tsv 内容写入 blastp 的结果文件。"""

    def _fake_run(cmd, check, capture_output, text):
        if cmd[1] == "blastp":
            result_tsv = cmd[cmd.index("-o") + 1]
            with open(result_tsv, "w", encoding="utf-8") as fh:
                fh.write("\n".join(tsv_lines))
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    return _fake_run


def test_align_all_vs_all_parses_result_and_filters_self_hit(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")
    tsv_lines = [
        "g1|1\tg1|1\t100.0\t100.0\t100.0\t0.0\t200.0\t100",  # 自身比对, 应被过滤
        "g1|1\tg1|2\t85.5\t90.0\t95.0\t1e-30\t150.0\t90",
    ]
    monkeypatch.setattr(subprocess, "run", _fake_run_factory(tsv_lines))

    aligner = DiamondAligner()
    genes = [_gene("g1|1", 100), _gene("g1|2", 100)]
    results = aligner.align_all_vs_all(genes)

    assert len(results) == 1
    [result] = results
    assert result.query_id == "g1|1"
    assert result.subject_id == "g1|2"
    assert result.identity == pytest.approx(0.855)
    assert result.query_cover == pytest.approx(0.9)
    assert result.subject_cover == pytest.approx(0.95)
    assert result.coverage == pytest.approx(0.9)
    assert result.evalue == pytest.approx(1e-30)
    assert result.bitscore == pytest.approx(150.0)
    assert result.alignment_length == 90


def test_align_query_vs_subject_uses_separate_fasta(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")
    tsv_lines = ["q1\ts1\t50.0\t50.0\t50.0\t1e-5\t20.0\t10"]
    monkeypatch.setattr(subprocess, "run", _fake_run_factory(tsv_lines))

    aligner = DiamondAligner()
    query = [_gene("q1", 10)]
    subject = [_gene("s1", 20)]
    results = aligner.align_query_vs_subject(query, subject)

    assert len(results) == 1
    assert results[0].identity == pytest.approx(0.5)


def test_empty_query_or_subject_returns_empty_without_calling_subprocess(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")

    def _should_not_be_called(*args, **kwargs):
        raise AssertionError("subprocess.run 不应被调用")

    monkeypatch.setattr(subprocess, "run", _should_not_be_called)

    aligner = DiamondAligner()
    assert aligner.align_all_vs_all([]) == []
    assert aligner.align_query_vs_subject([_gene("q1", 10)], []) == []


def test_malformed_output_line_is_skipped(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")
    tsv_lines = ["g1|1\tg1|2\t85.5\t90.0\t95.0\t1e-30\t150.0"]  # 缺一列
    monkeypatch.setattr(subprocess, "run", _fake_run_factory(tsv_lines))

    aligner = DiamondAligner()
    results = aligner.align_all_vs_all([_gene("g1|1", 100), _gene("g1|2", 100)])
    assert results == []


def test_non_numeric_field_line_is_skipped_not_raised(monkeypatch):
    # 列数正确但某字段是无法转为 float 的非数值字符串(如工具异常输出),
    # 应跳过该行而不抛异常中断流水线(注意 "nan"/"inf" 本身能被 float()
    # 合法解析, 这里用一个真正非数值的字符串 "N/A" 来触发 ValueError)。
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")
    tsv_lines = [
        "g1|1\tg1|2\tN/A\t90.0\t95.0\t1e-30\t150.0\t90",
        "g1|1\tg1|3\t85.5\t90.0\t95.0\t1e-30\t150.0\t90",
    ]
    monkeypatch.setattr(subprocess, "run", _fake_run_factory(tsv_lines))

    aligner = DiamondAligner()
    results = aligner.align_all_vs_all(
        [_gene("g1|1", 100), _gene("g1|2", 100), _gene("g1|3", 100)]
    )

    assert len(results) == 1
    assert results[0].subject_id == "g1|3"


def test_subprocess_failure_raises_runtime_error(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")

    def _fake_run_fail(cmd, check, capture_output, text):
        raise subprocess.CalledProcessError(
            returncode=1, cmd=cmd, stderr="diamond: command failed"
        )

    monkeypatch.setattr(subprocess, "run", _fake_run_fail)

    aligner = DiamondAligner()
    with pytest.raises(RuntimeError, match="diamond 命令执行失败"):
        aligner.align_all_vs_all([_gene("g1|1", 100), _gene("g1|2", 100)])


def test_custom_config_values_are_used_in_command(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")
    captured_cmd = {}

    def _fake_run(cmd, check, capture_output, text):
        if cmd[1] == "blastp":
            captured_cmd["cmd"] = cmd
            result_tsv = cmd[cmd.index("-o") + 1]
            with open(result_tsv, "w", encoding="utf-8") as fh:
                fh.write("")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = DiamondAligner(
        DiamondConfig(sensitivity="ultra-sensitive", evalue=1e-10, threads=4)
    )
    aligner.align_all_vs_all([_gene("g1|1", 100), _gene("g1|2", 100)])

    cmd = captured_cmd["cmd"]
    assert "--ultra-sensitive" in cmd
    assert "1e-10" in cmd
    assert "--threads" in cmd
    assert "4" in cmd


def test_default_threads_zero_omits_threads_flag(monkeypatch):
    """threads<=0(默认)时不应传递 --threads。"""
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")
    captured_cmd = {}

    def _fake_run(cmd, check, capture_output, text):
        if cmd[1] == "blastp":
            captured_cmd["cmd"] = cmd
            result_tsv = cmd[cmd.index("-o") + 1]
            with open(result_tsv, "w", encoding="utf-8") as fh:
                fh.write("")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = DiamondAligner()  # 默认 DiamondConfig(threads=0)
    aligner.align_all_vs_all([_gene("g1|1", 100), _gene("g1|2", 100)])

    cmd = captured_cmd["cmd"]
    assert "--threads" not in cmd


def test_build_index_runs_makedb_and_returns_persistent_db_path(monkeypatch, tmp_path):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")
    captured = {}

    def _fake_run(cmd, check, capture_output, text):
        captured["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = DiamondAligner(DiamondConfig(work_dir=_tmp_work_dir(tmp_path)))
    index = aligner.build_index([_gene("g1|1", 100)], index_id="genome_x")

    assert "makedb" in captured["cmd"]
    assert index is not None
    assert "genome_x" in index


def test_build_index_empty_genes_returns_none_without_calling_subprocess(monkeypatch, tmp_path):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")

    def _should_not_be_called(*args, **kwargs):
        raise AssertionError("subprocess.run 不应被调用")

    monkeypatch.setattr(subprocess, "run", _should_not_be_called)

    aligner = DiamondAligner(DiamondConfig(work_dir=_tmp_work_dir(tmp_path)))
    assert aligner.build_index([]) is None


def test_search_runs_blastp_against_persisted_index_without_query_makedb(monkeypatch, tmp_path):
    """diamond 的 search() 不需要像 mmseqs 那样对 query 单独建库, 只调用一次 blastp。"""
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")
    calls = []

    def _fake_run(cmd, check, capture_output, text):
        calls.append(cmd)
        result_tsv = cmd[cmd.index("-o") + 1]
        with open(result_tsv, "w", encoding="utf-8") as fh:
            fh.write("q1\ts1\t80.0\t0.8\t0.8\t1e-10\t100.0\t50")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = DiamondAligner(DiamondConfig(work_dir=_tmp_work_dir(tmp_path)))
    persisted_index = os.path.join(str(tmp_path), "persisted", "db")
    results = aligner.search([_gene("q1", 10)], persisted_index)

    # search() 只应调用一次 blastp(不需要像 mmseqs 那样先 createdb query)
    assert len(calls) == 1
    assert calls[0][1] == "blastp"
    assert calls[0][calls[0].index("-d") + 1] == persisted_index
    # 比对中间文件(query.fasta/result.tsv)应持久化保存在 index 所在目录的
    # search/ 子目录下, 不使用临时目录(与 mmseqs_adapter 保持一致)。
    assert os.path.exists(os.path.join(str(tmp_path), "persisted", "search", "query.fasta"))
    assert os.path.exists(os.path.join(str(tmp_path), "persisted", "search", "result.tsv"))
    assert len(results) == 1
    assert results[0].identity == pytest.approx(0.8)


def test_search_overwrites_previous_results_for_same_index_dir(monkeypatch, tmp_path):
    """重跑同一 index(如同一 genome_id)时应覆盖旧的 search/ 目录内容。"""
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")

    def _fake_run(cmd, check, capture_output, text):
        result_tsv = cmd[cmd.index("-o") + 1]
        with open(result_tsv, "w", encoding="utf-8") as fh:
            fh.write("q1\ts1\t80.0\t0.8\t0.8\t1e-10\t100.0\t50")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = DiamondAligner(DiamondConfig(work_dir=_tmp_work_dir(tmp_path)))
    persisted_index = os.path.join(str(tmp_path), "persisted", "db")

    search_dir = os.path.join(str(tmp_path), "persisted", "search")
    os.makedirs(search_dir, exist_ok=True)
    stale_file = os.path.join(search_dir, "stale.txt")
    with open(stale_file, "w", encoding="utf-8") as fh:
        fh.write("stale")

    aligner.search([_gene("q1", 10)], persisted_index)

    assert not os.path.exists(stale_file)
    assert os.path.exists(os.path.join(search_dir, "result.tsv"))


def test_search_empty_query_or_index_returns_empty_without_calling_subprocess(monkeypatch, tmp_path):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/diamond")

    def _should_not_be_called(*args, **kwargs):
        raise AssertionError("subprocess.run 不应被调用")

    monkeypatch.setattr(subprocess, "run", _should_not_be_called)

    aligner = DiamondAligner(DiamondConfig(work_dir=_tmp_work_dir(tmp_path)))
    assert aligner.search([], "/fake/db") == []
    assert aligner.search([_gene("q1", 10)], None) == []
