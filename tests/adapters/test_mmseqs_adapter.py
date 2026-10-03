"""MMseqsAligner 单元测试。

通过 monkeypatch 替换 shutil.which/subprocess.run, 不依赖真实安装的 mmseqs2
二进制程序, 可在任意环境下运行(包括未安装 mmseqs 的 CI)。
"""
import os
import subprocess

import pytest

from linc.adapters.aligners.mmseqs_adapter import MMseqsAligner, MMseqsConfig
from linc.core.entities import Gene
from linc.ports.aligner_port import AlignerPort


def _tmp_work_dir(tmp_path):
    return str(tmp_path / "linc_work")


def _gene(gene_id: str, length: int) -> Gene:
    return Gene(gene_id=gene_id, genome_id="g1", length=length, sequence="M" * length)


def test_missing_binary_raises_clear_error(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(RuntimeError, match="未找到 mmseqs"):
        MMseqsAligner()


def test_satisfies_aligner_port(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    assert isinstance(MMseqsAligner(), AlignerPort)


def _fake_run_factory(tsv_lines):
    """构造一个假的 subprocess.run, 把预设的 tsv 内容写入结果文件。"""

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        result_tsv = cmd[4]  # easy-search 的第5个位置参数是结果文件路径
        with open(result_tsv, "w", encoding="utf-8") as fh:
            fh.write("\n".join(tsv_lines))
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    return _fake_run


def test_align_all_vs_all_parses_result_and_filters_self_hit(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    tsv_lines = [
        "g1|1\tg1|1\t100.0\t1.0\t1.0\t0.0\t200.0\t100\t100\t100",  # 自身比对, 应被过滤
        "g1|1\tg1|2\t85.5\t0.9\t0.95\t1e-30\t150.0\t90\t100\t100",
    ]
    monkeypatch.setattr(subprocess, "run", _fake_run_factory(tsv_lines))

    aligner = MMseqsAligner()
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
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    tsv_lines = ["q1\ts1\t50.0\t0.5\t0.5\t1e-5\t20.0\t10\t10\t20"]
    monkeypatch.setattr(subprocess, "run", _fake_run_factory(tsv_lines))

    aligner = MMseqsAligner()
    query = [_gene("q1", 10)]
    subject = [_gene("s1", 20)]
    results = aligner.align_query_vs_subject(query, subject)

    assert len(results) == 1
    assert results[0].identity == pytest.approx(0.5)


def test_align_all_vs_all_parses_query_subject_length(monkeypatch):
    """qlen/tlen 列应解析为 query_length/subject_length, 为未来 MCL 权重
    归一化(如 OrthoFinder 式长度校正)预留(2026-10-02 与用户确认)。
    """
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    tsv_lines = ["g1|1\tg1|2\t85.5\t0.9\t0.95\t1e-30\t150.0\t90\t120\t130"]
    monkeypatch.setattr(subprocess, "run", _fake_run_factory(tsv_lines))

    aligner = MMseqsAligner()
    results = aligner.align_all_vs_all([_gene("g1|1", 120), _gene("g1|2", 130)])

    assert len(results) == 1
    assert results[0].query_length == 120
    assert results[0].subject_length == 130


def test_empty_query_or_subject_returns_empty_without_calling_subprocess(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")

    def _should_not_be_called(*args, **kwargs):
        raise AssertionError("subprocess.run 不应被调用")

    monkeypatch.setattr(subprocess, "run", _should_not_be_called)

    aligner = MMseqsAligner()
    assert aligner.align_all_vs_all([]) == []
    assert aligner.align_query_vs_subject([_gene("q1", 10)], []) == []


def test_malformed_output_line_is_skipped(monkeypatch, caplog):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    tsv_lines = ["g1|1\tg1|2\t85.5\t0.9\t0.95\t1e-30\t150.0"]  # 缺一列
    monkeypatch.setattr(subprocess, "run", _fake_run_factory(tsv_lines))

    aligner = MMseqsAligner()
    results = aligner.align_all_vs_all([_gene("g1|1", 100), _gene("g1|2", 100)])
    assert results == []


def test_non_numeric_field_line_is_skipped_not_raised(monkeypatch):
    # 列数正确但某字段是无法转为 float 的非数值字符串(如工具异常输出),
    # 应跳过该行而不抛异常中断流水线。
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    tsv_lines = [
        "g1|1\tg1|2\tN/A\t0.9\t0.95\t1e-30\t150.0\t90\t100\t100",
        "g1|1\tg1|3\t85.5\t0.9\t0.95\t1e-30\t150.0\t90\t100\t100",
    ]
    monkeypatch.setattr(subprocess, "run", _fake_run_factory(tsv_lines))

    aligner = MMseqsAligner()
    results = aligner.align_all_vs_all(
        [_gene("g1|1", 100), _gene("g1|2", 100), _gene("g1|3", 100)]
    )

    assert len(results) == 1
    assert results[0].subject_id == "g1|3"


def test_subprocess_failure_raises_runtime_error(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")

    def _fake_run_fail(cmd, check, capture_output, text, **kwargs):
        raise subprocess.CalledProcessError(
            returncode=1, cmd=cmd, stderr="mmseqs: command failed"
        )

    monkeypatch.setattr(subprocess, "run", _fake_run_fail)

    aligner = MMseqsAligner()
    with pytest.raises(RuntimeError, match="mmseqs easy-search 执行失败"):
        aligner.align_all_vs_all([_gene("g1|1", 100), _gene("g1|2", 100)])


def test_custom_config_values_are_used_in_command(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    captured_cmd = {}

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        captured_cmd["cmd"] = cmd
        result_tsv = cmd[4]
        with open(result_tsv, "w", encoding="utf-8") as fh:
            fh.write("")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner(MMseqsConfig(sensitivity=6.0, evalue=1e-10, threads=4))
    aligner.align_all_vs_all([_gene("g1|1", 100), _gene("g1|2", 100)])

    cmd = captured_cmd["cmd"]
    assert "6.0" in cmd
    assert "1e-10" in cmd
    assert "--threads" in cmd
    assert "4" in cmd


def test_easy_search_passes_max_seqs_to_avoid_default_truncation(monkeypatch):
    # 第二轮漏洞4: mmseqs 默认 --max-seqs=300 会静默截断大旁系家族的命中数,
    # 必须显式传递。
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    captured_cmd = {}

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        captured_cmd["cmd"] = cmd
        result_tsv = cmd[4]
        with open(result_tsv, "w", encoding="utf-8") as fh:
            fh.write("")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner(MMseqsConfig(max_seqs=50000))
    aligner.align_all_vs_all([_gene("g1|1", 100), _gene("g1|2", 100)])

    cmd = captured_cmd["cmd"]
    assert "--max-seqs" in cmd
    assert "50000" in cmd


def test_easy_search_default_max_seqs_is_not_mmseqs_default_300(monkeypatch):
    # 默认配置下也应显式传入一个远大于 mmseqs 自带默认值(300)的 --max-seqs,
    # 贴合"必须全量比对"的设计访对要求。
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    captured_cmd = {}

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        captured_cmd["cmd"] = cmd
        result_tsv = cmd[4]
        with open(result_tsv, "w", encoding="utf-8") as fh:
            fh.write("")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner()
    aligner.align_all_vs_all([_gene("g1|1", 100), _gene("g1|2", 100)])

    cmd = captured_cmd["cmd"]
    max_seqs_idx = cmd.index("--max-seqs")
    assert int(cmd[max_seqs_idx + 1]) > 300


def test_default_threads_zero_omits_threads_flag(monkeypatch):
    """threads<=0(默认)时不应传递 --threads, 避免 mmseqs 报 "Error in argument --threads"。"""
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    captured_cmd = {}

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        captured_cmd["cmd"] = cmd
        result_tsv = cmd[4]
        with open(result_tsv, "w", encoding="utf-8") as fh:
            fh.write("")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner()  # 默认 MMseqsConfig(threads=0)
    aligner.align_all_vs_all([_gene("g1|1", 100), _gene("g1|2", 100)])

    cmd = captured_cmd["cmd"]
    assert "--threads" not in cmd


def test_build_index_runs_createdb_and_returns_persistent_db_path(monkeypatch, tmp_path):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    captured = {}

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        captured["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner(MMseqsConfig(work_dir=_tmp_work_dir(tmp_path)))
    index = aligner.build_index([_gene("g1|1", 100)], index_id="genome_x")

    assert "createdb" in captured["cmd"]
    assert index is not None
    assert "genome_x" in index


def test_build_index_reuses_existing_index_without_calling_subprocess(monkeypatch, tmp_path):
    """幂等性回归测试(工程审查 force 风险修复, 2026-10-03 与用户确认):
    若 index_id 对应的索引文件已存在且有效(index_exists 为 True), build_index
    应直接复用返回现有路径, 不重新调用 createdb 子进程——避免 --force 全量
    重算等场景下对序列未变、仅算法参数变化的历史基因组做冗余建库。
    """
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")

    def _should_not_be_called(*args, **kwargs):
        raise AssertionError("subprocess.run 不应被调用(应复用已存在的索引)")

    monkeypatch.setattr(subprocess, "run", _should_not_be_called)

    work_dir = _tmp_work_dir(tmp_path)
    aligner = MMseqsAligner(MMseqsConfig(work_dir=work_dir))
    index_dir = os.path.join(work_dir, "genome_x")
    os.makedirs(index_dir, exist_ok=True)
    db_path = os.path.join(index_dir, "db")
    with open(db_path, "w", encoding="utf-8") as fh:
        fh.write("fake db")
    with open(f"{db_path}.dbtype", "w", encoding="utf-8") as fh:
        fh.write("fake dbtype")

    index = aligner.build_index([_gene("g1|1", 100)], index_id="genome_x")
    assert index == db_path


def test_build_index_empty_genes_returns_none_without_calling_subprocess(monkeypatch, tmp_path):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")

    def _should_not_be_called(*args, **kwargs):
        raise AssertionError("subprocess.run 不应被调用")

    monkeypatch.setattr(subprocess, "run", _should_not_be_called)

    aligner = MMseqsAligner(MMseqsConfig(work_dir=_tmp_work_dir(tmp_path)))
    assert aligner.build_index([]) is None


def test_index_exists_true_when_db_and_dbtype_files_present(tmp_path):
    # 第二轮漏洞5回归测试: index_exists 应据此判断 --resume 时索引是否可用。
    db_path = os.path.join(str(tmp_path), "db")
    with open(db_path, "w", encoding="utf-8") as fh:
        fh.write("fake db")
    with open(f"{db_path}.dbtype", "w", encoding="utf-8") as fh:
        fh.write("fake dbtype")

    aligner = MMseqsAligner(MMseqsConfig(work_dir=_tmp_work_dir(tmp_path)))
    assert aligner.index_exists(db_path) is True


def test_index_exists_false_when_index_files_missing(tmp_path):
    # 索引目录被手动删除(或从未建立)时应返回 False, 而不是让 search 子进程崩溃。
    db_path = os.path.join(str(tmp_path), "does_not_exist", "db")
    aligner = MMseqsAligner(MMseqsConfig(work_dir=_tmp_work_dir(tmp_path)))
    assert aligner.index_exists(db_path) is False
    assert aligner.index_exists(None) is False


def test_search_runs_createdb_search_convertalis_against_persisted_index(monkeypatch, tmp_path):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    calls = []

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        calls.append(cmd)
        if cmd[1] == "convertalis":
            result_tsv = cmd[5]
            with open(result_tsv, "w", encoding="utf-8") as fh:
                fh.write("q1\ts1\t80.0\t0.8\t0.8\t1e-10\t100.0\t50\t10\t10")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner(MMseqsConfig(work_dir=_tmp_work_dir(tmp_path)))
    persisted_index = os.path.join(str(tmp_path), "persisted", "db")
    results = aligner.search([_gene("q1", 10)], persisted_index)

    sub_commands = [c[1] for c in calls]
    assert sub_commands == ["createdb", "search", "convertalis"]
    # search/convertalis 都应直接使用传入的持久化 index 路径作为 target, 不重新 createdb
    assert calls[1][3] == persisted_index
    assert calls[2][3] == persisted_index
    # 第二轮漏洞4: search 子命令也应显式传递 --max-seqs, 避免默认值(300)截断。
    search_cmd = calls[1]
    assert "--max-seqs" in search_cmd
    assert int(search_cmd[search_cmd.index("--max-seqs") + 1]) > 300
    # 比对中间文件(query.fasta/result.tsv)应持久化保存在 index 所在目录的
    # search/ 子目录下, 不使用临时目录(2026-10-01, 与用户确认)。
    assert os.path.exists(os.path.join(str(tmp_path), "persisted", "search", "query.fasta"))
    assert os.path.exists(os.path.join(str(tmp_path), "persisted", "search", "result.tsv"))
    assert len(results) == 1
    assert results[0].identity == pytest.approx(0.8)


def test_search_overwrites_previous_results_for_same_index_dir(monkeypatch, tmp_path):
    """重跑同一 index(如同一 genome_id)时应覆盖旧的 search/ 目录内容。"""
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        if cmd[1] == "convertalis":
            result_tsv = cmd[5]
            with open(result_tsv, "w", encoding="utf-8") as fh:
                fh.write("q1\ts1\t80.0\t0.8\t0.8\t1e-10\t100.0\t50\t10\t10")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner(MMseqsConfig(work_dir=_tmp_work_dir(tmp_path)))
    persisted_index = os.path.join(str(tmp_path), "persisted", "db")

    search_dir = os.path.join(str(tmp_path), "persisted", "search")
    os.makedirs(search_dir, exist_ok=True)
    stale_file = os.path.join(search_dir, "stale.txt")
    with open(stale_file, "w", encoding="utf-8") as fh:
        fh.write("stale")

    aligner.search([_gene("q1", 10)], persisted_index)

    assert not os.path.exists(stale_file)
    assert os.path.exists(os.path.join(search_dir, "result.tsv"))


def test_search_tag_uses_separate_subdir_and_does_not_overwrite_other_tags(monkeypatch, tmp_path):
    """不同 tag 的 search 调用应写入各自独立的子目录, 互不覆盖
    (2026-10-02 与用户确认: 修复 intra自比对/join/full补边三类调用共用
    固定 `search/` 子目录而互相覆盖彼此中间结果的问题)。
    """
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        if cmd[1] == "convertalis":
            result_tsv = cmd[5]
            with open(result_tsv, "w", encoding="utf-8") as fh:
                fh.write("q1\ts1\t80.0\t0.8\t0.8\t1e-10\t100.0\t50\t10\t10")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner(MMseqsConfig(work_dir=_tmp_work_dir(tmp_path)))
    persisted_index = os.path.join(str(tmp_path), "persisted", "db")

    aligner.search([_gene("q1", 10)], persisted_index, tag="intra_self")
    aligner.search([_gene("q2", 10)], persisted_index, tag="join")

    intra_self_result = os.path.join(
        str(tmp_path), "persisted", "search_intra_self", "result.tsv"
    )
    join_result = os.path.join(str(tmp_path), "persisted", "search_join", "result.tsv")
    assert os.path.exists(intra_self_result)
    assert os.path.exists(join_result)


def test_search_self_search_skips_query_fasta_and_createdb(monkeypatch, tmp_path):
    """self_search=True(2026-10-02, 方案3, 与用户确认): query 与 index 为同一份
    序列集合(IntraGenomeStage 自比对场景)时, 应跳过为 query 重新写 fasta/
    createdb, 直接复用 index 本身作为 query_db, 避免 genes.fasta 与
    query.fasta 内容100%重复落盘。
    """
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    calls = []

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        calls.append(cmd)
        if cmd[1] == "convertalis":
            result_tsv = cmd[5]
            with open(result_tsv, "w", encoding="utf-8") as fh:
                fh.write("q1\tq1\t100.0\t1.0\t1.0\t0.0\t100.0\t50\t10\t10")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner(MMseqsConfig(work_dir=_tmp_work_dir(tmp_path)))
    persisted_index = os.path.join(str(tmp_path), "persisted", "db")

    aligner.search(
        [_gene("q1", 10)], persisted_index, tag="intra_self", self_search=True
    )

    sub_commands = [c[1] for c in calls]
    # 不应再调用 createdb(query 复用 index, 不重新建库)。
    assert sub_commands == ["search", "convertalis"]
    # search/convertalis 的 query_db 参数应直接是 index 路径本身。
    assert calls[0][2] == persisted_index
    assert calls[1][2] == persisted_index
    # 不应写入一份内容重复的 query.fasta。
    assert not os.path.exists(
        os.path.join(str(tmp_path), "persisted", "search_intra_self", "query.fasta")
    )
    assert os.path.exists(
        os.path.join(str(tmp_path), "persisted", "search_intra_self", "result.tsv")
    )


def test_cleanup_search_removes_tagged_subdir_but_not_index(monkeypatch, tmp_path):
    """数据管理优化第3项(2026-10-03 与用户确认): `cleanup_search(index, tag)`
    应删除 `<index_dir>/search_<tag>/` 中间文件子目录, 不影响索引本身。
    """
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        if cmd[1] == "convertalis":
            result_tsv = cmd[5]
            with open(result_tsv, "w", encoding="utf-8") as fh:
                fh.write("q1\ts1\t80.0\t0.8\t0.8\t1e-10\t100.0\t50\t10\t10")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner(MMseqsConfig(work_dir=_tmp_work_dir(tmp_path)))
    persisted_index = os.path.join(str(tmp_path), "persisted", "db")
    os.makedirs(os.path.dirname(persisted_index), exist_ok=True)
    with open(persisted_index, "w", encoding="utf-8") as fh:
        fh.write("fake db")

    aligner.search([_gene("q1", 10)], persisted_index, tag="join")
    search_dir = os.path.join(str(tmp_path), "persisted", "search_join")
    assert os.path.exists(search_dir)

    aligner.cleanup_search(persisted_index, tag="join")

    assert not os.path.exists(search_dir)
    # 索引本身不受影响。
    assert os.path.exists(persisted_index)


def test_cleanup_search_with_falsy_index_is_noop(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    aligner = MMseqsAligner()
    # 不应抛出任何异常。
    aligner.cleanup_search(None, tag="join")
    aligner.cleanup_search("", tag="join")


def test_search_empty_query_or_index_returns_empty_without_calling_subprocess(monkeypatch, tmp_path):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")

    def _should_not_be_called(*args, **kwargs):
        raise AssertionError("subprocess.run 不应被调用")

    monkeypatch.setattr(subprocess, "run", _should_not_be_called)

    aligner = MMseqsAligner(MMseqsConfig(work_dir=_tmp_work_dir(tmp_path)))
    assert aligner.search([], "/fake/db") == []
    assert aligner.search([_gene("q1", 10)], None) == []


def test_easy_search_passes_default_timeout_to_subprocess_run(monkeypatch, tmp_path):
    """工程审查问题8(2026-10-03 与用户确认): subprocess.run 应携带
    MMseqsConfig.timeout_seconds(默认3600秒), 避免外部工具卡死时进程
    永久挂起。"""
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")
    captured = {}

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        captured["timeout"] = kwargs.get("timeout")
        result_tsv = cmd[4]
        with open(result_tsv, "w", encoding="utf-8") as fh:
            fh.write("")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner(MMseqsConfig(work_dir=_tmp_work_dir(tmp_path)))
    aligner.align_all_vs_all([_gene("g1", 10), _gene("g2", 10)])

    assert captured["timeout"] == 3600


def test_easy_search_timeout_expired_raises_clear_runtime_error(monkeypatch, tmp_path):
    """工程审查问题8(2026-10-03 与用户确认): mmseqs 命令执行超时时,
    应被转换为清晰的 RuntimeError, 而不是让 subprocess.TimeoutExpired
    直接向上抛出。"""
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mmseqs")

    def _fake_run(cmd, check, capture_output, text, **kwargs):
        raise subprocess.TimeoutExpired(cmd=cmd, timeout=kwargs.get("timeout"))

    monkeypatch.setattr(subprocess, "run", _fake_run)

    aligner = MMseqsAligner(
        MMseqsConfig(work_dir=_tmp_work_dir(tmp_path), timeout_seconds=5)
    )
    with pytest.raises(RuntimeError, match="执行超时"):
        aligner.align_all_vs_all([_gene("g1", 10), _gene("g2", 10)])
