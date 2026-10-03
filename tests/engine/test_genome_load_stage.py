"""GenomeLoadStage 单元测试。"""
import os

import pytest

from linc.core.entities import GenomeStatus
from linc.engine.genome_load_stage import (
    GenomeIdConflictError,
    GenomeLoadStage,
    build_input_signature,
    discover_fasta_files,
)


def _write_fasta(path: str, records: list) -> None:
    """records: list of (header, sequence) tuples。"""
    with open(path, "w", encoding="utf-8") as fh:
        for header, sequence in records:
            fh.write(f">{header}\n{sequence}\n")


def test_run_parses_single_valid_genome(tmp_path):
    file_path = tmp_path / "genome01.faa"
    _write_fasta(
        str(file_path),
        [("gene1 desc", "MKV"), ("gene2 desc", "MKVLA")],
    )

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])

    assert set(result.keys()) == {"G00001"}
    genome = result["G00001"]
    assert genome.status == GenomeStatus.LOADED
    assert genome.genome_id == "G00001"
    assert genome.source_file_name == "genome01"
    assert len(genome.genes) == 2


def test_gene_id_format_and_sequential_local_id(tmp_path):
    file_path = tmp_path / "genomeA.fasta"
    _write_fasta(str(file_path), [("x", "MK"), ("y", "MKL"), ("z", "MKLA")])

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])
    genome = result["G00001"]

    gene_ids = [g.gene_id for g in genome.genes]
    assert gene_ids == ["G00001|1", "G00001|2", "G00001|3"]
    for gene in genome.genes:
        assert gene.genome_id == "G00001"


def test_original_header_preserved_on_gene(tmp_path):
    # 核心修复: Gene 必须保留原始 FASTA header(取首 token), 不再因内部
    # genome_id/local_id 重编码而丢失原始序列 ID(2026-10-03 与用户确认)。
    file_path = tmp_path / "genomeA.fasta"
    _write_fasta(
        str(file_path),
        [("protXYZ some description", "MK"), ("protABC", "MKL")],
    )

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])
    genome = result["G00001"]

    assert [g.original_header for g in genome.genes] == ["protXYZ", "protABC"]


def test_duplicate_original_header_logs_warning_but_not_fail(tmp_path, caplog):
    file_path = tmp_path / "dup.faa"
    _write_fasta(str(file_path), [("same", "MK"), ("same", "MKL")])

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])
    genome = result["G00001"]

    assert genome.status == GenomeStatus.LOADED
    assert [g.original_header for g in genome.genes] == ["same", "same"]


def test_gene_length_and_sequence_match_fasta_body(tmp_path):
    file_path = tmp_path / "genomeB.fa"
    # 序列跨多行,解析时应拼接为一整条
    with open(file_path, "w", encoding="utf-8") as fh:
        fh.write(">g1\nMKVL\nAQST\n")

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])
    gene = result["G00001"].genes[0]

    assert gene.sequence == "MKVLAQST"
    assert gene.length == 8


def test_genome_id_is_internal_code_not_filename(tmp_path):
    # 2026-10-03 与用户确认"方案B": genome_id 不再取自文件名(可能含任意
    # 字符), 而是内部顺序短编号; 原始文件名保存在 source_file_name。
    file_path = tmp_path / "my_genome_007.faa"
    _write_fasta(str(file_path), [("g", "M")])

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])

    assert set(result.keys()) == {"G00001"}
    assert result["G00001"].source_file_name == "my_genome_007"


def test_multiple_genomes_get_sequential_codes_sorted_by_filename(tmp_path):
    path_b = tmp_path / "b.faa"
    path_a = tmp_path / "a.faa"
    _write_fasta(str(path_a), [("g1", "MK")])
    _write_fasta(str(path_b), [("g1", "MK")])

    stage = GenomeLoadStage()
    # 故意乱序传入, 编号分配应与传入顺序无关, 按文件名排序后分配
    result = stage.run([str(path_b), str(path_a)])

    assert result["G00001"].source_file_name == "a"
    assert result["G00002"].source_file_name == "b"


def test_nonexistent_file_marks_failed_without_raising(tmp_path):
    missing_path = str(tmp_path / "does_not_exist.faa")

    stage = GenomeLoadStage()
    result = stage.run([missing_path])

    genome = result["G00001"]
    assert genome.status == GenomeStatus.FAILED
    assert genome.source_file_name == "does_not_exist"
    assert genome.genes == []


def test_empty_file_marks_failed(tmp_path):
    file_path = tmp_path / "empty.faa"
    file_path.write_text("", encoding="utf-8")

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])

    genome = result["G00001"]
    assert genome.status == GenomeStatus.FAILED
    assert genome.genes == []


def test_existing_genomes_filename_conflict_raises(tmp_path):
    # --append 场景: 新文件名若与历史 genomes 的 source_file_name 冲突,
    # 应同样视为无法自行消歧的歧义, 直接中断。
    from linc.core.entities import Genome

    file_path = tmp_path / "genome01.faa"
    _write_fasta(str(file_path), [("g1", "MK")])

    historical = {"G00001": Genome(genome_id="G00001", file_path="hist.faa", source_file_name="genome01")}

    stage = GenomeLoadStage()
    with pytest.raises(GenomeIdConflictError) as exc_info:
        stage.run([str(file_path)], existing_genomes=historical)

    assert "genome01" in exc_info.value.conflicts


def test_existing_genomes_no_conflict_succeeds_and_continues_numbering(tmp_path):
    from linc.core.entities import Genome

    file_path = tmp_path / "genome02.faa"
    _write_fasta(str(file_path), [("g1", "MK")])

    historical = {"G00001": Genome(genome_id="G00001", file_path="hist.faa", source_file_name="genome01")}

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)], existing_genomes=historical)

    # 新 genome_id 应接续历史最大编号分配, 不与历史冲突
    assert set(result.keys()) == {"G00002"}
    assert result["G00002"].source_file_name == "genome02"


def test_one_failed_genome_does_not_block_others(tmp_path):
    good_path = tmp_path / "good.faa"
    _write_fasta(str(good_path), [("g1", "MK")])
    missing_path = str(tmp_path / "missing.faa")

    stage = GenomeLoadStage()
    # 按文件名排序后分配编号: good < missing
    result = stage.run([str(good_path), missing_path])

    good = next(g for g in result.values() if g.source_file_name == "good")
    missing = next(g for g in result.values() if g.source_file_name == "missing")
    assert good.status == GenomeStatus.LOADED
    assert missing.status == GenomeStatus.FAILED
    assert len(good.genes) == 1


def test_discover_fasta_files_finds_all_supported_extensions(tmp_path):
    (tmp_path / "a.fasta").write_text(">x\nM\n", encoding="utf-8")
    (tmp_path / "b.fa").write_text(">x\nM\n", encoding="utf-8")
    (tmp_path / "c.faa").write_text(">x\nM\n", encoding="utf-8")
    (tmp_path / "d.txt").write_text("not a fasta", encoding="utf-8")

    found = discover_fasta_files(str(tmp_path))

    found_names = sorted(os.path.basename(p) for p in found)
    assert found_names == ["a.fasta", "b.fa", "c.faa"]


def test_discover_fasta_files_returns_sorted_list(tmp_path):
    (tmp_path / "z.faa").write_text(">x\nM\n", encoding="utf-8")
    (tmp_path / "a.faa").write_text(">x\nM\n", encoding="utf-8")

    found = discover_fasta_files(str(tmp_path))

    assert found == sorted(found)


def test_discover_fasta_files_empty_directory_returns_empty_list(tmp_path):
    found = discover_fasta_files(str(tmp_path))
    assert found == []


def test_duplicate_genome_filename_raises_and_aborts(tmp_path):
    # foo.fa 与 foo.faa 的文件名 stem 相同, 应检测到冲突并直接中断整个
    # run(), 而不是静默覆盖或只保留其中一个。
    fa_path = tmp_path / "foo.fa"
    faa_path = tmp_path / "foo.faa"
    _write_fasta(str(fa_path), [("g1", "MK")])
    _write_fasta(str(faa_path), [("g1", "MK")])

    stage = GenomeLoadStage()
    with pytest.raises(GenomeIdConflictError) as exc_info:
        stage.run([str(fa_path), str(faa_path)])

    assert "foo" in exc_info.value.conflicts
    assert set(exc_info.value.conflicts["foo"]) == {str(fa_path), str(faa_path)}


def test_no_conflict_when_filenames_are_distinct(tmp_path):
    # 对照: 文件名 stem 不同时不应触发冲突检测, run() 正常返回。
    foo_path = tmp_path / "foo.faa"
    bar_path = tmp_path / "bar.faa"
    _write_fasta(str(foo_path), [("g1", "MK")])
    _write_fasta(str(bar_path), [("g1", "MK")])

    stage = GenomeLoadStage()
    result = stage.run([str(foo_path), str(bar_path)])

    source_names = {g.source_file_name for g in result.values()}
    assert source_names == {"foo", "bar"}


def test_leading_sequence_before_first_header_is_ignored(tmp_path):
    # 文件开头在第一个 ">" 之前出现的序列行应被忽略, 不计入任何基因,
    # 也不应与后续正常记录的序列产生拼接(2026-10-02 健壮性修复)。
    file_path = tmp_path / "leading.faa"
    file_path.write_text("MKVLSTRAY\n>g1\nMKV\n", encoding="utf-8")

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])
    genome = result["G00001"]

    assert genome.status == GenomeStatus.LOADED
    assert len(genome.genes) == 1
    assert genome.genes[0].sequence == "MKV"


def test_consecutive_headers_produce_empty_record_is_skipped(tmp_path):
    # 连续两个 ">" 行(header 之间无序列)会产生空序列记录, 应被跳过,
    # 不生成 length=0 的 Gene, 也不打乱后续基因的 local_id 连续性。
    file_path = tmp_path / "empty_record.faa"
    file_path.write_text(">g1\n>g2\nMKV\n", encoding="utf-8")

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])
    genome = result["G00001"]

    assert genome.status == GenomeStatus.LOADED
    assert len(genome.genes) == 1
    assert genome.genes[0].sequence == "MKV"
    assert genome.genes[0].gene_id == "G00001|1"
    assert genome.genes[0].original_header == "g2"


def test_file_with_only_empty_records_marks_failed(tmp_path):
    # 若文件里所有记录序列都为空(全部被跳过), 最终应没有任何基因,
    # 与"空文件"一致地标记为 FAILED, 而不是产出 0 个 Gene 的 LOADED 状态。
    file_path = tmp_path / "all_empty.faa"
    file_path.write_text(">g1\n>g2\n", encoding="utf-8")

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])
    genome = result["G00001"]

    assert genome.status == GenomeStatus.FAILED
    assert genome.genes == []


def test_build_input_signature_returns_sorted_basename_size_mtime_tuples(tmp_path):
    # 签名元素为 (basename, size, mtime) 三元组, 按 basename 排序
    # (2026-10-03, 第二轮漏洞"签名只比 basename"修复)。
    path_b = tmp_path / "b.faa"
    path_a = tmp_path / "a.faa"
    path_b.write_text(">x\nMK\n", encoding="utf-8")
    path_a.write_text(">x\nMKV\n", encoding="utf-8")

    signature = build_input_signature([str(path_b), str(path_a)])

    assert [item[0] for item in signature] == ["a.faa", "b.faa"]
    assert all(len(item) == 3 for item in signature)


def test_build_input_signature_changes_when_file_content_and_size_change(tmp_path):
    # 内容被替换且大小发生变化时, 签名应跟着变化(能拦住该误用场景, 与
    # 仅比 basename 的旧逻辑形成对比)。
    path = tmp_path / "genome01.faa"
    path.write_text(">x\nMK\n", encoding="utf-8")
    sig_before = build_input_signature([str(path)])

    path.write_text(">x\nMKVLSTRAYQWERTY\n", encoding="utf-8")
    sig_after = build_input_signature([str(path)])

    assert sig_before != sig_after
    assert sig_before[0][0] == sig_after[0][0] == "genome01.faa"
