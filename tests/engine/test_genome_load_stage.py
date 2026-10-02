"""GenomeLoadStage 单元测试。"""
import os

import pytest

from linc.core.entities import GenomeStatus
from linc.engine.genome_load_stage import (
    GenomeIdConflictError,
    GenomeLoadStage,
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

    assert set(result.keys()) == {"genome01"}
    genome = result["genome01"]
    assert genome.status == GenomeStatus.LOADED
    assert genome.genome_id == "genome01"
    assert len(genome.genes) == 2


def test_gene_id_format_and_sequential_local_id(tmp_path):
    file_path = tmp_path / "genomeA.fasta"
    _write_fasta(str(file_path), [("x", "MK"), ("y", "MKL"), ("z", "MKLA")])

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])
    genome = result["genomeA"]

    gene_ids = [g.gene_id for g in genome.genes]
    assert gene_ids == ["genomeA|1", "genomeA|2", "genomeA|3"]
    for gene in genome.genes:
        assert gene.genome_id == "genomeA"


def test_gene_length_and_sequence_match_fasta_body(tmp_path):
    file_path = tmp_path / "genomeB.fa"
    # 序列跨多行,解析时应拼接为一整条
    with open(file_path, "w", encoding="utf-8") as fh:
        fh.write(">g1\nMKVL\nAQST\n")

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])
    gene = result["genomeB"].genes[0]

    assert gene.sequence == "MKVLAQST"
    assert gene.length == 8


def test_genome_id_extracted_from_filename_without_extension(tmp_path):
    file_path = tmp_path / "my_genome_007.faa"
    _write_fasta(str(file_path), [("g", "M")])

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])

    assert "my_genome_007" in result


def test_nonexistent_file_marks_failed_without_raising(tmp_path):
    missing_path = str(tmp_path / "does_not_exist.faa")

    stage = GenomeLoadStage()
    result = stage.run([missing_path])

    genome = result["does_not_exist"]
    assert genome.status == GenomeStatus.FAILED
    assert genome.genes == []


def test_empty_file_marks_failed(tmp_path):
    file_path = tmp_path / "empty.faa"
    file_path.write_text("", encoding="utf-8")

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])

    genome = result["empty"]
    assert genome.status == GenomeStatus.FAILED
    assert genome.genes == []


def test_one_failed_genome_does_not_block_others(tmp_path):
    good_path = tmp_path / "good.faa"
    _write_fasta(str(good_path), [("g1", "MK")])
    missing_path = str(tmp_path / "missing.faa")

    stage = GenomeLoadStage()
    result = stage.run([str(good_path), missing_path])

    assert result["good"].status == GenomeStatus.LOADED
    assert result["missing"].status == GenomeStatus.FAILED
    assert len(result["good"].genes) == 1


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


def test_duplicate_genome_id_raises_and_aborts(tmp_path):
    # foo.fa 与 foo.faa 的文件名 stem 相同, 均解析为 genome_id="foo",
    # 应检测到冲突并直接中断整个 run(), 而不是静默覆盖或只保留其中一个。
    fa_path = tmp_path / "foo.fa"
    faa_path = tmp_path / "foo.faa"
    _write_fasta(str(fa_path), [("g1", "MK")])
    _write_fasta(str(faa_path), [("g1", "MK")])

    stage = GenomeLoadStage()
    with pytest.raises(GenomeIdConflictError) as exc_info:
        stage.run([str(fa_path), str(faa_path)])

    assert "foo" in exc_info.value.conflicts
    assert set(exc_info.value.conflicts["foo"]) == {str(fa_path), str(faa_path)}


def test_no_conflict_when_genome_ids_are_distinct(tmp_path):
    # 对照: 文件名 stem 不同时不应触发冲突检测, run() 正常返回。
    foo_path = tmp_path / "foo.faa"
    bar_path = tmp_path / "bar.faa"
    _write_fasta(str(foo_path), [("g1", "MK")])
    _write_fasta(str(bar_path), [("g1", "MK")])

    stage = GenomeLoadStage()
    result = stage.run([str(foo_path), str(bar_path)])

    assert set(result.keys()) == {"foo", "bar"}


def test_leading_sequence_before_first_header_is_ignored(tmp_path):
    # 文件开头在第一个 ">" 之前出现的序列行应被忽略, 不计入任何基因,
    # 也不应与后续正常记录的序列产生拼接(2026-10-02 健壮性修复)。
    file_path = tmp_path / "leading.faa"
    file_path.write_text("MKVLSTRAY\n>g1\nMKV\n", encoding="utf-8")

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])
    genome = result["leading"]

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
    genome = result["empty_record"]

    assert genome.status == GenomeStatus.LOADED
    assert len(genome.genes) == 1
    assert genome.genes[0].sequence == "MKV"
    assert genome.genes[0].gene_id == "empty_record|1"


def test_file_with_only_empty_records_marks_failed(tmp_path):
    # 若文件里所有记录序列都为空(全部被跳过), 最终应没有任何基因,
    # 与"空文件"一致地标记为 FAILED, 而不是产出 0 个 Gene 的 LOADED 状态。
    file_path = tmp_path / "all_empty.faa"
    file_path.write_text(">g1\n>g2\n", encoding="utf-8")

    stage = GenomeLoadStage()
    result = stage.run([str(file_path)])
    genome = result["all_empty"]

    assert genome.status == GenomeStatus.FAILED
    assert genome.genes == []
