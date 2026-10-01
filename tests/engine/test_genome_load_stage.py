"""GenomeLoadStage 单元测试。"""
import os

from linc.core.entities import GenomeStatus
from linc.engine.genome_load_stage import GenomeLoadStage, discover_fasta_files


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
