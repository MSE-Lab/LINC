"""linc.cli 的 `pan` 子命令集成测试(使用 typer.testing.CliRunner)。

直接构造 Orthogroups.txt + gene_genome.tsv 两个输入文件(不依赖完整
`linc run` 流水线), 验证 presence/absence 矩阵、core/soft-core/shell/
cloud/singleton 分类与各基因组分解表的计算逻辑。
"""
from typer.testing import CliRunner

from linc.cli import app

runner = CliRunner()


def _write_lines(path: str, lines: list) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for line in lines:
            fh.write(line + "\n")


def _setup_20_genome_fixture(tmp_path):
    """构造 20 个基因组(g1..g20)的 gene_genome.tsv, 以及覆盖
    core(ratio=100%)/soft-core(ratio=95%)/shell(ratio=50%)/
    cloud(ratio=10%)/singleton(1个基因, 无同源)五类的 Orthogroups.txt。
    用 20 个基因组是为了让 95%/99% 阈值能落在整数基因组个数上(19/20=95%)。
    """
    genomes = [f"g{i}" for i in range(1, 21)]

    gene_genome_lines = []
    for genome in genomes:
        gene_genome_lines.append(f"{genome}_core\t{genome}")
    for genome in genomes[:19]:
        gene_genome_lines.append(f"{genome}_softcore\t{genome}")
    for genome in genomes[:10]:
        gene_genome_lines.append(f"{genome}_shell\t{genome}")
    for genome in genomes[:2]:
        gene_genome_lines.append(f"{genome}_cloud\t{genome}")
    gene_genome_lines.append(f"{genomes[0]}_singleton\t{genomes[0]}")

    gene_genome_path = tmp_path / "gene_genome.tsv"
    _write_lines(str(gene_genome_path), gene_genome_lines)

    orthogroups_lines = [
        "OG0000000\t" + " ".join(f"{g}_core" for g in genomes),
        "OG0000001\t" + " ".join(f"{g}_softcore" for g in genomes[:19]),
        "OG0000002\t" + " ".join(f"{g}_shell" for g in genomes[:10]),
        "OG0000003\t" + " ".join(f"{g}_cloud" for g in genomes[:2]),
        f"OG0000004\t{genomes[0]}_singleton",
    ]
    orthogroups_path = tmp_path / "Orthogroups.txt"
    _write_lines(str(orthogroups_path), orthogroups_lines)

    return str(orthogroups_path), str(gene_genome_path)


def test_pan_classifies_core_softcore_shell_cloud_singleton(tmp_path):
    orthogroups_path, gene_genome_path = _setup_20_genome_fixture(tmp_path)
    output_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        [
            "pan",
            "-g", orthogroups_path,
            "-m", gene_genome_path,
            "-o", str(output_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout

    summary_path = output_dir / "pan_genome_summary.tsv"
    assert summary_path.exists()
    summary_text = summary_path.read_text(encoding="utf-8")
    rows = {
        line.split("\t")[0]: line.split("\t")[1]
        for line in summary_text.strip().splitlines()[1:]
    }
    assert rows["core"] == "1"
    assert rows["soft-core"] == "1"
    assert rows["shell"] == "1"
    assert rows["cloud"] == "1"
    assert rows["singleton"] == "1"
    # OG0000000(core) 在每个基因组内都只有 1 条成员基因, 应计入 single-copy-core。
    assert rows["single-copy-core"] == "1"


def test_pan_single_copy_core_excludes_multi_copy_core_family(tmp_path):
    """core 家族若存在任意基因组内拷贝数 > 1(多拷贝), 不应计入 single-copy-core,
    但仍应正常计入 core 分类。"""
    genomes = [f"g{i}" for i in range(1, 3)]
    gene_genome_lines = [f"{g}_core_a\t{g}" for g in genomes]
    gene_genome_lines.append(f"{genomes[0]}_core_b\t{genomes[0]}")

    gene_genome_path = tmp_path / "gene_genome.tsv"
    _write_lines(str(gene_genome_path), gene_genome_lines)

    orthogroups_lines = [
        "OG0000000\t"
        + " ".join([f"{g}_core_a" for g in genomes] + [f"{genomes[0]}_core_b"]),
    ]
    orthogroups_path = tmp_path / "Orthogroups.txt"
    _write_lines(str(orthogroups_path), orthogroups_lines)

    output_dir = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            "pan",
            "-g", str(orthogroups_path),
            "-m", str(gene_genome_path),
            "-o", str(output_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout
    summary_text = (output_dir / "pan_genome_summary.tsv").read_text(encoding="utf-8")
    rows = {
        line.split("\t")[0]: line.split("\t")[1]
        for line in summary_text.strip().splitlines()[1:]
    }
    assert rows["core"] == "1"
    assert rows["single-copy-core"] == "0"


def test_pan_writes_presence_absence_matrix_with_copy_numbers(tmp_path):
    orthogroups_path = tmp_path / "Orthogroups.txt"
    gene_genome_path = tmp_path / "gene_genome.tsv"

    # g1 对 OG0000000 有 2 个拷贝(同基因组内多条成员基因), g2 只有 1 个。
    _write_lines(
        str(gene_genome_path),
        ["g1_a\tg1", "g1_b\tg1", "g2_a\tg2"],
    )
    _write_lines(
        str(orthogroups_path),
        ["OG0000000\tg1_a g1_b g2_a"],
    )

    output_dir = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            "pan",
            "-g", str(orthogroups_path),
            "-m", str(gene_genome_path),
            "-o", str(output_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout

    matrix_path = output_dir / "presence_absence_matrix.tsv"
    lines = matrix_path.read_text(encoding="utf-8").strip().splitlines()
    header = lines[0].split("\t")
    assert header == ["Orthogroup", "g1", "g2"]
    row = lines[1].split("\t")
    assert row == ["OG0000000", "2", "1"]


def test_pan_per_genome_breakdown_counts_genes_not_families(tmp_path):
    orthogroups_path = tmp_path / "Orthogroups.txt"
    gene_genome_path = tmp_path / "gene_genome.tsv"

    # 单一基因组 g1, 单一家族(=100% core, 因为 n_total_genomes=1), 但
    # 该家族在 g1 内有 3 个拷贝, per_genome_breakdown 应统计为 3 个基因
    # 而不是 1 个家族。
    _write_lines(
        str(gene_genome_path),
        ["g1_a\tg1", "g1_b\tg1", "g1_c\tg1"],
    )
    _write_lines(
        str(orthogroups_path),
        ["OG0000000\tg1_a g1_b g1_c"],
    )

    output_dir = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            "pan",
            "-g", str(orthogroups_path),
            "-m", str(gene_genome_path),
            "-o", str(output_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout

    breakdown_path = output_dir / "per_genome_breakdown.tsv"
    lines = breakdown_path.read_text(encoding="utf-8").strip().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert row["Genome"] == "g1"
    assert row["core"] == "3"
    assert row["singleton"] == "0"


def test_pan_singleton_takes_priority_over_ratio(tmp_path):
    orthogroups_path = tmp_path / "Orthogroups.txt"
    gene_genome_path = tmp_path / "gene_genome.tsv"

    # 只有 1 个基因组时, 任何存在于该基因组的家族比例都是 100%(本应是
    # core), 但若该家族总成员数恰好为 1(无任何同源序列), 应优先分类
    # 为 singleton 而不是 core。
    _write_lines(str(gene_genome_path), ["g1_only\tg1"])
    _write_lines(str(orthogroups_path), ["OG0000000\tg1_only"])

    output_dir = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            "pan",
            "-g", str(orthogroups_path),
            "-m", str(gene_genome_path),
            "-o", str(output_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout
    summary_text = (output_dir / "pan_genome_summary.tsv").read_text(encoding="utf-8")
    rows = {
        line.split("\t")[0]: line.split("\t")[1]
        for line in summary_text.strip().splitlines()[1:]
    }
    assert rows["singleton"] == "1"
    assert rows["core"] == "0"


def test_pan_missing_orthogroups_file_exits_with_error(tmp_path):
    gene_genome_path = tmp_path / "gene_genome.tsv"
    _write_lines(str(gene_genome_path), ["g1_a\tg1"])

    result = runner.invoke(
        app,
        [
            "pan",
            "-g", str(tmp_path / "does_not_exist.txt"),
            "-m", str(gene_genome_path),
            "-o", str(tmp_path / "out"),
        ],
    )

    assert result.exit_code != 0


def test_pan_missing_gene_genome_file_exits_with_error(tmp_path):
    orthogroups_path = tmp_path / "Orthogroups.txt"
    _write_lines(str(orthogroups_path), ["OG0000000\tg1_a"])

    result = runner.invoke(
        app,
        [
            "pan",
            "-g", str(orthogroups_path),
            "-m", str(tmp_path / "does_not_exist.tsv"),
            "-o", str(tmp_path / "out"),
        ],
    )

    assert result.exit_code != 0


def test_pan_unknown_gene_id_in_orthogroups_is_warned_and_ignored(tmp_path):
    orthogroups_path = tmp_path / "Orthogroups.txt"
    gene_genome_path = tmp_path / "gene_genome.tsv"

    _write_lines(str(gene_genome_path), ["g1_a\tg1"])
    _write_lines(
        str(orthogroups_path),
        ["OG0000000\tg1_a unknown_gene_id"],
    )

    output_dir = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            "pan",
            "-g", str(orthogroups_path),
            "-m", str(gene_genome_path),
            "-o", str(output_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert "unknown_gene_id" in result.stdout
    matrix_text = (output_dir / "presence_absence_matrix.tsv").read_text(encoding="utf-8")
    assert "g1" in matrix_text
