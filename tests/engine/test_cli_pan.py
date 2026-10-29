"""`linc.cli` 的 `pan` 子命令集成测试(使用 typer.testing.CliRunner)。

`pan` 现在读取 `linc run` 的 state(而非 Orthogroups.txt + gene_genome.tsv),
本文件只验证 CLI 接线与错误处理; 泛基因组分类/统计逻辑由 tests/pan/ 下的
单元测试覆盖。
"""
from typer.testing import CliRunner

from linc.adapters.storage.file_repo import FileRepository
from linc.cli import app
from linc.core.entities import CompleteSSN, Gene, Genome, OrthologFamily, SSNEdge

runner = CliRunner()


def _save_state(out_dir):
    repo = FileRepository(str(out_dir / ".linc_work" / "state"))
    genomes = {}
    for gid, name in (("G00001", "g1"), ("G00002", "g2"), ("G00003", "g3")):
        genes = [Gene(gene_id=f"{gid}|g{i}", genome_id=gid, length=100,
                      original_header=f"{name}|g{i}") for i in range(2)]
        genomes[gid] = Genome(genome_id=gid, file_path=f"/x/{name}.faa",
                              genes=genes, source_file_name=name)
    families = [
        OrthologFamily("OG0000000", 0, {"G00001|g0", "G00002|g0", "G00003|g0"}),
        OrthologFamily("OG0000001", 0, {"G00001|g1", "G00002|g1"}),
    ]
    edges = {
        SSNEdge("G00001|g0", "G00002|g0", 0.9, 0.9, 100),
        SSNEdge("G00001|g1", "G00002|g1", 0.9, 0.9, 100),
    }
    repo.save_genomes(genomes)
    repo.save_families(families)
    repo.save_complete_ssn(CompleteSSN(edges=edges, edge_completion_mode="full"))
    repo.save_config({"family_mode": "full_graph", "inflation": 1.5,
                      "edge_weight_mode": "identity_x_coverage"})


def test_pan_generates_report_from_state(tmp_path):
    out = tmp_path / "out"
    _save_state(out)

    result = runner.invoke(app, ["pan", "-o", str(out)])

    assert result.exit_code == 0, result.stdout
    assert (out / "report" / "report.html").exists()
    assert (out / "report" / "run_summary.json").exists()
    assert (out / "run_summary.json").exists()
    assert (out / "report" / "pan_genome_summary.tsv").exists()
    assert (out / "report" / "ssn.html").exists()


def test_pan_missing_state_exits_with_error(tmp_path):
    result = runner.invoke(app, ["pan", "-o", str(tmp_path / "no_such_run")])
    assert result.exit_code != 0


def test_pan_with_family_and_graphml(tmp_path):
    out = tmp_path / "out"
    _save_state(out)

    result = runner.invoke(
        app, ["pan", "-o", str(out), "--family", "OG0000000", "--graphml"],
    )
    assert result.exit_code == 0, result.stdout
    assert (out / "report" / "ssn_view.graphml").exists()
