"""`linc/pan/report.py` 的单元测试: run_summary.json + report.html + 编排。"""
import json

from linc.adapters.storage.file_repo import FileRepository
from linc.core.entities import CompleteSSN, Gene, Genome, OrthologFamily, SSNEdge
from linc.pan.analysis import analyze
from linc.pan.loader import LoadedState
from linc.pan.report import build_report_html, build_run_summary, generate_report


def _make_state() -> LoadedState:
    genomes = {}
    gene_to_original = {}
    genome_to_name = {}
    for gid, name in (("G00001", "g1"), ("G00002", "g2"), ("G00003", "g3")):
        genes = []
        for i in range(3):
            gene_id = f"{gid}|g{i}"
            genes.append(Gene(gene_id=gene_id, genome_id=gid, length=100,
                              original_header=f"{name}|g{i}"))
            gene_to_original[gene_id] = f"{name}|g{i}"
        genomes[gid] = Genome(genome_id=gid, file_path=f"/x/{name}.faa",
                              genes=genes, source_file_name=name)
        genome_to_name[gid] = name

    families = [
        OrthologFamily("OG0000000", 0, {"G00001|g0", "G00002|g0", "G00003|g0"}),
        OrthologFamily("OG0000001", 0, {"G00001|g1", "G00002|g1"}),
        OrthologFamily("OG0000002", 1, {"G00003|g2"}),
    ]
    edges = {
        SSNEdge("G00001|g0", "G00002|g0", 0.9, 0.9, 100),
        SSNEdge("G00002|g0", "G00003|g0", 0.9, 0.9, 100),
        SSNEdge("G00001|g1", "G00002|g1", 0.9, 0.9, 100),
    }
    ssn = CompleteSSN(edges=edges, edge_completion_mode="full")
    return LoadedState(
        output_dir="/x", genomes=genomes, families=families, complete_ssn=ssn,
        intra_results={}, join_result=None,
        config={"family_mode": "full_graph", "inflation": 1.5,
                "edge_weight_mode": "identity_x_coverage"},
        gene_id_to_original=gene_to_original, genome_id_to_name=genome_to_name,
    )


def test_build_run_summary_structure():
    state = _make_state()
    data = analyze(state)
    summary = build_run_summary(state, data, visualization={"tsv": [], "charts": [], "ssn": []})

    assert summary["schema"] == "linc-run-summary/1"
    assert summary["counts"]["n_genomes"] == 3
    assert summary["counts"]["n_families"] == 3
    assert summary["pan_genome"]["categories"]["core"]["n_families"] == 1
    assert "generated_at" in summary
    assert "parameters" in summary


def test_build_report_html_contains_sections():
    state = _make_state()
    data = analyze(state)
    html = build_report_html(state, data)

    assert html.startswith("<!DOCTYPE html>")
    assert "关键指标" in html
    assert "泛基因组分析" in html
    assert "ssn.html" in html
    assert "Top families" in html
    assert "运行参数" in html


def _save_state_to_disk(out_dir):
    state = _make_state()
    repo = FileRepository(str(out_dir / ".linc_work" / "state"))
    repo.save_genomes(state.genomes)
    repo.save_families(state.families)
    repo.save_complete_ssn(state.complete_ssn)
    repo.save_config(state.config)
    return state


def test_generate_report_end_to_end(tmp_path):
    out = tmp_path / "out"
    _save_state_to_disk(out)

    report_dir = generate_report(str(out))

    assert report_dir == str(out / "report")
    # run_summary.json 写到两处。
    assert (out / "run_summary.json").exists()
    assert (out / "report" / "run_summary.json").exists()
    # report.html + 5 SVG + ssn.html + 9 TSV。
    assert (out / "report" / "report.html").exists()
    assert (out / "report" / "ssn.html").exists()
    assert (out / "report" / "pan_genome_summary.tsv").exists()
    assert (out / "report" / "pan_genome_families.svg").exists()

    summary = json.loads((out / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["counts"]["n_families"] == 3
    assert summary["visualization"]["report_dir"].endswith("report")
