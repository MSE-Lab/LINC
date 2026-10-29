"""`linc/pan/ssn_html.py` 的单元测试: 家族子图 + GraphML。"""
import xml.etree.ElementTree as ET

from linc.core.entities import CompleteSSN, Gene, Genome, OrthologFamily, SSNEdge
from linc.pan.analysis import analyze
from linc.pan.loader import LoadedState
from linc.pan.ssn_html import render_ssn_html, select_families


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
        config={"family_mode": "full_graph", "edge_weight_mode": "identity_x_coverage"},
        gene_id_to_original=gene_to_original, genome_id_to_name=genome_to_name,
    )


def test_select_families_explicit_and_top_n():
    state = _make_state()
    fams = select_families(state, family_ids=["OG0000001"])
    assert [f.family_id for f in fams] == ["OG0000001"]

    # 不指定时取 size 最大的 top_n。
    fams = select_families(state, top_n=2)
    assert [f.family_id for f in fams] == ["OG0000000", "OG0000001"]


def test_render_ssn_html_for_single_family(tmp_path):
    state = _make_state()
    data = analyze(state)
    names = render_ssn_html(state, data, str(tmp_path), family_ids=["OG0000000"])

    assert names == ["ssn.html"]
    html = (tmp_path / "ssn.html").read_text(encoding="utf-8")
    # 该家族的三个原始基因 ID 都应出现在节点里。
    assert "g1|g0" in html and "g2|g0" in html and "g3|g0" in html
    # 不属于该家族的基因不应作为节点出现(其原 ID 不出现)。
    assert "g1|g1" not in html
    assert html.startswith("<!DOCTYPE html>")


def test_render_ssn_html_with_graphml(tmp_path):
    state = _make_state()
    data = analyze(state)
    names = render_ssn_html(state, data, str(tmp_path),
                            family_ids=["OG0000001"], graphml=True)

    assert names == ["ssn.html", "ssn_view.graphml"]
    # GraphML 是合法 XML。
    ET.fromstring((tmp_path / "ssn_view.graphml").read_text(encoding="utf-8"))
