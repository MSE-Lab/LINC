"""`linc/pan/charts.py` 的单元测试: 5 张 SVG 生成 + 合法 XML。"""
import xml.etree.ElementTree as ET

from linc.core.entities import CompleteSSN, Gene, Genome, OrthologFamily, SSNEdge
from linc.pan.analysis import analyze
from linc.pan.charts import render_charts
from linc.pan.loader import LoadedState


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
    config = {"family_mode": "full_graph", "inflation": 1.5,
              "edge_weight_mode": "identity_x_coverage"}
    return LoadedState(
        output_dir="/x", genomes=genomes, families=families, complete_ssn=ssn,
        intra_results={}, join_result=None, config=config,
        gene_id_to_original=gene_to_original, genome_id_to_name=genome_to_name,
    )


def test_render_charts_writes_five_valid_svgs(tmp_path):
    state = _make_state()
    data = analyze(state)
    names = render_charts(state, data, str(tmp_path))

    assert names == [
        "cluster_size_distribution.svg",
        "mcl_before_after.svg",
        "mcl_top_components.svg",
        "pan_genome_families.svg",
        "pan_genome_genes.svg",
    ]
    for name in names:
        content = (tmp_path / name).read_text(encoding="utf-8")
        # 每张图都是合法 XML。
        ET.fromstring(content)
        assert content.startswith("<svg")
        assert content.endswith("</svg>\n")


def test_charts_handle_empty_data(tmp_path):
    state = LoadedState(
        output_dir="/x", genomes={}, families=[], complete_ssn=None,
        intra_results={}, join_result=None, config={},
        gene_id_to_original={}, genome_id_to_name={},
    )
    data = analyze(state)
    names = render_charts(state, data, str(tmp_path))
    for name in names:
        ET.fromstring((tmp_path / name).read_text(encoding="utf-8"))
