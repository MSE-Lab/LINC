"""`linc/pan/analysis.py` 的单元测试: 泛基因组分类 + 统计 + TSV 导出。"""
from linc.core.entities import CompleteSSN, Gene, Genome, OrthologFamily, SSNEdge
from linc.pan.analysis import (
    PAN_CATEGORIES,
    ReportData,
    analyze,
    build_gene_to_genome,
    build_pan_genome_data,
    classify_family,
    connected_components,
    write_tsvs,
)
from linc.pan.loader import LoadedState


# ---------------------------------------------------------------------------
# classify_family
# ---------------------------------------------------------------------------
def test_classify_family_five_categories():
    assert classify_family(20, 20, 20) == "core"       # 100%
    assert classify_family(19, 20, 19) == "soft-core"  # 95%
    assert classify_family(3, 20, 3) == "shell"        # 15%
    assert classify_family(2, 20, 2) == "cloud"        # 10%
    assert classify_family(1, 20, 1) == "singleton"    # 1 成员


def test_classify_family_singleton_takes_priority_over_ratio():
    # 仅 1 个基因组时比例是 100%(本应 core), 但成员数为 1 应优先判 singleton。
    assert classify_family(1, 1, 1) == "singleton"


# ---------------------------------------------------------------------------
# build_pan_genome_data
# ---------------------------------------------------------------------------
def test_build_pan_genome_data_and_single_copy_core():
    fam_core = OrthologFamily("OG0", 0, {"A|g0", "B|g0", "C|g0"})
    # OG1 覆盖 3/3 基因组(core), 但 A 内有两个拷贝 -> 不应计入 single-copy-core。
    fam_multi = OrthologFamily("OG1", 0, {"A|g1", "A|g2", "B|g1", "C|g1"})
    gene_to_genome = {"A|g0": "A", "B|g0": "B", "C|g0": "C",
                      "A|g1": "A", "A|g2": "A", "B|g1": "B", "C|g1": "C"}
    data = build_pan_genome_data([fam_core, fam_multi], gene_to_genome,
                                 ["A", "B", "C"])

    assert data["family_category"] == {"OG0": "core", "OG1": "core"}
    # OG0 每基因组 1 拷贝 -> single-copy-core; OG1 在 A 内 2 拷贝 -> 不计。
    assert data["n_single_copy_core"] == 1
    assert data["presence_absence"]["OG0"] == {"A": 1, "B": 1, "C": 1}
    assert data["presence_absence"]["OG1"] == {"A": 2, "B": 1, "C": 1}
    # genome_breakdown: core 基因数 = 3(OG0) + 4(OG1) = 7
    assert data["genome_breakdown"]["A"]["core"] == 3
    assert data["genome_breakdown"]["B"]["core"] == 2
    assert data["genome_breakdown"]["C"]["core"] == 2


# ---------------------------------------------------------------------------
# connected_components
# ---------------------------------------------------------------------------
def test_connected_components():
    edges = {
        SSNEdge("a", "b", 0.9, 0.9),
        SSNEdge("b", "c", 0.9, 0.9),
        SSNEdge("d", "e", 0.9, 0.9),
    }
    comp_of_gene, components = connected_components(edges)
    assert len(components) == 2
    assert comp_of_gene["a"] == comp_of_gene["b"] == comp_of_gene["c"]
    assert comp_of_gene["d"] == comp_of_gene["e"]
    assert comp_of_gene["a"] != comp_of_gene["d"]


# ---------------------------------------------------------------------------
# analyze + write_tsvs 集成
# ---------------------------------------------------------------------------
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

    fam0 = OrthologFamily("OG0000000", 0,
                          {"G00001|g0", "G00002|g0", "G00003|g0"})
    fam1 = OrthologFamily("OG0000001", 0, {"G00001|g1", "G00002|g1"})
    fam2 = OrthologFamily("OG0000002", 1, {"G00003|g2"})
    families = [fam0, fam1, fam2]

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


def test_analyze_counts_and_ssn():
    state = _make_state()
    data = analyze(state)

    assert isinstance(data, ReportData)
    assert data.counts["n_genomes"] == 3
    assert data.counts["n_genes"] == 9
    assert data.counts["n_families"] == 3
    assert data.counts["n_ssn_edges"] == 3
    assert data.counts["n_genes_with_family"] == 6
    assert data.counts["n_genes_without_family"] == 3  # g2 三个基因组各剩 1 个

    assert data.ssn["n_nodes"] == 5  # 5 个基因出现在边里
    assert data.ssn["n_components"] == 2
    assert data.ssn["family_mode"] == "full_graph"

    # pan-genome: OG0 = 3/3 core; OG1 = 2/3 shell(67%); OG2 = singleton。
    assert data.family_category == {
        "OG0000000": "core", "OG0000001": "shell", "OG0000002": "singleton",
    }
    assert data.pan_genome["categories"]["core"]["n_families"] == 1
    assert data.pan_genome["categories"]["shell"]["n_families"] == 1
    assert data.pan_genome["categories"]["singleton"]["n_families"] == 1
    assert data.n_single_copy_core == 1

    # MCL 前后: 2 个连通分量, 但 component 0 被拆成 OG0/OG1 两个家族。
    assert data.mcl["n_components_before"] == 2
    assert data.mcl["n_families_after"] == 3
    assert data.mcl["n_components_split"] == 1


def test_write_tsvs_writes_all_nine_files(tmp_path):
    state = _make_state()
    data = analyze(state)
    names = write_tsvs(state, data, str(tmp_path))

    assert len(names) == 9
    for name in names:
        assert (tmp_path / name).exists()

    # pan_genome_summary.tsv 应含五类 + single-copy-core。
    summary = (tmp_path / "pan_genome_summary.tsv").read_text(encoding="utf-8")
    rows = {line.split("\t")[0]: line.split("\t")[1]
            for line in summary.strip().splitlines()[1:]}
    assert rows == {"core": "1", "soft-core": "0", "shell": "1",
                    "cloud": "0", "singleton": "1", "single-copy-core": "1",
                    "cloud+singleton": "1"}


def test_analyze_cloud_including_singleton_and_averages():
    data = analyze(_make_state())
    pan = data.pan_genome

    ci = pan["cloud_including_singleton"]
    assert ci["n_families"] == 1  # cloud(0) + singleton(1)
    assert ci["n_genes"] == 1
    # 基础 5 类仍互斥, 合并视图不改变它们。
    assert pan["categories"]["cloud"]["n_families"] == 0
    assert pan["categories"]["singleton"]["n_families"] == 1

    core = pan["categories"]["core"]
    assert core["avg_genes_per_genome"] == 1.0  # 3 基因 / 3 基因组
    assert core["avg_families_per_genome"] == 1.0
    assert pan["single_copy_core_pct_of_core"] == 1.0


def test_analyze_top_n_limits_top_lists():
    data = analyze(_make_state(), top_n=1)
    assert len(data.top_families) == 1
    assert len(data.top_components) == 1
