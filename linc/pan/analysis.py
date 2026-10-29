"""统计与聚合: 泛基因组分类(迁移自 cli.py)+ SSN/家族/MCL 统计 + TSV 导出。

设计依据: docs/LINC_结果描述增强开发计划.md 第 3/4 节。本模块只做"读入
`LoadedState` -> 计算 -> 写 TSV", 不产生任何比对或聚类副作用。
"""
import logging
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

import networkx as nx

from ..core.entities import OrthologFamily, SSNEdge, VALID_EDGE_WEIGHT_MODES
from .loader import LoadedState

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# 泛基因组分类(迁移自 linc/cli.py, 与用户确认的 5 类互斥约定)
# ---------------------------------------------------------------------------
PAN_CATEGORIES = ["core", "soft-core", "shell", "cloud", "singleton"]


def classify_family(n_genomes_present: int, n_total_genomes: int, n_members: int) -> str:
    """按用户确认的规则对一个直系同源家族(OG)分类:
    - 家族总成员数为 1(即只出现在 1 个基因组且该基因组内只有 1 条基因,
      无任何同源序列) -> singleton, 优先于下面的比例判断。
    - 否则按"覆盖的基因组数 / 总基因组数"的比例:
      >=99% -> core | [95%, 99%) -> soft-core |
      [15%, 95%) -> shell | <15% -> cloud。
    """
    if n_members == 1:
        return "singleton"
    ratio = n_genomes_present / n_total_genomes
    if ratio >= 0.99:
        return "core"
    if ratio >= 0.95:
        return "soft-core"
    if ratio >= 0.15:
        return "shell"
    return "cloud"


def build_pan_genome_data(
    families: List[OrthologFamily],
    gene_to_genome: Dict[str, str],
    genome_ids: List[str],
) -> Dict[str, Any]:
    """基于 state 里的家族与 gene->genome 映射, 一次性计算:
    presence_absence / family_category / genome_breakdown / n_single_copy_core。

    gene_to_genome: 内部 gene_id -> 内部 genome_id。
    """
    n_total_genomes = len(genome_ids)
    presence_absence: Dict[str, Dict[str, int]] = {}
    family_category: Dict[str, str] = {}
    n_single_copy_core = 0
    genome_breakdown = {
        gid: {category: 0 for category in PAN_CATEGORIES} for gid in genome_ids
    }

    for fam in families:
        copy_counts: Dict[str, int] = {}
        for gene_id in fam.member_genes:
            gid = gene_to_genome.get(gene_id)
            if gid is None:
                # state 应当自洽; 防御性跳过(不影响分类主干)。
                logger.warning("家族 %s 的基因 %s 未在 genomes 中找到, 已忽略",
                               fam.family_id, gene_id)
                continue
            copy_counts[gid] = copy_counts.get(gid, 0) + 1

        presence_absence[fam.family_id] = copy_counts
        n_genomes_present = len(copy_counts)
        n_members = sum(copy_counts.values())
        category = classify_family(n_genomes_present, n_total_genomes, n_members)
        family_category[fam.family_id] = category
        if category == "core" and all(count == 1 for count in copy_counts.values()):
            n_single_copy_core += 1

        for gid, count in copy_counts.items():
            genome_breakdown[gid][category] += count

    return {
        "presence_absence": presence_absence,
        "family_category": family_category,
        "n_single_copy_core": n_single_copy_core,
        "genome_breakdown": genome_breakdown,
    }


# ---------------------------------------------------------------------------
# 辅助: 建映射 / 连通分量
# ---------------------------------------------------------------------------
def build_gene_to_genome(genomes) -> Dict[str, str]:
    """内部 gene_id -> 内部 genome_id。"""
    return {
        gene.gene_id: genome.genome_id
        for genome in genomes.values()
        for gene in genome.genes
    }


def build_gene_to_family(families: List[OrthologFamily]) -> Dict[str, str]:
    """内部 gene_id -> family_id。"""
    return {
        gene_id: fam.family_id
        for fam in families
        for gene_id in fam.member_genes
    }


def connected_components(edges: Set[SSNEdge]) -> Tuple[Dict[str, int], Dict[int, Set[str]]]:
    """对无向 SSN 边集求连通分量。

    返回 (gene_id -> component_id, component_id -> set(gene_id))。孤立节点
    (不出现在任何边中的基因)不在此返回, 由调用方按 n_genes 另行推导。
    """
    graph = nx.Graph()
    for edge in edges:
        graph.add_edge(edge.source, edge.target)
    comp_of_gene: Dict[str, int] = {}
    components: Dict[int, Set[str]] = {}
    for idx, comp in enumerate(nx.connected_components(graph)):
        for node in comp:
            comp_of_gene[node] = idx
        components[idx] = set(comp)
    return comp_of_gene, components


# ---------------------------------------------------------------------------
# 报告数据
# ---------------------------------------------------------------------------
@dataclass
class ReportData:
    """一次报告所需的全部聚合统计。字段名与 run_summary.json 顶层对齐。"""

    counts: Dict[str, Any] = field(default_factory=dict)
    derived: Dict[str, Any] = field(default_factory=dict)
    ssn: Dict[str, Any] = field(default_factory=dict)
    families: Dict[str, Any] = field(default_factory=dict)
    pan_genome: Dict[str, Any] = field(default_factory=dict)
    mcl: Dict[str, Any] = field(default_factory=dict)
    cluster_size_distribution: Dict[str, Any] = field(default_factory=dict)
    genomes: Dict[str, Any] = field(default_factory=dict)
    top_families: List[Dict[str, Any]] = field(default_factory=list)
    top_components: List[Dict[str, Any]] = field(default_factory=list)
    # 供 TSV / charts / report 复用的中间数据(不直接进 run_summary.json)。
    family_category: Dict[str, str] = field(default_factory=dict)
    presence_absence: Dict[str, Dict[str, int]] = field(default_factory=dict)
    genome_breakdown: Dict[str, Dict[str, int]] = field(default_factory=dict)
    n_single_copy_core: int = 0
    comp_of_gene: Dict[str, int] = field(default_factory=dict)
    components: Dict[int, Set[str]] = field(default_factory=dict)
    gene_to_genome: Dict[str, str] = field(default_factory=dict)
    gene_to_family: Dict[str, str] = field(default_factory=dict)
    degree: Dict[str, int] = field(default_factory=dict)
    representative_genes: Set[str] = field(default_factory=set)


def _edge_weight_mode(config: Dict[str, Any]) -> str:
    mode = config.get("edge_weight_mode", "identity_x_coverage")
    return mode if mode in VALID_EDGE_WEIGHT_MODES else "identity_x_coverage"


def analyze(state: LoadedState) -> ReportData:
    """从 LoadedState 计算全部统计, 返回 ReportData。"""
    genomes = state.genomes
    families = state.families
    ssn = state.complete_ssn
    intra_results = state.intra_results
    join_result = state.join_result
    config = state.config

    genome_ids = sorted(genomes.keys())
    n_genomes = len(genome_ids)
    n_genes = sum(len(g.genes) for g in genomes.values())

    gene_to_genome = build_gene_to_genome(genomes)
    gene_to_family = build_gene_to_family(families)
    families_size_1 = {fam.family_id: fam.size == 1 for fam in families}

    # ---- counts ----
    n_local_clusters = sum(len(r.local_clusters) for r in intra_results.values())
    n_local_representatives = sum(r.n_representatives for r in intra_results.values())
    global_clusters = 0
    global_representatives = 0
    representative_genes: Set[str] = set()
    if join_result is not None:
        global_clusters = join_result.global_lib.n_active_clusters
        global_representatives = join_result.global_lib.n_representatives
        for cluster in join_result.global_lib.clusters:
            if cluster.active:
                representative_genes |= cluster.representative_ids

    edges: Set[SSNEdge] = ssn.edges if ssn is not None else set()
    n_ssn_edges = len(edges)

    n_genes_with_family = sum(fam.size for fam in families)
    n_genes_without_family = n_genes - n_genes_with_family
    n_mcl_components = len({fam.component_id for fam in families})

    counts = {
        "n_genomes": n_genomes,
        "n_genes": n_genes,
        "n_local_clusters": n_local_clusters,
        "n_local_representatives": n_local_representatives,
        "n_global_clusters": global_clusters,
        "n_global_representatives": global_representatives,
        "n_ssn_edges": n_ssn_edges,
        "n_mcl_components": n_mcl_components,
        "n_families": len(families),
        "n_genes_with_family": n_genes_with_family,
        "n_genes_without_family": n_genes_without_family,
    }

    # ---- 连通分量 + degree ----
    comp_of_gene, components = connected_components(edges)
    ssn_nodes = len(comp_of_gene)
    n_isolated = n_genes - ssn_nodes

    degree: Dict[str, int] = {node: 0 for node in comp_of_gene}
    for edge in edges:
        if edge.source in degree:
            degree[edge.source] += 1
        if edge.target in degree:
            degree[edge.target] += 1
    degree_values = sorted(degree.values())
    median_degree = _median(degree_values)

    def _cross_genome(edge: SSNEdge) -> bool:
        return gene_to_genome.get(edge.source) != gene_to_genome.get(edge.target)

    n_cross_genome_edges = sum(1 for e in edges if _cross_genome(e))
    n_same_family_edges = sum(
        1 for e in edges if gene_to_family.get(e.source) == gene_to_family.get(e.target)
    )

    weight_mode = _edge_weight_mode(config)
    identities = [e.identity for e in edges]
    coverages = [e.coverage for e in edges]
    weights = [e.weight(weight_mode) for e in edges]

    ssn_stats = {
        "n_nodes": ssn_nodes,
        "n_edges": n_ssn_edges,
        "n_components": len(components),
        "density": _graph_density(ssn_nodes, n_ssn_edges),
        "n_isolated_nodes": n_isolated,
        "cross_genome_edge_fraction": _ratio(n_cross_genome_edges, n_ssn_edges),
        "same_family_edge_fraction": _ratio(n_same_family_edges, n_ssn_edges),
        "degree": {
            "min": degree_values[0] if degree_values else 0,
            "median": median_degree,
            "max": degree_values[-1] if degree_values else 0,
        },
        "identity": _min_median_max(identities),
        "coverage": _min_median_max(coverages),
        "weight": _min_median_max(weights),
        "edge_completion_mode": ssn.edge_completion_mode if ssn is not None else None,
        "family_mode": config.get("family_mode", "full_graph"),
    }

    # ---- families 聚合 ----
    family_sizes = sorted((fam.size for fam in families), reverse=True)
    singleton_families = sum(1 for s in family_sizes if s == 1)
    families_stats = {
        "n_families": len(families),
        "mean_size": _mean(family_sizes),
        "median_size": _median(family_sizes),
        "max_size": family_sizes[0] if family_sizes else 0,
        "n_singleton": singleton_families,
    }

    # ---- MCL 前后对比 ----
    component_sizes = sorted((len(c) for c in components.values()), reverse=True)
    comp_to_families: Dict[int, Set[str]] = {}
    for fam in families:
        comp_to_families.setdefault(fam.component_id, set()).add(fam.family_id)
    n_split = sum(1 for fams in comp_to_families.values() if len(fams) > 1)
    max_subfamilies = max((len(fams) for fams in comp_to_families.values()), default=0)

    mcl_stats = {
        "inflation": config.get("inflation", 1.5),
        "edge_weight_mode": weight_mode,
        "n_components_before": n_mcl_components,
        "component_size": {
            "mean": _mean(component_sizes),
            "median": _median(component_sizes),
            "max": component_sizes[0] if component_sizes else 0,
        },
        "n_families_after": len(families),
        "family_size": {
            "mean": _mean(family_sizes),
            "median": _median(family_sizes),
            "max": family_sizes[0] if family_sizes else 0,
        },
        "n_components_split": n_split,
        "max_subfamilies_per_component": max_subfamilies,
    }

    # ---- 泛基因组 ----
    pan_data = build_pan_genome_data(families, gene_to_genome, genome_ids)
    family_category = pan_data["family_category"]
    presence_absence = pan_data["presence_absence"]
    genome_breakdown = pan_data["genome_breakdown"]
    n_single_copy_core = pan_data["n_single_copy_core"]

    pan_counts = {c: 0 for c in PAN_CATEGORIES}
    for cat in family_category.values():
        pan_counts[cat] += 1
    total_families = len(families)
    pan_genome = {
        "categories": {
            cat: {
                "n_families": pan_counts[cat],
                "family_pct": _ratio(pan_counts[cat], total_families),
                "n_genes": sum(genome_breakdown[gid][cat] for gid in genome_ids),
            }
            for cat in PAN_CATEGORIES
        },
        "n_single_copy_core": n_single_copy_core,
    }
    for cat in PAN_CATEGORIES:
        pan_genome["categories"][cat]["gene_pct"] = _ratio(
            pan_genome["categories"][cat]["n_genes"], n_genes
        )

    # ---- 每个基因组概况 ----
    genome_stats: Dict[str, Any] = {}
    for gid in genome_ids:
        g_genes = genomes[gid].genes
        fam_ids = {
            gene_to_family[g.gene_id]
            for g in g_genes
            if g.gene_id in gene_to_family
        }
        n_singleton_genes = sum(
            1 for g in g_genes if families_size_1.get(gene_to_family.get(g.gene_id))
        )
        genome_stats[gid] = {
            "n_genes": len(g_genes),
            "n_families": len(fam_ids),
            "n_singleton_genes": n_singleton_genes,
        }

    # ---- top families / top components ----
    fam_by_size = sorted(families, key=lambda f: (-f.size, f.family_id))
    top_families = []
    for fam in fam_by_size[:30]:
        gids = {gene_to_genome.get(m) for m in fam.member_genes} - {None}
        members = sorted(fam.member_genes)
        top_families.append({
            "family_id": fam.family_id,
            "size": fam.size,
            "n_genomes": len(gids),
            "edge_density": fam.edge_density,
            "category": family_category.get(fam.family_id, "cloud"),
            "representative": members[0] if members else None,
        })

    comp_by_size = sorted(components.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    top_components = []
    for cid, members in comp_by_size[:30]:
        fam_ids = {gene_to_family.get(m) for m in members} - {None}
        top_components.append({
            "component_id": cid,
            "n_genes": len(members),
            "n_families": len(fam_ids),
        })

    # ---- 簇大小分布(分箱) ----
    cluster_size_distribution = _size_distribution(family_sizes)

    # ---- derived ----
    largest_component_size = component_sizes[0] if component_sizes else 0
    derived = {
        "avg_genes_per_genome": _ratio(n_genes, n_genomes),
        "avg_local_reps_per_genome": _ratio(n_local_representatives, n_genomes),
        "avg_global_reps_per_cluster": _ratio(global_representatives, global_clusters),
        "avg_family_size": _ratio(n_genes_with_family, len(families)),
        "largest_component_size": largest_component_size,
        "largest_component_fraction": _ratio(largest_component_size, n_genes),
        "largest_family_size": family_sizes[0] if family_sizes else 0,
        "median_degree": median_degree,
        "isolated_gene_fraction": _ratio(n_isolated, n_genes),
        "cross_genome_edge_fraction": _ratio(n_cross_genome_edges, n_ssn_edges),
    }

    return ReportData(
        counts=counts,
        derived=derived,
        ssn=ssn_stats,
        families=families_stats,
        pan_genome=pan_genome,
        mcl=mcl_stats,
        cluster_size_distribution=cluster_size_distribution,
        genomes=genome_stats,
        top_families=top_families,
        top_components=top_components,
        family_category=family_category,
        presence_absence=presence_absence,
        genome_breakdown=genome_breakdown,
        n_single_copy_core=n_single_copy_core,
        comp_of_gene=comp_of_gene,
        components=components,
        gene_to_genome=gene_to_genome,
        gene_to_family=gene_to_family,
        degree=degree,
        representative_genes=representative_genes,
    )


# ---------------------------------------------------------------------------
# 数值辅助
# ---------------------------------------------------------------------------
def _ratio(numerator: float, denominator: float) -> float:
    return round(numerator / denominator, 6) if denominator else 0.0


def _mean(values: List[float]) -> float:
    return round(sum(values) / len(values), 6) if values else 0.0


def _median(values: List[float]) -> float:
    if not values:
        return 0.0
    n = len(values)
    if n % 2 == 1:
        return values[n // 2]
    return (values[n // 2 - 1] + values[n // 2]) / 2


def _min_median_max(values: List[float]) -> Dict[str, float]:
    if not values:
        return {"min": 0.0, "median": 0.0, "max": 0.0}
    return {"min": min(values), "median": _median(sorted(values)), "max": max(values)}


def _graph_density(n_nodes: int, n_edges: int) -> float:
    if n_nodes < 2:
        return 0.0
    return round(2 * n_edges / (n_nodes * (n_nodes - 1)), 6)


def _size_distribution(sizes: List[int]) -> Dict[str, Any]:
    """按 [1,2,3-5,6-10,11-20,21-50,51-100,101-500,501+] 分箱计数。"""
    bins = [
        ("1", 1, 1), ("2", 2, 2), ("3-5", 3, 5), ("6-10", 6, 10),
        ("11-20", 11, 20), ("21-50", 21, 50), ("51-100", 51, 100),
        ("101-500", 101, 500), ("501+", 501, None),
    ]
    result: Dict[str, Any] = {"bins": []}
    for label, lo, hi in bins:
        n = sum(1 for s in sizes if s >= lo and (hi is None or s <= hi))
        result["bins"].append({"label": label, "n": n})
    result["total"] = len(sizes)
    return result


# ---------------------------------------------------------------------------
# TSV 导出
# ---------------------------------------------------------------------------
def _write_tsv(path: str, header: List[str], rows: List[List[Any]]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(str(x) for x in header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")


def _original(gene_id: str, state: LoadedState) -> str:
    return state.gene_id_to_original.get(gene_id, gene_id)


def write_tsvs(state: LoadedState, data: ReportData, report_dir: str) -> List[str]:
    """把所有 TSV 写入 report_dir, 返回写入的文件名列表(按写入顺序)。"""
    os.makedirs(report_dir, exist_ok=True)
    written: List[str] = []

    genome_ids = sorted(state.genomes.keys())
    genome_names = [state.genome_id_to_name.get(gid, gid) for gid in genome_ids]

    # 1. presence_absence_matrix.tsv(行=family, 列=基因组, 值=拷贝数)
    rows = []
    for fam in sorted(state.families, key=lambda f: f.family_id):
        counts = data.presence_absence.get(fam.family_id, {})
        rows.append([fam.family_id] + [counts.get(gid, 0) for gid in genome_ids])
    _write_tsv(os.path.join(report_dir, "presence_absence_matrix.tsv"),
               ["Orthogroup"] + genome_names, rows)
    written.append("presence_absence_matrix.tsv")

    # 2. per_genome_breakdown.tsv(每行一个基因组, 列为五类基因数)
    rows = []
    for gid in genome_ids:
        rows.append([state.genome_id_to_name.get(gid, gid)]
                    + [data.genome_breakdown.get(gid, {}).get(cat, 0) for cat in PAN_CATEGORIES])
    _write_tsv(os.path.join(report_dir, "per_genome_breakdown.tsv"),
               ["Genome"] + PAN_CATEGORIES, rows)
    written.append("per_genome_breakdown.tsv")

    # 3. pan_genome_summary.tsv(五类 + single-copy-core)
    total = len(state.families)
    rows = []
    for cat in PAN_CATEGORIES:
        n = data.pan_genome["categories"][cat]["n_families"]
        pct = (n / total * 100) if total else 0.0
        rows.append([cat, n, f"{pct:.2f}%"])
    n_core = data.pan_genome["categories"]["core"]["n_families"]
    sc_pct = (data.n_single_copy_core / n_core * 100) if n_core else 0.0
    rows.append(["single-copy-core", data.n_single_copy_core, f"{sc_pct:.2f}%"])
    _write_tsv(os.path.join(report_dir, "pan_genome_summary.tsv"),
               ["Category", "N_orthogroups", "Percentage"], rows)
    written.append("pan_genome_summary.tsv")

    # 4. ssn_nodes.tsv(全部基因节点, 含 degree/代表标记/长度)
    rows = []
    for genome in sorted(state.genomes.values(), key=lambda g: g.genome_id):
        for gene in sorted(genome.genes, key=lambda g: g.gene_id):
            gid = gene.gene_id
            rows.append([
                _original(gid, state),
                data.gene_to_family.get(gid, ""),
                data.comp_of_gene.get(gid, ""),
                data.gene_to_genome.get(gid, ""),
                data.degree.get(gid, 0),
                "yes" if gid in data.representative_genes else "no",
                gene.length,
            ])
    _write_tsv(os.path.join(report_dir, "ssn_nodes.tsv"),
               ["gene_id", "family_id", "component_id", "genome_id",
                "degree", "is_representative", "length"], rows)
    written.append("ssn_nodes.tsv")

    # 5. ssn_edges.tsv
    rows = []
    if state.complete_ssn is not None:
        weight_mode = _edge_weight_mode(state.config)
        for e in sorted(state.complete_ssn.edges, key=lambda e: (e.source, e.target)):
            same = data.gene_to_family.get(e.source) == data.gene_to_family.get(e.target)
            rows.append([
                _original(e.source, state), _original(e.target, state),
                e.identity, e.coverage, e.bitscore, e.weight(weight_mode),
                "yes" if same else "no",
            ])
    _write_tsv(os.path.join(report_dir, "ssn_edges.tsv"),
               ["source", "target", "identity", "coverage", "bitscore",
                "weight", "same_family"], rows)
    written.append("ssn_edges.tsv")

    # 6. family_members.tsv
    rows = []
    for fam in sorted(state.families, key=lambda f: f.family_id):
        for gid in sorted(fam.member_genes):
            rows.append([fam.family_id, _original(gid, state),
                         data.gene_to_genome.get(gid, "")])
    _write_tsv(os.path.join(report_dir, "family_members.tsv"),
               ["family_id", "gene_id", "genome_id"], rows)
    written.append("family_members.tsv")

    # 7. family_stats.tsv
    rows = []
    for fam in sorted(state.families, key=lambda f: f.family_id):
        gids = {data.gene_to_genome.get(m) for m in fam.member_genes} - {None}
        rows.append([fam.family_id, fam.size, len(gids), fam.edge_density,
                     data.family_category.get(fam.family_id, "cloud")])
    _write_tsv(os.path.join(report_dir, "family_stats.tsv"),
               ["family_id", "size", "n_genomes", "edge_density", "category"], rows)
    written.append("family_stats.tsv")

    # 8. mcl_components.tsv(MCL 前连通分量及亚家族拆分)
    rows = []
    comp_to_fams: Dict[int, Set[str]] = {}
    for fam in state.families:
        comp_to_fams.setdefault(fam.component_id, set()).add(fam.family_id)
    for cid in sorted(data.components.keys()):
        fams = comp_to_fams.get(cid, set())
        rows.append([cid, len(data.components[cid]), len(fams),
                     ",".join(sorted(fams))])
    _write_tsv(os.path.join(report_dir, "mcl_components.tsv"),
               ["component_id", "n_genes", "n_families", "family_ids"], rows)
    written.append("mcl_components.tsv")

    # 9. cluster_size_distribution.tsv
    rows = [[b["label"], b["n"]] for b in data.cluster_size_distribution["bins"]]
    _write_tsv(os.path.join(report_dir, "cluster_size_distribution.tsv"),
               ["bin", "n_families"], rows)
    written.append("cluster_size_distribution.tsv")

    return written
