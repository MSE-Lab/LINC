"""折叠式家族聚类(family_mode=collapsed)的折叠映射计算。

设计依据: docs/LINC_折叠式家族聚类设计方案.md(Q1~Q4 已由用户确认)。

折叠对象:
    1. 全局成员: GlobalJoin 中以"仅成员"身份并入全局簇的局部代表, 且在
       `ssn_edges` 中与同簇某全局代表存在 identity >= cutoff_ii 且
       coverage >= cutoff_c 的边(高置信依附)。依附到 identity 最高的
       代表(并列时取 gene_id 字典序最小者, 保证确定性)。
       仅因代表数达上限而成为成员、但无高置信边的基因不折叠。
    2. 局部成员: 其所在局部簇的全部局部代表都已被折叠时, 该局部成员随
       局部代表一起折叠, 依附到其"最相似的已折叠局部代表"的锚点代表
       (无 pairwise 记录时取字典序最小的局部代表)。

返回 `{被折叠基因 -> 锚点全局代表}`。锚点必为未被折叠的全局代表,
因此不存在依附链。
"""
from typing import Dict, Optional

from ..core.entities import GenomeIntraResult, IncrementalResult


def compute_collapse_map(
    join_result: IncrementalResult,
    intra_results: Dict[str, GenomeIntraResult],
    cutoff_ii: float,
    cutoff_c: float,
) -> Dict[str, str]:
    gene_to_active_cluster: Dict[str, str] = {}
    global_reps = set()
    members_only = set()
    for cluster in join_result.global_lib.clusters:
        if not cluster.active:
            continue
        global_reps |= cluster.representative_ids
        members_only |= cluster.member_ids - cluster.representative_ids
        for gene_id in cluster.member_ids:
            gene_to_active_cluster[gene_id] = cluster.cluster_id

    # 全局成员 -> (最佳 identity, 锚点代表)
    best: Dict[str, tuple] = {}
    for edge in join_result.ssn_edges:
        if edge.identity < cutoff_ii or edge.coverage < cutoff_c:
            continue
        for member, rep in ((edge.source, edge.target), (edge.target, edge.source)):
            if member not in members_only or rep not in global_reps:
                continue
            if gene_to_active_cluster.get(member) != gene_to_active_cluster.get(rep):
                continue
            cand = (-edge.identity, rep)
            if member not in best or cand < best[member]:
                best[member] = cand
    collapse_map: Dict[str, str] = {member: cand[1] for member, cand in best.items()}

    # 局部成员: 所有局部代表均被折叠时才折叠
    for intra in intra_results.values():
        pair_identity: Optional[Dict[frozenset, float]] = None
        for local_cluster in intra.local_clusters:
            if local_cluster.is_orphan:
                continue
            reps = sorted(local_cluster.representative_ids)
            if not reps or any(r not in collapse_map for r in reps):
                continue
            if pair_identity is None:
                pair_identity = {}
                for hit in intra.pairwise_hits:
                    key = frozenset((hit.query_id, hit.subject_id))
                    if hit.identity > pair_identity.get(key, -1.0):
                        pair_identity[key] = hit.identity
            for member in sorted(local_cluster.member_ids - set(reps)):
                best_rep = max(
                    reps,
                    key=lambda r: (pair_identity.get(frozenset((r, member)), -1.0),
                                   [-ord(c) for c in r]),
                )
                collapse_map[member] = collapse_map[best_rep]
    return collapse_map
