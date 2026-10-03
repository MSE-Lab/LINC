"""FamilyRefineStage: MCL 亚家族细分阶段(对应 v1 Step6)。

设计依据: docs/归档/STEP6_MCL亚家族细分.md + 2026-10-02 与用户确认的方案:
    1. 输入是 SSNBuilderStage 产出的 `CompleteSSN`(无论
       `edge_completion_mode` 是什么都要经过该阶段, 即使 `none` 模式下
       `CompleteSSN` 里也只有代表-代表边, 仍然是合法的、本阶段的输入)。
    2. 先把 `CompleteSSN.edges` 构建成一张 networkx 图(边权重用
       `SSNEdge.weight(edge_weight_mode, edge_weight_alpha)` 惰性计算,
       延续"存储原始值、按需计算权重"的既有设计), 按连通分量拆分
       (一个连通分量 = 一个"超家族", 对应 `OrthologFamily.component_id`)。
    3. 对每个连通分量分别调用外部 `mcl` 命令做细分(而不是整图一次性调用):
       - 连通分量只有 1 个基因(孤立节点, 无边): 直接作为一个 size=1 的
         独立家族, 不调用 mcl(mcl 对无边输入没有意义, 与用户确认)。
       - 其余连通分量: 调用 `ClusteringPort.cluster(subgraph, inflation=...)`
         (生产环境传入 `MCLClustering`), 得到的每个子簇生成一个
         `OrthologFamily`。
       - 注意(2026-10-02 修复): 图的节点不能只靠 `CompleteSSN.edges` 隐式
         添加 —— 若某基因在局部簇内是孤儿(无其他成员)且全局比对中也没有
         找到任何可连的代表, 它在 edges 里完全不出现, 若不显式补充节点,
         会被 `nx.connected_components()` 彻底忽略、从最终结果中静默消失
         (而不是被当作一个合法的 size=1 独立家族保留)。`run()` 因此新增
         `all_gene_ids` 参数, 用于在构建图后把全部基因显式 `add_node()`
         补齐。
    4. `edge_density` 按标准图密度公式计算(对每个家族的诱导子图):
       `2 * E / (N * (N-1))`, N<=1 时为 0.0(与用户确认, 2026-10-02)。
    5. `family_id` 采用 OrthoFinder 风格编号 `OG0000000`(与用户确认),
       按 `(component_id, 成员排序)` 排序后依次分配, 以保证跨运行确定性
       (mcl/连通分量遍历顺序本身不保证跨运行一致, 排序后消除该非确定性,
       延续项目一贯的"结果排序以消除非确定性来源"原则)。
"""
import logging
from typing import Iterable, List, Optional

import networkx as nx

from ..core.entities import CompleteSSN, OrthologFamily
from ..ports.clustering_port import ClusteringPort

logger = logging.getLogger(__name__)


class FamilyRefineStage:
    """对应 v1 Step6 的 MCL 亚家族细分实现。"""

    def __init__(
        self,
        clustering: ClusteringPort,
        inflation: float = 1.5,
        edge_weight_mode: str = "identity_x_coverage",
        edge_weight_alpha: float = 0.5,
        cleanup_intermediate: bool = True,
        threads: int = 0,
    ) -> None:
        self._clustering = clustering
        self._inflation = inflation
        self._edge_weight_mode = edge_weight_mode
        self._edge_weight_alpha = edge_weight_alpha
        # 数据管理优化第3项(2026-10-03 与用户确认): mcl 落盘的中间文件
        # (component_*/) 结果已转存进 OrthologFamily, 默认用完即删。
        self._cleanup_intermediate = cleanup_intermediate
        # 工程审查问题9(2026-10-03 与用户确认): 配合"全图一次mcl调用"的
        # 优化, 通过 mcl 原生 `-te` 参数利用多核(0 表示不显式指定, 与
        # CLI `--threads 0` 既有语义对齐)。
        self._threads = threads

    def run(
        self,
        complete_ssn: CompleteSSN,
        all_gene_ids: Optional[Iterable[str]] = None,
    ) -> List[OrthologFamily]:
        graph = self._build_graph(complete_ssn, all_gene_ids)

        # 按最小节点排序连通分量, 保证遍历顺序跨运行确定(graph 的节点/边
        # 插入顺序已经按排序后的边构建, 这里再显式排序进一步保证稳妥)。
        components = sorted(
            (sorted(component) for component in nx.connected_components(graph)),
            key=lambda nodes: nodes[0],
        )

        # 工程审查问题9(2026-10-03 与用户确认): 把"每个连通分量一次mcl
        # 子进程"合并为"全图一次mcl调用", 调用次数从 O(分量数) 降到 O(1)。
        # mcl 对不连通的子图天然不会产生交叉簇(markov转移概率在分量间为0),
        # 因此合并调用的结果与逐分量调用在数学上是等价的, 只是子进程启动
        # 开销被摊薄, 并可通过 `-te` 利用多核(由 MCLClustering.cluster
        # 的 threads 参数传入)。

        family_specs = []  # List[Tuple[component_id, member_ids: set]]
        gene_to_component_id = {}
        non_trivial_nodes = set()
        for component_id, component_nodes in enumerate(components, start=1):
            if len(component_nodes) == 1:
                # 孤立节点(无边的连通分量): 不调用 mcl, 直接作为独立家族
                # (与用户确认, 2026-10-02)。
                family_specs.append((component_id, set(component_nodes)))
                continue
            for gene_id in component_nodes:
                gene_to_component_id[gene_id] = component_id
            non_trivial_nodes.update(component_nodes)

        if non_trivial_nodes:
            combined_subgraph = graph.subgraph(non_trivial_nodes)
            index_id = "all_components"
            sub_clusters = self._clustering.cluster(
                combined_subgraph,
                inflation=self._inflation,
                index_id=index_id,
                threads=self._threads,
            )
            if self._cleanup_intermediate:
                self._clustering.cleanup(index_id=index_id)
            for member_ids in sub_clusters:
                member_ids = set(member_ids)
                # 每个 sub_cluster 必然完全落在某一个原始连通分量内部
                # (mcl 不会产生跨分量簇), 任取一个成员反查其 component_id
                # 即可, 以保持 OrthologFamily.component_id 的既有语义不变。
                sample_gene = next(iter(member_ids))
                component_id = gene_to_component_id[sample_gene]
                family_specs.append((component_id, member_ids))

            # 工程审查问题5(2026-10-03 与用户确认): 防御性校验 mcl 返回的
            # 分区是否覆盖了全部非孤立节点。外部工具异常(退化输入、版本
            # 差异等)可能返回不完整的分区, 若不校验, 未被任何子簇包含的
            # 基因会在后续生成 family 时被静默丢弃, 从最终结果中彻底消失
            # 而不报任何错误。这里采取"兜底补齐"而非直接抛异常中断整个
            # pipeline(因为这是外部工具输出异常, 不是用户输入错误); 缺失
            # 的基因各自作为独立家族(size=1)补回, 并打印 warning 便于
            # 事后排查(按其原始 component_id 归属, 而非笼统归为一个分量)。
            covered = set().union(*sub_clusters) if sub_clusters else set()
            missing = non_trivial_nodes - covered
            if missing:
                logger.warning(
                    "FamilyRefineStage: mcl 对全图(共 %d 个非孤立节点)返回的"
                    "分区未覆盖全部节点, 缺失 %d 个基因: %s; 已将这些基因各自"
                    "作为独立家族补回, 避免静默丢失(可能是 mcl 异常退出或"
                    "退化输入导致)。",
                    len(non_trivial_nodes), len(missing), sorted(missing),
                )
                for gene_id in sorted(missing):
                    family_specs.append((gene_to_component_id[gene_id], {gene_id}))

        # 按 (component_id, 排序后的成员列表) 排序, 保证 family_id 分配顺序
        # 跨运行确定(mcl 输出簇的顺序/ dict 遍历顺序不保证跨运行一致)。
        family_specs.sort(key=lambda spec: (spec[0], sorted(spec[1])))

        families: List[OrthologFamily] = []
        for idx, (component_id, member_ids) in enumerate(family_specs):
            families.append(
                OrthologFamily(
                    family_id=f"OG{idx:07d}",
                    component_id=component_id,
                    member_genes=member_ids,
                    edge_density=self._edge_density(member_ids, graph),
                )
            )

        logger.info(
            "FamilyRefineStage 完成: inflation=%.2f, edge_weight_mode=%s, "
            "共 %d 个连通分量, %d 个直系同源家族",
            self._inflation, self._edge_weight_mode, len(components), len(families),
        )
        return families

    def _build_graph(
        self,
        complete_ssn: CompleteSSN,
        all_gene_ids: Optional[Iterable[str]] = None,
    ) -> nx.Graph:
        # 按 (source, target) 排序后再插入, 消除 CompleteSSN.edges(set)
        # 遍历顺序的非确定性(与 mmseqs/diamond adapter 结果排序的既有原则
        # 一致), 保证后续连通分量/节点插入顺序跨运行确定。
        sorted_edges = sorted(
            complete_ssn.edges,
            key=lambda e: tuple(sorted((e.source, e.target))),
        )
        graph = nx.Graph()
        for edge in sorted_edges:
            graph.add_edge(
                edge.source, edge.target,
                weight=edge.weight(self._edge_weight_mode, self._edge_weight_alpha),
            )
        # 显式补全全部基因节点(2026-10-02 修复): 零边的孤立基因(局部无其他
        # 成员、全局也未找到任何可连代表)不会通过 add_edge 被加入图中,
        # 若不补齐会在 nx.connected_components() 中彻底消失。按排序后的
        # gene_id 依次 add_node, 保证节点插入顺序跨运行确定。
        if all_gene_ids is not None:
            for gene_id in sorted(all_gene_ids):
                if gene_id not in graph:
                    graph.add_node(gene_id)
        return graph

    @staticmethod
    def _edge_density(member_ids, graph: nx.Graph) -> float:
        n = len(member_ids)
        if n <= 1:
            return 0.0
        subgraph = graph.subgraph(member_ids)
        e = subgraph.number_of_edges()
        return 2 * e / (n * (n - 1))
