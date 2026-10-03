"""领域实体定义。

对应 docs/LINC_技术框架设计_v2.md 第 2.1 节,是 v1 `models/data_structures.py`
的现代化重写版本:字段语义保持一致,但引入了 frozen=True(不可变值对象)、
property(派生字段如 coverage/size)等现代 dataclass 写法。
"""
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Set


class GenomeStatus(str, Enum):
    """基因组加载/处理状态。"""

    PENDING = "pending"
    LOADED = "loaded"
    FAILED = "failed"


@dataclass(frozen=True)
class Gene:
    """基因实体。

    gene_id 格式约定为 "{genome_id}|{local_id}",保证跨基因组全局唯一。
    """

    gene_id: str
    genome_id: str
    length: int
    sequence: Optional[str] = None   # 允许延迟加载,大规模场景下不常驻内存


@dataclass
class Genome:
    """基因组实体,持有该基因组下的全部基因列表。"""

    genome_id: str
    file_path: str
    genes: List[Gene] = field(default_factory=list)
    status: Optional[GenomeStatus] = None


@dataclass
class AlignmentResult:
    """一次两两比对的结果记录(来自 MMseqs2/DIAMOND 等比对引擎)。

    接口契约(各 AlignerPort 实现必须遵守): identity/query_cover/subject_cover
    统一为 [0, 1] 范围的小数, 而非百分数(0-100)。比对工具原始输出若为百分数
    (如 diamond 的 pident/qcovhsp/scovhsp),adapter 必须在解析阶段自行 ÷100
    换算后再构造本类实例,不允许把百分数原样传出,以保证不同后端(mmseqs/
    diamond/mock/未来的 blastp 等)之间的数值可比、且与 cutoff_i/cutoff_ii
    等按 [0, 1] 定义的阈值语义一致。
    """

    query_id: str
    subject_id: str
    identity: float
    query_cover: float
    subject_cover: float
    evalue: float
    bitscore: float
    alignment_length: int
    # query/subject 序列全长(氨基酸数), 直接取自比对工具输出(mmseqs qlen/tlen,
    # diamond qlen/slen), 而非重新查 Gene.length, 保证与该次比对实际使用的
    # 序列长度口径一致。默认 0 表示未知(如历史数据/未升级的测试构造),
    # 为未来 MCL 权重归一化(如 OrthoFinder 式长度校正)预留(2026-10-02 与用户确认)。
    query_length: int = 0
    subject_length: int = 0

    @property
    def coverage(self) -> float:
        """取 query/subject 覆盖度中较小的一个,作为保守覆盖度估计。

        等价于 `coverage_for_mode(self, cov_mode=0)`,保留此 property 是为了
        向后兼容既有调用方;新代码如需按 cov_mode 可配置地取 coverage,
        请使用 `coverage_for_mode`。
        """
        return min(self.query_cover, self.subject_cover)


# coverage 合成方式(参考 mmseqs2 的 --cov-mode 语义,与用户确认,2026-10-02):
# 0 = 双向(取 query_cover/subject_cover 中较小者,两边都要达标才算数, 默认);
# 1 = 仅看 subject_cover(target coverage);
# 2 = 仅看 query_cover(query coverage)。
COV_MODE_BIDIRECTIONAL = 0
COV_MODE_SUBJECT = 1
COV_MODE_QUERY = 2
VALID_COV_MODES = (COV_MODE_BIDIRECTIONAL, COV_MODE_SUBJECT, COV_MODE_QUERY)


def coverage_for_mode(hit: "AlignmentResult", cov_mode: int = COV_MODE_BIDIRECTIONAL) -> float:
    """按 cov_mode 从一次比对结果中取出用于与 cutoff_c 比较的覆盖度数值。"""
    if cov_mode == COV_MODE_SUBJECT:
        return hit.subject_cover
    if cov_mode == COV_MODE_QUERY:
        return hit.query_cover
    return min(hit.query_cover, hit.subject_cover)


def passes_inclusion_cutoff(
    hit: "AlignmentResult",
    cutoff_i: float,
    cutoff_c: float,
    cov_mode: int = COV_MODE_BIDIRECTIONAL,
) -> bool:
    """统一的"纳入判据": identity >= cutoff_i 且 coverage(按cov_mode) >= cutoff_c。

    IntraGenomeStage/GlobalJoinStage 均应调用本函数判定一次比对是否成簇,
    避免两阶段各自实现而出现判据不一致(2026-10-02 与用户确认修复高危2)。
    """
    return hit.identity >= cutoff_i and coverage_for_mode(hit, cov_mode) >= cutoff_c


@dataclass
class Cluster:
    """基因组内或全局代表库中的一个簇。"""

    cluster_id: str
    representative_ids: Set[str] = field(default_factory=set)
    member_ids: Set[str] = field(default_factory=set)
    active: bool = True   # 被合并后的旧簇标记为 False,而非物理删除(便于追溯)


@dataclass(eq=False)
class SSNEdge:
    """序列相似性网络(SSN)中的一条边。

    SSN 在生物学语义上是无向图(A-B 相似 = B-A 相似),且补边阶段需要按
    "基因对"去重,而 identity/coverage/bitscore 是可能因重新计算产生浮点误差的
    派生值,不应参与相等性/哈希判断。因此本类关闭 dataclass 默认生成的 __eq__,
    改为自定义 __eq__/__hash__: 仅依据 frozenset({source, target}) 判定两条边
    是否代表同一对基因,与字段顺序、identity/coverage/bitscore 的具体数值无关。
    这样 SSNEdge 实例既可以安全地放入 set/用作 dict key(满足
    RepositoryPort.save_ssn_edges(edges: set) 的使用场景),又保留了
    identity/coverage/bitscore 可变、可原地更新的能力(未设为 frozen)。

    权重采用"存储原始度量值 + 按需惰性计算"的设计(2026-10-02 与用户确认):
    identity/coverage/bitscore 三个原始值随边一起存储,不固化成单一权重数字,
    便于后续(如 MCL 细分阶段)按需切换权重公式而无需重新构建 SSN。
    """

    source: str
    target: str
    identity: float
    coverage: float
    bitscore: float = 0.0
    # 与 AlignmentResult.query_length/subject_length 对应, 从产生该边的比对结果
    # 转存而来, 为未来 MCL 权重归一化(如 OrthoFinder 式长度校正)预留,
    # 默认 0 表示未知(2026-10-02 与用户确认)。source/target 对应 query/subject
    # 的先后关系由构造方保证一致, 不在本类内部做强制校验。
    #
    # bitscore 方向约定说明(第二轮漏洞7, 2026-10-03 与用户确认, 文档澄清、
    # 不改动逻辑): bitscore 是有方向性的量(query→subject 与 subject→query
    # 的比对结果数值可能不同)。SSNBuilderStage 的 full 模式在
    # _complete_cross_genome_member_edges() 中跨基因组补边时, 只按
    # genome_id 字典序固定一个方向发起比对(genome_a < genome_b 时以
    # genome_a 的基因作 query、genome_b 的基因作 subject), 不会反向再比对
    # 一次。因此 source/target 对应的 query/subject 身份(以及由此而来的
    # bitscore、source_length/target_length)由"谁的 genome_id 字典序更小"
    # 决定, 是一个人为约定而非语义上的"真实"query/subject 关系。使用
    # --edge-weight-mode bitscore 时, 边权重会按此约定确定性产生(同一输入
    # 多次运行结果一致), 但数值本身受"谁排前"影响, 不代表双向比对的对称值。
    source_length: int = 0
    target_length: int = 0

    def weight(self, mode: str = "identity_x_coverage", alpha: float = 0.5) -> float:
        """根据指定模式计算边权重。

        Args:
            mode: "identity"(含兼容别名 "identity_only") |
                "coverage"(含兼容别名 "coverage_only") |
                "identity_x_coverage"(默认, identity*coverage) |
                "identity_coverage_weighted"(alpha*identity + (1-alpha)*coverage) |
                "bitscore"(原始值, 不做归一化, 调用方需自行决定是否/如何归一化)。
            alpha: 仅 mode="identity_coverage_weighted" 时生效, 默认 0.5。
        """
        if mode in ("identity", "identity_only"):
            return self.identity
        if mode in ("coverage", "coverage_only"):
            return self.coverage
        if mode == "identity_coverage_weighted":
            return alpha * self.identity + (1 - alpha) * self.coverage
        if mode == "bitscore":
            return self.bitscore
        return self.identity * self.coverage

    def _pair_key(self) -> frozenset:
        return frozenset((self.source, self.target))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SSNEdge):
            return NotImplemented
        return self._pair_key() == other._pair_key()

    def __hash__(self) -> int:
        return hash(self._pair_key())


# SSNEdge.weight() 支持的合法 mode 取值(用于 CLI --edge-weight-mode 的参数
# 校验, 2026-10-02 与用户确认, 不包含向后兼容别名 identity_only/coverage_only,
# 对外只暴露这5个正式名称)。
VALID_EDGE_WEIGHT_MODES = (
    "identity",
    "coverage",
    "identity_x_coverage",
    "identity_coverage_weighted",
    "bitscore",
)


@dataclass
class LocalCluster:
    """基因组内自聚类(Step2)产生的一个局部簇。

    单成员簇(孤儿基因)的 representative_ids 与 member_ids 相同,is_orphan=True。
    多成员簇根据簇内最小 pairwise identity 是否达到 Cutoff-II(即任意两个
    成员是否都高度相似),representative_ids 可能只有 1 个(整体均质)或
    多个(存在亚型分化,由 multi_representative_clustering 贪心覆盖选出,
    数量不超过 max_per_cluster)。
    """

    cluster_id: str
    member_ids: Set[str] = field(default_factory=set)
    representative_ids: List[str] = field(default_factory=list)
    is_orphan: bool = False


@dataclass
class GenomeIntraResult:
    """单个基因组完成 Step2(基因组内自聚类)后的结果。"""

    genome_id: str
    local_clusters: List[LocalCluster] = field(default_factory=list)
    local_representatives: List[str] = field(default_factory=list)
    orphan_genes: List[str] = field(default_factory=list)
    n_genes: int = 0
    n_clusters: int = 0
    n_representatives: int = 0
    genome_index: Any = None   # AlignerPort.build_index 返回的持久化索引 handle,
                                 # 供 GlobalJoinStage 复用, 避免重复建库(见 v2 设计文档 2.4节)
    # 本基因组全对全比对的原始结果(不限于代表-成员, 也包括成员-成员),
    # 供 SSNBuilderStage 补边时直接复用, 避免对同一基因组内已经比对过的基因对
    # 重复发起比对(2026-10-02 与用户确认)。
    pairwise_hits: List["AlignmentResult"] = field(default_factory=list)


@dataclass
class GlobalRepresentativeLibrary:
    """全局代表库(Step3 增量并入阶段持有/更新的核心状态)。

    持有全局 Cluster 列表 + "代表基因 -> 簇" 的映射(仅代表基因参与映射,
    普通成员基因不在此追踪, 其归属关系体现在 Cluster.member_ids 里)。
    簇合并时旧簇标记 active=False 而非物理删除(复用 Cluster.active 字段,
    便于追溯, 与 IntraGenomeStage 的既有设计风格保持一致)。
    """

    clusters: List[Cluster] = field(default_factory=list)
    gene_to_cluster: Dict[str, str] = field(default_factory=dict)  # 代表基因id -> cluster_id
    _genes_by_id: Dict[str, "Gene"] = field(default_factory=dict, repr=False)
    # cluster_id -> Cluster 的索引, 与 self.clusters 保持同步(新增簇时一并写入),
    # 使 get_cluster 从 O(n) 线性扫描降为 O(1) 字典查找(2026-10-02 修复低危1,
    # global_join_stage._decide_and_apply 对每个命中都会调用一次 get_cluster,
    # 大样本下原实现会退化为 O(命中数 × 簇数))。
    _cluster_index: Dict[str, Cluster] = field(default_factory=dict, repr=False)

    def get_cluster(self, cluster_id: str) -> Optional[Cluster]:
        return self._cluster_index.get(cluster_id)

    def representative_genes(self) -> List["Gene"]:
        """返回当前所有存活簇的代表基因列表(作为下一次 query 使用)。"""
        # 注意: representative_ids 是 Set[str], Python 字符串哈希随机化
        # (PYTHONHASHSEED 每次进程启动随机)会导致 set 遍历顺序在不同运行间
        # 不一致,进而影响喂给 mmseqs 的 query 序列顺序,产生跨进程的非确定性
        # 结果差异;此处显式按 rep_id 字典序排序遍历,消除该非确定性来源。
        reps: List["Gene"] = []
        for cluster in self.clusters:
            if not cluster.active:
                continue
            for rep_id in sorted(cluster.representative_ids):
                gene = self._genes_by_id.get(rep_id)
                if gene is not None:
                    reps.append(gene)
        return reps

    def add_new_cluster(self, cluster_id: str, gene: "Gene") -> None:
        """命中 0 个簇: 以该代表为初始代表创建新簇。"""
        self._genes_by_id[gene.gene_id] = gene
        cluster = Cluster(
            cluster_id=cluster_id,
            representative_ids={gene.gene_id},
            member_ids={gene.gene_id},
        )
        self.clusters.append(cluster)
        self._cluster_index[cluster_id] = cluster
        self.gene_to_cluster[gene.gene_id] = cluster_id

    def add_representative(self, cluster_id: str, gene: "Gene") -> None:
        """命中 1 个簇、仅达 Cutoff-I: 增补为该簇的新代表。"""
        self._genes_by_id[gene.gene_id] = gene
        cluster = self.get_cluster(cluster_id)
        cluster.representative_ids.add(gene.gene_id)
        cluster.member_ids.add(gene.gene_id)
        self.gene_to_cluster[gene.gene_id] = cluster_id

    def add_member_only(self, cluster_id: str, gene: "Gene") -> None:
        """命中 1 个簇、达到 Cutoff-II: 沿用已有代表, 仅记录成员, 不新增代表。"""
        cluster = self.get_cluster(cluster_id)
        cluster.member_ids.add(gene.gene_id)

    def merge_clusters(
        self,
        cluster_ids: List[str],
        merged_id: str,
        new_representative: Optional["Gene"] = None,
    ) -> Cluster:
        """合并多个簇为一个新簇(旧簇标记 active=False, 不物理删除)。

        new_representative: 路径③"均仅达 Cutoff-I"时, 连接远缘家族的该代表
        本身也需要新增为合并后簇的代表, 一并传入即可。
        """
        representative_ids: Set[str] = set()
        member_ids: Set[str] = set()
        for cid in cluster_ids:
            cluster = self.get_cluster(cid)
            if cluster is None or not cluster.active:
                continue
            cluster.active = False
            representative_ids |= cluster.representative_ids
            member_ids |= cluster.member_ids
            for rep_id in cluster.representative_ids:
                self.gene_to_cluster[rep_id] = merged_id
        if new_representative is not None:
            self._genes_by_id[new_representative.gene_id] = new_representative
            representative_ids.add(new_representative.gene_id)
            member_ids.add(new_representative.gene_id)
            self.gene_to_cluster[new_representative.gene_id] = merged_id
        merged = Cluster(
            cluster_id=merged_id,
            representative_ids=representative_ids,
            member_ids=member_ids,
        )
        self.clusters.append(merged)
        self._cluster_index[merged_id] = merged
        return merged

    @property
    def n_active_clusters(self) -> int:
        return sum(1 for c in self.clusters if c.active)

    @property
    def n_representatives(self) -> int:
        return sum(len(c.representative_ids) for c in self.clusters if c.active)


@dataclass
class IncrementalResult:
    """Step3(增量并入全局代表库)的最终输出。

    `cluster_counter`(2026-10-02, 第8步持久化/--append 设计): 记录本次运行
    结束时 `GlobalJoinStage._cluster_counter` 的当前值, 供下一次 `--append`
    时作为 `initial_cluster_counter` 传入, 使新基因组产生的 cluster_id 能
    接续历史编号, 不与历史簇冲突。
    """

    global_lib: GlobalRepresentativeLibrary = field(default_factory=GlobalRepresentativeLibrary)
    ssn_edges: Set["SSNEdge"] = field(default_factory=set)
    processing_order: List[str] = field(default_factory=list)
    n_clusters_final: int = 0
    n_representatives_final: int = 0
    cluster_counter: int = 0


# 补边模式(Step4+5 合并, SSNBuilderStage, 与用户确认, 2026-10-02):
# none = 不补边, 仅保留代表-代表边(稀疏SSN);
# representative_member = 默认, 补"代表-成员"边, 复用 IntraGenomeStage 已有
#   结果, 不产生任何新比对;
# full = 在 representative_member 基础上, 补全簇内全部基因对(含成员-成员)
#   关系, 同基因组内部分依旧复用已有结果, 仅跨基因组的成员对才真正发起新比对。
EDGE_COMPLETION_NONE = "none"
EDGE_COMPLETION_REPRESENTATIVE_MEMBER = "representative_member"
EDGE_COMPLETION_FULL = "full"
VALID_EDGE_COMPLETION_MODES = (
    EDGE_COMPLETION_NONE,
    EDGE_COMPLETION_REPRESENTATIVE_MEMBER,
    EDGE_COMPLETION_FULL,
)


@dataclass
class CompleteSSN:
    """SSNBuilderStage(合并 Step4 稀疏SSN构建 + Step5 最终补边)的最终输出。

    与用户确认(2026-10-02): 对于 full 模式下规模过大的簇("巨簇"), 本阶段仅
    做信息提醒(日志 warning), 不做自动降级或拆分 —— 避免悄悄产出有损/近似的
    SSN, 把"是否要为此付出更长运行时间"的决定权交还给用户。
    """

    edges: Set["SSNEdge"] = field(default_factory=set)
    edge_completion_mode: str = EDGE_COMPLETION_REPRESENTATIVE_MEMBER


@dataclass
class OrthologFamily:
    """最终输出的直系同源基因家族(由 MCL 细分 SSN 连通分量得到)。"""

    family_id: str
    component_id: int
    member_genes: Set[str] = field(default_factory=set)
    edge_density: float = 0.0

    @property
    def size(self) -> int:
        return len(self.member_genes)
