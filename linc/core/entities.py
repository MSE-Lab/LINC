"""领域实体定义。

对应 docs/LINC_技术框架设计_v2.md 第 2.1 节,是 v1 `models/data_structures.py`
的现代化重写版本:字段语义保持一致,但引入了 frozen=True(不可变值对象)、
property(派生字段如 coverage/size)等现代 dataclass 写法。
"""
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Set


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
    """一次两两比对的结果记录(来自 MMseqs2/DIAMOND 等比对引擎)。"""

    query_id: str
    subject_id: str
    identity: float
    query_cover: float
    subject_cover: float
    evalue: float
    bitscore: float
    alignment_length: int

    @property
    def coverage(self) -> float:
        """取 query/subject 覆盖度中较小的一个,作为保守覆盖度估计。"""
        return min(self.query_cover, self.subject_cover)


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
    "基因对"去重,而 identity/coverage 是可能因重新计算产生浮点误差的派生值,
    不应参与相等性/哈希判断。因此本类关闭 dataclass 默认生成的 __eq__,
    改为自定义 __eq__/__hash__: 仅依据 frozenset({source, target}) 判定两条边
    是否代表同一对基因,与字段顺序、identity/coverage 的具体数值无关。
    这样 SSNEdge 实例既可以安全地放入 set/用作 dict key(满足
    RepositoryPort.save_ssn_edges(edges: set) 的使用场景),又保留了
    identity/coverage 可变、可原地更新的能力(未设为 frozen)。
    """

    source: str
    target: str
    identity: float
    coverage: float

    def weight(self, mode: str = "identity_x_coverage") -> float:
        """根据指定模式计算边权重。

        Args:
            mode: "identity_only" | "coverage_only" | "identity_x_coverage"(默认)。
        """
        if mode == "identity_only":
            return self.identity
        if mode == "coverage_only":
            return self.coverage
        return self.identity * self.coverage

    def _pair_key(self) -> frozenset:
        return frozenset((self.source, self.target))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SSNEdge):
            return NotImplemented
        return self._pair_key() == other._pair_key()

    def __hash__(self) -> int:
        return hash(self._pair_key())


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
