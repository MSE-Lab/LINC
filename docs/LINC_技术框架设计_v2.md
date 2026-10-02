# ============================================================
# LINC 技术框架设计文档 v2 (现代化架构)
# ============================================================

> 本文档是 LINC 项目的**权威工程设计依据**。每次重新打开项目(含重启 PyCharm 后)应首先阅读本文档,
> 按文档中的模块划分、接口契约、开发顺序继续推进开发工作。
> 配套进度记录见 `docs/开发进度.md`,请在每次完成一个模块后同步更新进度文档。

---

## 0. 背景与目标

LINC 的核心目标:**以增量方式比较多个基因组,构建序列相似性网络(SSN),并通过 MCL 算法
实现蛋白质家族的快速、可复现聚类**。

v1 设计(`LINC_算法设计文档.md`)按"Step1~Step6"六个流程步骤平铺组织模块,存在以下问题:

- 比对引擎(MMseqs2/DIAMOND)与核心算法逻辑耦合,难以替换/测试
- 没有容错/断点续跑能力,大规模增量场景下一次失败需要全部重跑
- 缺少可观测性(除文本日志外无法追溯"某个基因为什么被合并")
- 配置缺乏类型安全校验

v2 设计采用**六边形架构(Ports & Adapters)**解决上述问题,详见下文。

---

## 1. 架构总览

```
linc/
├── core/            # 领域核心: 实体、事件、策略接口 (不依赖任何外部工具)
├── ports/            # 抽象接口 (Protocol): Aligner / Storage / Clustering
├── adapters/         # Port 的具体实现: mmseqs / diamond / mock / sqlite / mcl
├── engine/           # 增量引擎: 五个 Stage 串联核心算法
├── pipeline.py       # 声明式 Pipeline,支持 checkpoint/resume
├── config/           # 配置 Schema 与加载
└── cli.py            # 命令行入口
```

详细模块职责见第 2 节。

---

## 2. 模块详细设计

### 2.1 `core/` — 领域核心层

**原则**: 本层代码**禁止** import 任何第三方比对工具或存储实现,只能依赖标准库 + `numpy`/`networkx`
这类纯计算库。这是六边形架构的"内核",保证核心算法可以脱离具体基础设施被单独测试。

#### `core/entities.py` — 领域实体

```python
from dataclasses import dataclass, field
from typing import Set, Dict, List, Optional
from enum import Enum

@dataclass(frozen=True)
class Gene:
    gene_id: str          # 格式: "{genome_id}|{local_id}"
    genome_id: str
    length: int
    sequence: Optional[str] = None   # 允许延迟加载,大规模场景下不常驻内存

class GenomeStatus(str, Enum):
    PENDING = "pending"
    LOADED = "loaded"
    FAILED = "failed"

@dataclass
class Genome:
    genome_id: str
    file_path: str
    genes: List[Gene]
    status: Optional[GenomeStatus] = None

@dataclass
class AlignmentResult:
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
        return min(self.query_cover, self.subject_cover)

@dataclass
class Cluster:
    cluster_id: str
    representative_ids: Set[str]
    member_ids: Set[str]
    active: bool = True   # 被合并后的旧簇标记为 False,而非物理删除(便于追溯)

@dataclass
class SSNEdge:
    source: str
    target: str
    identity: float
    coverage: float

    def weight(self, mode: str = "identity_x_coverage") -> float:
        if mode == "identity_only":
            return self.identity
        if mode == "coverage_only":
            return self.coverage
        return self.identity * self.coverage

@dataclass
class OrthologFamily:
    family_id: str
    component_id: int
    member_genes: Set[str]
    edge_density: float

    @property
    def size(self) -> int:
        return len(self.member_genes)

@dataclass
class LocalCluster:
    """Step2 基因组内自聚类产生的局部簇。"""
    cluster_id: str
    member_ids: Set[str] = field(default_factory=set)
    representative_ids: List[str] = field(default_factory=list)
    is_orphan: bool = False

@dataclass
class GenomeIntraResult:
    """单个基因组完成 Step2 后的结果。"""
    genome_id: str
    local_clusters: List[LocalCluster] = field(default_factory=list)
    local_representatives: List[str] = field(default_factory=list)
    orphan_genes: List[str] = field(default_factory=list)
    n_genes: int = 0
    n_clusters: int = 0
    n_representatives: int = 0
```

#### `core/events.py` — 领域事件(可观测性基础)

```python
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Dict

@dataclass(frozen=True)
class DomainEvent:
    timestamp: datetime
    payload: Dict[str, Any]

class ClusterCreated(DomainEvent): ...
class ClusterMerged(DomainEvent): ...
class RepresentativeAdded(DomainEvent): ...
class GenomeJoined(DomainEvent): ...
class EdgeCompleted(DomainEvent): ...

class EventBus:
    """极简内存事件总线,订阅者可以是日志记录器/进度条/未来的可视化前端"""
    def __init__(self):
        self._subscribers: list = []

    def subscribe(self, handler) -> None:
        self._subscribers.append(handler)

    def publish(self, event: DomainEvent) -> None:
        for h in self._subscribers:
            h(event)
```

#### `core/policies.py` — 可替换的决策策略

```python
from typing import Protocol, List
from .entities import Gene, AlignmentResult, Cluster

class OrderingPolicy(Protocol):
    """Order-D 排序策略接口,默认实现见 adapters 或直接内置在 engine 中"""
    def sort(self, genomes: list) -> list: ...

class ThresholdDecision(Protocol):
    action: str            # "create_new" | "attach_no_rep" | "attach_add_rep" | "merge_no_rep" | "merge_add_rep"
    target_clusters: list

class ThresholdPolicy(Protocol):
    cutoff_i: float
    cutoff_ii: float
    def decide(self, new_rep: Gene, hits: List[AlignmentResult], lib) -> ThresholdDecision: ...
```

---

### 2.2 `ports/` — 抽象接口层

```python
# ports/aligner_port.py
from typing import Protocol, Sequence, List
from ..core.entities import Gene, AlignmentResult

class AlignerPort(Protocol):
    def align_all_vs_all(self, genes: Sequence[Gene]) -> List[AlignmentResult]:
        """基因组内全对全比对 (Step2 用途)"""
        ...

    def align_query_vs_subject(self, query: Sequence[Gene], subject: Sequence[Gene]) -> List[AlignmentResult]:
        """新代表 vs 全局代表库 / 补边批量比对 (Step3/Step5 用途)"""
        ...
```

```python
# ports/storage_port.py
from typing import Protocol, Optional
from ..core.entities import Cluster, SSNEdge

class RepositoryPort(Protocol):
    def save_checkpoint(self, stage: str, state: dict) -> None: ...
    def load_checkpoint(self, stage: str) -> Optional[dict]: ...
    def save_global_library(self, clusters: dict) -> None: ...
    def load_global_library(self) -> dict: ...
    def save_ssn_edges(self, edges: set) -> None: ...
```

```python
# ports/clustering_port.py
from typing import Protocol
import networkx as nx

class ClusteringPort(Protocol):
    def cluster(self, graph: nx.Graph, **params) -> list[set[str]]:
        """输入图,输出聚类结果(每个 set 为一个家族的节点集合)"""
        ...
```

> **设计说明**: 使用 `typing.Protocol`(结构化子类型),而非强制继承的抽象基类(ABC)。
> 这意味着任何实现了对应方法签名的类都自动满足接口,不需要显式继承,
> 降低了 adapters 和第三方库之间的耦合,也便于给已有的第三方库(如 `markov_clustering` 模块函数)
> 直接写一个薄包装即可符合接口,无需侵入式改造。

---

### 2.3 `adapters/` — 具体实现层

| Adapter | 实现的 Port | 说明 |
|---------|------------|------|
| `adapters/aligners/mmseqs_adapter.py` | `AlignerPort` | 已实现。统一为 `build_index`(createdb,持久化索引)+ `search`(手动拆分 createdb(仅query)→search→convertalis 三步,复用已建索引)两步原语,`align_all_vs_all`/`align_query_vs_subject` 作为一把梭的便捷方法保留(走 `easy-search`) |
| `adapters/aligners/diamond_adapter.py` | `AlignerPort` | 已实现。接口/组织方式与 `mmseqs_adapter.py` 对齐(`build_index`/`search` + 便捷方法),因 diamond 的 `blastp` 可直接以 fasta 作为 query,`search()` 比 mmseqs 更简单(query 侧无需单独建库) |
| `adapters/aligners/mock_adapter.py` | `AlignerPort` | 测试用,返回预设/随机比对结果,不依赖外部二进制;同样实现了 `build_index`/`search`(索引 handle 直接是基因列表本身,无实际开销) |
| `adapters/storage/memory_repo.py` | `RepositoryPort` | 纯内存实现,小规模/测试场景,不支持跨进程 resume |
| `adapters/storage/file_repo.py` | `RepositoryPort` | **已完成(第8步,2026-10-02)**。与文档原计划的 `sqlite_repo.py` 不同,经与用户确认改为最简单的基于文件(pickle)的持久化方案:每个 Stage 产出序列化为独立文件,存放于 `<output_dir>/.linc_work/state/`,与现有 `mmseqs_db`/`diamond_db` 等索引目录风格一致;支持跨进程 `--resume`/`--append`,`sqlite_repo.py` 本身未实现 |
| `adapters/clustering/union_find_adapter.py` | `ClusteringPort` | Step2 基因组内聚类用 |
| `adapters/clustering/mcl_adapter.py` | `ClusteringPort` | **已完成**。与用户确认(2026-10-02)不做 `markov_clustering` Python 库 fallback,直接调用用户已安装的系统 `mcl` 命令行(`.abc` 格式输入,`-I` 指定 inflation) |

**Adapter 选择与装配**: 统一在 `config/schema.py` 中声明 `backend` 字段,由一个简单的工厂函数
(`adapters/factory.py`)根据配置创建对应实例,`engine/` 层只依赖 Port 类型注解,不关心具体实现:

```python
# adapters/factory.py
def create_aligner(settings) -> AlignerPort:
    if settings.aligner.backend == "mmseqs":
        return MMseqsAdapter(settings.aligner.mmseqs)
    if settings.aligner.backend == "diamond":
        return DiamondAdapter(settings.aligner.diamond)
    if settings.aligner.backend == "mock":
        return MockAligner()
    raise ValueError(f"Unknown aligner backend: {settings.aligner.backend}")
```

---

### 2.4 `engine/` — 增量引擎(核心)

`engine/` 是唯一包含"业务逻辑"的地方,所有依赖都通过构造函数注入 Port 接口,**不直接 import 任何
adapters 子模块**(避免反向依赖)。

```
engine/
├── genome_load_stage.py    # 对应 v1 Step1
├── intra_genome_stage.py   # 对应 v1 Step2 (Union-Find 聚类 + 代表选取)
├── global_join_stage.py    # 对应 v1 Step3 (Order-D 排序 + 双阈值决策树,核心)
├── ssn_builder_stage.py    # 对应 v1 Step4+5 合并 (稀疏图构建 + 延迟补边,对外一个阶段)
└── family_refine_stage.py  # 对应 v1 Step6 (MCL 细分)
```

#### 对外统一入口:`pipeline.Pipeline`(取代原计划的 `IncrementalSSNEngine`)

> **与用户确认(2026-10-02,第5步)**: 不单独新建 `engine/incremental_engine.py`。
> "增量并入"的核心算法已完整实现在 `GlobalJoinStage` 里,"Stage 编排顺序"这件事
> 直接合并进 `pipeline.py` 的 `Pipeline` 类(见 2.5 节),不再额外设计一层同义的
> `IncrementalSSNEngine` 外壳。

```python
class Pipeline:
    def __init__(
        self,
        aligner: AlignerPort,
        clustering_for_intra: ClusteringPort,
        clustering_for_mcl: ClusteringPort,
        repository: RepositoryPort,
        config: "PipelineConfig",
    ):
        ...

    def run(self, genome_files: list[str]) -> "PipelineResult":
        """首次完整运行(无历史状态), 依次跑完 5 个 Stage 并逐一持久化"""
        ...

    def resume(self, genome_files: list[str]) -> "PipelineResult":
        """从上次中断处接着跑完剩余 Stage(Stage 级别恢复粒度, 见 2.5 节)"""
        ...

    def append(self, new_genome_files: list[str]) -> "PipelineResult":
        """把新基因组追加到历史完整结果上(第8步, 2026-10-02)"""
        ...
```

#### 各 Stage 的输入输出契约(与 v1 保持算法一致,仅接口形式现代化)

| Stage | 输入 | 输出 | 核心算法(详见 `LINC_算法设计文档.md` 第3节,算法本身不变) |
|-------|------|------|------------------------------------------------------|
| `GenomeLoadStage` | `list[str]`(文件路径) | `dict[str, Genome]` | FASTA 解析 + 统计 |
| `IntraGenomeStage` | `dict[str, Genome]` | `dict[str, GenomeIntraResult]` | 全对全比对 + Union-Find + 代表选取 |
| `GlobalJoinStage` | `dict[str, GenomeIntraResult]` | `IncrementalResult` | Order-D 排序 + 双阈值决策树(**不可并行**,见 v1 文档 6.3) |
| `SSNBuilderStage` | `IncrementalResult` | `CompleteSSN` | 稀疏图构建 + 长度过滤/去重补边 |
| `FamilyRefineStage` | `CompleteSSN` | `MCLResult` | 自适应 Inflation + MCL 聚类 |

> **算法细节不在本文档重复**,请始终以 `docs/归档/LINC_算法设计文档.md` 第 3 节(Order-D 排序、
> Union-Find、双阈值决策树、补边剪枝、MCL Inflation)作为算法真理来源(Single Source of Truth),
> 本文档只定义"这些算法应该被放进哪个类、以什么接口对外暴露"。
>
> **已知偏离**: `IntraGenomeStage`(Step2)代表选取的判断依据已修正为"簇内最小 pairwise
> identity 是否 >= Cutoff-II"(而非归档文档字面的"最大"值——该依据在语义上不成立,会导致
> 多代表分支的贪心覆盖永远无法提前停止),并已实现 `multi_representative_clustering`
> 为贪心覆盖(k-center 近似)算法。该项算法以 `linc/engine/intra_genome_stage.py` 的实现
> 及其 docstring 为准,归档文档此处内容已过时。

#### `GlobalJoinStage` 性能方案(与用户讨论确认,2026-10-01)

**问题背景**: `AlignerPort` 的 mmseqs 实现每次调用都对 subject 重新 `createdb`+建索引
(无状态、无缓存)。若严格按 STEP3 文档字面"逐个局部代表单独与全局库比对"的粒度实现,
全局库(作为 subject)会在整个 Stage 运行期间从 0 增长到上千条,每次都重建索引的开销
是 O(总局部代表数 × 平均库规模),可能成为实际性能瓶颈。

**最终方案(三轮讨论确认)**:

1. **query/target 方向反转**: `IntraGenomeStage` 的 all-vs-all 本质是"基因组整体(query)
   vs 该基因组自己的 DB(target)"。类比过来,`GlobalJoinStage` 应该是"全局代表库整体
   (query) vs 新进来的基因组 DB(target)"——**复用 `IntraGenomeStage` 阶段已经建好的
   基因组索引**,避免 target 侧重复建库(这是最贵的部分)。
2. **`AlignerPort` 统一为 build_index(建库) + search(查询) 两步原语**,`IntraGenomeStage`
   / `GlobalJoinStage` 均遵循同一操作模式:先 `build_index` 得到可复用的索引 handle,
   再用该 handle 反复 `search`。`IntraGenomeStage.run()` 为基因组建一次索引,存入
   `GenomeIntraResult.genome_index`,供 `GlobalJoinStage` 直接取用,不重新建库。
3. **持久化工作目录**: 运行期间维护一个不自动清理的工作目录(如 `output/.linc_work/`)
   保存所有基因组的 mmseqs DB 文件,贴合 LINC 本身的断点续跑/增量添加数据的定位,
   避免过程文件被提前清理导致无法复用。
4. **`GlobalJoinStage` 批处理粒度**: 按 Order-D 顺序逐个基因组处理;对每个基因组,
   用当前全局代表库整体作为 query,对该基因组持久化的索引做**一次**比对,
   从结果中只取"局部代表"对应的行,按确定顺序(字典序)逐个走双阈值决策树并串行应用决策。
   **正确性论证**: `IntraGenomeStage` 的 `multi_representative_clustering` 本身是贪心
   覆盖算法,保证同一基因组内所有被选出的局部代表彼此 pairwise identity 必然 < Cutoff-II,
   因此"基因组级批处理快照"不会引入同基因组内本应共享代表却被误判为各自新建簇的问题,
   与"严格逐代表单独 query"正确性等价,但比对调用次数降至 O(基因组数)。

该方案已在 `开发进度.md` 变更日志中记录,是 `GlobalJoinStage` 实现的既定设计依据。

#### 跨进程结果非确定性问题与修复(与用户多轮排查确认,2026-10-01)

实测中发现同样输入多次 `linc run --backend mmseqs` 的结果存在微小差异(簇数/代表数/边数
相差个位数到百余个成员归属),经排查定位到两个独立根因并均已修复:

1. **`representative_ids` 遍历顺序受 `PYTHONHASHSEED` 影响**: `Cluster.representative_ids`
   字段类型为 `Set[str]`,`GlobalRepresentativeLibrary.representative_genes()` 原直接遍历该
   set 构造下一次 query 的代表基因列表;Python 默认开启字符串哈希随机化,导致内容相同的
   set 在不同进程间遍历顺序不同,使喂给比对工具的 query 序列顺序产生跨进程差异。
   **修复**: 改为 `sorted(cluster.representative_ids)` 按字典序遍历,字段类型本身不变。
2. **mmseqs2/diamond 比对结果内容确定,但输出行序不保证跨运行一致**: 通过最小复现实验
   (两基因组反复 `createdb+search+convertalis`)证实,相同输入两次比对的结果集合完全相同,
   但原始输出行序不同;而 `GlobalJoinStage` 的决策逻辑(如 `hit_by_cluster` 的平局判定)
   依赖"先遇到的 hit"顺序,行序抖动被放大为簇归属差异。
   **修复**: `mmseqs_adapter.py`/`diamond_adapter.py` 的结果解析方法均在返回前按
   `(query_id, subject_id)` 字典序排序,使上层拿到的 `List[AlignmentResult]` 顺序恒定。
3. **mmseqs 侧额外加固**: `search`/`easy-search` 命令均显式传 `--split 1`,禁用其基于运行时
   可用内存的自动数据库分片机制(该机制会导致边界候选命中结果的微小抖动,与线程数无关);
   CLI 新增 `--threads` 参数暴露给用户,默认 0(不显式指定,使用全部核心)。

以上三处修复后,经用户在单线程(`--threads 1`)及多线程(默认全部核心)场景下各重跑两次
验证,`join_summary.json` 的聚合指标与全部簇明细均完全一致,非确定性问题已彻底解决。

#### 持久化工作目录的文件组织(`output/.linc_work/<backend>_db/`)

`build_index`/`search` 两步原语约定的持久化目录结构(mmseqs/diamond 一致):

```
<work_dir>/<genome_id>/            # build_index() 建的索引目录(长期保留,不自动清理)
├── genes.fasta                    # 该基因组全部基因的 fasta
├── db*                            # createdb/makedb 产出的索引文件
└── search/                        # search() 产出的比对中间/结果文件(同样不自动清理)
    ├── query.fasta                # 本次调用的 query(通常是当前全局代表库)
    ├── result.tsv                 # 比对结果(已按 (query_id, subject_id) 排序)
    └── ...                        # mmseqs 额外产出 query_db*/result_db*/tmp/
```

`search/` 子目录每次调用前会先整体清空再重建,即"同一个 genome_id 重跑直接覆盖旧结果"
(与用户确认的约定),索引本身(`genes.fasta`/`db*`)则长期保留、不随 `search` 调用清理。

---

### 2.5 `pipeline.py` — `Pipeline` + Checkpoint/Resume/Append(已完成,第5/8步)

> **实际实现与本文档最初草图(通用 `load_checkpoint("pipeline")`)不同**,经与用户
> 多轮确认(2026-10-02)改为`RepositoryPort` 按 5 个 Stage 产出类型逐一提供专用
> `save_*`/`load_*`/`has_*` 方法(而非存一个笼统的 `dict` checkpoint),`Pipeline`
> 提供 `run`/`resume`/`append` 三个方法,而非 `run(resume: bool)` 这种参数形式:

```python
class Pipeline:
    def __init__(
        self,
        aligner: AlignerPort,
        clustering_for_intra: ClusteringPort,
        clustering_for_mcl: ClusteringPort,
        repository: RepositoryPort,
        config: "PipelineConfig",
    ):
        ...

    def run(self, genome_files: list[str]) -> "PipelineResult":
        """首次完整运行: 依次跑完 5 个 Stage, 每个 Stage 跑完后整体持久化一次"""
        ...

    def resume(self, genome_files: list[str]) -> "PipelineResult":
        """按 genomes->intra_results->join_result->complete_ssn->families 顺序,
        依次用 has_*() 判断某 Stage 是否已持久化过: 已持久化的直接 load_*() 复用,
        未持久化的才调用对应 Stage 重新跑"""
        ...

    def append(self, new_genome_files: list[str]) -> "PipelineResult":
        """新基因组只跑 GenomeLoadStage+IntraGenomeStage, 再与历史 join_result
        接续运行 GlobalJoinStage; SSNBuilderStage/FamilyRefineStage 对合并后的
        完整数据全量重跑(不做增量优化)"""
        ...
```

**Resume/Checkpoint 粒度(已与用户核实确认, 2026-10-02)**: 以"**Stage 级别**"为
恢复粒度,且是 Stage **内部循环全部基因组跑完后才整体持久化一次**,不是"基因组
级别"续跑。具体表现为: `_run_intra()`/`_run_join()` 内部各是一个对全部基因组
逐一处理的循环, 只有循环整体结束才调用一次 `save_intra_results()`/
`save_join_result()`; 若失败发生在 `IntraGenomeStage`(或 `GlobalJoinStage`)
处理第 N 个基因组(N 小于总数)时, resume 后会**重新执行该 Stage 涉及的全部
基因组**, 而不是只补做剩下未完成的那些 —— 这是为了保证正确性优先、避免细粒度
续跑引入的额外状态管理复杂度, 但 Stage 与 Stage 之间的断点是有效的(已完整持久化
的 Stage 不会被重新执行)。细粒度到基因组级别的 checkpoint 仍是一项未落地的
待办(见第7节开放问题), 如需要应改造各 Stage 内部循环 + `RepositoryPort`/
`FileRepository` 接口粒度。

---

### 2.6 `config/schema.py` — 类型安全配置

**取舍说明**: 经过与你确认,考虑到项目依赖应保持极简(当前仅 `numpy`/`networkx`),
**本期不引入 Pydantic**,改用标准库 `dataclasses` + 手写校验函数实现类型安全配置,
未来如果配置复杂度显著上升,可以无痛升级为 Pydantic(字段定义方式高度相似)。

```python
from dataclasses import dataclass, field

@dataclass
class ThresholdConfig:
    cutoff_i: float = 0.50
    cutoff_ii: float = 0.80
    max_per_cluster: int = 5

    def validate(self) -> None:
        assert 0.0 <= self.cutoff_i <= self.cutoff_ii <= 1.0, "必须满足 0 <= cutoff_i <= cutoff_ii <= 1"
        assert self.max_per_cluster >= 1

@dataclass
class MMseqsConfig:
    sensitivity: float = 4.5
    evalue: float = 1e-5
    max_accept: int = 500
    min_alignment_score: int = 20

@dataclass
class AlignerConfig:
    backend: str = "mmseqs"   # mmseqs | diamond | mock
    mmseqs: MMseqsConfig = field(default_factory=MMseqsConfig)

@dataclass
class MCLConfig:
    inflation: float = 2.0
    auto_inflation: bool = True
    edge_weight_mode: str = "identity_x_coverage"

@dataclass
class PruningConfig:
    length_ratio_threshold: float = 10.0
    dedup_enabled: bool = True

@dataclass
class ParallelConfig:
    n_threads: int = 0
    inter_genome_parallel: bool = True
    finalize_parallel: bool = True

@dataclass
class StorageConfig:
    backend: str = "memory"   # memory | sqlite
    sqlite_path: str = "output/linc_state.db"

@dataclass
class LINCSettings:
    threshold: ThresholdConfig = field(default_factory=ThresholdConfig)
    aligner: AlignerConfig = field(default_factory=AlignerConfig)
    mcl: MCLConfig = field(default_factory=MCLConfig)
    pruning: PruningConfig = field(default_factory=PruningConfig)
    parallel: ParallelConfig = field(default_factory=ParallelConfig)
    storage: StorageConfig = field(default_factory=StorageConfig)
    genome_list: list[str] = field(default_factory=list)

    def validate(self) -> None:
        self.threshold.validate()
        assert self.aligner.backend in {"mmseqs", "diamond", "mock"}
        assert self.storage.backend in {"memory", "sqlite"}
```

---

### 2.7 `cli.py` — 命令行入口(已实现,较文档初始草图有所扩展)

实际实现提前采用了 `typer`(非占位),并在 `linc run` 单命令内直接依次接通了
`GenomeLoadStage`→`IntraGenomeStage`→`GlobalJoinStage` 三个 Stage(后续 Stage 完成后按
同样方式扩展,命令行接口形态不再变动):

```python
@app.callback()
def _callback() -> None:
    """空 callback, 强制 typer 保留子命令结构(否则单命令会被折叠,无法通过
    `linc run ...` 调用)。"""


@app.command()
def run(
    input_dir: str = typer.Option(..., "-i", "--input-dir"),
    output_dir: str = typer.Option("output", "-o", "--output-dir"),
    cutoff_i: float = typer.Option(0.5, "--cutoff-i"),
    cutoff_ii: float = typer.Option(0.8, "--cutoff-ii"),
    cutoff_c: float = typer.Option(0.0, "--cutoff-c"),
    cov_mode: int = typer.Option(0, "--cov-mode"),      # 0双向|1subject|2query
    max_per_cluster: int = typer.Option(5, "--max-per-cluster"),
    aligner: str = typer.Option("mmseqs", "--aligner"),  # mmseqs | diamond(mock不暴露, 仅测试用)
    threads: int = typer.Option(0, "--threads"),       # 0 = 不显式指定, 使用工具自身默认值
    edge_completion_mode: str = typer.Option("representative_member", "--edge-completion-mode"),
    mega_cluster_warning_threshold: int = typer.Option(0, "--mega-cluster-warning-threshold"),
    inflation: float = typer.Option(1.5, "--inflation"),
    edge_weight_mode: str = typer.Option("identity_x_coverage", "--edge-weight-mode"),
    edge_weight_alpha: float = typer.Option(0.5, "--edge-weight-alpha"),
    resume: bool = typer.Option(False, "--resume"),     # 与 -a/--append 互斥
    append: str = typer.Option(None, "-a", "--append"), # 追加新基因组目录(支持链式多次追加)
    verbose: bool = typer.Option(False, "--verbose"),
) -> None:
    """构造 Pipeline 并调用 run()/resume()/append(),
    输出 load/intra/join/ssn/family_summary.json + Orthogroups.txt。"""
    ...
```

> **实际实现较本文档最初草图大幅扩展(第6/8步, 2026-10-02)**: `--resume`(开关)
> 与 `-a`/`--append <目录>`(接收路径, 非开关)均已实现并接入, 二者互斥; 两者均
> 要求 `-i` 与历史 `input_signature`(discover 出的文件名集合, 非路径字符串)
> 一致, 否则报错退出; `-a` 额外要求历史 5 个 Stage 产出均已完整持久化(否则
> 提示先用 `--resume` 补全)。未采用 `config.yaml` 配置文件驱动方式, 全部参数
> 通过 CLI Option 直接传入(`config/schema.py` 的 dataclass+validate 仍未独立
> 抽取, 校验逻辑现状直接写在 `cli.py` 的 `run` 函数体内)。
> `--backend` 已按用户要求重命名为 `--aligner`, 默认值由 `mock` 改为 `mmseqs`,
> 取值收窄为 `mmseqs|diamond`(`MockAligner` 保留仅供测试内部构造夹具使用, 不再
> 作为 CLI 可选项暴露)。
> `--aligner mmseqs|diamond` 时,`work_dir` 固定为 `<output_dir>/.linc_work/<aligner>_db/`,
> 随 `-o` 参数变化,运行期间不自动清理(贴合断点续跑/增量添加数据的定位); 持久化
> 状态(`FileRepository`)存放于 `<output_dir>/.linc_work/state/`。

---

## 3. 数据流与调用关系总览

```
cli.py
  │
  ▼
config/schema.py (LINCSettings.validate())
  │
  ▼
adapters/factory.py  → 创建 AlignerPort / ClusteringPort / RepositoryPort 具体实例
  │
  ▼
pipeline.py (Pipeline, 注入上述三个 Port, 提供 run/resume/append)
  │
  ├─ GenomeLoadStage
  ├─ IntraGenomeStage        (可并行: 基因组间)
  ├─ GlobalJoinStage         (不可并行: Order-D 顺序依赖)
  ├─ SSNBuilderStage         (可并行: 簇间)
  └─ FamilyRefineStage       (可并行: 连通分量间)
  │
  ▼
PipelineResult { genomes, intra_results, join_result, complete_ssn, families }
```

> `engine/incremental_engine.py` 未单独新建, "Stage 编排顺序"这件事已直接合并进
> `pipeline.py` 的 `Pipeline` 类(见 2.4/2.5 节, 与用户确认, 2026-10-02)。

---

## 4. 与 v1 文档的关系

- **算法正确性的唯一权威来源仍是 `LINC_算法设计文档.md` 第 3 节**(Order-D、Union-Find、
  双阈值决策树、补边剪枝、MCL Inflation),本文档不重复算法细节,只定义工程组织方式。
- **数据结构**: v2 的 `core/entities.py` 是 v1 `models/data_structures.py` 的现代化重写版本,
  字段语义保持一致,但引入了 `frozen=True`(不可变值对象)、`property`(派生字段如 `coverage`/`size`)
  等现代 dataclass 写法,减少冗余字段和潜在的状态不一致风险。
- **v1 文档中"待实现/待确认事项清单"(第7节)在 v2 中原样有效**,不因架构升级而自动解决,
  需要在对应 adapter/engine 实现时逐项确认。

---

## 5. 推荐开发顺序(供 `docs/开发进度.md` 初始化参考)

1. `core/entities.py` + `core/events.py`(无外部依赖,最先可写,也最容易单测)
2. `ports/*`(接口定义,工作量小但统领后续所有实现)
3. `adapters/aligners/mock_adapter.py` + `adapters/clustering/union_find_adapter.py`(让 engine 尽早可测)
4. `engine/genome_load_stage.py` → `intra_genome_stage.py` → `global_join_stage.py`
   → `ssn_builder_stage.py` → `family_refine_stage.py`(严格按流水线顺序实现并逐个写单测)
5. `pipeline.py`(串联已实现的 Stage, 不单独建 `incremental_engine.py`)— **已完成**
6. `adapters/storage/memory_repo.py`(先用内存版跑通全流程)— **已完成**
7. `config/schema.py` + `cli.py`(让整个系统可以通过命令行驱动)— `cli.py` 已完成, `config/schema.py` 未独立抽取(校验逻辑现状在 `cli.py` 内)
8. `adapters/aligners/mmseqs_adapter.py`(接入真实比对引擎,替换 Mock)— **已完成**
9. `adapters/storage/file_repo.py`(基于pickle文件, 取代原计划的 `sqlite_repo.py`)+ `--resume`/`--append` 完整验证 — **已完成**(2026-10-02)
10. `adapters/aligners/diamond_adapter.py` + `adapters/clustering/mcl_adapter.py`(调用系统 `mcl` 命令, 不做 `markov_clustering` 库 fallback)— **均已完成**

---

## 6. 性能优化建议(第1轮,待后续实现时逐项落地)

> 以下建议针对上层算法执行逻辑,按预期收益排序。优先级最高的是 6.1、6.3(直接降低复杂度),
> 其余可作为后续性能调优阶段再补充,不阻塞第一版实现(先 correctness,后 performance)。

### 6.1 `GlobalJoinStage` 的比对/决策解耦 + 全局库预过滤

- **比对与决策解耦**:决策逻辑本身必须严格按 Order-D 顺序串行执行,但"新代表 vs 全局代表库"的比对
  可以提前批量/异步发起(生产者),决策环节只负责按顺序消费比对结果(消费者),避免"比对一个、决策一个"
  的同步阻塞链。
- **全局代表库预过滤**:随着代表库增长(M 可能达到成千上万),避免每次都用全量代表库做精确比对。
  建议引入预过滤(如长度/k-mer 粗筛),先筛出候选子集,再对子集做精确比对,进一步压缩 M² 复杂度中的 M。

### 6.2 `IntraGenomeStage` 并行粒度

- 基因组间并行是对的,但单个基因组若基因数很大,会成为长尾瓶颈(木桶效应)。
- 建议直接利用 MMseqs2/DIAMOND 自带多线程参数(`--threads`),而不是在 Python 层对基因做切分比对。
- Union-Find 本身近 O(n),瓶颈几乎全在比对阶段,无需额外优化。

### 6.3 `SSNBuilderStage` 内存优化(整数 ID + 延迟构图)

- `compared_pairs` 等去重集合建议用**整数 ID 替代字符串 ID**:在 `core/entities.py` 的 `Gene`/`Cluster`
  增加内部自增整数索引,仅对外展示时映射回字符串 `gene_id`,降低哈希/比较开销与内存占用。
- 节点规模较大(如超过 10 万)时,`networkx.Graph` 性能会明显下降,建议 `union_find_adapter`/
  `mcl_adapter` 内部优先使用轻量数组/字典邻接表,仅在需要调用 networkx 算法(如连通分量分析)时
  才临时构建 `nx.Graph`,避免全程依赖 networkx 对象作为主数据结构。

### 6.4 `FamilyRefineStage` 任务调度均衡

- 连通分量间并行是对的,但分量大小分布往往极不均匀。
- 建议按节点数对连通分量排序,采用"大任务优先调度"策略(类似 LPT, Longest Processing Time First),
  让最耗时的任务尽早开始执行,避免拖慢整体收尾时间。

### 6.5 `ParallelConfig.n_threads` 语义明确

- `n_threads = 0` 约定为"自动检测 CPU 核数",该约定需在 `config/schema.py` 加载/校验逻辑中统一解析
  (如调用 `os.cpu_count()`),避免各 Stage 各自解读不一致。

### 6.6 `EventBus` 预留性能观测能力

- 建议为 `DomainEvent`(或新增 `StageTimingEvent`)增加耗时字段,使性能调优无需额外埋点,
  直接订阅事件即可获取"每个 Stage/每个基因组处理耗时"数据,便于后续火焰图分析与瓶颈定位。

---

## 7. 开放问题(需后续迭代中明确)

- Checkpoint 的"Stage 级别"粒度是否足够?若未来基因组规模达到数千级,是否需要升级为更细粒度的恢复策略?
  (2026-10-02 与用户核实确认: 当前实现确实是 Stage 内部全部基因组跑完才整体持久化一次,
  中途中断会导致该 Stage 涉及的全部基因组重新处理一遍, 细粒度到基因组级别的 checkpoint 仍待后续决定是否实现)
- `EventBus` 是否需要支持异步/多进程场景(当前默认假设单进程内存事件总线)?
- `typer` vs `argparse` 的最终选择待定,取决于是否接受新增依赖。
