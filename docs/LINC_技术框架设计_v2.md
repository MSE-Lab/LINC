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

@dataclass
class Genome:
    genome_id: str
    file_path: str
    genes: List[Gene]
    status: "GenomeStatus" = None

class GenomeStatus(str, Enum):
    PENDING = "pending"
    LOADED = "loaded"
    FAILED = "failed"

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
| `adapters/aligners/mmseqs_adapter.py` | `AlignerPort` | 封装 `mmseqs createdb/prefilter/align` 三步调用 |
| `adapters/aligners/diamond_adapter.py` | `AlignerPort` | 封装 DIAMOND 调用(可选) |
| `adapters/aligners/mock_adapter.py` | `AlignerPort` | 测试用,返回预设/随机比对结果,不依赖外部二进制 |
| `adapters/storage/memory_repo.py` | `RepositoryPort` | 纯内存实现,小规模/测试场景,不支持跨进程 resume |
| `adapters/storage/sqlite_repo.py` | `RepositoryPort` | SQLite 持久化,支持 checkpoint/resume,适合大规模增量场景 |
| `adapters/clustering/union_find_adapter.py` | `ClusteringPort` | Step2 基因组内聚类用 |
| `adapters/clustering/mcl_adapter.py` | `ClusteringPort` | 优先尝试 `markov_clustering` Python 库,找不到则 fallback 到系统 `mcl` 命令行 |

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
├── incremental_engine.py   # IncrementalSSNEngine: 对外统一入口,持有 EventBus 和各 Stage
├── genome_load_stage.py    # 对应 v1 Step1
├── intra_genome_stage.py   # 对应 v1 Step2 (Union-Find 聚类 + 代表选取)
├── global_join_stage.py    # 对应 v1 Step3 (Order-D 排序 + 双阈值决策树,核心)
├── ssn_builder_stage.py    # 对应 v1 Step4+5 合并 (稀疏图构建 + 延迟补边,对外一个阶段)
└── family_refine_stage.py  # 对应 v1 Step6 (MCL 细分)
```

#### `IncrementalSSNEngine` 对外接口

```python
class IncrementalSSNEngine:
    def __init__(
        self,
        aligner: AlignerPort,
        clustering: ClusteringPort,
        repository: RepositoryPort,
        threshold_policy: ThresholdPolicy,
        event_bus: Optional[EventBus] = None,
    ):
        ...

    def run(self, genome_files: list[str], resume: bool = False) -> "EngineResult":
        """串联五个 Stage,返回最终 SSN + 家族划分结果"""
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

> **算法细节不在本文档重复**,请始终以 `LINC_算法设计文档.md` 第 3 节(Order-D 排序、
> Union-Find、双阈值决策树、补边剪枝、MCL Inflation)作为算法真理来源(Single Source of Truth),
> 本文档只定义"这些算法应该被放进哪个类、以什么接口对外暴露"。

---

### 2.5 `pipeline.py` — 声明式 Pipeline + Checkpoint/Resume

```python
@dataclass
class StageResult:
    stage_name: str
    output: Any
    event_count: int

class Pipeline:
    def __init__(self, engine: IncrementalSSNEngine, repository: RepositoryPort):
        self.engine = engine
        self.repository = repository

    def run(self, genome_files: list[str], resume: bool = False) -> "EngineResult":
        checkpoint = self.repository.load_checkpoint("pipeline") if resume else None
        start_stage = checkpoint["completed_stage"] + 1 if checkpoint else 0
        # 依次执行 engine 内部 5 个 stage,每完成一个 stage 调用
        # self.repository.save_checkpoint("pipeline", {...}) 落盘进度
        ...
```

**Resume 语义**: 以"Stage 级别"为最小恢复粒度(不做到"基因组级别"的细粒度续跑,
后者复杂度过高且收益有限)。即:如果失败发生在 `GlobalJoinStage` 处理第 50 个基因组时,
resume 后会**重新执行整个 `GlobalJoinStage`**,而不是从第 50 个基因组继续
—— 这是为了保证"一次比对全局复用"等不变量不被断点逻辑破坏,优先保证正确性而非恢复效率。

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

### 2.7 `cli.py` — 命令行入口

```python
import typer   # 新增轻量依赖,替代裸脚本调用,体现"现代化" CLI 体验

app = typer.Typer()

@app.command()
def run(config: str = "config.yaml", resume: bool = False):
    """运行完整 LINC 流程"""
    ...

@app.command()
def resume(checkpoint_dir: str = "output/"):
    """从上次 checkpoint 继续运行"""
    ...

if __name__ == "__main__":
    app()
```

> 若不希望引入 `typer`,可退化为标准库 `argparse` 实现,功能等价但交互体验略逊。
> 本期先用 `argparse` 占位实现,`typer` 作为后续可选升级项(见开发进度文档)。

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
engine/incremental_engine.py (IncrementalSSNEngine, 注入上述三个 Port)
  │
  ▼
pipeline.py (Pipeline.run, 支持 resume)
  │
  ├─ GenomeLoadStage
  ├─ IntraGenomeStage        (可并行: 基因组间)
  ├─ GlobalJoinStage         (不可并行: Order-D 顺序依赖)
  ├─ SSNBuilderStage         (可并行: 簇间)
  └─ FamilyRefineStage       (可并行: 连通分量间)
  │
  ▼
EngineResult { ortholog_families, ssn_edges, global_library, run_summary }
```

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
5. `engine/incremental_engine.py` + `pipeline.py`(串联已实现的 Stage)
6. `adapters/storage/memory_repo.py`(先用内存版跑通全流程)
7. `config/schema.py` + `cli.py`(让整个系统可以通过命令行驱动)
8. `adapters/aligners/mmseqs_adapter.py`(接入真实比对引擎,替换 Mock)
9. `adapters/storage/sqlite_repo.py` + Checkpoint/Resume 完整验证
10. `adapters/aligners/diamond_adapter.py` + `adapters/clustering/mcl_adapter.py`(补齐可选后端)

---

## 6. 开放问题(需后续迭代中明确)

- Checkpoint 的"Stage 级别"粒度是否足够?若未来基因组规模达到数千级,是否需要升级为更细粒度的恢复策略?
- `EventBus` 是否需要支持异步/多进程场景(当前默认假设单进程内存事件总线)?
- `typer` vs `argparse` 的最终选择待定,取决于是否接受新增依赖。
