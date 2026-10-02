"""RepositoryPort: 存储/持久化的抽象接口 (Protocol)。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.2 节, 并于第5步(pipeline打通)
与用户确认重新设计: 接口按 engine/ 五个 Stage 的实际产出类型逐一对应提供
save/load 方法(而非早期版本里粗粒度的 checkpoint/global_library/ssn_edges),
使 `Pipeline` 在每个 Stage 跑完后可以直接、类型安全地做存取, 也为未来
adapters/storage/sqlite_repo.py(第8步持久化+断点续跑)实现同一接口打好基础。

具体实现见 adapters/storage/memory_repo.py(纯内存, 不支持跨进程恢复)。
"""
from typing import Dict, List, Optional, Protocol, runtime_checkable

from ..core.entities import CompleteSSN, Genome, GenomeIntraResult, IncrementalResult, OrthologFamily


@runtime_checkable
class RepositoryPort(Protocol):
    """存储/持久化抽象接口, 按 Stage 产出类型划分存取方法。"""

    def save_genomes(self, genomes: Dict[str, Genome]) -> None:
        """保存 GenomeLoadStage 产出(genome_id -> Genome)。"""
        ...

    def load_genomes(self) -> Dict[str, Genome]:
        """加载 GenomeLoadStage 产出, 不存在则返回空字典。"""
        ...

    def save_intra_results(self, results: Dict[str, GenomeIntraResult]) -> None:
        """保存 IntraGenomeStage 产出(genome_id -> GenomeIntraResult)。"""
        ...

    def load_intra_results(self) -> Dict[str, GenomeIntraResult]:
        """加载 IntraGenomeStage 产出, 不存在则返回空字典。"""
        ...

    def save_join_result(self, result: IncrementalResult) -> None:
        """保存 GlobalJoinStage 产出。"""
        ...

    def load_join_result(self) -> Optional[IncrementalResult]:
        """加载 GlobalJoinStage 产出, 不存在则返回 None。"""
        ...

    def save_complete_ssn(self, ssn: CompleteSSN) -> None:
        """保存 SSNBuilderStage 产出。"""
        ...

    def load_complete_ssn(self) -> Optional[CompleteSSN]:
        """加载 SSNBuilderStage 产出, 不存在则返回 None。"""
        ...

    def save_families(self, families: List[OrthologFamily]) -> None:
        """保存 FamilyRefineStage 产出。"""
        ...

    def load_families(self) -> List[OrthologFamily]:
        """加载 FamilyRefineStage 产出, 不存在则返回空列表。"""
        ...
