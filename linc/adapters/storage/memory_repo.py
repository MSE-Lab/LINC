"""MemoryRepository: RepositoryPort 的纯内存实现。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.3 节: "纯内存实现, 小规模/测试场景,
不支持跨进程 resume"。所有状态只保存在进程内存中, 进程结束即丢失, 因此
`load_checkpoint` 在同一进程内可用于"同次运行中途重试", 但不能用于真正的
跨进程断点续跑(该需求由 adapters/storage/sqlite_repo.py 负责, 见开发进度第8步)。
"""
import copy
from typing import Dict, Optional, Set

from ...core.entities import Cluster, SSNEdge


class MemoryRepository:
    """RepositoryPort 的内存实现。"""

    def __init__(self) -> None:
        self._checkpoints: Dict[str, dict] = {}
        self._global_library: Dict[str, Cluster] = {}
        self._ssn_edges: Set[SSNEdge] = set()

    def save_checkpoint(self, stage: str, state: dict) -> None:
        # 深拷贝, 避免调用方后续修改原 dict 影响已保存的快照
        self._checkpoints[stage] = copy.deepcopy(state)

    def load_checkpoint(self, stage: str) -> Optional[dict]:
        state = self._checkpoints.get(stage)
        return copy.deepcopy(state) if state is not None else None

    def save_global_library(self, clusters: dict) -> None:
        self._global_library = dict(clusters)

    def load_global_library(self) -> dict:
        return dict(self._global_library)

    def save_ssn_edges(self, edges: set) -> None:
        self._ssn_edges = set(edges)
