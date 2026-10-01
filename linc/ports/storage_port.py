"""RepositoryPort: 存储/持久化的抽象接口 (Protocol)。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.2 节。
负责 checkpoint 读写、全局代表库读写、SSN 边集合保存, 具体实现见
adapters/storage/memory_repo.py (内存) 与 adapters/storage/sqlite_repo.py (持久化)。
"""
from typing import Optional, Protocol, runtime_checkable

from ..core.entities import Cluster, SSNEdge


@runtime_checkable
class RepositoryPort(Protocol):
    """存储/持久化抽象接口。"""

    def save_checkpoint(self, stage: str, state: dict) -> None:
        """保存指定 Stage 的断点状态。"""
        ...

    def load_checkpoint(self, stage: str) -> Optional[dict]:
        """加载指定 Stage 的断点状态, 不存在则返回 None。"""
        ...

    def save_global_library(self, clusters: dict) -> None:
        """保存全局代表库 (cluster_id -> Cluster)。"""
        ...

    def load_global_library(self) -> dict:
        """加载全局代表库 (cluster_id -> Cluster)。"""
        ...

    def save_ssn_edges(self, edges: set) -> None:
        """保存 SSN 边集合。"""
        ...
