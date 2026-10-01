"""领域事件(可观测性基础)。

对应 docs/LINC_技术框架设计_v2.md 第 2.1 节。事件可用于日志记录、进度展示、
未来的实时可视化或审计追溯(例如回答"某个基因为什么被合并进了这个家族")。
"""
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Callable, Dict, List


@dataclass(frozen=True)
class DomainEvent:
    """所有领域事件的基类。"""

    timestamp: datetime
    payload: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ClusterCreated(DomainEvent):
    """基因组内自聚类阶段创建了一个新簇。"""


@dataclass(frozen=True)
class ClusterMerged(DomainEvent):
    """增量并入阶段发生了簇合并。"""


@dataclass(frozen=True)
class RepresentativeAdded(DomainEvent):
    """全局代表库新增了一个代表序列。"""


@dataclass(frozen=True)
class GenomeJoined(DomainEvent):
    """一个基因组完成了并入全局代表库的处理。"""


@dataclass(frozen=True)
class EdgeCompleted(DomainEvent):
    """补边阶段确认/新增了一条 SSN 边。"""


EventHandler = Callable[[DomainEvent], None]


class EventBus:
    """极简内存事件总线,订阅者可以是日志记录器/进度条/未来的可视化前端。"""

    def __init__(self) -> None:
        self._subscribers: List[EventHandler] = []

    def subscribe(self, handler: EventHandler) -> None:
        self._subscribers.append(handler)

    def publish(self, event: DomainEvent) -> None:
        for handler in self._subscribers:
            handler(event)
