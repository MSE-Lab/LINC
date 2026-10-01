"""可替换的决策策略接口。

对应 docs/LINC_技术框架设计_v2.md 第 2.1 节。使用 typing.Protocol(结构化子类型)
定义接口,具体算法实现(Order-D 排序、双阈值决策树等)由 engine/ 或 adapters/ 提供,
算法正确性的唯一权威来源为 docs/归档/LINC_算法设计文档.md 第 3 节。
"""
from dataclasses import dataclass, field
from typing import List, Protocol, Sequence

from .entities import AlignmentResult, Gene, Genome


class OrderingPolicy(Protocol):
    """Order-D 排序策略接口,默认实现见 adapters 或直接内置在 engine 中。"""

    def sort(self, genomes: Sequence[Genome]) -> List[Genome]:
        """按 Order-D 规则(基因总数→局部代表数→文件名字典序)对基因组排序。"""
        ...


@dataclass(frozen=True)
class ThresholdDecision:
    """双阈值决策树的输出结果。"""

    action: str            # "create_new" | "attach_no_rep" | "attach_add_rep" | "merge_no_rep" | "merge_add_rep"
    target_clusters: List[str] = field(default_factory=list)


class ThresholdPolicy(Protocol):
    """双阈值(Cutoff-I/Cutoff-II)决策策略接口。"""

    cutoff_i: float
    cutoff_ii: float

    def decide(self, new_rep: Gene, hits: List[AlignmentResult], lib: object) -> ThresholdDecision:
        """根据新代表序列与全局代表库的比对结果,给出决策动作。"""
        ...
