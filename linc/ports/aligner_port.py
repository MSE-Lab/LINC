"""AlignerPort: 序列比对引擎的抽象接口 (Protocol)。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.2 节。
使用 typing.Protocol (结构化子类型), 任何实现了对应方法签名的类都自动满足该接口,
不需要显式继承, 便于 adapters (mmseqs/diamond/mock) 解耦接入。
"""
from typing import List, Protocol, Sequence, runtime_checkable

from ..core.entities import AlignmentResult, Gene


@runtime_checkable
class AlignerPort(Protocol):
    """序列比对引擎抽象接口。"""

    def align_all_vs_all(self, genes: Sequence[Gene]) -> List[AlignmentResult]:
        """基因组内全对全比对 (Step2 用途)。"""
        ...

    def align_query_vs_subject(
        self, query: Sequence[Gene], subject: Sequence[Gene]
    ) -> List[AlignmentResult]:
        """新代表 vs 全局代表库 / 补边批量比对 (Step3/Step5 用途)。"""
        ...
