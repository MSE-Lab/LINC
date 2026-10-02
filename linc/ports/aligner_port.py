"""AlignerPort: 序列比对引擎的抽象接口 (Protocol)。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.2 节。
使用 typing.Protocol (结构化子类型), 任何实现了对应方法签名的类都自动满足该接口,
不需要显式继承, 便于 adapters (mmseqs/diamond/mock) 解耦接入。
"""
from typing import Any, List, Optional, Protocol, Sequence, runtime_checkable

from ..core.entities import AlignmentResult, Gene


@runtime_checkable
class AlignerPort(Protocol):
    """序列比对引擎抽象接口。

    数值口径契约: 所有方法返回的 AlignmentResult.identity/query_cover/
    subject_cover 必须统一为 [0, 1] 范围的小数(详见 core.entities.
    AlignmentResult 的说明)。接入新的比对工具(如 blastp)时,若该工具原始
    输出是百分数(0-100),必须在对应 adapter 内部换算后再返回,不得依赖
    上层代码再次转换。

    build_index/search 是统一的"建库-查询"两步原语(2026-10-01 新增,
    详见 docs/LINC_技术框架设计_v2.md 第2.4节"GlobalJoinStage 性能方案"):
    - build_index: 为给定序列集合建立一个可复用的索引(builder 持久化到工作目录,
      不随调用结束而清理),返回不透明 handle,后续反复 search 时不再重新建库。
    - search: 用已建好的 index 作为 target, 对 query 做比对。
    IntraGenomeStage/GlobalJoinStage 均遵循该模式: 为基因组 build_index 一次,
    handle 随 GenomeIntraResult 向后传递复用, 避免反复重建目标库的开销。
    align_all_vs_all/align_query_vs_subject 作为无状态的便捷方法保留(内部等价于
    build_index+search 的组合), 供一次性/不需要复用索引的场景(如补边)使用。
    """

    def build_index(
        self, genes: Sequence[Gene], index_id: Optional[str] = None
    ) -> Any:
        """为 genes 建立可复用的比对索引(建库), 返回不透明 handle。

        index_id 用于指定持久化存储的标识(如 genome_id), 省略时使用临时标识。
        """
        ...

    def search(
        self,
        query: Sequence[Gene],
        index: Any,
        tag: Optional[str] = None,
        self_search: bool = False,
    ) -> List[AlignmentResult]:
        """用已建好的 index(作为 target/subject)对 query 做比对, 不重新建库。

        tag: 可选, 标识本次调用的性质/来源(如 "intra_self"/"join"/
        "full_from_<genome_id>"), 用于持久化中间文件时区分子目录,
        避免同一个 index 被不同阶段/不同来源反复调用时互相覆盖彼此的
        中间结果(2026-10-02 与用户确认)。省略时各 adapter 回退到旧的
        固定子目录名, 保持向后兼容。

        self_search: 可选, 当调用方明确告知"本次 query 就是 index 对应的
        同一份序列集合"(典型场景: IntraGenomeStage 的基因组内自比对)时
        置为 True(2026-10-02 与用户确认, 方案3)。此时实现应直接复用
        build_index 已经持久化的序列/数据库, 不再为 query 重新写入
        fasta/重新建库, 既节省磁盘(避免 genes.fasta 与 query.fasta 内容
        100%重复)也节省一次建库子进程调用。默认 False 保持向后兼容
        (query 与 index 来自不同序列集合的常规场景)。
        """
        ...

    def align_all_vs_all(self, genes: Sequence[Gene]) -> List[AlignmentResult]:
        """基因组内全对全比对 (Step2 用途)。"""
        ...

    def align_query_vs_subject(
        self, query: Sequence[Gene], subject: Sequence[Gene]
    ) -> List[AlignmentResult]:
        """新代表 vs 全局代表库 / 补边批量比对 (Step3/Step5 用途)。"""
        ...
