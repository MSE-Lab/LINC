"""MockAligner: AlignerPort 的测试用假实现,不依赖任何外部二进制程序。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.3 节(adapters 一览表)。

两种工作模式:
1. **查表模式**: 调用方预先通过 `set_result` / 构造函数注入一张
   `(query_id, subject_id) -> AlignmentResult` 的结果表, 需要精确控制某对基因的
   比对结果时使用(例如测试双阈值决策树在 Cutoff 边界值附近的行为)。
2. **默认规则模式**: 对于未预设的基因对, 按"长度相似度"生成一个确定性的
   AlignmentResult(两条序列长度越接近, identity/coverage 越高), 保证同一对输入
   无论调用多少次结果都一致, 便于回归测试和端到端流程验证。
"""
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ...core.entities import AlignmentResult, Gene

_ResultTable = Dict[Tuple[str, str], AlignmentResult]


class MockAligner:
    """AlignerPort 的 Mock 实现。"""

    def __init__(self, preset_results: Optional[_ResultTable] = None) -> None:
        # 以 (query_id, subject_id) 为 key, 双向都会查表(比对结果视为无向)
        self._table: _ResultTable = dict(preset_results) if preset_results else {}

    def set_result(self, query_id: str, subject_id: str, result: AlignmentResult) -> None:
        """预设某一对基因的比对结果, 查表时对双向 key 均生效。"""
        self._table[(query_id, subject_id)] = result
        self._table[(subject_id, query_id)] = result

    def _lookup_or_compute(self, query: Gene, subject: Gene) -> AlignmentResult:
        key = (query.gene_id, subject.gene_id)
        if key in self._table:
            return self._table[key]
        return self._default_result(query, subject)

    @staticmethod
    def _default_result(query: Gene, subject: Gene) -> AlignmentResult:
        """根据序列长度差异生成一个确定性的"伪相似度"结果。

        长度完全相同 -> identity/coverage = 1.0;
        长度差异越大 -> identity/coverage 越低(线性衰减, 下限为 0.0)。
        """
        shorter = min(query.length, subject.length)
        longer = max(query.length, subject.length)
        ratio = shorter / longer if longer > 0 else 0.0
        alignment_length = shorter
        return AlignmentResult(
            query_id=query.gene_id,
            subject_id=subject.gene_id,
            identity=ratio,
            query_cover=ratio,
            subject_cover=ratio,
            evalue=0.0 if ratio == 1.0 else 1e-3 * (1.0 - ratio + 1e-9),
            bitscore=ratio * 100.0,
            alignment_length=alignment_length,
            # Mock 场景下无真实比对工具可取 qlen/tlen, 直接用 Gene.length
            # (2026-10-02 与用户确认, 为未来 MCL 权重归一化预留数据)。
            query_length=query.length,
            subject_length=subject.length,
        )

    def align_all_vs_all(self, genes: Sequence[Gene]) -> List[AlignmentResult]:
        """基因组内全对全比对(不含自身, 不重复生成 i,j 与 j,i)。"""
        results: List[AlignmentResult] = []
        for i in range(len(genes)):
            for j in range(i + 1, len(genes)):
                results.append(self._lookup_or_compute(genes[i], genes[j]))
        return results

    def align_query_vs_subject(
        self, query: Sequence[Gene], subject: Sequence[Gene]
    ) -> List[AlignmentResult]:
        """新代表 vs 全局代表库 / 补边批量比对。"""
        return [self._lookup_or_compute(q, s) for q in query for s in subject]

    def build_index(self, genes: Sequence[Gene], index_id: Optional[str] = None) -> Any:
        """Mock 场景无需真正建库, 直接以基因列表本身作为"索引" handle。"""
        return list(genes)

    def index_exists(self, index: Any) -> bool:
        """Mock 场景的"索引"只是内存中的基因列表, 不涉及真实文件, 只要该
        handle 非空即视为"存在"(第二轮漏洞5, 2026-10-03 与用户确认)。
        """
        return bool(index)

    def search(
        self,
        query: Sequence[Gene],
        index: Any,
        tag: Optional[str] = None,
        self_search: bool = False,
    ) -> List[AlignmentResult]:
        """用 build_index 返回的基因列表作为 target, 对 query 做笛卡尔积比对。

        tag/self_search 仅为满足 AlignerPort 接口签名而接受, Mock 实现不落盘
        任何中间文件, 因此直接忽略(2026-10-02 与用户确认)。
        """
        return [self._lookup_or_compute(q, s) for q in query for s in index]
