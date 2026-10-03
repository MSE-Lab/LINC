"""MemoryRepository: RepositoryPort 的纯内存实现。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.3 节: "纯内存实现, 小规模/测试场景,
不支持跨进程 resume"。所有状态只保存在进程内存中, 进程结束即丢失, 因此
`load_*` 在同一进程内可用于"同次运行中途重试", 但不能用于真正的跨进程断点
续跑/增量追加(该需求由 adapters/storage/sqlite_repo.py 负责, 见开发进度第8步)。

第5步(pipeline打通)与用户确认, 接口按 engine/ 五个 Stage 的实际产出类型
逐一对应(而非早期版本里粗粒度的 checkpoint/global_library/ssn_edges)。
"""
import copy
from typing import Any, Dict, List, Optional

from ...core.entities import CompleteSSN, Genome, GenomeIntraResult, IncrementalResult, OrthologFamily


class MemoryRepository:
    """RepositoryPort 的内存实现。"""

    def __init__(self) -> None:
        self._genomes: Dict[str, Genome] = {}
        self._intra_results: Dict[str, GenomeIntraResult] = {}
        self._join_result: Optional[IncrementalResult] = None
        self._complete_ssn: Optional[CompleteSSN] = None
        self._families: List[OrthologFamily] = []
        self._has_genomes = False
        self._has_intra_results = False
        self._has_join_result = False
        self._has_complete_ssn = False
        self._has_families = False
        self._input_signature: Optional[List[Any]] = None
        self._append_history: List[Dict[str, Any]] = []

    def save_genomes(self, genomes: Dict[str, Genome]) -> None:
        self._genomes = copy.deepcopy(genomes)
        self._has_genomes = True

    def load_genomes(self) -> Dict[str, Genome]:
        return copy.deepcopy(self._genomes)

    def save_intra_results(self, results: Dict[str, GenomeIntraResult]) -> None:
        self._intra_results = copy.deepcopy(results)
        self._has_intra_results = True

    def load_intra_results(self) -> Dict[str, GenomeIntraResult]:
        return copy.deepcopy(self._intra_results)

    def save_join_result(self, result: IncrementalResult) -> None:
        self._join_result = copy.deepcopy(result)
        self._has_join_result = True

    def load_join_result(self) -> Optional[IncrementalResult]:
        return copy.deepcopy(self._join_result) if self._join_result is not None else None

    def save_complete_ssn(self, ssn: CompleteSSN) -> None:
        self._complete_ssn = copy.deepcopy(ssn)
        self._has_complete_ssn = True

    def load_complete_ssn(self) -> Optional[CompleteSSN]:
        return copy.deepcopy(self._complete_ssn) if self._complete_ssn is not None else None

    def save_families(self, families: List[OrthologFamily]) -> None:
        self._families = copy.deepcopy(families)
        self._has_families = True

    def load_families(self) -> List[OrthologFamily]:
        return copy.deepcopy(self._families)

    def has_genomes(self) -> bool:
        return self._has_genomes

    def has_intra_results(self) -> bool:
        return self._has_intra_results

    def has_join_result(self) -> bool:
        return self._has_join_result

    def has_complete_ssn(self) -> bool:
        return self._has_complete_ssn

    def has_families(self) -> bool:
        return self._has_families

    def save_input_signature(self, signature: List[Any]) -> None:
        self._input_signature = list(signature)

    def load_input_signature(self) -> Optional[List[Any]]:
        return list(self._input_signature) if self._input_signature is not None else None

    def append_audit_log(self, entry: Dict[str, Any]) -> None:
        self._append_history.append(copy.deepcopy(entry))

    def load_append_history(self) -> List[Dict[str, Any]]:
        return copy.deepcopy(self._append_history)
