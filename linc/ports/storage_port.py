"""RepositoryPort: 存储/持久化的抽象接口 (Protocol)。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.2 节, 并于第5步(pipeline打通)
与用户确认重新设计: 接口按 engine/ 五个 Stage 的实际产出类型逐一对应提供
save/load 方法(而非早期版本里粗粒度的 checkpoint/global_library/ssn_edges),
使 `Pipeline` 在每个 Stage 跑完后可以直接、类型安全地做存取。

第8步(2026-10-02, 持久化存储 + Checkpoint/Resume + --append)与用户确认:
不采用 sqlite, 改用基于文件(pickle)的持久化, 数据统一存放在 `.linc_work/
state/` 下(与现有 `.linc_work/mmseqs_db` 等索引目录风格一致)。为支持
`--resume`/`--append` 的核验逻辑, 本接口新增:
- `has_*()` 系列: 精确判断某个 Stage 的产出是否已被持久化保存过(不能仅靠
  `load_*()` 返回空集合/None 来判断, 因为空结果本身也可能是合法产出);
- `save_input_signature`/`load_input_signature`: 保存/加载本次(或历史)
  `-i` 参数对应目录下 discover 出的基因组文件签名列表, 供 `--resume`/
  `--append` 校验 `-i` 是否与历史一致。签名元素为 `(basename, size, mtime)`
  三元组(2026-10-03 起, 第二轮漏洞"签名只比 basename"修复, 与用户权衡性能后
  采用 大小+mtime 而非完整内容 MD5 的折中方案), 能拦住"同名但内容已变"的
  误用场景; 历史遗留数据可能仍是纯 basename 字符串列表, 读取时原样返回,
  由调用方(cli.py)识别格式并做兼容处理。

具体实现见 adapters/storage/memory_repo.py(纯内存, 不支持跨进程恢复)与
adapters/storage/file_repo.py(基于 pickle 文件, 支持跨进程 resume/append)。
"""
from typing import Any, Dict, List, Optional, Protocol, runtime_checkable

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

    def has_genomes(self) -> bool:
        """GenomeLoadStage 产出是否已被保存过(区别于"保存过但为空")。"""
        ...

    def has_intra_results(self) -> bool:
        """IntraGenomeStage 产出是否已被保存过。"""
        ...

    def has_join_result(self) -> bool:
        """GlobalJoinStage 产出是否已被保存过。"""
        ...

    def has_complete_ssn(self) -> bool:
        """SSNBuilderStage 产出是否已被保存过。"""
        ...

    def has_families(self) -> bool:
        """FamilyRefineStage 产出是否已被保存过。"""
        ...

    def save_input_signature(self, signature: List[Any]) -> None:
        """保存 `-i` 对应目录 discover 出的基因组文件签名列表(已按 basename 排序)。"""
        ...

    def load_input_signature(self) -> Optional[List[Any]]:
        """加载历史保存的 `-i` 文件签名列表, 不存在则返回 None。"""
        ...
