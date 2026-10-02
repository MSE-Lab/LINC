"""FileRepository: RepositoryPort 的基于文件(pickle)的持久化实现。

设计依据: 第8步(2026-10-02, 持久化存储 + Checkpoint/Resume + --append)
与用户确认: 不采用 sqlite, 改用最简单的文件持久化方案——把每个 Stage 的
产出对象用 pickle 序列化, 统一存放在 `<output_dir>/.linc_work/state/` 下,
与现有 `.linc_work/mmseqs_db`/`diamond_db` 等索引目录的持久化风格保持一致。

这些对象只是程序内部中间状态, 不要求人类可读, 因此选用 pickle(而非 json)
以避免为 dataclass/Set 等结构额外编写编解码器。
"""
import logging
import os
import pickle
from typing import Dict, List, Optional

from ...core.entities import CompleteSSN, Genome, GenomeIntraResult, IncrementalResult, OrthologFamily

logger = logging.getLogger(__name__)


class StateCorruptedError(RuntimeError):
    """持久化状态文件存在但无法正常反序列化(文件损坏/格式不兼容)。

    与"文件不存在"严格区分: 不存在 -> 视为"该 Stage 尚未跑过", 按既有策略
    继续(resume 时从该 Stage 开始跑); 文件存在但损坏 -> 绝不能静默当作
    "不存在"去重跑, 否则会悄悄丢失已完成的工作或产生不一致的状态
    (2026-10-02 与用户确认的校验原则)。
    """


class FileRepository:
    """RepositoryPort 的文件(pickle)实现, 支持跨进程 resume/append。"""

    _FILES = {
        "genomes": "genomes.pkl",
        "intra_results": "intra_results.pkl",
        "join_result": "join_result.pkl",
        "complete_ssn": "complete_ssn.pkl",
        "families": "families.pkl",
        "input_signature": "input_signature.pkl",
    }

    def __init__(self, work_dir: str) -> None:
        """`work_dir` 通常为 `<output_dir>/.linc_work/state`。"""
        self._state_dir = work_dir
        os.makedirs(self._state_dir, exist_ok=True)

    def _path(self, key: str) -> str:
        return os.path.join(self._state_dir, self._FILES[key])

    def _save(self, key: str, obj) -> None:
        path = self._path(key)
        tmp_path = path + ".tmp"
        with open(tmp_path, "wb") as fh:
            pickle.dump(obj, fh)
        os.replace(tmp_path, path)

    def _load(self, key: str):
        path = self._path(key)
        if not os.path.exists(path):
            return None
        try:
            with open(path, "rb") as fh:
                return pickle.load(fh)
        except Exception as exc:  # noqa: BLE001 - 任何反序列化失败都视为损坏
            logger.error("持久化状态文件损坏, 无法反序列化: %s", path, exc_info=True)
            raise StateCorruptedError(
                f"持久化状态文件损坏, 无法反序列化: {path}; 请检查 .linc_work/state/ "
                f"目录是否被意外修改/截断, 不建议直接删除后重跑(会丢失已完成的工作)。"
            ) from exc

    def _exists(self, key: str) -> bool:
        return os.path.exists(self._path(key))

    # ---- GenomeLoadStage ----
    def save_genomes(self, genomes: Dict[str, Genome]) -> None:
        self._save("genomes", genomes)

    def load_genomes(self) -> Dict[str, Genome]:
        result = self._load("genomes")
        return result if result is not None else {}

    def has_genomes(self) -> bool:
        return self._exists("genomes")

    # ---- IntraGenomeStage ----
    def save_intra_results(self, results: Dict[str, GenomeIntraResult]) -> None:
        self._save("intra_results", results)

    def load_intra_results(self) -> Dict[str, GenomeIntraResult]:
        result = self._load("intra_results")
        return result if result is not None else {}

    def has_intra_results(self) -> bool:
        return self._exists("intra_results")

    # ---- GlobalJoinStage ----
    def save_join_result(self, result: IncrementalResult) -> None:
        self._save("join_result", result)

    def load_join_result(self) -> Optional[IncrementalResult]:
        return self._load("join_result")

    def has_join_result(self) -> bool:
        return self._exists("join_result")

    # ---- SSNBuilderStage ----
    def save_complete_ssn(self, ssn: CompleteSSN) -> None:
        self._save("complete_ssn", ssn)

    def load_complete_ssn(self) -> Optional[CompleteSSN]:
        return self._load("complete_ssn")

    def has_complete_ssn(self) -> bool:
        return self._exists("complete_ssn")

    # ---- FamilyRefineStage ----
    def save_families(self, families: List[OrthologFamily]) -> None:
        self._save("families", families)

    def load_families(self) -> List[OrthologFamily]:
        result = self._load("families")
        return result if result is not None else []

    def has_families(self) -> bool:
        return self._exists("families")

    # ---- -i 输入签名(供 --resume/--append 校验一致性) ----
    def save_input_signature(self, file_names: List[str]) -> None:
        self._save("input_signature", list(file_names))

    def load_input_signature(self) -> Optional[List[str]]:
        return self._load("input_signature")
