"""FileRepository: RepositoryPort 的基于文件(pickle)的持久化实现。

设计依据: 第8步(2026-10-02, 持久化存储 + Checkpoint/Resume + --append)
与用户确认: 不采用 sqlite, 改用最简单的文件持久化方案——把每个 Stage 的
产出对象用 pickle 序列化, 统一存放在 `<output_dir>/.linc_work/state/` 下,
与现有 `.linc_work/mmseqs_db`/`diamond_db` 等索引目录的持久化风格保持一致。

这些对象只是程序内部中间状态, 不要求人类可读, 因此选用 pickle(而非 json)
以避免为 dataclass/Set 等结构额外编写编解码器。

存储体积优化(2026-10-02, 与用户确认方案A): 大规模(数百至上千基因组)场景下
`genomes.pkl`/`intra_results.pkl` 体积可能较大(前者存有全部基因序列, 后者
存有 IntraGenomeStage 全对全比对的全部原始结果以供下游免费复用), 为降低
磁盘占用, 所有状态文件统一改为 gzip 压缩写入(序列/比对结果等高重复度文本
数据压缩比通常可达 3-4 倍), 不改变任何读取逻辑或数据结构本身。读取时优先
按 gzip 解压, 若文件是本次改动之前写入的旧版非压缩 pickle 文件, 自动回退为
直接 pickle 读取, 保持向后兼容(不强制用户手动迁移已有的 `.linc_work/state/`)。
"""
import gzip
import json
import logging
import os
import pickle
from typing import Any, Dict, List, Optional

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
        "config": "config.pkl",
    }

    def __init__(self, work_dir: str) -> None:
        """`work_dir` 通常为 `<output_dir>/.linc_work/state`。"""
        self._state_dir = work_dir
        os.makedirs(self._state_dir, exist_ok=True)
        self._append_history_path = os.path.join(self._state_dir, "append_history.jsonl")

    def _path(self, key: str) -> str:
        return os.path.join(self._state_dir, self._FILES[key])

    def _save(self, key: str, obj) -> None:
        # gzip 压缩写入(2026-10-02, 方案A): 大规模场景下显著降低磁盘占用,
        # 不改变整体"写临时文件再原子替换"的落盘语义。
        path = self._path(key)
        tmp_path = path + ".tmp"
        with gzip.open(tmp_path, "wb") as fh:
            pickle.dump(obj, fh)
        os.replace(tmp_path, path)

    def _load(self, key: str):
        path = self._path(key)
        if not os.path.exists(path):
            return None
        try:
            with gzip.open(path, "rb") as fh:
                return pickle.load(fh)
        except gzip.BadGzipFile:
            # 向后兼容(2026-10-02引入gzip压缩之前写入的旧版非压缩pickle文件):
            # 不视为损坏, 直接按原始pickle格式读取。
            try:
                with open(path, "rb") as fh:
                    return pickle.load(fh)
            except Exception as exc:  # noqa: BLE001 - 任何反序列化失败都视为损坏
                logger.error("持久化状态文件损坏, 无法反序列化: %s", path, exc_info=True)
                raise StateCorruptedError(
                    f"持久化状态文件损坏, 无法反序列化: {path}; 请检查 .linc_work/state/ "
                    f"目录是否被意外修改/截断, 不建议直接删除后重跑(会丢失已完成的工作)。"
                ) from exc
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
    # 签名元素类型为 `(basename, size, mtime)` 三元组(2026-10-03 起, 第二轮
    # 漏洞"签名只比 basename"修复); 历史遗留数据可能仍是纯 basename 字符串,
    # 读取时原样返回, 由调用方(cli.py)识别格式并做兼容处理。
    def save_input_signature(self, signature: List[Any]) -> None:
        self._save("input_signature", list(signature))

    def load_input_signature(self) -> Optional[List[Any]]:
        return self._load("input_signature")

    # ---- 算法参数一致性校验(工程审查问题1, 2026-10-03 与用户确认) ----
    def save_config(self, config: Dict[str, Any]) -> None:
        self._save("config", dict(config))

    def load_config(self) -> Optional[Dict[str, Any]]:
        return self._load("config")

    # ---- --append 审计轨迹(数据管理优化第1项, 2026-10-03 与用户确认) ----
    # 用 JSON Lines(每行一条独立的 JSON 对象)而非 pickle 存储: 人类可读,
    # 且天然支持"只追加、不需要整体反序列化重写"的写入模式。
    def append_audit_log(self, entry: Dict[str, Any]) -> None:
        with open(self._append_history_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(entry, ensure_ascii=False) + "\n")

    def load_append_history(self) -> List[Dict[str, Any]]:
        if not os.path.exists(self._append_history_path):
            return []
        history = []
        with open(self._append_history_path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    history.append(json.loads(line))
        return history
