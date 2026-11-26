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

多文件原子批量提交(2026-11-26, 代码深度检查问题3, 与用户确认): `--append`
需要一次性更新 genomes/intra_results/join_result/complete_ssn/families/config
多个状态文件, 单个文件的写入虽然原子, 但文件之间没有整体原子性, 中途崩溃会
留下"新旧混杂"的 state(例如 genomes 已含新基因组, intra/join 仍是旧的, 之后
--resume 会静默漏掉新基因组)。`commit_state_batch()` 采用两阶段提交:
  1. 阶段一: 每个对象写成 `<name>.pkl.new`, 旧状态完全未动;
  2. 提交点: 原子写入标记文件 `commit_pending.json`(列出本批涉及的键);
  3. 阶段二: 逐个 `os.replace(<name>.pkl.new, <name>.pkl)`, 全部完成后删除标记。
`recover_pending_commit()` 在 resume/append/CLI 启动时调用: 标记存在 -> 说明
所有 .new 都已写完, 把剩余的 .new 改名前滚并删除标记(可重复执行/幂等); 标记
不存在但有零散 .new -> 崩溃发生在提交点之前, 旧状态完整, 直接删除这些 .new。
最终要么整批生效, 要么整批不生效。
"""
import gzip
import json
import logging
import os
import pickle
import uuid
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


class SchemaVersionMismatchError(StateCorruptedError):
    """持久化状态文件的 schema 版本与当前代码不兼容(工程审查问题2,
    2026-10-03 与用户确认)。

    代码升级后给某个 Stage 的产出 dataclass 新增/删除字段时, 旧版本的
    pickle 文件仍能"成功"反序列化, 但后续代码访问新字段会抛出语义不清的
    `AttributeError`。引入 `SCHEMA_VERSION` 标记后, 版本不匹配时改为在
    反序列化阶段就抛出这个明确的异常, 提示用户版本不兼容, 而不是让错误
    延迟到下游某个随机位置才暴露。
    """


class FileRepository:
    """RepositoryPort 的文件(pickle)实现, 支持跨进程 resume/append。"""

    # schema 版本控制(工程审查问题2, 2026-10-03 与用户确认): 每次给任意
    # 持久化 dataclass(GenomeIntraResult/IncrementalResult/CompleteSSN 等)
    # 新增/删除/改变语义的字段时, 必须手动递增此常量。`_save()` 会把该版本
    # 号与实际数据一并写入每个 pickle 文件; `_load()` 读取时校验版本号,
    # 不匹配则抛出 `SchemaVersionMismatchError`, 而不是让程序在下游某个
    # 随机位置因缺字段抛出语义不清的 `AttributeError`。
    # 2026-10-03 递增至2(与用户确认): "方案B"改动中 Genome 新增
    # source_file_name 字段、Gene 新增 original_header 字段, 且 genome_id
    # 的语义从"文件名 stem"变为"内部短编号(G00001...)", 属于本注释所述
    # "新增字段"与"改变语义"两种触发条件, 故递增版本号。
    # 2026-10-04 递增至3(接触判据/纳入判据解耦设计方案): IncrementalResult
    # 新增 touching_map 字段, 命中"新增字段"触发条件, 故递增版本号。
    SCHEMA_VERSION = 3

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
        self._run_id_path = os.path.join(self._state_dir, "run_id.pkl")
        self._run_id: Optional[str] = None

    _COMMIT_MARKER = "commit_pending.json"
    _NEW_SUFFIX = ".new"

    def _path(self, key: str) -> str:
        return os.path.join(self._state_dir, self._FILES[key])

    def _get_or_create_run_id(self) -> str:
        """获取本次(逻辑意义上的)运行标识(工程审查问题3, 2026-10-03 与用户
        确认)。首次调用(目录下尚无 `run_id.pkl`)时生成一个新 UUID 并落盘;
        resume/append 场景下该文件已存在, 直接复用, 代表与上次是"同一次
        逻辑运行"的延续。
        """
        if self._run_id is not None:
            return self._run_id
        if os.path.exists(self._run_id_path):
            with open(self._run_id_path, "rb") as fh:
                self._run_id = pickle.load(fh)
        else:
            self._run_id = str(uuid.uuid4())
            with open(self._run_id_path, "wb") as fh:
                pickle.dump(self._run_id, fh)
        return self._run_id

    def _save(self, key: str, obj) -> None:
        # gzip 压缩写入(2026-10-02, 方案A): 大规模场景下显著降低磁盘占用,
        # 不改变整体"写临时文件再原子替换"的落盘语义。
        # schema 版本标记(2026-10-03, 工程审查问题2): 实际数据包一层
        # {"__schema_version__": N, "payload": obj}, 供 `_load()` 校验。
        path = self._path(key)
        tmp_path = path + ".tmp"
        self._write_wrapped(tmp_path, obj)
        os.replace(tmp_path, path)

    def _write_wrapped(self, target_path: str, obj) -> None:
        """把 obj 包上 schema/run_id 标记后 gzip+pickle 写到 target_path。"""
        wrapped = {
            "__schema_version__": self.SCHEMA_VERSION,
            "__run_id__": self._get_or_create_run_id(),
            "payload": obj,
        }
        with gzip.open(target_path, "wb") as fh:
            pickle.dump(wrapped, fh)

    # ---- 多文件原子批量提交(2026-11-26, 代码深度检查问题3) ----
    def commit_state_batch(self, entries: Dict[str, Any]) -> None:
        """把多个状态对象作为一个整体原子提交(两阶段, 见模块 docstring)。

        `entries`: 状态键(genomes/intra_results/join_result/complete_ssn/
        families/input_signature/config)-> 待保存对象。键按字典顺序处理。
        """
        unknown = [key for key in entries if key not in self._FILES]
        if unknown:
            raise KeyError(f"未知的状态键: {unknown}")
        if not entries:
            return
        # 先处理上一次遗留的未完成提交, 避免与本批的 .new 文件混淆。
        self.recover_pending_commit()

        keys = list(entries)
        # 阶段一: 全部写成 .new, 旧状态完全未动。
        for key in keys:
            self._write_wrapped(self._path(key) + self._NEW_SUFFIX, entries[key])

        # 提交点: 原子写入标记文件。
        marker_path = os.path.join(self._state_dir, self._COMMIT_MARKER)
        marker_tmp = marker_path + ".tmp"
        with open(marker_tmp, "w", encoding="utf-8") as fh:
            json.dump({"keys": keys}, fh)
        os.replace(marker_tmp, marker_path)

        # 阶段二: 逐个改名生效, 全部完成后删除标记。
        self._roll_forward(keys)
        os.remove(marker_path)

    def _roll_forward(self, keys: List[str]) -> None:
        for key in keys:
            new_path = self._path(key) + self._NEW_SUFFIX
            if os.path.exists(new_path):
                os.replace(new_path, self._path(key))

    def recover_pending_commit(self) -> bool:
        """检查并处理上一次崩溃遗留的批量提交; 可重复调用(幂等)。

        返回 True 表示发现并处理了遗留(前滚或清理), False 表示无事可做。
        - 标记存在: 前滚——把剩余 .new 改名, 再删除标记;
        - 标记不存在但有零散 .new: 崩溃在提交点之前, 旧状态完整, 删除 .new。
        """
        marker_path = os.path.join(self._state_dir, self._COMMIT_MARKER)
        handled = False

        stray_marker_tmp = marker_path + ".tmp"
        if os.path.exists(stray_marker_tmp):
            os.remove(stray_marker_tmp)
            handled = True

        if os.path.exists(marker_path):
            try:
                with open(marker_path, "r", encoding="utf-8") as fh:
                    keys = list(json.load(fh)["keys"])
            except Exception as exc:  # noqa: BLE001
                raise StateCorruptedError(
                    f"批量提交标记文件损坏, 无法判断上次提交的范围: {marker_path}; "
                    f"请检查 .linc_work/state/ 目录, 不要直接删除后重跑。"
                ) from exc
            logger.warning("检测到上次未完成的批量提交, 正在前滚恢复: keys=%s", keys)
            self._roll_forward([key for key in keys if key in self._FILES])
            os.remove(marker_path)
            handled = True

        for key in self._FILES:
            new_path = self._path(key) + self._NEW_SUFFIX
            if os.path.exists(new_path):
                logger.warning("清理提交点之前遗留的临时文件: %s", new_path)
                os.remove(new_path)
                handled = True
        return handled

    def _load(self, key: str):
        path = self._path(key)
        if not os.path.exists(path):
            return None
        try:
            with gzip.open(path, "rb") as fh:
                raw = pickle.load(fh)
        except gzip.BadGzipFile:
            # 向后兼容(2026-10-02引入gzip压缩之前写入的旧版非压缩pickle文件):
            # 不视为损坏, 直接按原始pickle格式读取。
            try:
                with open(path, "rb") as fh:
                    raw = pickle.load(fh)
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
        return self._unwrap_schema(raw, path)

    def _unwrap_schema(self, raw, path: str):
        """校验并剥离 schema 版本标记(工程审查问题2, 2026-10-03 与用户确认)。

        - 带有 `__schema_version__` 标记(本次改动之后写入的新格式): 版本号
          必须与当前 `SCHEMA_VERSION` 一致, 否则抛出 `SchemaVersionMismatchError`。
        - 不带标记(本次改动之前写入的旧版文件, 无版本概念): 向后兼容, 仅打印
          warning 并原样返回数据, 不报错(与 `input_signature` 旧格式兼容的
          处理原则一致)。
        """
        if isinstance(raw, dict) and "__schema_version__" in raw and "payload" in raw:
            version = raw["__schema_version__"]
            if version != self.SCHEMA_VERSION:
                raise SchemaVersionMismatchError(
                    f"持久化状态文件版本不兼容: {path} 的 schema_version={version}, "
                    f"当前代码要求 schema_version={self.SCHEMA_VERSION}; 这通常发生在"
                    f"代码升级后旧的 .linc_work/state/ 目录未重新生成, 请删除该目录"
                    f"后重新完整运行(会丢失已完成的进度, 需要重新开始)。"
                )
            file_run_id = raw.get("__run_id__")
            if file_run_id is not None and file_run_id != self._get_or_create_run_id():
                raise StateCorruptedError(
                    f"检测到状态文件来自不同运行: {path} 的 run_id={file_run_id}, "
                    f"与当前 .linc_work/state/ 目录的 run_id={self._run_id} 不一致; "
                    f"这通常意味着 state 目录中的文件被手工替换/从别处拷贝拼接过, "
                    f"请删除 .linc_work 后重新完整运行, 不要混用不同运行的状态文件。"
                )
            return raw["payload"]
        logger.warning(
            "持久化状态文件 %s 不带 schema_version 标记(本次引入版本控制之前写入的"
            "旧版文件), 将按当前版本原样使用, 不做兼容性校验。",
            path,
        )
        return raw

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
