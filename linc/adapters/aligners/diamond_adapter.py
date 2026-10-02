"""DiamondAligner: AlignerPort 的另一个真实实现,调用系统级 diamond 命令行工具。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.3 节(adapters 一览表)。
接口/组织方式与 `mmseqs_adapter.py` 保持一致(统一的"建库-查询"两步原语 +
一次性便捷方法),便于 cli.py 通过 `--backend` 切换且不影响调用方代码。

与 mmseqs 实现的差异说明:
- diamond 的 `blastp` 命令可以直接以 fasta 文件作为 query 输入, 不需要像
  mmseqs 那样先对 query 单独 createdb, 因此 `search()` 比 mmseqs 版本更简单
  (少了 query 侧建库这一步)。
- `makedb`(建库) 只需对 target/subject 执行一次, 产出 `<name>.dmnd` 文件,
  这正是 `build_index` 要持久化复用的对象。
"""
import logging
import os
import shutil
import subprocess
import tempfile
import uuid
from dataclasses import dataclass
from typing import Any, List, Optional, Sequence

from ...core.entities import AlignmentResult, Gene

logger = logging.getLogger(__name__)

# diamond --outfmt 6 自定义列顺序, 与 AlignmentResult 字段逐一对应
# (pident 为百分数, qcovhsp/scovhsp 为百分数, 需与 mmseqs 侧保持同口径)。
# qlen/slen(query/subject 全长)取自 diamond 自身输出, 为未来 MCL 权重归一化
# (如 OrthoFinder 式长度校正)预留, 而非回查 Gene.length(2026-10-02 与用户确认)。
_OUTFMT_COLUMNS = "qseqid sseqid pident qcovhsp scovhsp evalue bitscore length qlen slen"

# 持久化索引(build_index/search 原语)的默认工作目录,运行期间不自动清理,
# 贴合 LINC 断点续跑/增量添加数据的定位(与 mmseqs_adapter 保持一致)。
_DEFAULT_WORK_DIR = os.path.join("output", ".linc_work", "diamond_db")


@dataclass
class DiamondConfig:
    """diamond_adapter 的可配置参数。"""

    sensitivity: str = "sensitive"  # diamond 的敏感度档位: fast/sensitive/more-sensitive/very-sensitive/ultra-sensitive
    evalue: float = 1e-5
    threads: int = 0  # 0 = 不显式指定 --threads, 让 diamond 使用默认值
    work_dir: str = _DEFAULT_WORK_DIR  # build_index 持久化索引的存放目录


class DiamondAligner:
    """AlignerPort 的 diamond 实现。"""

    def __init__(self, config: "DiamondConfig" = None) -> None:
        self._config = config or DiamondConfig()
        if shutil.which("diamond") is None:
            raise RuntimeError(
                "未找到 diamond 可执行文件, 请确认已安装 diamond 并添加到 PATH"
                "(推荐: conda install -c bioconda diamond), 或确认当前终端"
                "已激活正确的 conda 环境。"
            )

    @staticmethod
    def _write_fasta(genes: Sequence[Gene], path: str) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for gene in genes:
                sequence = gene.sequence or ""
                fh.write(f">{gene.gene_id}\n{sequence}\n")

    def _sensitivity_flag(self) -> str:
        return f"--{self._config.sensitivity}"

    def _run_diamond_search(
        self, query_genes: Sequence[Gene], subject_genes: Sequence[Gene]
    ) -> List[AlignmentResult]:
        """一次性建库+比对(临时目录, 用完即删), 风格对齐 `_run_mmseqs_search`。

        未来若需要对"反复作为 subject 的全局代表库"做 DB 复用优化,
        只需改动此方法内部实现, 不影响外部调用方和 AlignerPort 接口。
        """
        if not query_genes or not subject_genes:
            return []

        with tempfile.TemporaryDirectory(prefix="linc_diamond_") as tmp_dir:
            query_fasta = os.path.join(tmp_dir, "query.fasta")
            subject_fasta = os.path.join(tmp_dir, "subject.fasta")
            db_path = os.path.join(tmp_dir, "db")
            result_tsv = os.path.join(tmp_dir, "result.tsv")

            self._write_fasta(query_genes, query_fasta)
            self._write_fasta(subject_genes, subject_fasta)

            makedb_cmd = ["diamond", "makedb", "--in", subject_fasta, "-d", db_path]
            search_cmd = self._build_search_cmd(query_fasta, db_path, result_tsv)

            for cmd in (makedb_cmd, search_cmd):
                logger.debug("调用 diamond: %s", " ".join(cmd))
                try:
                    subprocess.run(cmd, check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as exc:
                    raise RuntimeError(
                        f"diamond 命令执行失败(returncode={exc.returncode}): "
                        f"{exc.stderr}"
                    ) from exc

            return self._parse_result_tsv(result_tsv)

    def _build_search_cmd(self, query_fasta: str, db_path: str, result_tsv: str) -> List[str]:
        cmd = [
            "diamond", "blastp",
            "-q", query_fasta,
            "-d", db_path,
            "-o", result_tsv,
            "--outfmt", "6", *_OUTFMT_COLUMNS.split(),
            "-e", str(self._config.evalue),
            self._sensitivity_flag(),
            "--quiet",
        ]
        if self._config.threads > 0:
            cmd += ["--threads", str(self._config.threads)]
        return cmd

    @staticmethod
    def _parse_result_tsv(path: str) -> List[AlignmentResult]:
        results: List[AlignmentResult] = []
        if not os.path.exists(path):
            return results
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) != 10:
                    logger.warning("忽略格式异常的 diamond 输出行: %r", line)
                    continue
                query_id, subject_id, pident, qcov, tcov, evalue, bits, alnlen, qlen, slen = parts
                if query_id == subject_id:
                    # 过滤自身比对(全对全场景下 query/subject 来自同一份 fasta)。
                    continue
                try:
                    result = AlignmentResult(
                        query_id=query_id,
                        subject_id=subject_id,
                        identity=float(pident) / 100.0,
                        query_cover=float(qcov) / 100.0,
                        subject_cover=float(tcov) / 100.0,
                        evalue=float(evalue),
                        bitscore=float(bits),
                        alignment_length=int(alnlen),
                        query_length=int(qlen),
                        subject_length=int(slen),
                    )
                except ValueError:
                    # 列数正确但某个字段非数值(如 diamond 输出 nan/inf),
                    # 跳过该行而不中断整条流水线(2026-10-02 健壮性修复)。
                    logger.warning("忽略含非数值字段的 diamond 输出行: %r", line)
                    continue
                results.append(result)
        # diamond 与 mmseqs 类似, 不保证跨进程/跨运行的输出行序一致, 而上层代码
        # 多处依赖"先遇到的hit"做平局判定(如 GlobalJoinStage 的 hit_by_cluster),
        # 此处按 (query_id, subject_id) 字典序排序, 消除这个非确定性来源。
        results.sort(key=lambda r: (r.query_id, r.subject_id))
        return results

    def align_all_vs_all(self, genes: Sequence[Gene]) -> List[AlignmentResult]:
        """基因组内全对全比对(Step2 用途): query/subject 为同一份序列集合。"""
        return self._run_diamond_search(genes, genes)

    def align_query_vs_subject(
        self, query: Sequence[Gene], subject: Sequence[Gene]
    ) -> List[AlignmentResult]:
        """新代表 vs 全局代表库 / 补边批量比对(Step3/Step5 用途)。"""
        return self._run_diamond_search(query, subject)

    def build_index(
        self, genes: Sequence[Gene], index_id: Optional[str] = None
    ) -> Optional[str]:
        """为 genes 建立持久化的 diamond DB(makedb), 返回 `.dmnd` 路径(不含后缀)作为 handle。

        与 `_run_diamond_search`(一把梭, 每次都重建+清理)不同, 本方法把 DB 建在
        `DiamondConfig.work_dir` 下、以 index_id 命名的子目录, 调用结束后不清理,
        供后续 `search` 反复复用(不重新 makedb)。
        """
        if not genes:
            return None
        index_id = index_id or uuid.uuid4().hex
        index_dir = os.path.join(self._config.work_dir, index_id)
        os.makedirs(index_dir, exist_ok=True)
        fasta_path = os.path.join(index_dir, "genes.fasta")
        db_path = os.path.join(index_dir, "db")
        self._write_fasta(genes, fasta_path)
        cmd = ["diamond", "makedb", "--in", fasta_path, "-d", db_path]
        logger.debug("调用 diamond: %s", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"diamond makedb 执行失败(returncode={exc.returncode}): "
                f"{exc.stderr}"
            ) from exc
        return db_path

    def search(
        self,
        query: Sequence[Gene],
        index: Any,
        tag: Optional[str] = None,
        self_search: bool = False,
    ) -> List[AlignmentResult]:
        """用 build_index 建好的持久化 DB(index, target/subject) 对 query 做比对。

        diamond 的 query 可直接是 fasta 文件, 不需要像 mmseqs 那样先对 query
        单独建库, 因此这里只需写 query.fasta + 调用一次 blastp。

        比对的中间文件(query.fasta/result.tsv)持久化保存在 index 所在的
        工作目录下(`<index_dir>/search_<tag>/`, tag 省略时回退为
        `<index_dir>/search/`), 不自动清理, 与 build_index 的索引一并保留,
        贴合断点续跑/增量添加数据的定位(与 mmseqs_adapter 保持一致)。

        tag 用于区分"同一个 index 被不同阶段/不同来源反复调用"的场景
        (2026-10-02 与用户确认, 修复 intra自比对/global join/full模式补边
        三类调用因共用固定的 `search/` 子目录而互相覆盖彼此中间结果的问题):
        调用方应传入能唯一标识本次调用性质的 tag(如 "intra_self"/"join"/
        "full_from_<genome_id>"), 同一 index_dir 下不同 tag 各自独立保存,
        同一 tag 重跑时直接覆盖旧结果(贴合"同一标识覆盖旧结果"的既有语义)。

        self_search: 当 query 与 index 对应同一份序列集合时(IntraGenomeStage
        基因组内自比对)置为 True(2026-10-02 与用户确认, 方案3): 跳过重新
        写 query.fasta, 直接复用 build_index 时已经写在 index 同目录下的
        `genes.fasta` 作为 blastp 的查询输入, 避免与 genes.fasta 内容
        100%重复的 query.fasta 再落盘一份。
        """
        if not query or not index:
            return []

        index_dir = os.path.dirname(str(index))
        subdir_name = f"search_{tag}" if tag else "search"
        if index_dir:
            tmp_dir = os.path.join(index_dir, subdir_name)
        else:
            tmp_dir = tempfile.mkdtemp(prefix="linc_diamond_search_")
        result_tsv = os.path.join(tmp_dir, "result.tsv")
        # 每次重跑先清空旧的 search_<tag>/ 子目录, 实现"同一标识(tag + index
        # 所在目录名, 通常是 genome_id)覆盖旧结果"的语义, 与 mmseqs_adapter
        # .search() 一致。
        shutil.rmtree(tmp_dir, ignore_errors=True)
        os.makedirs(tmp_dir, exist_ok=True)

        if self_search and index_dir:
            # query 与 index 是同一份序列集合: 直接复用 build_index 已写好的
            # genes.fasta, 不再重新写一份内容完全相同的 query.fasta
            # (2026-10-02, 方案3, 与用户确认)。
            query_fasta = os.path.join(index_dir, "genes.fasta")
        else:
            query_fasta = os.path.join(tmp_dir, "query.fasta")
            self._write_fasta(query, query_fasta)

        search_cmd = self._build_search_cmd(query_fasta, index, result_tsv)
        logger.debug("调用 diamond: %s", " ".join(search_cmd))
        try:
            subprocess.run(search_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"diamond 命令执行失败(returncode={exc.returncode}): "
                f"{exc.stderr}"
            ) from exc

        return self._parse_result_tsv(result_tsv)
