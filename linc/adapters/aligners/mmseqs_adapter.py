"""MMseqsAligner: AlignerPort 的真实实现,调用系统级 mmseqs2 命令行工具。

设计依据: docs/LINC_技术框架设计_v2.md 第 2.3 节(adapters 一览表) +
docs/归档/MMseqs2_Aligner.md(参数/输出字段映射)。

实现方式说明(已与用户确认):
- 统一使用 `mmseqs easy-search`(内部已封装 createdb/prefilter/align/convertalis
  四步调用,单次调用成本与手动三步调用完全一致,只是少写代码)。
- `align_all_vs_all` / `align_query_vs_subject` 均委托给私有方法
  `_run_mmseqs_search(query_fasta, subject_fasta)`——这是未来做"全局代表库
  DB 复用/分库优化"时的唯一改造点,不影响 AlignerPort 接口和调用方代码。
  当前版本不做 DB 复用,每次调用都对 subject 重新建库,属于已知的、留待后续
  遇到真实性能瓶颈时再优化的事项(与设计文档第6节"全局代表库查找结构要注意
  复杂度"的既定节奏一致)。
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

# easy-search 输出列的顺序, 与 AlignmentResult 字段逐一对应(pident 为百分数, 需转为 0-1)。
_FORMAT_OUTPUT = "query,target,pident,qcov,tcov,evalue,bits,alnlen"

# 持久化索引(build_index/search 原语)的默认工作目录,运行期间不自动清理,
# 贴合 LINC 断点续跑/增量添加数据的定位(与用户确认,2026-10-01)。
_DEFAULT_WORK_DIR = os.path.join("output", ".linc_work", "mmseqs_db")


@dataclass
class MMseqsConfig:
    """mmseqs_adapter 的可配置参数。"""

    sensitivity: float = 4.5
    evalue: float = 1e-5
    threads: int = 0  # 0 = mmseqs 自动检测
    work_dir: str = _DEFAULT_WORK_DIR  # build_index 持久化索引的存放目录


class MMseqsAligner:
    """AlignerPort 的 mmseqs2 实现。"""

    def __init__(self, config: "MMseqsConfig" = None) -> None:
        self._config = config or MMseqsConfig()
        if shutil.which("mmseqs") is None:
            raise RuntimeError(
                "未找到 mmseqs 可执行文件, 请确认已安装 mmseqs2 并添加到 PATH"
                "(推荐: conda install -c bioconda mmseqs2), 或确认当前终端"
                "已激活正确的 conda 环境。"
            )

    @staticmethod
    def _write_fasta(genes: Sequence[Gene], path: str) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for gene in genes:
                sequence = gene.sequence or ""
                fh.write(f">{gene.gene_id}\n{sequence}\n")

    def _run_mmseqs_search(
        self, query_genes: Sequence[Gene], subject_genes: Sequence[Gene]
    ) -> List[AlignmentResult]:
        """调用一次 `mmseqs easy-search`, 把结果解析为 AlignmentResult 列表。

        未来若需要对"反复作为 subject 的全局代表库"做 DB 复用/分库优化,
        只需要改动此方法内部实现, 不影响外部调用方和 AlignerPort 接口。
        """
        if not query_genes or not subject_genes:
            return []

        with tempfile.TemporaryDirectory(prefix="linc_mmseqs_") as tmp_dir:
            query_fasta = os.path.join(tmp_dir, "query.fasta")
            subject_fasta = os.path.join(tmp_dir, "subject.fasta")
            result_tsv = os.path.join(tmp_dir, "result.tsv")
            work_dir = os.path.join(tmp_dir, "tmp")
            os.makedirs(work_dir, exist_ok=True)

            self._write_fasta(query_genes, query_fasta)
            self._write_fasta(subject_genes, subject_fasta)

            cmd = [
                "mmseqs",
                "easy-search",
                query_fasta,
                subject_fasta,
                result_tsv,
                work_dir,
                "-s",
                str(self._config.sensitivity),
                "-e",
                str(self._config.evalue),
                "--format-output",
                _FORMAT_OUTPUT,
                "--split",
                "1",
                "-v",
                "1",
            ]
            if self._config.threads > 0:
                # mmseqs 要求 --threads >= 1, 0 表示"不显式指定",
                # 让 mmseqs 使用其自身默认的全部核心数。
                cmd += ["--threads", str(self._config.threads)]
            logger.debug("调用 mmseqs: %s", " ".join(cmd))
            try:
                subprocess.run(
                    cmd, check=True, capture_output=True, text=True
                )
            except subprocess.CalledProcessError as exc:
                raise RuntimeError(
                    f"mmseqs easy-search 执行失败(returncode={exc.returncode}): "
                    f"{exc.stderr}"
                ) from exc

            return self._parse_result_tsv(result_tsv)

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
                if len(parts) != 8:
                    logger.warning("忽略格式异常的 mmseqs 输出行: %r", line)
                    continue
                query_id, subject_id, pident, qcov, tcov, evalue, bits, alnlen = parts
                if query_id == subject_id:
                    # 过滤自身比对(全对全场景下 query/subject 来自同一份 fasta)。
                    continue
                try:
                    result = AlignmentResult(
                        query_id=query_id,
                        subject_id=subject_id,
                        identity=float(pident) / 100.0,
                        query_cover=float(qcov),
                        subject_cover=float(tcov),
                        evalue=float(evalue),
                        bitscore=float(bits),
                        alignment_length=int(alnlen),
                    )
                except ValueError:
                    # 列数正确但某个字段非数值(如 mmseqs 输出 nan/inf),
                    # 跳过该行而不中断整条流水线(2026-10-02 健壮性修复)。
                    logger.warning("忽略含非数值字段的 mmseqs 输出行: %r", line)
                    continue
                results.append(result)
        # mmseqs 输出的行序不保证跨进程/跨运行一致(已通过最小复现实验验证:
        # 相同输入两次运行的结果集合完全一致, 但行序不同), 而上层代码多处
        # 依赖"先遇到的hit"做平局判定(如 GlobalJoinStage 的 hit_by_cluster)。
        # 此处按 (query_id, subject_id) 字典序排序, 消除这个非确定性来源。
        results.sort(key=lambda r: (r.query_id, r.subject_id))
        return results

    def align_all_vs_all(self, genes: Sequence[Gene]) -> List[AlignmentResult]:
        """基因组内全对全比对(Step2 用途): query/subject 为同一份序列集合。"""
        return self._run_mmseqs_search(genes, genes)

    def align_query_vs_subject(
        self, query: Sequence[Gene], subject: Sequence[Gene]
    ) -> List[AlignmentResult]:
        """新代表 vs 全局代表库 / 补边批量比对(Step3/Step5 用途)。"""
        return self._run_mmseqs_search(query, subject)

    def build_index(
        self, genes: Sequence[Gene], index_id: Optional[str] = None
    ) -> Optional[str]:
        """为 genes 建立持久化的 mmseqs DB(createdb), 返回 DB 路径作为 handle。

        与 `_run_mmseqs_search`(easy-search 一把梭, 每次都重建+清理)不同,
        本方法把 DB 建在 `MMseqsConfig.work_dir` 下、以 index_id 命名的子目录,
        调用结束后不清理, 供后续 `search` 反复复用(不重新 createdb)。
        """
        if not genes:
            return None
        index_id = index_id or uuid.uuid4().hex
        index_dir = os.path.join(self._config.work_dir, index_id)
        os.makedirs(index_dir, exist_ok=True)
        fasta_path = os.path.join(index_dir, "genes.fasta")
        db_path = os.path.join(index_dir, "db")
        self._write_fasta(genes, fasta_path)
        cmd = ["mmseqs", "createdb", fasta_path, db_path, "-v", "1"]
        logger.debug("调用 mmseqs: %s", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"mmseqs createdb 执行失败(returncode={exc.returncode}): "
                f"{exc.stderr}"
            ) from exc
        return db_path

    def search(self, query: Sequence[Gene], index: Any) -> List[AlignmentResult]:
        """用 build_index 建好的持久化 DB(index, target/subject) 对 query 做比对。

        内部手动拆分 createdb(仅 query) → search → convertalis 三步,
        而非 easy-search 一把梭, 避免对已持久化的 index 重复建库。

        比对的中间文件(query.fasta/result.tsv 等)持久化保存在 index 所在的
        工作目录下(`<index_dir>/search/`), 不自动清理, 与 build_index 的
        索引一并保留, 贴合断点续跑/增量添加数据的定位(与用户确认,
        2026-10-01); 以 index 所在目录名(通常是 genome_id)为标识, 重跑
        同一 genome_id 时直接覆盖旧结果。
        """
        if not query or not index:
            return []

        index_dir = os.path.dirname(str(index))
        if index_dir:
            tmp_dir = os.path.join(index_dir, "search")
        else:
            tmp_dir = tempfile.mkdtemp(prefix="linc_mmseqs_search_")
        query_fasta = os.path.join(tmp_dir, "query.fasta")
        query_db = os.path.join(tmp_dir, "query_db")
        result_db = os.path.join(tmp_dir, "result_db")
        result_tsv = os.path.join(tmp_dir, "result.tsv")
        work_dir = os.path.join(tmp_dir, "tmp")
        # 每次重跑先清空旧的 query_db/result_db/work 子目录(mmseqs createdb/
        # search 遇到已存在的同名 DB 文件会报错), 同时实现"同一标识(index
        # 所在目录名, 通常是 genome_id)覆盖旧结果"的语义。
        shutil.rmtree(tmp_dir, ignore_errors=True)
        os.makedirs(work_dir, exist_ok=True)

        self._write_fasta(query, query_fasta)

        create_cmd = ["mmseqs", "createdb", query_fasta, query_db, "-v", "1"]
        search_cmd = [
            "mmseqs", "search", query_db, index, result_db, work_dir,
            "-s", str(self._config.sensitivity),
            "-e", str(self._config.evalue),
            "--split", "1",
            "-v", "1",
        ]
        if self._config.threads > 0:
            search_cmd += ["--threads", str(self._config.threads)]
        convert_cmd = [
            "mmseqs", "convertalis", query_db, index, result_db, result_tsv,
            "--format-output", _FORMAT_OUTPUT,
            "-v", "1",
        ]

        for cmd in (create_cmd, search_cmd, convert_cmd):
            logger.debug("调用 mmseqs: %s", " ".join(cmd))
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as exc:
                raise RuntimeError(
                    f"mmseqs 命令执行失败(returncode={exc.returncode}): "
                    f"{exc.stderr}"
                ) from exc

        return self._parse_result_tsv(result_tsv)
