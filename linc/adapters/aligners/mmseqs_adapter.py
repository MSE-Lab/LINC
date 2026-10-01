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
from dataclasses import dataclass
from typing import List, Sequence

from ...core.entities import AlignmentResult, Gene

logger = logging.getLogger(__name__)

# easy-search 输出列的顺序, 与 AlignmentResult 字段逐一对应(pident 为百分数, 需转为 0-1)。
_FORMAT_OUTPUT = "query,target,pident,qcov,tcov,evalue,bits,alnlen"


@dataclass
class MMseqsConfig:
    """mmseqs_adapter 的可配置参数。"""

    sensitivity: float = 4.5
    evalue: float = 1e-5
    threads: int = 0  # 0 = mmseqs 自动检测


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
                results.append(
                    AlignmentResult(
                        query_id=query_id,
                        subject_id=subject_id,
                        identity=float(pident) / 100.0,
                        query_cover=float(qcov),
                        subject_cover=float(tcov),
                        evalue=float(evalue),
                        bitscore=float(bits),
                        alignment_length=int(alnlen),
                    )
                )
        return results

    def align_all_vs_all(self, genes: Sequence[Gene]) -> List[AlignmentResult]:
        """基因组内全对全比对(Step2 用途): query/subject 为同一份序列集合。"""
        return self._run_mmseqs_search(genes, genes)

    def align_query_vs_subject(
        self, query: Sequence[Gene], subject: Sequence[Gene]
    ) -> List[AlignmentResult]:
        """新代表 vs 全局代表库 / 补边批量比对(Step3/Step5 用途)。"""
        return self._run_mmseqs_search(query, subject)
