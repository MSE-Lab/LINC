"""GenomeLoadStage: 多基因组输入阶段。

对应 v1 Step1(docs/归档/STEP1_多基因组输入.md)与设计文档 2.4 节的接口契约:
    输入: list[str](FASTA 文件路径)
    输出: dict[str, Genome]

本阶段职责:
    1. 逐个解析 FASTA 文件,提取基因序列,构造 Genome/Gene 实体。
    2. 单个文件解析失败(不存在/为空/格式异常)不中断整体流程,
       该基因组标记为 GenomeStatus.FAILED 并记录明确的日志提示,
       继续处理其余文件(增量场景下一个文件坏了不该拖垮全局)。
    3. 若多个文件解析出相同的 genome_id(文件名 stem 相同, 如 foo.fa 与
       foo.faa), 则视为无法自行消歧的歧义, 直接抛出 GenomeIdConflictError
       中断整个流程, 要求用户手动解决目录中的文件冲突后重新运行。

不引入 biopython 等第三方依赖,使用标准库手写一个极简 FASTA 解析器。
"""
import glob
import logging
import os
from collections import defaultdict
from typing import Dict, Iterator, List, Optional, Set, Tuple

from ..core.entities import Gene, Genome, GenomeStatus

logger = logging.getLogger(__name__)

_FASTA_EXTENSIONS = ("*.fasta", "*.fa", "*.faa")


class GenomeIdConflictError(ValueError):
    """多个 FASTA 文件解析出相同的 genome_id(文件名 stem 相同)。

    例如同目录下存在 foo.fa 与 foo.faa: 二者会被视为"同一个基因组",
    这通常意味着目录里混入了不该共存的文件(同一基因组的不同格式拷贝、
    新旧版本等),程序无法自行判断该保留哪一个,因此直接中断整个流程,
    要求用户先手动解决目录中的文件冲突后再重新运行(2026-10-02 与用户确认)。
    """

    def __init__(self, conflicts: Dict[str, List[str]]) -> None:
        self.conflicts = conflicts
        detail = "; ".join(
            f"{genome_id} -> {paths}" for genome_id, paths in conflicts.items()
        )
        super().__init__(f"检测到重名的 genome_id, 已中断加载: {detail}")


def discover_fasta_files(directory: str) -> List[str]:
    """扫描目录下的 FASTA 文件(*.fasta/*.fa/*.faa),按路径排序后返回。

    供未来 cli.py 的 -i 目录参数调用:先用本函数把目录展开成文件列表,
    再传给 GenomeLoadStage.run()。GenomeLoadStage 本身的接口不因此改变。
    """
    paths: List[str] = []
    for pattern in _FASTA_EXTENSIONS:
        paths.extend(glob.glob(os.path.join(directory, pattern)))
    return sorted(set(paths))


def build_input_signature(file_paths: List[str]) -> List[Tuple[str, int, float]]:
    """构造 `-i` 输入文件集合的签名, 用于 `--resume`/`-a` 的一致性校验。

    签名内容为 `(basename, size, mtime)` 三元组列表(按 basename 排序),
    而不是单纯的 basename 集合(2026-10-03, 第二轮漏洞"签名只比 basename"
    修复)。选用 "文件大小+mtime" 而非完整内容 hash(如 MD5), 是与用户
    权衡性能后的折中方案: `os.stat` 几乎零开销, 能拦住"同名但内容已被
    替换"的绝大多数真实误用场景(文件大小或 mtime 通常会随内容变化而变),
    但不保证能拦住"刻意伪造相同大小+相同mtime"的极端场景——这类场景
    需要用户蓄意构造, 现实中概率极低, 为此引入全文件 MD5 校验的运行时
    开销(resume/append 每次启动都要重新读一遍所有历史文件内容)并不划算。
    """
    signature: List[Tuple[str, int, float]] = []
    for file_path in file_paths:
        stat = os.stat(file_path)
        signature.append((os.path.basename(file_path), stat.st_size, stat.st_mtime))
    return sorted(signature, key=lambda item: item[0])


def _parse_fasta(file_path: str) -> Iterator[Tuple[str, str]]:
    """极简 FASTA 解析器,逐条产出 (header, sequence)。

    header 为 ">" 之后到行尾的原始文本(本阶段不保留，仅解析时使用);
    sequence 为去除换行/空白后拼接的序列字符串。

    健壮性边界处理(2026-10-02, 与"不静默"原则保持一致):
    - 若文件开头在第一个 ">" 之前出现非空序列行(无合法 header 归属),
      记录 warning 后丢弃该片段, 不计入任何基因;
    - 若某条记录的序列为空(如连续两个 ">" 行, header 之间无序列),
      记录 warning 后跳过该记录, 不产出 length=0 的空序列。
    """
    header: str = ""
    seen_header: bool = False
    seq_chunks: List[str] = []
    with open(file_path, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header:
                    sequence = "".join(seq_chunks)
                    if sequence:
                        yield header, sequence
                    else:
                        logger.warning(
                            "FASTA 记录序列为空,已跳过: file_path=%s header=%r",
                            file_path, header,
                        )
                header = line[1:].strip()
                seen_header = True
                seq_chunks = []
            else:
                if not seen_header:
                    logger.warning(
                        "FASTA 文件在首个 '>' 之前出现序列行,已忽略: file_path=%s",
                        file_path,
                    )
                    continue
                seq_chunks.append(line)
        if seen_header:
            sequence = "".join(seq_chunks)
            if sequence:
                yield header, sequence
            else:
                logger.warning(
                    "FASTA 记录序列为空,已跳过: file_path=%s header=%r",
                    file_path, header,
                )


_GENOME_CODE_PREFIX = "G"
_GENOME_CODE_WIDTH = 5


class GenomeLoadStage:
    """对应 v1 Step1 的现代化实现。"""

    def run(
        self,
        genome_files: List[str],
        existing_genomes: Optional[Dict[str, Genome]] = None,
    ) -> Dict[str, Genome]:
        """解析一批 FASTA 文件,返回 genome_id -> Genome 映射。

        genome_id 为 LINC 内部生成的短编号(如 "G00001"),与原始文件名完全
        解耦(2026-10-03 与用户确认"方案B": 原文件名可能包含任意字符,不
        适合直接参与 gene_id 拼接,因此改为内部顺序编号); 原始文件名
        (去掉扩展名)保存在 `Genome.source_file_name`,仅用于展示、重名
        检测与 --resume/--append 一致性校验。

        gene_id 格式为 "{genome_id}|{local_id}",local_id 为基因在文件中的
        序号(从 1 开始); 原始 FASTA header 保存在 `Gene.original_header`,
        用于最终结果回填为用户可读的原始 ID。

        `existing_genomes`(2026-10-02, 第8步 `--append` 设计; 2026-10-03
        为配合 genome_id 短编号化改为传入完整的历史 genomes 字典): 若传入,
        表示历史(此前已完整运行过)的 genome_id -> Genome 映射。本次用于:
        ①按历史 `source_file_name` 检测新文件是否与历史重名(无法自行消歧
        的歧义, 直接中断,与本次新文件之间互相重名的处理方式一致);
        ②新 genome_id 的编号从历史最大编号之后接续分配, 保证不与历史
        genome_id 冲突、且历史已分配的 genome_id 保持不变。
        """
        existing_genomes = existing_genomes or {}
        existing_source_file_names = {
            genome.source_file_name for genome in existing_genomes.values()
        }
        self._check_filename_conflicts(genome_files, existing_source_file_names)

        next_code_num = self._next_code_number(existing_genomes)
        result: Dict[str, Genome] = {}
        for file_path in sorted(genome_files, key=self._extract_source_file_name):
            source_file_name = self._extract_source_file_name(file_path)
            genome_id = self._format_genome_code(next_code_num)
            next_code_num += 1
            genome = self._load_single_genome(genome_id, source_file_name, file_path)
            result[genome_id] = genome
        self._log_summary(result)
        return result

    @staticmethod
    def _format_genome_code(num: int) -> str:
        return f"{_GENOME_CODE_PREFIX}{num:0{_GENOME_CODE_WIDTH}d}"

    @classmethod
    def _next_code_number(cls, existing_genomes: Dict[str, Genome]) -> int:
        """取历史 genome_id 编号中的最大值 + 1, 作为本次新编号的起点;
        历史为空则从 1 开始。
        """
        max_num = 0
        prefix_len = len(_GENOME_CODE_PREFIX)
        for genome_id in existing_genomes:
            suffix = genome_id[prefix_len:]
            if genome_id.startswith(_GENOME_CODE_PREFIX) and suffix.isdigit():
                max_num = max(max_num, int(suffix))
        return max_num + 1

    def _check_filename_conflicts(
        self,
        genome_files: List[str],
        existing_source_file_names: Optional[Set[str]] = None,
    ) -> None:
        """检测原始文件名(stem)重名, 一旦发现立即中断整个流程。

        与单文件解析失败(不存在/为空/格式异常)不同, 文件名冲突不是
        单个文件自身的问题, 而是程序无法自行消歧的歧义, 因此不适用
        "单文件失败不中断"的既有策略, 必须直接报错停止。
        """
        files_by_name: Dict[str, List[str]] = defaultdict(list)
        for file_path in genome_files:
            source_file_name = self._extract_source_file_name(file_path)
            files_by_name[source_file_name].append(file_path)
        conflicts = {
            name: paths
            for name, paths in files_by_name.items()
            if len(paths) > 1
        }
        if existing_source_file_names:
            for name, paths in files_by_name.items():
                if name in existing_source_file_names:
                    conflicts[name] = conflicts.get(name, []) + paths + ["<历史已存在>"]
        if conflicts:
            logger.error("检测到重名的基因组文件名, 中断加载: %s", conflicts)
            raise GenomeIdConflictError(conflicts)

    @staticmethod
    def _extract_source_file_name(file_path: str) -> str:
        base_name = os.path.basename(file_path)
        source_file_name, _ = os.path.splitext(base_name)
        return source_file_name

    def _load_single_genome(self, genome_id: str, source_file_name: str, file_path: str) -> Genome:
        if not os.path.isfile(file_path):
            logger.warning(
                "基因组文件不存在,跳过并标记为 FAILED: genome_id=%s file_path=%s",
                genome_id, file_path,
            )
            return Genome(
                genome_id=genome_id,
                file_path=file_path,
                source_file_name=source_file_name,
                genes=[],
                status=GenomeStatus.FAILED,
            )

        genes: List[Gene] = []
        try:
            for local_id, (header, sequence) in enumerate(_parse_fasta(file_path), start=1):
                gene_id = f"{genome_id}|{local_id}"
                # original_header 取 header 第一个空白前的 token, 与 OrthoFinder
                # 习惯一致(2026-10-03 与用户确认"方案B"), 用于最终结果回填为
                # 用户可读的原始序列 ID。
                original_header = header.split()[0] if header else None
                genes.append(
                    Gene(
                        gene_id=gene_id,
                        genome_id=genome_id,
                        length=len(sequence),
                        sequence=sequence,
                        original_header=original_header,
                    )
                )
        except OSError:
            logger.warning(
                "基因组文件读取失败,标记为 FAILED: genome_id=%s file_path=%s",
                genome_id, file_path, exc_info=True,
            )
            return Genome(
                genome_id=genome_id,
                file_path=file_path,
                source_file_name=source_file_name,
                genes=[],
                status=GenomeStatus.FAILED,
            )

        if not genes:
            logger.warning(
                "基因组文件为空(未解析出任何序列),标记为 FAILED: genome_id=%s file_path=%s",
                genome_id, file_path,
            )
            return Genome(
                genome_id=genome_id,
                file_path=file_path,
                source_file_name=source_file_name,
                genes=[],
                status=GenomeStatus.FAILED,
            )

        headers_seen = [g.original_header for g in genes if g.original_header is not None]
        if len(headers_seen) != len(set(headers_seen)):
            logger.warning(
                "基因组文件内存在重复的原始序列 ID(original_header), "
                "回填 Orthogroups.txt 时将自动加 #序号 消歧: genome_id=%s file_path=%s",
                genome_id, file_path,
            )

        return Genome(
            genome_id=genome_id,
            file_path=file_path,
            source_file_name=source_file_name,
            genes=genes,
            status=GenomeStatus.LOADED,
        )

    @staticmethod
    def _log_summary(genomes: Dict[str, Genome]) -> None:
        total = len(genomes)
        loaded = [g for g in genomes.values() if g.status == GenomeStatus.LOADED]
        failed = [g for g in genomes.values() if g.status == GenomeStatus.FAILED]
        total_genes = sum(len(g.genes) for g in loaded)
        logger.info(
            "GenomeLoadStage 完成: 共 %d 个基因组, 成功 %d 个, 失败 %d 个, 成功基因组共 %d 条基因",
            total, len(loaded), len(failed), total_genes,
        )
        for g in failed:
            logger.warning("基因组加载失败: genome_id=%s file_path=%s", g.genome_id, g.file_path)
