"""MCLClustering: ClusteringPort 的外部 mcl 命令行实现。

设计依据: docs/归档/STEP6_MCL亚家族细分.md + 2026-10-02 与用户确认的方案:
- 直接调用系统已安装的外部 `mcl` 二进制(用户已在 conda `linc` 环境中装好),
  不做 `markov_clustering` Python 库的 fallback(与归档文档"优先
  markov_clustering 库"的原始设想不同,按用户本次的明确指示,以系统命令
  为唯一实现)。
- 输入: networkx.Graph(节点=基因id, 边带 weight 属性), 转换为 mcl 的
  `.abc` 格式(每行 "source\\ttarget\\tweight")。
- 输出: mcl 聚类结果(每行一个簇, tab 分隔的成员id列表), 解析为
  List[Set[str]], 与 ClusteringPort.cluster() 的返回类型一致。
- 中间文件(.abc 输入/输出)持久化到 `work_dir` 下以 index_id 命名的子目录,
  不自动清理, 与 mmseqs/diamond adapter 的既有风格保持一致(2026-10-02
  与用户确认)。
"""
import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from typing import List, Optional, Set

import networkx as nx

logger = logging.getLogger(__name__)

# 持久化中间文件(.abc 输入/mcl 输出)的默认工作目录,运行期间不自动清理,
# 与 mmseqs_db/diamond_db 的既有约定保持一致(2026-10-02 与用户确认)。
_DEFAULT_WORK_DIR = os.path.join("output", ".linc_work", "mcl")


@dataclass
class MCLConfig:
    """mcl_adapter 的可配置参数。"""

    work_dir: str = _DEFAULT_WORK_DIR


class MCLClustering:
    """ClusteringPort 的外部 mcl 命令行实现(Step6 亚家族细分用途)。"""

    def __init__(self, config: Optional["MCLConfig"] = None) -> None:
        self._config = config or MCLConfig()
        if shutil.which("mcl") is None:
            raise RuntimeError(
                "未找到 mcl 可执行文件, 请确认已安装 mcl 并添加到 PATH"
                "(推荐: conda install -c bioconda mcl), 或确认当前终端"
                "已激活正确的 conda 环境。"
            )

    def cluster(
        self,
        graph: nx.Graph,
        inflation: float = 1.5,
        index_id: Optional[str] = None,
        **params,
    ) -> List[Set[str]]:
        """用外部 mcl 命令对 graph 做聚类, 返回每个簇的节点集合。

        调用方(FamilyRefineStage)负责在调用前过滤掉无边的连通分量(单个
        孤立节点), 直接作为独立家族处理, 不经过本方法(与用户确认,
        2026-10-02: mcl 对无边输入没有意义)。这里仍保留一个兜底分支,
        避免误传空边图时直接崩溃。

        index_id 用于区分不同连通分量各自的持久化子目录(如
        "component_1"), 同一 index_id 重跑时覆盖旧结果, 与
        mmseqs/diamond adapter 的 `search(tag=...)` 风格保持一致。
        """
        if graph.number_of_edges() == 0:
            return [{node} for node in graph.nodes]

        index_id = index_id or "default"
        work_dir = os.path.join(self._config.work_dir, index_id)
        shutil.rmtree(work_dir, ignore_errors=True)
        os.makedirs(work_dir, exist_ok=True)

        abc_path = os.path.join(work_dir, "input.abc")
        output_path = os.path.join(work_dir, "output.mcl")
        self._write_abc(graph, abc_path)

        cmd = [
            "mcl", abc_path, "--abc",
            "-I", str(inflation),
            "-o", output_path,
        ]
        logger.debug("调用 mcl: %s", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"mcl 执行失败(returncode={exc.returncode}): {exc.stderr}"
            ) from exc

        return self._parse_output(output_path)

    @staticmethod
    def _write_abc(graph: nx.Graph, path: str) -> None:
        # 按 (u, v) 字典序排序后再写入, 消除 networkx 图边遍历顺序(取决于
        # 节点/边插入顺序)带来的非确定性, 与项目中 mmseqs/diamond 结果排序
        # 的既有原则保持一致。
        edges = sorted(
            ((min(u, v), max(u, v), data.get("weight", 1.0))
             for u, v, data in graph.edges(data=True)),
        )
        with open(path, "w", encoding="utf-8") as fh:
            for u, v, weight in edges:
                fh.write(f"{u}\t{v}\t{weight}\n")

    @staticmethod
    def _parse_output(path: str) -> List[Set[str]]:
        clusters: List[Set[str]] = []
        if not os.path.exists(path):
            return clusters
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                members = set(line.split("\t"))
                clusters.append(members)
        return clusters
