"""MCLClustering 单元测试。

通过 monkeypatch 替换 shutil.which/subprocess.run, 不依赖真实安装的 mcl
二进制程序, 可在任意环境下运行(包括未安装 mcl 的 CI)。
"""
import os
import subprocess

import networkx as nx
import pytest

from linc.adapters.clustering.mcl_adapter import MCLClustering, MCLConfig
from linc.ports.clustering_port import ClusteringPort


def test_missing_binary_raises_clear_error(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(RuntimeError, match="未找到 mcl"):
        MCLClustering()


def test_satisfies_clustering_port(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mcl")
    assert isinstance(MCLClustering(), ClusteringPort)


def _fake_run_factory(cluster_lines):
    """构造一个假的 subprocess.run, 把预设的 mcl 输出内容写入结果文件。"""

    def _fake_run(cmd, check, capture_output, text):
        output_path = cmd[cmd.index("-o") + 1]
        with open(output_path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(cluster_lines))
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    return _fake_run


def test_cluster_writes_abc_and_parses_output(tmp_path, monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mcl")
    cluster_lines = ["A\tB", "C"]
    monkeypatch.setattr(subprocess, "run", _fake_run_factory(cluster_lines))

    graph = nx.Graph()
    graph.add_edge("A", "B", weight=0.9)
    graph.add_edge("B", "C", weight=0.3)

    clustering = MCLClustering(MCLConfig(work_dir=str(tmp_path / "mcl_work")))
    clusters = clustering.cluster(graph, inflation=1.5, index_id="component_1")

    assert clusters == [{"A", "B"}, {"C"}]

    abc_path = tmp_path / "mcl_work" / "component_1" / "input.abc"
    assert abc_path.exists()
    content = abc_path.read_text(encoding="utf-8")
    # 按 (min(u,v), max(u,v)) 排序写入, 消除边遍历顺序的非确定性。
    assert content.splitlines() == ["A\tB\t0.9", "B\tC\t0.3"]


def test_cluster_passes_inflation_to_command(tmp_path, monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mcl")
    recorded_cmds = []

    def _fake_run(cmd, check, capture_output, text):
        recorded_cmds.append(cmd)
        output_path = cmd[cmd.index("-o") + 1]
        with open(output_path, "w", encoding="utf-8") as fh:
            fh.write("A\tB")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_run)

    graph = nx.Graph()
    graph.add_edge("A", "B", weight=1.0)

    clustering = MCLClustering(MCLConfig(work_dir=str(tmp_path / "mcl_work")))
    clustering.cluster(graph, inflation=4.0, index_id="component_1")

    [cmd] = recorded_cmds
    assert cmd[cmd.index("-I") + 1] == "4.0"


def test_cluster_overwrites_previous_results_for_same_index_id(tmp_path, monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mcl")

    work_dir = str(tmp_path / "mcl_work")
    component_dir = os.path.join(work_dir, "component_1")
    os.makedirs(component_dir, exist_ok=True)
    stale_path = os.path.join(component_dir, "stale.txt")
    with open(stale_path, "w", encoding="utf-8") as fh:
        fh.write("stale")

    monkeypatch.setattr(subprocess, "run", _fake_run_factory(["A\tB"]))

    graph = nx.Graph()
    graph.add_edge("A", "B", weight=1.0)

    clustering = MCLClustering(MCLConfig(work_dir=work_dir))
    clustering.cluster(graph, inflation=1.5, index_id="component_1")

    assert not os.path.exists(stale_path)


def test_empty_edge_graph_returns_singleton_clusters_without_calling_subprocess(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/mcl")

    def _should_not_be_called(*args, **kwargs):
        raise AssertionError("subprocess.run 不应被调用")

    monkeypatch.setattr(subprocess, "run", _should_not_be_called)

    graph = nx.Graph()
    graph.add_node("A")
    graph.add_node("B")

    clustering = MCLClustering()
    clusters = clustering.cluster(graph, inflation=1.5)

    assert clusters == [{"A"}, {"B"}]
