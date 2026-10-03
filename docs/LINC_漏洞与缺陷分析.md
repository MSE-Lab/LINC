# LINC 漏洞与缺陷分析报告

> 分析日期:2026-10-01
> 分析方式:静态代码审查(通读 `linc/` 全部源码、`tests/` 测试、`docs/` 设计文档,并核对 git 历史)
> 分析范围:已实现的 core / ports / adapters / engine 三阶段 / cli

---

## 一、项目当前状态

架构清晰、分层规范(六边形架构 Ports & Adapters),是一个增量比较基因组、构建
序列相似性网络(SSN)并用 MCL 细分蛋白家族的工具。进度如下:

| 模块 | 状态 |
|---|---|
| `linc/core/`(实体/事件/策略) | ✅ 完成 |
| `linc/ports/`(aligner/clustering/storage 三接口) | ✅ 完成 |
| `linc/adapters/` mock / union-find / memory_repo / mmseqs / diamond | ✅ 完成(mcl 未做) |
| `linc/engine/` GenomeLoadStage / IntraGenomeStage / GlobalJoinStage | ✅ 完成 |
| `linc/cli.py`(`linc run`,支持 `--backend mock/mmseqs/diamond`) | ✅ 完成 |
| SSN 构建、MCL 细分、sqlite 持久化、pipeline 串联、config schema | ⬜ 未开始 |

此前反复排查的 **mmseqs 跨进程非确定性问题已被彻底修复**(`representative_genes()`
排序遍历 + `_parse_result_tsv()` 稳定排序 + `--split 1`),git log 与开发进度文档一致,
这部分没有问题。

---

## 二、发现的漏洞 / 缺陷(按严重度排序)

### 🔴 高危 1:Diamond 后端的 coverage 单位错误

**位置**:`linc/adapters/aligners/diamond_adapter.py` `_parse_result_tsv()`

```python
identity=float(pident) / 100.0,   # pident 是百分数,÷100 正确
query_cover=float(qcov),          # ❌ qcovhsp 也是百分数(0-100),漏了 ÷100
subject_cover=float(tcov),        # ❌ scovhsp 同理
```

diamond 的 `qcovhsp`/`scovhsp` 是**百分数(0–100)**,与 `pident` 同口径;而 mmseqs 的
`qcov`/`tcov` 是**小数(0–1)**。代码只对 `pident` 除了 100,导致 diamond 后端的 coverage
比 mmseqs/mock 大整整 100 倍。

**后果**:

- `global_join_stage.py` 的 `hit.coverage >= cutoff_i` 过滤形同虚设(任何有值的 coverage 都 ≥ 0.5);
- `AlignmentResult.coverage`、SSN 边权重(`identity × coverage`)全部偏大 100 倍;
- `--backend diamond` 与 `mmseqs` 结果不可比。

**为何测试没抓到**:`tests/adapters/test_diamond_adapter.py` 的假数据写的是 `0.9/0.95`
(小数),这本身就是对 diamond 输出格式的错误假设——真实 diamond 输出的是 `90.0/95.0`。

**建议修复**:改为 `query_cover=float(qcov) / 100.0`、`subject_cover=float(tcov) / 100.0`,
并同步把测试假数据改成百分数。

---

### 🔴 高危 2:Intra 与 Global 两个阶段的"纳入判据"不一致

**位置**:

- 基因组内:`linc/engine/intra_genome_stage.py` 只判 `hit.identity >= cutoff_i`;
- 跨基因组:`linc/engine/global_join_stage.py` 判 `hit.identity >= cutoff_i and hit.coverage >= cutoff_i`。

同一个相似度关系,在"基因组内聚类"和"并入全局库"用了不同标准:一对 80% identity、
5% coverage 的基因在基因组内会被聚成同簇,但跨基因组时却连不上。

**后果**:局部聚类与全局并入结果不自洽,与"Cutoff-I 是判定是否成簇的统一阈值"的设计
语义矛盾。

**建议修复**:明确 coverage 是否参与 Cutoff-I,并让两个阶段保持一致。

---

### 🟡 中危 3:genome_id 冲突导致静默数据丢失

**位置**:`linc/engine/genome_load_stage.py` `run()`

用文件名 stem 作 dict key,但同目录可同时存在 `foo.fa`、`foo.faa`、`foo.fasta`,它们
stem 都是 `foo`;后解析的文件会**无告警地覆盖**前者。`discover_fasta_files` 会把这些
文件都列出来,导致一个基因组被静默吞掉。

**建议修复**:检测 `genome_id` 重复并告警/重命名。

---

### 🟡 中危 4:CLI 参数完全无校验

**位置**:`linc/cli.py` 的 `cutoff_i` / `cutoff_ii` / `max_per_cluster` / `threads`

- `--cutoff-i 0.9 --cutoff-ii 0.5`(I > II)会产出语义错乱的结果;
- `--max-per-cluster 0` 或负数会使 `multi_representative_clustering` 的循环直接跳过
  (退化为恒 1 个代表);
- `--threads -1` 静默被忽略。

**建议修复**:在入口做参数校验。

---

### 🟢 低危 / 健壮性

1. **`get_cluster` 是 O(n) 线性扫描**:`linc/core/entities.py` 的 `GlobalRepresentativeLibrary.get_cluster`
   在 `global_join_stage.py` 的 `_decide_and_apply` 对每个命中都调用一次,整体退化为
   O(命中数 × 簇数),大样本下是性能瓶颈(设计文档自己也标注过"全局代表库查找结构要注意
   复杂度",但尚未落地)。建议改为 `dict` 索引。

2. **`_cluster_counter` 不复位**:同一 `GlobalJoinStage` 实例重复调用 `run`,cluster_id 会跨
   运行递增(不影响唯一性,但语义易误解)。

3. **解析容错不足**:`mmseqs_adapter.py` / `diamond_adapter.py` 的 `_parse_result_tsv` 只对
   "列数 ≠ 8"容错,一旦工具输出 `nan`/非数值列,`float()` 会抛 `ValueError` 直接中断整条
   流水线,而非跳过该行。

4. **FASTA 边界**:`genome_load_stage.py` 的 `_parse_fasta` 若文件开头无 `>` 直接是序列行,
   该段序列会被静默丢弃;连续两个 `>` 会生成 length=0 的 Gene(MockAligner 下 ratio=0)。

---

## 三、安全面评估

整体安全面很小,未发现可被利用的漏洞:

- `subprocess.run` 全部用 **list 参数 + `shell=False`**(无命令注入);
- 输入仅本地文件路径、无网络、无反序列化;
- git 跟踪文件中无硬编码密钥/凭证。

真正的风险集中在**正确性/一致性**而非安全。

---

## 四、优先级建议

| 优先级 | 缺陷 | 影响 |
|---|---|---|
| P0 | Diamond coverage 单位错误 | `--backend diamond` 结果错误、与 mmseqs 不可比 |
| P0 | Intra/Global 纳入判据不一致 | 局部聚类与全局并入不自洽 |
| P1 | genome_id 冲突静默覆盖 | 数据静默丢失 |
| P1 | CLI 参数无校验 | 无意义/退化结果 |
| P2 | `get_cluster` O(n) | 大样本性能 |
| P2 | 解析容错不足 | 异常输出导致流水线崩溃 |
| P2 | `_cluster_counter` 不复位 / FASTA 边界 | 语义/健壮性 |
