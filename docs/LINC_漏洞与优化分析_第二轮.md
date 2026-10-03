# LINC 漏洞与优化分析(第二轮)

> 分析日期:2026-10-02
> 分析背景:本轮基于项目完成"5 阶段完整流水线 + 持久化 + resume/append + 全后端"后的最新源码
> 说明:第一轮缺陷(见 `docs/LINC_漏洞与缺陷分析.md`)已全部闭环,本文仅记录**本轮新发现**的漏洞与优化空间

---

## 一、算法 / 逻辑漏洞

### 🔴 漏洞 1:full 模式补边有缺口(同基因组、不同局部簇、同属一个全局簇的基因对不会被补边)

**位置**:`linc/engine/ssn_builder_stage.py`

full 模式声称"补全簇内(全局簇)全部基因对",但两条补边路径之间存在空白:

- `_add_known_intra_genome_edges()` 只遍历 `intra.local_clusters`,补的是**同一局部簇内**的
  rep-member / member-member 边;
- `_complete_cross_genome_member_edges()` 只处理 `genome_a != genome_b` 的**跨基因组**对。

于是:**来自同一基因组、但落在不同局部簇里、后来又被全局合并进同一个全局簇的两个基因,
两边都不覆盖**——它们的边会永久缺失。

**后果**:这种情形在真实数据中并不罕见(全局合并路径③经常把同基因组的不同局部簇拉进同一
全局簇),而这些基因对在 `pairwise_hits` 里其实已经比对过了(数据在,只是没被取用)。导致
full 模式得到的 SSN 并不"full",进而影响 MCL 细分结果。

**建议修复**:在 `_complete_cross_genome_member_edges()` 中增加"同基因组、跨局部簇"的复用补边
分支(直接从 `pairwise_hits` 取,无需新比对)。

---

### 🟠 漏洞 2:SSN 边的纳入判据前后不一致

**位置**:`linc/engine/global_join_stage.py` vs `linc/engine/ssn_builder_stage.py`

- 稀疏边(join 产生)经过 `passes_inclusion_cutoff()`(identity ≥ cutoff_i 且 coverage ≥ cutoff_c);
- `_add_known_intra_genome_edges()` 复用 `pairwise_hits` 时**完全不过滤**,只要存在命中(高于
  evalue)就加边。

**后果**:SSN 里会混入 `identity < cutoff_i` 的弱边,与"cutoff_i 定义相似性"的语义、以及与跨
基因组边的过滤口径不一致,会直接改变 MCL 看到的权重分布。

**建议修复**:确认补边是否也应统一走 `passes_inclusion_cutoff()`,保持 SSN 边纳入判据一致。

---

### 🟠 漏洞 3:cutoff_ii 高置信判定只看 identity,不看 coverage

**位置**:`linc/engine/global_join_stage.py` `_decide_and_apply()`

代码用 `best_hit.identity >= cutoff_ii` 判"高置信 → 只做成员、不当代表"。但 v1 设计文档里
Step3 的 cutoff_ii 是 `identity AND coverage >= 0.8` 双条件。

**后果**:默认 `cutoff_c=0.0`(覆盖度过滤关闭)下,"identity 0.9 但 coverage 0.1"的碎片匹配会被
误判为高置信,该基因只当 member、不再新增代表,丢失家族内亚型分辨率。

**建议**:确认这是有意简化还是重构时的回归;若需严格,cutoff_ii 应同时考虑 coverage。

---

### 🟡 漏洞 4:mmseqs 命中数被默认上限截断,违背"必须全量比对"

**位置**:`linc/adapters/aligners/mmseqs_adapter.py` `search()` 的 `search_cmd`

未显式设置 `--max-accept` / `--max-seqs`,mmseqs 两者默认均为 300,意味着单个 query 最多报
300 条命中。

**后果**:含大旁系家族(如转座子、>300 个近缘拷贝)的基因组,Step2"全对全"会被**静默截断**,
超过 300 的旁系边直接丢失。

**建议**:按"全量"目标显式调大或设 `--max-accept 0`(不限)并评估内存。

---

### 🟡 漏洞 5:`--resume` 只校验 state、不校验索引

**位置**:`linc/cli.py`(resume 签名校验)+ `linc/pipeline.py`(resume)

resume 校验的是 `.linc_work/state/*.pkl`,但加载的 `intra_results.genome_index` 指向
`.linc_work/mmseqs_db/<genome_id>/db`。

**后果**:若用户删了索引而 state 还在,resume 会在 GlobalJoinStage 的 `search` 处以 mmseqs
子进程报错崩溃,错误信息不友好。

**建议**:resume 前检查索引目录存在性并给出明确指引。

---

### 🟢 漏洞 6:MCL 自环(self-loop)需实测确认

**位置**:`linc/adapters/clustering/mcl_adapter.py` `_write_abc()`

`.abc` 输入只写 `u\tv\tweight`,没有自环。MCL 算法本身需要对角线自环;`mcl` 命令通常会自动补,
但依赖版本/参数。

**建议**:用一个小连通分量实测确认,或在文档中明确"已验证 mcl 自动补自环"。

---

### 🟢 漏洞 7:bitscore 的方向约定

**位置**:`linc/engine/ssn_builder_stage.py` `_complete_cross_genome_member_edges()`

full 模式跨基因组补边只搜 `genome_a→genome_b` 一个方向(靠字典序保证确定性),而 bitscore 是
有方向性的量。

**后果**:用 `--edge-weight-mode bitscore` 时,边权重会受"谁排前"影响,虽确定但不直观。

**建议**:在文档说明该方向约定。

---

## 二、性能优化空间

1. **full 模式是 O(G²) 次子进程调用**:`_complete_cross_genome_member_edges()` 每个基因组对都
   走一次 `search`(mmseqs 内部是 createdb + search + convertalis 三次子进程)。G=100 → ~5000 次
   子进程,G=500 → 百万级。可考虑把多个基因组对的 query 合并成更少的批量调用。

2. **`pairwise_hits` 全量落盘 + 常驻内存**:`intra_results.pkl` 存 O(命中数) 条原始比对结果,
   是 `genomes.pkl` 之外最大的状态文件。可**按 `edge_completion_mode` 裁剪**:`none` 模式根本
   不需要,`representative_member` 只需 rep-member 子集,只有 `full` 才需要全量。

3. **`family_refine_stage._build_graph` 对全图边做 `sorted()`**:对每条边构造
   `tuple(sorted(...))` 再全量排序只为确定性,百万级边时是 O(E log E) 且带大量临时对象分配。
   可考虑用规范化整数 ID 降低排序常数。

4. **序列全程常驻内存**:`genomes` dict 在整个 pipeline 生命周期持有所有 `Gene.sequence`,
   而 IntraGenomeStage 之后序列不再被使用。大样本下可在 intra 完成后把 `sequence` 置空释放。

---

## 三、数据管理优化空间

1. **append 没有审计轨迹**:`pipeline.append` 不记录"追加了哪些目录/哪些 genome_id/顺序",与
   项目"可追溯/可复现"目标相悖。建议持久化一份 append 历史(追加目录 + 时间戳 + 新增
   genome_id 列表)。

2. **resume/append 签名只比 basename**:`cli.py` 用 `os.path.basename` 集合做一致性校验,
   同名不同内容的文件会通过校验、静默错配。建议改用文件内容 hash(如 MD5)做签名。

3. **`.linc_work` 无 GC 机制**:索引、`search_*` 中间产物、mcl 中间文件、state pkl 全部无限
   累积,大样本下可达数十 GB。建议提供分层清理策略(索引保留以支持 append,`search_*` 中间
   产物可安全清理)。

4. **`FileRepository` 无 fsync**:`os.replace` 是原子的,但无 `fsync` 就换名,掉电仍可能留下
   损坏的 pkl,影响 resume 可靠性。

5. **`genomes.pkl` 仍是全序列冗余**:gzip(方案A)只缓解磁盘,本质还是把 FASTA 内容重复存了
   一份。更彻底是只存长度+路径,序列在需要时从原 FASTA 重读。

---

## 四、优先级建议

| 优先级 | 项 | 类型 |
|---|---|---|
| P0 | full 模式补边缺口(同基因组跨局部簇) | 正确性漏洞 |
| P0 | SSN 边纳入判据不一致 | 正确性漏洞 |
| P1 | cutoff_ii 只看 identity(默认 coverage 关闭下碎片匹配误判) | 正确性/设计 |
| P1 | mmseqs `--max-accept` 截断 | 正确性/性能 |
| P1 | full 模式 O(G²) 子进程 | 性能 |
| P2 | resume 索引缺失、签名只比 basename | 健壮性 |
| P2 | pairwise_hits 按模式裁剪、序列延迟释放、.linc_work GC | 数据管理 |

> 漏洞 1 和漏洞 2 建议优先处理,因为它们直接影响 full 模式产出的 SSN 完整性和 MCL 细分结果的正确性。
