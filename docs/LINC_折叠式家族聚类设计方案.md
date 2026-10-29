# LINC 折叠式家族聚类(collapsed 模式)设计方案

> 状态: **草案, 待用户确认**(2026-10-29)。本文档尚未实施任何代码改动。
> 目的: 在保留现有"敏感模式"(补边 + 全 SSN 上 MCL)的基础上, 增加一条可选的快速路径,
> 缩小进入 MCL 的图规模, 解决上万基因组规模下 MCL 成为主要瓶颈的问题。

---

## 1. 背景与动机

- LINC 增量比对的最终目的是让 MCL 聚类结果尽可能接近全量 all-vs-all 的结果。
- 要达到 OrthoFinder 式的高质量聚类, 需要对边权做完整归一化, 但增量架构拿不到全量 hit 数据, 难以实现(`docs/LINC_归一化bitscore设计方案.md` 只能做长度层面的归一化)。
- LINC 真正的优势在于速度: 比对次数 O(基因组数) + 支持 `--append`。
- 在 10000+ 基因组规模下, 进入 MCL 的节点/边数随基因数增长, MCL 将成为瓶颈。
  同一家族的高相似成员(identity ≥ Cutoff-II)高度冗余, 完全可以折叠。

## 2. 现有算法回顾(据代码整理)

### 2.1 GenomeLoad
读取 FASTA, 基因 ID 格式 `{genome_id}|{local_id}`。

### 2.2 IntraGenomeStage(每个基因组独立)
1. 基因组内全对全比对(`build_index` + `search(self_search=True)`)。
2. 用 Cutoff-I(默认 0.5, 并要求 coverage ≥ Cutoff-C 0.5)过滤, Union-Find 得到局部簇。
3. 每个局部簇选代表:
   - 孤儿基因: 自己为代表。
   - 簇内任意两成员 identity ≥ Cutoff-II(0.8)且覆盖度达标(整体均质): 选 1 个中心代表。
   - 否则: 贪心覆盖(k-center 近似)选多个代表, 数量 ≤ `max_per_cluster`(默认 5);
     同基因组内局部代表彼此 identity 必然 < Cutoff-II。
4. 持久化索引供 GlobalJoin 复用; `pairwise_hits` 按 `edge_completion_mode` 裁剪后保存。

### 2.3 GlobalJoinStage(增量, 核心)
- Order-D 排序: (基因数 ↓, 局部代表数 ↓, genome_id ↑)。
- 每个基因组一次比对: 全局代表库整体为 query, 对该基因组持久化索引搜索
  (基因组级批处理快照, 同基因组内所有局部代表基于同一份快照决策)。
- 对每个局部代表, 按命中(过 Cutoff-I)的全局簇数分三路:
  1. 零簇命中: 新建簇, 自己为代表。
  2. 单簇命中: identity ≥ Cutoff-II 且 coverage ≥ Cutoff-C → `add_member_only`(仅成员);
     否则簇代表数 < `max_per_cluster` 时作为新代表加入, 达上限则仅作成员。
  3. ≥2 簇命中: 合并这些簇; 任一命中 ≥ Cutoff-II 则新基因仅作成员, 否则自己成为新代表。
- 每个命中的簇产生一条"局部代表 – 全局代表"的 SSNEdge(稀疏 SSN)。
- 接触判据(evalue)只记录基因组接触了哪些簇(`touching_map`), 不参与合并。

### 2.4 SSNBuilderStage(补边, 三档 `--edge-completion-mode`)
- `none`: 只有稀疏 SSN(代表–代表边)。
- `representative_member`(默认): 补各基因组内"局部代表–局部成员"边, 复用 `pairwise_hits`, 无新比对。
- `full`: 再补簇内全部基因对; 同基因组内复用, 跨基因组成员对才发起新比对
  (按 target 基因组聚合, 每个 target 基因组仅一次 `search`)。
- 所有补边统一过 `passes_inclusion_cutoff`。

### 2.5 FamilyRefineStage(MCL)
1. 把 `CompleteSSN.edges` 建成 networkx 图(权重由 `edge_weight_mode`, 默认 identity×coverage), 补入孤立基因节点。
2. 按连通分量拆分; 孤立节点直接成为 size=1 家族。
3. 对所有非平凡分量合并, **一次**调用外部 `mcl`(inflation 默认 1.5), 得到 `OG0000000…` 家族。

**现状瓶颈**: 进入 MCL 的节点 = 图中出现的全部基因, 边数随成员数增长(`full` 模式下簇内近似 O(N²))。

## 3. 新模式设计: 折叠 – 聚类 – 展开

建议开关: `--family-mode {full_graph, collapsed}`(默认 `full_graph`, 即现有行为)。

### 3.1 核心思路
1. 将高相似(identity ≥ Cutoff-II 且 coverage ≥ Cutoff-C)的成员折叠到其所依附的代表上, 不作为 MCL 节点, 也不为其补边。
2. MCL 只在代表(及按规则保留的非折叠基因)构成的小图上运行。
3. MCL 结束后, 按"成员 → 依附代表"的记录将成员展开到对应家族。

### 3.2 数据流
```
IntraGenome ─► GlobalJoin ─► SSNBuilder(collapsed 分支) ─► FamilyRefine ─► 展开
                │                         │                      │
        记录成员→代表+identity    只构造代表图, 不补折叠成员边   MCL 后展开成员
```

### 3.3 各阶段改动点(预计)
| 文件 | 改动 |
|------|------|
| `linc/core/entities.py` | 在 `IncrementalResult`/簇结构中增加"成员 → 依附代表 + identity"记录; 递增 `SCHEMA_VERSION` |
| `linc/engine/global_join_stage.py` | 在 `add_member_only` 及簇合并处记录依附关系; 合并后保证依附关系跟随最终簇归属 |
| `linc/engine/ssn_builder_stage.py` | 新增折叠分支: 仅构造代表图, 跳过被折叠成员的补边 |
| `linc/engine/family_refine_stage.py` | MCL 后增加展开步骤; `edge_density` 在折叠模式下的口径需说明 |
| `linc/cli.py` / `linc/pipeline.py` | 新增参数与组合校验(与 `--edge-completion-mode` 的关系) |
| `tests/` | Mock 比对小例子: 多代表簇、孤立代表、跨簇合并后展开、两种模式结果对照 |
| `docs/` | 更新 `LINC_技术框架设计_v2.md` 与 `开发进度.md` |

### 3.4 展开规则(草案)
- 全局成员(GlobalJoin 中 `add_member_only` 的局部代表)按其依附代表所属家族展开;
  多代表簇下依附到 identity 最高的代表。
- 局部成员(IntraGenome 局部簇内, Cutoff-I ≤ identity < Cutoff-II)随其局部代表一起展开。
- 簇合并后, 依附关系以最终簇内的代表为准, 展开时需用最终归属。

## 4. 风险与权衡
- 成员只依附单个代表, 可能丢失其与其他家族的相似信息, 旁系同源可能被合并。
- 簇合并与代表变更会影响已记录的依附关系, 需保证最终一致。
- 折叠模式下 `edge_density` 仅反映代表图, 与敏感模式不可直接比较。
- 速度收益须实测, 质量损失须量化(与 `full_graph` 模式及 OrthoFinder 对比 Rand index/ARI、家族数、单拷贝核心基因数)。

## 5. 评估计划
1. 只读统计脚本(由用户在终端运行): 统计 ≥ Cutoff-II 成员占比、MCL 输入节点数/边数(现有 benchmark 输出)。
2. 两种模式在同一数据集上对比: 耗时、峰值内存、MCL 输入规模、家族划分差异。
3. 与 OrthoFinder 结果对比, 确定质量损失是否可接受。

## 6. 待确认问题
1. "折叠"对象仅限 GlobalJoin 中 ≥ Cutoff-II 并入的成员, 还是也包括局部簇内(Cutoff-I ≤ identity < Cutoff-II)的成员?
2. 多代表簇里, 成员依附规则采用"identity 最高的代表"是否可以?
3. 新模式使用独立开关(`--family-mode`), 还是作为 `--edge-completion-mode` 的第四档?
4. 是否先编写只读统计脚本量化收益(由用户运行)?

> 以上问题确认后, 再按"大方案层面"流程给出最终实现清单并开始实现。
