# LINC vs OrthoFinder gap 根因分析报告

> 日期:2026-10-03
> 数据:Streptomyces(15 个链霉菌物种,117,405 基因)、E. coli(15 株,70,451 基因)
> 结论:gap 的根因是**构图判据 + 搜索完备性**两层——OrthoFinder 用 RBB(互惠最佳命中)+ 全对全搜索,
> LINC 用 identity+coverage 阈值 + 增量搜索。

---

## 摘要

在 E. coli(同物种)上 LINC 与 OrthoFinder 的 family F1 为 0.959,但在 Streptomyces(跨物种)上
骤降到 0.864~0.888。经系统排查(权重归一化 → 阈值 → 搜索工具 → 搜索 evalue → 构图判据),
最终用**直接法**锁定并**决定性证实**:

> **OrthoFinder 的图 = RBB(reciprocal best hit,互惠最佳命中)稀疏化;LINC 的图 = identity+coverage 阈值。**
> **用 OrthoFinder 的全量 RBB 图 + MCL 可 100% 复现 OrthoFinder(F1 = 1.0000),而 LINC 的 identity+coverage SSN + MCL 只有 0.888。**
> **但进一步验证(第 7 节)发现:RBB 在 LINC 的增量搜索产物上只能算到全量的 3.3%——因为增量搜索
> 没比对足够多的基因对,RBB 判据拿不到完整的 hit 集。故 gap 的完整根因是两层:判据(RBB vs identity)
> + 搜索完备性(全对全 vs 增量)。**

---

## 1. 问题背景

LINC 与 OrthoFinder 在同口径(`Orthogroups.txt`,均为 MCL 聚类结果,不涉及基因树推断)下对比:

| 数据集 | 分歧度 | LINC 家族数 / OF 家族数 | family F1 |
|---|---|---|---|
| E. coli | 低(同物种) | 8,028 / 7,287(+10%) | 0.959 |
| Streptomyces | 高(跨物种) | 31,121 / 20,497(+52%) | 0.888(diamond 0.4) |

跨物种时 gap 显著放大,LINC 表现为**过度切分**(家族数远多于 OrthoFinder),且伴随大量孤立基因
(LINC 22% vs 正常水平)。

---

## 2. 实验数据(Streptomyces,全部 15 个物种)

| config | 搜索 | cutoff-i/c | 家族数 | F1 | precision | recall |
|---|---|---|---|---|---|---|
| 0.5 | mmseqs | 0.5/0.5 | 39,382 | 0.8635 | 0.982 | 0.770 |
| **0.4** | mmseqs | 0.4/0.4 | 31,121 | **0.8878** | 0.943 | 0.839 |
| 0.3 | mmseqs | 0.3/0.3 | 24,327 | 0.8467 | 0.835 | 0.859 |
| diamond 0.4 | diamond | 0.4/0.4 | 28,758 | 0.8830 | 0.915 | 0.853 |
| diamond 0.25 | diamond | 0.25/0.4 | 21,023 | 0.7519 | 0.699 | 0.814 |
| OrthoFinder(参照) | diamond | — | 20,497 | — | — | — |

要点:降 identity 阈值能收窄家族数,但 F1 出现"峰"(0.4 最优 0.888),继续降到 0.25 时家族数
虽撞上 20k,precision 却崩到 0.70——这是"撞巧"而非"做对"。

---

## 3. 假设排查(排除法)

### 3.1 边权重归一化(OrthoFinder 式长度校正)→ ✗ 排除

用 `offline_normalize.py` 对 3 个场景做 5 模式消融(identity_x_coverage / fixed / fit_a /
fit_ab / fit_ab_bin),结果:

| 场景 | 基线 F1 | 归一化 F1 范围 |
|---|---|---|
| E. coli 0.5 | 0.9592 | 0.958~0.959 |
| Streptomyces 0.5 | 0.8635 | 0.863~0.864 |
| Streptomyces 0.4 | 0.8878 | 0.886~0.888 |

**归一化 bitscore 在所有场景下都无效**,连零拟合的 `fixed`(a=0.5,b=0)模式都不动结果。
原因:`identity×coverage` 本身已长度无关,长度校正是冗余的。

### 3.2 搜索 evalue(1e-5 → 1e-3)→ ✗ 排除

对"被 LINC 漏掉的 254,013 条边"(only-OF)查 evalue 分布:

| evalue | 边数 | 占比 |
|---|---|---|
| ≤ 1e-5(搜索已返回,被 identity 阈值丢) | 247,718 | **97.5%** |
| 1e-5 ~ 1e-3(提 evalue 才能搜到) | 6,295 | 2.5% |
| > 1e-3 | 0 | 0% |

**97.5% 的漏边其实早已被 LINC 的搜索返回,只是被 identity≥0.4 丢弃。** 提高搜索 evalue 只能
多找回 2.5%,几乎无用——病灶不在搜索灵敏度。

### 3.3 搜索工具(mmseqs → diamond)→ ✗ 排除

diamond 0.4 找到更多边(家族数 28,758 vs mmseqs 的 31,121),recall 略升(0.853 vs 0.839),
但 precision 降(0.915 vs 0.943),**F1 反而略降(0.883 vs 0.888)**。多出的边真假混杂,identity
判据分不清。搜索工具不是根因。

### 3.4 identity 阈值高低 → ✗ 只能部分解释

0.5→0.4 提升 2.4pp,但 0.25 翻车(precision 崩)。identity 是一维裸百分比,无法区分"低一致但
显著(真同源)"与"低一致不显著(噪声)",因此**在这条轴上找不到两全的甜点**。

---

## 4. 根因:边判据(RBB vs identity+coverage)

### 4.1 直接法证据:边集对比(OrthoFinder vs LINC diamond 0.4,同一搜索工具)

| 边集 | 数量 | 占比 |
|---|---|---|
| OrthoFinder 总边 | 718,441 | — |
| LINC 总边 | 620,707 | — |
| 共同边 | 464,428 | 64.6% of OF |
| **仅 OrthoFinder(被 LINC 漏掉)** | **254,013** | **35.4% of OF** |
| 仅 LINC | 156,279 | 25.2% of LINC |

**仅-OF 边的 identity(pident)分布:**

| identity | 边数 | 占比 |
|---|---|---|
| < 30% | 30,310 | 11.9% |
| **30~40%** | **132,468** | **52.2%** |
| 40~50% | 28,153 | 11.1% |
| 50~60% | 19,998 | 7.9% |
| ≥ 60% | 42,876 | 16.9% |

**64.1% 的被漏边 identity < 40%,被 LINC 的 identity 阈值直接挡掉;其中 58.9%
(149,683 条)是"identity < 40% 但 bitscore ≥ 50"的远缘同源。**

### 4.2 决定性证据:RBB 图 + MCL = OrthoFinder 完全一致

读 OrthoFinder 源码 `__main__.py:474-479` 的构图逻辑:

```python
m1     = LoadMatrix("connect", iSpec, jSpec)       # i→j 的最佳命中连接
m2tr   = transpose(LoadMatrix("connect", jSpec, iSpec))
connect2 = m1 + m2tr                                # 2-way 连接 = 互惠最佳命中(RBB)
B_connect = MatricesAnd_s(connect2, B)              # RBB ∧ 归一化 bitscore
```

OrthoFinder 的 718,441 条边是 **RBB(互惠最佳命中)稀疏化**的结果,把 514 万条原始 hit 压缩到
71.8 万条(14%)。

用 `offline_rbb.py` 验证:RBB 图 + MCL 与 `Orthogroups.txt` 对比:

| 对比 | family F1 | precision | recall |
|---|---|---|---|
| **RBB 图 + MCL** | **1.0000** | **1.0000** | **1.0000** |
| LINC identity+coverage SSN + MCL | 0.8880 | 0.9425 | 0.8391 |

> **把 LINC 的构图判据换成 RBB,结果与 OrthoFinder 完全一致(F1=1.0)。**
> **0.888 → 1.0 的整个 gap,100% 由"边判据"这一件事解释。**

---

## 5. 机制解释

远缘同源(跨物种)的特征是:**identity 低(常 30~40%),但比对显著、且互为对方的最佳命中。**

- **RBB 判据**不 care identity 高低,只问"A 的最佳命中是不是 B,且 B 的最佳命中是不是 A"——
  所以能抓住低 identity 的远缘同源;
- **identity+coverage 阈值**是"裸百分比"这一维,看不见"最佳命中"和"显著性"这些信息,于是
  把低 identity 的远缘同源一刀切掉 → 图变稀疏 → MCL 过度切分 → 家族数膨胀。

这解释了同物种/跨物种的差异:同物种时 identity 高、两个判据基本一致(gap 小);跨物种时
identity 低、两个判据分歧大(gap 大)。

---

## 6. 对 LINC 的启示:修复有一个"增量张力"

RBB 是**全局判据**,与 LINC 的增量/append 范式存在真实冲突:

- RBB 要求"知道每个基因在所有基因组中的最佳命中",而 `--append` 加入新基因组时,
  **历史基因的最佳命中可能被新基因取代**,历史 RBB 边可能失效,需要回溯重判;
- 而 identity+coverage 是 per-hit 的**局部**判定,天然可以增量做。

因此 LINC 引入 RBB 不是"加一个参数"那么轻,而是需要在以下之间做架构级取舍:

1. **完全引入 RBB**:最忠实 OrthoFinder,但 append 可能需要回溯重判、牺牲增量纯粹性;
2. **近似 RBB**:如"仅保留每个基因的 top-k 命中"、"对历史基因仅在 append 时增量重判受影响部分"等折中;
3. **保持 identity+coverage 但接受与 OrthoFinder 的固有 gap**(0.888 已是 identity 判据的天花板)。

---

## 7. RBB 增量实现验证(负结果)

### 7.1 验证思路

第 4 节证明"全量 RBB + MCL = F1 1.0"用的是 **OrthoFinder 自己算好的 RBB 图**(基于全对全)。
为了回答"LINC 能否在自己的增量搜索产物上实现 RBB",加了 `--dump-raw-hits` 调试出口,把
join/full 阶段 cutoff 过滤**之前**的全部原始跨基因组 hit 落盘,再离线算 RBB。

### 7.2 实验(diamond, cutoff-i 0.25 / cutoff-c 0.4)

| 指标 | 值 |
|---|---|
| 落盘原始 hit 数 | 1,655,499 |
| 从中算出的 RBB 边 | **23,992** |
| OrthoFinder 全量 RBB 边 | 718,441 |
| **LINC RBB / 全量 RBB** | **3.3%** |

| 模式 | F1 | precision | recall | n_fam |
|---|---|---|---|---|
| SSN 基线(identity+coverage) | 0.7519 | 0.6985 | 0.8142 | 21,023 |
| SSN ∪ RBB | 0.7585 | 0.7026 | 0.8240 | 20,516 |

**RBB 只提升 0.66pp(F1 0.7519 → 0.7585),几乎无效。**

### 7.3 为什么:增量搜索拿不出"完整 hit 集"

"单向完整 hit 集能算 RBB"的原理没错,但前提是**完整**。LINC 的增量搜索:

- join 只搜 **代表 → B**;
- full 只搜 **"接触 B 的簇"的成员 → B**;
- **"不接触 B 的簇"的成员根本没被拿去和 B 比对。**

于是"B 基因的最佳命中"在一个残缺的 query 集合上计算,大量本该是"互惠最佳命中"的边,因其中
一方基因从未被搜到而判不出来——这就是 3.3% 的来源。

### 7.4 结论:RBB 增量不可行

**"在 SSNBuilderStage 补边阶段实现 RBB"这条路实测走不通**:不是实现问题,而是数据源不足——
LINC 增量搜索产出的 hit 集太残缺,RBB 只能算到全量的 3.3%,F1 收益可忽略。

---

## 8. 结论与建议

**结论(修正版):** LINC 与 OrthoFinder 的 gap 根因是**两层叠加**:

1. **判据差异**:OrthoFinder 用 RBB(互惠最佳命中)稀疏化构图,LINC 用 identity+coverage 阈值;
2. **搜索完备性差异**:RBB 需要完整 hit 集,而 LINC 的增量搜索只覆盖部分基因对,拿不出完整 RBB。

> **RBB 判据本身没错(全量 RBB 能 100% 复现 OrthoFinder),但它在 LINC 的增量 hit 集上算不出完整结果。
> 因此 gap 靠"换判据"救不回来——根源在"增量 vs 全对全"的搜索完备性鸿沟。**
> 要真正追平 OrthoFinder,只有回到全对全(等于放弃 LINC 的增量卖点),或接受该 gap 作为
> "增量范式的固有代价"。

**建议:**

1. 将"判据 + 搜索完备性"双层根因固化为项目知识,避免后续重复排查;
2. **不建议**为追平 OrthoFinder 而引入 RBB——增量搜索提供不了完整 hit 集,收益可忽略;
3. 若未来要做跨物种场景,优先评估"增量搜索覆盖率的代价 vs 全对全"这一根本权衡,
   而非在判据层面继续迭代(权重/阈值/判据均已被数据排除)。

---

## 附:本次调查用到的离线脚本(均不动生产代码)

| 脚本 | 用途 |
|---|---|
| `benchmark/offline_normalize.py` | 归一化 bitscore 消融(排除权重) |
| `benchmark/edge_set_compare.py` | 边集对比 + identity/evalue 分布 |
| `benchmark/offline_significance.py` | 显著性逃逸口探索(揭示 RBB 稀疏化) |
| `benchmark/offline_rbb.py` | 全量 RBB 图 + MCL 复现 OrthoFinder(F1=1.0) |
| `benchmark/offline_rbb_member.py` | member 级 RBB 验证(3.3%, 负结果) |
| `--dump-raw-hits`(生产代码调试出口) | 落盘 join/full 阶段 cutoff 前原始 hit, 供离线 RBB 验证 |
