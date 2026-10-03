# LINC 归一化 bitscore(OrthoFinder 式幂律拟合)设计方案

> 状态:设计稿(待评审)
> 日期:2026-10-03
> 关联:`benchmark/RESULTS.md`(OrthoFinder 对照)、`docs/LINC_技术框架设计_v2.md`(长度归一化待办)

---

## 1. 背景与目标

LINC 当前的边权重用 `identity × coverage`(或原始 `bitscore`),与 OrthoFinder 在
`Orthogroups.txt`(MCL 聚类结果)这一同口径下的 family F1 为 **0.959**,存在约 4% 的 gap。

经源码核实(`/opt/miniconda3/envs/orthofinder2/bin/scripts_of/__main__.py`),OrthoFinder 的
核心优势之一是其**按物种对拟合的幂律归一化 bitscore**。本方案目标:在 LINC 中实现一个
**增量兼容的、per-pair 幂律拟合的归一化 bitscore**,以期收窄与 OrthoFinder 的 gap,同时
不破坏 LINC 的增量/append 卖点。

---

## 2. OrthoFinder 归一化算法(源码还原)

对每一对物种 (i, j),设原始 bit score 矩阵 B,B[q,h] = query q 对 hit h 的原始 bit score:

1. 取每个非零 hit 的 query 长度 `Lq`、hit 长度 `Lh`、原始分数 `S`;
2. 定义**长度乘积** `Lf = Lq × Lh`;
3. 按 `Lf` 分 bin,每个 bin 内取 **top 95%** 的高分 hit(去噪声、拟合"分数-长度"的上包络);
4. 对 top hit 做 log-log 线性拟合(源码用 scipy `curve_fit`):

```
log10(S) = a · log10(Lf) + b        ⟺   S ≈ 10^b · (Lq·Lh)^a
```

5. 归一化(除以拟合出的长度幂律):

```
S_norm = 10^(-b) · Lq^(-a) · S · Lh^(-a)  =  S / (10^b · (Lq·Lh)^a)
```

**参数含义**:`a`(指数)度量长度偏差强度,全长比对下理论值 ≈ 0.5;`b`(截距)编码
"每单位长度乘积的典型 bit score",即比对质量/进化距离基线,按物种对拟合。

---

## 3. 关键结论:可增量 + LINC 数据充足

### 3.1 per-pair 的 (a, b) 天然增量兼容

`a_ij`、`b_ij` **只依赖基因组 i 与 j 之间的 hits**,与其它基因组无关:

- `--append` 新基因组 k 时,只需对**新产生的 pair (k, i)**(对所有历史 i)拟合 (a, b);
- 历史 pair (i, j) 的 (a, b) **保持不变**(i、j 之间的 hits 没变)。

因此"按物种对拟合"并非全量批处理,而是"每新增一个基因组,补拟合一批新 pair",与 LINC
的增量范式相容。

### 3.2 LINC 的增量 hits 足够拟合

OrthoFinder 拟合用的是"top 95% 高分 hit"(真同源)。LINC 的 SSN 实测:

- 边级 precision **0.9999**、recall **0.89**;
- 漏掉的 11% 绝大多数是**低分/冗余边**(家族 F1 仍 0.995),本就不会进入"top"拟合集合。

结论:**用 LINC 的 SSN 边拟合 (a, b) ≈ 用全量 hits 拟合 (a, b)**,误差来自被漏掉的低分边,
对幂律拟合影响可忽略。

---

## 4. 总体设计

新增一种边权重模式 `normalized_bitscore`,权重 = `S / (10^b · (Lq·Lh)^a)`,其中 (a, b)
为**按基因组对拟合**的参数。流程上在 `SSNBuilderStage`(已产出带 `bitscore`/长度 的 SSN)
之后、`FamilyRefineStage`(MCL)之前插入一步**"归一化参数拟合"**。

关键点:

- 拟合用 **numpy 一次线性回归**(`np.polyfit(..., 1)`),**不引入 scipy**;
- 拟合输入为 SSN 边(已带 `bitscore`、`source_length`、`target_length`),无需新增比对;
- 归一化后的权重只在 MCL 构图时惰性计算,不改动 SSN 边的原始存储。

---

## 5. 详细设计

### 5.1 拟合函数(纯 numpy)

```python
import numpy as np
from typing import Tuple

# 单个 genome pair 的拟合: 输入若干 (bitscore, Lq, Lh)
def fit_pair(bitscores, Lq_list, Lh_list) -> Tuple[float, float]:
    """返回 (a, b), 满足 log10(S) = a*log10(Lq*Lh) + b。"""
    S = np.asarray(bitscores, dtype=float)
    Lf = np.asarray(Lq_list, dtype=float) * np.asarray(Lh_list, dtype=float)
    mask = (S > 0) & (Lf > 0)          # 过滤非正数(理论上不应出现)
    S, Lf = S[mask], Lf[mask]
    if len(S) < _MIN_HITS_PER_PAIR:    # 样本太少, 回退
        return _FALLBACK_A, _FALLBACK_B
    a, b = np.polyfit(np.log10(Lf), np.log10(S), 1)   # 线性最小二乘
    if not np.isfinite(a) or not np.isfinite(b):
        return _FALLBACK_A, _FALLBACK_B
    return float(a), float(b)
```

> `np.polyfit(log10(Lf), log10(S), 1)` 与 OrthoFinder 的
> `curve_fit(scnorm.loglinear, Lf, np.log10(S))`(其中 `loglinear = a*log10(x)+b`)
> 数学等价,只是用 numpy 代替 scipy。

### 5.2 归一化权重公式

```python
def normalized_bitscore_weight(bitscore, Lq, Lh, a, b) -> float:
    return bitscore / (10.0 ** b * (Lq * Lh) ** a)
```

在 `SSNEdge.weight()` 中新增 `"normalized_bitscore"` 模式;由于该模式需要 (a, b) 参数,
不能只用 `mode: str` 表达,需扩展 `weight()` 签名(见 5.3)。

### 5.3 在 pipeline 中的位置与依赖

```
SSNBuilderStage(产出 CompleteSSN: 边带 bitscore/source_length/target_length)
        │
        ▼
【新增】ScoreNormalizer.fit(edges, gene_to_genome)  →  Dict[(genome_a, genome_b), (a, b)]
        │
        ▼
FamilyRefineStage(用归一化权重构图 → MCL)
```

- `ScoreNormalizer.fit` 是**纯函数模块**(建议 `linc/core/score_normalizer.py` 或
  `linc/engine/score_normalizer.py`),只依赖 numpy + 标准库;
- `gene_to_genome`:把每条 SSN 边的 source/target 映射到 genome_id。`SSNBuilderStage`
  已持有该映射(`genes_by_id` 的 `Gene.genome_id`),应作为 fit 的输入传入;
- `Pipeline` 在 `_run_ssn` 之后、`_run_family` 之前调用 fit,把结果传入 `FamilyRefineStage`。

### 5.4 数据流与持久化

- 拟合结果 `{(genome_a, genome_b): (a, b)}` 需**随 state 持久化**(供 `--append`/
  `--resume` 复用),建议 `RepositoryPort` 新增 `save_score_params`/`load_score_params`,
  与既有 state 文件风格一致;
- append 时只拟合"新 genome × 历史 genome"的新 pair,历史 pair 参数直接复用(见 3.1)。

### 5.5 边界情况与回退

| 情况 | 处理 |
|---|---|
| 某 pair 的 SSN 边数 < `_MIN_HITS_PER_PAIR`(建议 20) | 回退 `a=0.5, b=0`(即 `S/sqrt(Lq·Lh)`) |
| 拟合结果非有限(nan/inf) | 回退 `a=0.5, b=0` |
| 某 pair 长度多样性极低(斜率不可靠) | 可加"Lf 跨度"检查,过窄则回退(可选) |
| 基因组对无任何 SSN 边 | 不产生该 pair 参数;相关边本就不存在 |

### 5.6 append / resume 策略

- **append 新基因组 k**:仅对新 pair (k, i) 拟合;历史 pair 参数从 state 加载,不重拟合;
- **resume**:拟合结果作为 state 的一部分,若 `score_params` 已存在则复用,否则重拟合;
- **force(参数不一致)**:`_check_config_consistency` 会把新增归一化参数纳入比较;纯"重拟合"
  不改变算法参数语义,需评估是否加入 `non_semantic_fields`。

---

## 6. CLI 设计

`linc run` 的 `--edge-weight-mode` 新增一个取值:

```
--edge-weight-mode normalized_bitscore
```

可选暴露拟合细节:

```
--norm-min-hits 20          # 每 pair 至少多少条边才拟合, 否则回退 a=0.5,b=0
--norm-top-percentile 95    # 拟合前取每个长度 bin 的 top 百分位(与 OrthoFinder 对齐)
```

> 默认 `--edge-weight-mode` 仍为 `identity_x_coverage`,归一化 bitscore 作为 opt-in。

---

## 7. 验证方案

复用 `benchmark/compare_orthofinder.py`(无需改动):

1. `--edge-weight-mode normalized_bitscore` 重跑 `linc run`;
2. 对比 OrthoFinder `Orthogroups.txt`,看 family F1 相对基线(0.959)是否提升;
3. **消融**:分别验证 (a) 只做 a(固定 b=0)vs (b) a+b 全拟合 vs (c) 是否启用 top-95% 分 bin,
   定位增益来源;
4. **增量一致性**:验证 `append` 后新 pair 拟合与"全量重跑"的 family 结果接近;
5. 记录 (a, b) 分布(是否集中在 ~0.5、跨 pair 是否稳定),作为可解释性证据。

---

## 8. 开放问题与决策点

1. **拟合输入范围**:MVP 用"cutoff 过滤后的 SSN 边";是否扩展为"full 阶段存下的 raw hits
   (过 evalue 未过 cutoff)"一起拟合?——建议 MVP 先用 SSN 边,收益不够再扩。
2. **top-95% 分 bin 是否纳入 MVP**:它是 OrthoFinder 去噪关键,但增加复杂度;建议 MVP 先
   全量拟合、把 top-95% 作为可选项。
3. **a、b 都拟合还是只拟合 a**:numpy 一次给出两者,成本相同;建议都拟合,样本不足时回退。
4. **是否影响等价性结论**:归一化只改边权重、不改边集,不影响"增量 ≈ 全量"的等价性验证;
   但会改变 MCL 结果,需同步重跑等价性指标。

---

## 9. 实现清单(拆解)

| 步骤 | 内容 | 依赖 |
|---|---|---|
| 1 | `linc/core/score_normalizer.py`:`fit_pair` / `fit_all_pairs` / `normalized_bitscore_weight` | numpy |
| 2 | `SSNEdge.weight()` 支持 `normalized_bitscore`(传入 pair 参数) | 步骤 1 |
| 3 | `Pipeline` 在 `_run_ssn` 后调用 fit,结果传给 `FamilyRefineStage` | 步骤 1 |
| 4 | `RepositoryPort`/`FileRepository` 增 `save_score_params`/`load_score_params` | — |
| 5 | append/resume 接入拟合参数复用(3.1/5.6) | 步骤 4 |
| 6 | CLI:`--edge-weight-mode normalized_bitscore` + 可选参数 | 步骤 2 |
| 7 | 单测:拟合数值正确性、回退逻辑、append 复用 | 步骤 1-6 |
| 8 | 基准:`compare_orthofinder.py` 验证 F1 增益 + 消融 | 步骤 1-6 |

---

## 10. 推进建议(评审结论, 2026-10-03 追加)

### 10.1 核实结论

对照现有代码(`SSNEdge`/`CompleteSSN`/`FamilyRefineStage`/`Pipeline`/`FileRepository`)逐条核实后,
本方案**技术上可行**, 与现有代码结构吻合度较高, 未发现"方案在当前代码库下走不通"的硬伤:

- `SSNEdge` 已预留 `bitscore`/`source_length`/`target_length` 字段, 注释明确写着"为未来 MCL
  权重归一化(如 OrthoFinder 式长度校正)预留", 上游阶段无需改动;
- `SSNEdge.weight(mode, alpha)` + `VALID_EDGE_WEIGHT_MODES` 的既有扩展机制可直接复用, 新增
  `normalized_bitscore` 只是多加一个分支, 但需要扩展 `weight()` 签名以传入 per-pair 的 `(a, b)`,
  这是唯一一处需要修改既有接口签名的地方, 实现时要重点关注调用方
  (`FamilyRefineStage` 的 `edge.weight(self._edge_weight_mode, self._edge_weight_alpha)`)
  怎么拿到并传入 `(a, b)`;
- `Pipeline._run_ssn`/`_run_family` 在 `run()`/`resume()`/`append()` 中紧邻调用, 插入点位置准确,
  且 `loaded_genomes` 在该作用域内可直接取到, 构建 `gene_to_genome` 映射无需改动上游;
- `Genome.genome_id` 现为内部短编号(`G00001`...), `append` 场景下历史基因组 `genome_id` 不变,
  "按 `(genome_a, genome_b)` 存取 `(a, b)`, 历史 pair 不受影响"的增量设计站得住脚;
- `np.polyfit(log10(Lf), log10(S), 1)` 与 OrthoFinder `curve_fit` 拟合同一个线性模型数学等价,
  用 numpy 代替 scipy 没有问题;
- `_check_config_consistency` 已有 `non_semantic_fields = {"threads"}` 机制, "评估拟合参数是否归入
  `non_semantic_fields`"可直接复用该机制。

### 10.2 编码前需要明确/补强的点

1. **第8节两个开放问题需要先用数据验证, 不能假设结论**: top-95% 分 bin 去噪是 OrthoFinder 认为
   的"去噪关键", 若 MVP 省略, 有可能导致 F1 提升不如预期; 拟合用 SSN 边还是 full 阶段原始 hits
   同理需要实测对比, 不能直接假设"SSN 边已经足够"。
2. **`_MIN_HITS_PER_PAIR` 回退阈值需结合真实数据校验**: 若较多 genome pair 因边数不足触发回退
   (`a=0.5, b=0`), 方案会退化为"长度几何平均归一化", 实际收益打折; 验证阶段应统计"回退触发比例"
   作为是否真正生效的指标, 不能只看最终 F1。
3. **bitscore 方向性约定需在实现时保持一致**: `source_length`/`target_length` 的 query/subject
   身份按 `genome_id` 字典序单向固定(人为约定, 非双向对称值); 拟合 `fit_pair` 的 `Lq_list`/
   `Lh_list` 取值必须与归一化计算权重时的 `Lq`/`Lh` 角色严格对应, 建议单测专门覆盖这一致性。
4. **新增 CLI 参数(`--norm-min-hits`/`--norm-top-percentile`)需补充参数校验**, 并明确是否纳入
   `_check_config_consistency` 的历史一致性比较。

### 10.3 建议的推进顺序

**第一步: 先离线脚本验证收益, 再决定是否写入 pipeline。**
不改动 `linc/` 生产代码, 先写一个独立脚本(类似 `benchmark/compare_orthofinder.py` 旁边新增),
直接读取已有的 `CompleteSSN`/中间产出, 离线做 `(a, b)` 拟合 + 重算权重 + 重跑 MCL, 对比 family F1。
同时把第8节两个开放问题的消融一并做完(分 bin vs 不分 bin; 只拟合 a vs a+b 全拟合, 共4种组合),
确认哪种组合收益最大、复杂度是否值得, 并记录 `_MIN_HITS_PER_PAIR` 回退触发比例、`(a, b)` 分布。

**第二步: 根据验证数据决定最终方案形态。**
若"不分 bin"已能拿到大部分收益, 维持 MVP 从简思路; 若收益明显打折, 分 bin 逻辑需直接写入 MVP、
不能留作可选项, 相应需重新评估工作量。

**第三步: 验证达到预期后才进入正式编码(按第9节实现清单)。**
只有确认"归一化确实能收窄与 OrthoFinder 的 gap", 才值得改动生产代码(`SSNEdge.weight()` 签名、
新增 `score_normalizer.py`、接入 `Pipeline`/`FileRepository`/CLI, 可能再次递增 `SCHEMA_VERSION`)。
涉及接口签名变化、持久化字段变化等"大方案层面"决策, 仍需先说明方案并等待确认后再动手实现。
