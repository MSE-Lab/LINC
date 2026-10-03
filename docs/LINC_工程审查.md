# LINC 工程审查报告

> 审查日期:2026-10-02
> 审查视角:安全、错误处理、崩溃一致性、状态兼容性、并发、资源管理
> 说明:算法层漏洞见 `docs/LINC_漏洞与优化分析_第二轮.md`,本文聚焦**工程层**问题

---

## 一、🔴 状态一致性(最值得先修)

### 1. `--resume`/`--append` 不校验算法参数是否与历史运行一致

**位置**:`linc/adapters/storage/file_repo.py` + `linc/cli.py` + `linc/pipeline.py`

`FileRepository` 只持久化了 `input_signature`(文件名集合)和 5 个 Stage 的产出,**没有
持久化任何算法参数**(cutoff_i/ii/c、cov_mode、max_per_cluster、edge_completion_mode、
inflation、edge_weight_mode、edge_weight_alpha 等)。

**后果**:用户用**不同的阈值**去 `--resume`/`--append` 时,系统完全不报错——加载用旧阈值算好
的 `intra_results`(state 里),再用新阈值去跑后续 Stage,产出**新旧阈值混搭的不一致结果**。例如:

- append 时 `historical_intra`(旧 cutoff)与 `new_intra_results`(新 cutoff)直接合并,无任何校验;
- resume 时加载旧 `join_result`,再用新 inflation 重跑 MCL。

**建议**:把 `PipelineConfig`(或其哈希)一并持久化到 state,resume/append 前比对,不一致则
报错并要求显式覆盖。

---

### 2. pickle 状态无 schema 版本控制

**位置**:`linc/adapters/storage/file_repo.py`

`FileRepository` 直接 `pickle.dump` dataclass 对象,无任何版本号。代码升级(如给
`IncrementalResult` 加 `cluster_counter`、`SSNEdge` 加 `bitscore`/`source_length`)后,旧
state 文件能"成功"反序列化但**缺字段**,随后在 `historical_join.cluster_counter` 等处抛
`AttributeError`——不是干净的"版本不兼容"报错。`StateCorruptedError` 只能兜住反序列化
失败,兜不住"字段缺失"这类静默不兼容。

**建议**:state 里加 `SCHEMA_VERSION` 标记,load 时校验,不匹配给明确提示。

---

### 3. 6 个 state 文件之间无"同一次运行"的绑定

**位置**:`linc/adapters/storage/file_repo.py`

`genomes / intra_results / join_result / complete_ssn / families` 是 6 个独立 pickle,
没有任何 run-id 或校验和把它们绑在一起。崩溃/并发/手工删改都可能导致 resume 时加载到
"跨多次运行拼凑"的状态,而 `has_*()` + `load_*()` 的设计无法察觉。

**建议**:存一个 run-id + 各文件校验和(或至少记录参数哈希,见第 1 条)。

---

## 二、🟠 错误处理

### 4. 缺失外部二进制 → 未捕获的 traceback

**位置**:`linc/cli.py` `_build_aligner()` / `_build_clustering()`

`MMseqsAligner()`/`DiamondAligner()`/`MCLClustering()` 构造函数在找不到二进制时
`raise RuntimeError`(信息本身清晰),但:

- `_build_aligner()`(约 L376)和 `_build_clustering()`(约 L434)都在 `try/except`
  (约 L450 才开始的 try)之外。

所以用户没装 mmseqs/diamond/mcl 时,看到的是 Python traceback,而非干净的 `typer.Exit`。

**建议**:把这些构造纳入 try,或统一转成 `typer.secho(..., err=True) + Exit(1)`。

---

### 5. MCL 输出完整性未校验 —— 基因可能静默丢失

**位置**:`linc/engine/family_refine_stage.py` `run()`

对每个连通分量调用 MCL,然后 `for member_ids in sub_clusters:` 逐个生成 family。若某分量
MCL 返回**空/不完整的分区**(工具异常、退化输入),代码不检查"所有分量节点是否都被覆盖",
这部分基因会**静默从结果中消失**。

**建议**:加防御性断言——所有分量节点必须出现在 `sub_clusters` 的并集中。

---

### 6. 空输入 / 空 append 目录静默产出空结果

**位置**:`linc/cli.py` `run()`

`discover_fasta_files()` 没找到 FASTA 时只打黄色 warning,**继续跑**,最终产出空
Orthogroups.txt;`--append` 指向空目录时同理("成功"但实际未追加任何东西)。

**建议**:无输入时直接报错退出。

---

## 三、🟡 并发 / 资源

### 7. 无文件锁,并发跑同一 `-o` 会互相破坏

**位置**:`linc/adapters/aligners/*.py` `search()` + `linc/adapters/storage/file_repo.py`

`search()` 里 `shutil.rmtree(search_<tag>)` + 重建、`FileRepository._save` 的 `.tmp` +
`os.replace`,都没有加锁。两个 `linc run -o output` 同时跑同一输出目录会竞争、互相覆盖
中间文件/state。

**建议**:在 `-o` 目录放 lock 文件,或文档明确"禁止并发同 -o"。

---

### 8. `subprocess.run` 没有 `timeout`

**位置**:`linc/adapters/aligners/mmseqs_adapter.py` / `diamond_adapter.py` /
`linc/adapters/clustering/mcl_adapter.py`

所有 mmseqs/diamond/mcl 调用都 `check=True, capture_output=True` 但无 `timeout`。外部工具
一旦卡死(极端输入、磁盘满等),进程会永远挂起。

**建议**:加可配置超时(可留大默认值)。

---

### 9. 每个连通分量一次 `mcl` 子进程

**位置**:`linc/engine/family_refine_stage.py` + `linc/adapters/clustering/mcl_adapter.py`

`FamilyRefineStage` 对每个非孤立连通分量单独调一次 `mcl`(每次一次子进程启动)。数据有几千个
小分量时,就是几千次子进程 spawn,启动开销累积显著(与第二轮的 O(G²) 是不同维度)。

**建议**:考虑"合并多个分量进一次 mcl 调用"或改用进程内库。

---

## 四、🟢 小问题

- **summary echo 用 CLI 参数而非实际加载值**:resume 后 `cli.py` 打印
  `edge_completion_mode={edge_completion_mode}` 用的是本次命令行参数,而实际 SSN 可能来自
  旧 state 的另一种模式,显示会误导。
- **`--inflation` 只校验 `> 0`**:MCL 的 inflation 通常应 ≥1(常用 1.2~5),当前允许 0.5 这类
  退化解,建议收紧为 `>= 1`。
- **`--append` 不记录追加历史**:无"追加了哪些目录/哪些 genome_id"的审计,与项目"可追溯"
  目标相关(第二轮已提)。

---

## 五、正面结论(安全面)

- `subprocess` 全部用 **list 参数 + `shell=False`**,无命令注入;
- `genome_id` 经 `os.path.basename` + `splitext`,无法构造路径穿越;
- 无硬编码密钥、无网络、无反序列化不可信输入(唯一例外是第一类 pickle,但那是"本地 state
  被篡改"的低风险场景)。

---

## 六、优先级建议

| 优先级 | 项 | 类型 |
|---|---|---|
| P0 | resume/append 参数一致性校验(第 1 条) | 状态一致性 |
| P1 | pickle schema 版本控制(第 2 条) | 兼容性 |
| P1 | state 文件 run-id 绑定(第 3 条) | 崩溃一致性 |
| P1 | 缺失二进制 → 干净报错(第 4 条) | 错误处理 |
| P1 | MCL 输出完整性校验(第 5 条) | 正确性兜底 |
| P2 | 空输入/空 append 报错(第 6 条) | 错误处理 |
| P2 | 文件锁(第 7 条) | 并发 |
| P2 | subprocess 超时(第 8 条) | 资源 |
| P2 | MCL 子进程批量(第 9 条) | 性能 |

> 第 1 条(resume/append 参数一致性)应最高优先——它会在**不报错的情况下产出错误结果**,
> 是工程上最危险的类型。第 2/3 条与第 1 条一起构成"状态安全"三件套。
