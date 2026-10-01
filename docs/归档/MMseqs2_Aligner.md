# ============================================================
# MMseqs2 比对引擎 - LINC 集成说明
# ============================================================

## 概述

LINC 采用 **MMseqs2** (Multi-Mapping Sequence Search) 作为底层序列比对引擎。

MMseqs2 是一个超快速、高灵敏度的蛋白质序列搜索与聚类套件。

---

## MMseqs2 性能特性

| 特性 | MMseqs2 | BLAST | DIAMOND |
|------|--------|-------|---------|
| 速度 | 基准 | 基准 | 快 500-10000x |
| 灵敏度 | 接近 BLAST | 100% | ~95% |
| 内存需求 | 中等 | 低 | 低 |
| 适用场景 | 大规模数据 | 小规模数据 | 超大规模数据 |

---

## MMseqs2 在 LINC 中的三个关键任务

### 1. 基因组内全对全比对 (Step 2)

```python
# Step 2 使用 MMseqs2 进行基因组内自聚类
aligner.run_all_vs_all(fasta_file, identity_threshold=0.5)
```

**参数配置:**
```yaml
aligner:
  backend: "mmseqs"
  mmseqs:
    sensitivity: 4.5  # 平衡模式
    evalue: 1e-5
    max_accept: 500
```

### 2. 新局部代表与全局代表库比对 (Step 3)

```python
# Step 3 使用 MMseqs2 进行增量比对
aligner.run_query(
    query_file=new_local_reps.faa,
    subject_file=global_library.faa
)
```

### 3. 补边阶段批量比对 (Step 5)

```python
# Step 5 使用 MMseqs2 进行批量补边
aligner.run_query(
    query_file=batch_genes.faa,
    subject_file=cluster_members.faa
)
```

---

## MMseqs2 工作流程

```
┌─────────────────────────────────────────────────────────┐
│           MMseqs2 三步比对流程                           │
└─────────────────────────────────────────────────────────┘

输入: FASTA 文件
    │
    ▼
┌─────────────────────────────────────────────────────────┐
│  Step 1: createdb - 构建数据库                         │
│  ├─ 序列数据 (.db)                                    │
│  ├─ 索引 (.index)                                    │
│  └─ 元数据 (.dbtype)                                 │
└─────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────┐
│  Step 2: prefilter - 快速预过滤                       │
│  ├─ 使用段存储 (Segment Hash) 快速查找候选匹配        │
│  ├─ 极快的 k-mer/tuple 索引                          │
│  └─ 大幅减少需要精确比对的序列对数量                   │
└─────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────┐
│  Step 3: align - 精确动态规划比对                     │
│  ├─ Myers Bit-Vector DP 算法                         │
│  ├─ SIMD 并行化 (AVX2/AVX512)                       │
│  └─ 生成高质量的 Smith-Waterman 比对                  │
└─────────────────────────────────────────────────────────┘
    │
    ▼
输出: AlignmentResult 列表
```

---

## MMseqs2 敏感度参数

MMseqs2 的敏感度通过 `-s` 参数控制,范围 1.0-7.0:

| 敏感度 | 速度 | 适用场景 |
|--------|------|---------|
| 1.0-2.0 | 最快 | 快速筛选,近缘同源 |
| 3.0-4.0 | 快 | 中等灵敏度筛选 |
| **4.5** | 平衡 | **默认,推荐** |
| 5.0-6.0 | 较慢 | 高灵敏度,远缘同源 |
| 7.0 | 最慢 | 最高灵敏度 |

### LINC 中的推荐配置

```yaml
# 快速模式 (适合初步筛选)
aligner:
  backend: "mmseqs"
  mmseqs:
    sensitivity: 3.5
    evalue: 1e-3

# 平衡模式 (默认)
aligner:
  backend: "mmseqs"
  mmseqs:
    sensitivity: 4.5
    evalue: 1e-5

# 高灵敏度 (适合远缘同源检测)
aligner:
  backend: "mmseqs"
  mmseqs:
    sensitivity: 6.0
    evalue: 1e-5
```

---

## MMseqs2 vs DIAMOND

| 特性 | MMseqs2 | DIAMOND |
|------|---------|---------|
| 算法 | k-mer + DP | 压缩索引 + DP |
| 敏感度 | 可调 (1-7) | 可调 (fast-sensitive) |
| 多线程 | 原生支持 | 原生支持 |
| 增量搜索 | 支持 | 有限支持 |
| 聚类 | 内置 linclust | 需要外部工具 |
| 许可证 | BSD | GPL |

### LINC 选择 MMseqs2 的原因

1. **可调的灵敏度**: 可以根据需要调整灵敏度和速度的平衡
2. **内置聚类**: linclust 模式适合 LINC 的增量架构
3. **现代设计**: 更好的 API 和错误处理
4. **持续更新**: 活跃的开发社区

---

## 安装 MMseqs2

### 方法 1: conda (推荐)

```bash
conda install -c bioconda mmseqs2
```

### 方法 2: 编译安装

```bash
git clone https://github.com/soedinglab/MMseqs2.git
mkdir mmseqs2-build && cd mmseqs2-build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
make install
```

### 方法 3: 预编译二进制

```bash
# Linux
wget https://github.com/soedinglab/MMseqs2/releases/latest/MMseqs2-Linux.tar.gz
tar -xzf MMseqs2-Linux.tar.gz

# macOS
wget https://github.com/soedinglab/MMseqs2/releases/latest/MMseqs2-macOS.tar.gz
tar -xzf MMseqs2-macOS.tar.gz

# Windows (通过 WSL 或 Docker)
docker run -v $(pwd):/data -it soedinglab/mmseqs2
```

---

## 验证安装

```bash
mmseqs --version
```

输出应类似: `MMseqs2 Version: xxx`

---

## LINC 配置示例

### 完整 MMseqs2 配置

```yaml
aligner:
  backend: "mmseqs"

  mmseqs:
    # 敏感度 (1.0-7.0)
    sensitivity: 4.5

    # E-value 阈值
    evalue: 1e-5

    # 最大保留比对数
    max_accept: 500

    # 最小比对分数
    min_alignment_score: 20
```

### 并行配置

```yaml
parallel:
  # 使用所有可用核心
  n_threads: 0

  # 启用基因组间并行
  inter_genome_parallel: true

  # 启用补边并行
  finalize_parallel: true
```

---

## 输出格式

MMseqs2 输出转换为 LINC 的 `AlignmentResult`:

```python
@dataclass
class AlignmentResult:
    query_id: str           # 查询序列 ID
    subject_id: str         # 目标序列 ID
    identity: float         # 序列一致性 (0-1)
    coverage: float         # 覆盖度 (min query_cover, subject_cover)
    evalue: float          # E-value
    bitscore: float         # Bit-score
    query_cover: float      # 查询覆盖度
    subject_cover: float    # 目标覆盖度
    alignment_length: int   # 比对长度
```

---

## 故障排除

### 问题: MMseqs2 未找到

```
Error: MMseqs2 not found
```

**解决方案:**
1. 确保 MMseqs2 已安装并添加到 PATH
2. 使用 conda 安装: `conda install -c bioconda mmseqs2`
3. 验证: `which mmseqs` 或 `mmseqs --version`

### 问题: 内存不足

MMseqs2 需要足够的内存来处理大型数据库。

**解决方案:**
1. 减少 `max_accept` 参数
2. 使用较低的敏感度 (更快的预过滤)
3. 分批处理大型数据集

### 问题: 比对超时

**解决方案:**
1. 增加 `n_threads` 数量
2. 使用较低的敏感度
3. 增加超时设置

---

## 相关文档

- [MMseqs2 官方文档](https://mmseqs2.readthedocs.io/)
- [MMseqs2 GitHub](https://github.com/soedinglab/MMseqs2)
- [LINC Step 2 - 基因组内自聚类](./STEP2_基因组内自聚类.md)
- [LINC Step 3 - 增量并入全局代表库](./STEP3_增量并入全局代表库.md)
