# LINC

LINC 是一个以**增量方式比较多个基因组**、构建**序列相似性网络(SSN)**,
并通过 **MCL 算法**实现蛋白质家族快速聚类的框架。

## 项目状态

当前处于架构设计阶段,尚未开始正式编码。详细的技术方案与开发进度见 `docs/` 目录。

## 文档索引

- [docs/LINC_技术框架设计_v2.md](./docs/LINC_技术框架设计_v2.md) — 现代化技术框架设计(权威工程依据)
- [docs/开发进度.md](./docs/开发进度.md) — 开发进度记录
- [docs/归档/](./docs/归档) — 历史文档归档(算法设计文档 v1、STEP1~STEP6 系列、MMseqs2_Aligner 等)

## 环境

- Python 3.11(通过 conda 环境 `linc` 管理)
- 依赖: `numpy`、`networkx`
- 第三方比对/聚类工具: `diamond`、`mmseqs2`、`mcl`(通过 conda 安装)

## 运行

```bash
conda activate linc
python main.py
```

更多内容后续补充。
