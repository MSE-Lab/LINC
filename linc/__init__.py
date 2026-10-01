"""LINC: 增量式基因组比较与序列相似性网络(SSN)构建、MCL 蛋白质家族聚类框架。

包结构(六边形架构,详见 docs/LINC_技术框架设计_v2.md):
    core/       领域核心: 实体、事件、策略接口 (不依赖任何外部工具)
    ports/      抽象接口 (Protocol): Aligner / Storage / Clustering
    adapters/   Port 的具体实现: mmseqs / diamond / mock / sqlite / mcl
    engine/     增量引擎: 五个 Stage 串联核心算法
    config/     配置 Schema 与加载
"""

__version__ = "0.1.0"
