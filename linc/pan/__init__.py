"""linc pan 的结果描述/可视化增强实现包。

设计依据: docs/LINC_结果描述增强开发计划.md。本包只读 `.linc_work/state`,
不重跑任何 Stage,不修改任何 state 文件。入口为 `linc pan` 子命令(见
`linc/cli.py`),本包被其调用。
"""
