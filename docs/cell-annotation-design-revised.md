# 细胞注释模块修订设计

## 设计取舍

本版本以类器官为主要对象，保留“自动预测 → Marker 证据 → 冲突识别 → 人工复核”的主线。
不把 CellTypist、SingleR、Azimuth 或 scANVI 作为默认最终注释器：当前数据体系缺少稳定的同协议参考时，强行映射会把新状态压到最相近的已知标签。现在仅将已下载的人类 CellTypist 模型作为可选的独立参考证据，不覆盖 Marker 最终标签。

## 输出层级

细胞类型和细胞状态分开保存：

- `cell_lineage` / `cell_type_l1`：主要谱系
- `cell_type_l2`：主要细胞类型
- `cell_type_l3`：亚型或组织内细分类型
- `cell_state`：单一主状态，如 Cycling、Stress response
- `cell_state_flags`：可并存的多个状态
- `developmental_state`：类器官表达成熟度，如 `progenitor_like`、`transitional`、`mature_like`

未知或证据不足时保留 `Unknown`，不强制分类。

## 类器官成熟度

不要求用户主观选择“早期/成熟”。系统优先读取 `culture_day`、`day`、`timepoint`、`hour` 等 `adata.obs` 元数据，并独立计算：

- 前体模块分数
- 成熟模块分数
- 增殖模块分数
- `organoid_maturity_index`

成熟度结果只作为证据和复核信息，不参与细胞类型打分。没有时间元数据时，系统明确提示仅使用表达模块。

## 证据与冲突

- 正向 Marker：当前组织 panel 和数据驱动 Marker
- 负向 Marker：降低不匹配候选分数，不单独删除标签
- Doublet：两个独立 Marker 模块同时获得高分时标记 `suspect_doublet`
- 环境 RNA：高普遍性异源 Marker 的启发式提示，不能替代 empty-droplet 去污染
- Cluster 一致性、Top1/Top2、Marker 覆盖度和人工 review 共同写入 h5ad 与复核表

## 可追溯性

每次运行记录：

- `annotation_version`
- `annotation_comment`
- `celltypist_label`、`celltypist_confidence`、`celltypist_status`、`celltypist_comparison`（启用参考时）
- 使用的 marker_set、类器官类型和 cluster 列
- Marker 覆盖度、状态证据、Doublet/环境 RNA 统计

人工映射仍通过 `manual` 方法执行；新的输出不会覆盖旧 h5ad，后续版本可以使用不同的 `annotation_version` 区分。

## CellTypist 参考模式

界面中的“启用 CellTypist 参考交叉验证”只读取 `data/references/celltypist/` 下的本地 `.pkl`，不会在分析期间联网下载。`prob match` 是默认模式，可保留 `Unknown`、低置信度和多标签结果；CellTypist 与 Marker 的冲突写入 `celltypist_comparison`，并进入逐簇复核证据。缺少 `celltypist` 依赖或模型文件时，任务回退到 Marker 注释并记录 warning。
