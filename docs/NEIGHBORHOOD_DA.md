# 单细胞邻域差异丰度

`neighborhood_da` 面向连续而非完全离散的细胞状态，例如 IBD 类器官中的 stress、TA/cycling、absorptive 与 metabolic 上皮谱系。

## 当前实现与统计边界

模块在 `X_pca_harmony`、`X_scVI`、`X_scanorama`、`X_pca_combat` 或 `X_pca` 中（按优先顺序自动选择，也可显式指定）抽取可复现的 seed cell，并构建重叠 KNN 邻域。二维 UMAP 只用于显示，绝不用于定义统计邻域。每个邻域在一个样本中的 abundance 是：

\[
\frac{\text{该样本落在邻域内的细胞数}}{\text{该样本范围内的全部细胞数}}
\]

两条件以这些 sample-level proportions 做双侧 Mann–Whitney U 检验；多个邻域在每个指定条件比较内做 BH-FDR。低于 `min_cells_per_sample` 的样本只保留在表格中、不进入检验；任一条件的独立生物学样本少于 `min_samples_per_condition` 时直接停止。细胞数不是重复数。

这是适合当前内网小规模平台的 Milo-inspired 实现，保留“局部状态 + 样本重复 + 邻域 FDR”的关键语义；它不是 R Milo 的 negative-binomial 图模型，也不应在论文中称为 Milo 结果。

## 结果与建议

- `neighbourhood_differential_abundance.csv`：每个邻域、每个比较的样本数、mean proportion、log2 proportion fold-change、P、BH-FDR 与推断状态；
- `neighbourhood_sample_proportions.csv`：所有样本的零比例也保留，便于复核；
- `neighbourhood_metadata.csv`：seed、UMAP 位置和主导 celltype；
- `neighbourhood_da_<comparison>.png/svg`：灰色背景细胞，彩色点表示邻域 log2 proportion fold-change，FDR 达标邻域带外圈。

对于本项目，建议在 `scope_key=celltype` 中选择 epithelial 相关类型，再运行 `IBD-vs-Control`。将显著 DA 邻域与 per-celltype pseudobulk DEG 和 functional-state score 并列：前者回答局部状态是否富集，后两者回答同类细胞内部的表达/功能程序是否改变。
