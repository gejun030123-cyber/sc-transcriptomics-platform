# 单细胞转录组模块完善设计

**日期**: 2026-06-05
**参考**: Scanpy 标准教程 https://scanpy.scverse.org/en/stable/tutorials/basics/clustering.html

---

## 1. qc.py 增强

**新增参数**:
- `batch_key` (text, default 'batch'): Scrublet 批次感知的批次列名

**代码变更**:
- 在 `sc.pp.scrublet()` 添加 `batch_key` 参数
- 归一化前保存 `adata.layers["counts"] = adata.X.copy()`
- 新增散点图：`sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")` 风格的 Plotly 图

---

## 2. preprocess.py 增强

**新增参数**:
- `batch_key` (text, default ''): HVG 批次校正的批次列名

**代码变更**:
- 归一化前保存 `adata.layers["counts"] = adata.X.copy()`
- `sc.pp.highly_variable_genes()` 添加 `batch_key` 参数（非空时生效）

---

## 3. clustering.py 增强

**变更**:
- 使用 `flavor="igraph"` 和 `n_iterations=2`
- 生成多分辨率聚类后，添加 UMAP 并排比较图（每个分辨率一个 UMAP 子图）
- 聚类后自动生成 marker 基因 dotplot（使用内置 marker 基因集，参照 annotation.py 的 DEFAULT_TME_MARKERS）

---

## 4. 新增 qc_reassess.py — 聚类后 QC 重新评估

**模块名**: `qc_reassess`
**显示名**: QC 重新评估
**描述**: 聚类后检查 doublet 和 QC 指标，标记低质量簇

**参数**:
- `cluster_key` (text, default 'leiden'): 聚类列名
- `doublet_threshold` (number, default 0.3): doublet 比例高于此值的簇标记为低质量
- `mt_threshold` (number, default 15): MT 比例高于此值的簇标记为低质量

**功能**:
1. 计算每个簇的 doublet 比例和平均 MT 比例
2. 生成 UMAP 图：按 doublet_score 和 pct_counts_mt 着色
3. 生成 UMAP 图：按聚类着色 + 低质量簇高亮
4. 输出 low_quality_clusters.csv（簇ID、doublet比例、MT比例、是否低质量）
5. 在 summary 中报告低质量簇数量

---

## 5. annotation.py 增强

**新增参数**:
- `marker_set` (select, default 'TME'): marker 基因集选择（TME / Immune / Blood / Custom）

**变更**:
- 扩展 marker 基因集：添加 Immune 和 Blood marker sets（参照 Scanpy tutorial）
- 聚类注释后生成 dotplot（使用 Plotly heatmap 模拟，或保存为 base64 PNG）
- 支持用户自定义 marker 基因（通过 `custom_markers` 文本参数，格式：`CellType1:GENE1,GENE2;CellType2:GENE3,GENE4`）

---

## 6. deg.py 增强

**新增参数**:
- `show_dotplot` (checkbox, default True): 是否生成 DEG dotplot
- `plot_genes_umap` (text, default ''): 在 UMAP 上展示的基因名（逗号分隔）

**变更**:
- 新增 DEG dotplot：Top N 基因 × 聚类的点图（Plotly heatmap 模拟）
- 新增基因表达 UMAP 图：指定基因在 UMAP 上的表达分布
- 添加 `sc.get.rank_genes_groups_df()` 结果导出为 CSV

---

## 7. 路由/参数变更

- `routes/analysis.py`:
  - `SC_MODULE_LIST` 新增 `qc_reassess` 模块
  - `PARAM_SCHEMAS` 新增 `qc_reassess` 参数定义
  - 更新 `qc`、`preprocess`、`clustering`、`annotation`、`deg` 的参数定义
- `modules/__init__.py`: 注册 `qc_reassess` 模块
