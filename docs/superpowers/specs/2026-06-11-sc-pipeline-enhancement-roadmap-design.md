# 单细胞转录组全流程参数增强路线图

**日期**: 2026-06-11
**状态**: 已完成（P0-P6 全部实施，2026-06-18）
**方案**: 按分析流程顺序分 6 阶段渐进式增强

---

## 背景

当前平台有 10 个单细胞模块共 38 个参数，部分参数为死代码（声明但未实现）。本方案参照完整的单细胞分析流程规范，将参数扩展到 ~125 个，新增 1 个模块（细胞通讯），修复 4 个死代码参数。

### 约束

- **前端**：保持 Jinja2 服务端渲染，新增参数分"基础/高级"折叠组
- **后端**：纯 Python 生态（scanpy + omicverse + scvi-tools + squidpy + liana）
- **架构**：不改变 BaseAnalysis 接口，不改 worker/数据库结构
- **向后兼容**：所有新参数都有默认值，不填则行为与当前完全一致
- **不包含**：空间转录组、R 桥接（DESeq2/edgeR/Monocle3）

---

## 架构变更

### preprocess.py 拆分

将现有 `preprocess.py`（同时做标准化 + HVG）拆分为两个独立模块：

- **`normalize.py`**（模块名 `normalize`）：仅负责标准化
- **`hvg.py`**（模块名 `hvg`）：仅负责高变基因选择

PIPELINE_ORDER 变更：`qc → normalize → hvg → dimred → ...`

原 `preprocess.py` 中的标准化逻辑迁移到 `normalize.py`，HVG 逻辑迁移到 `hvg.py`。`preprocess.py` 保留为兼容别名或移除。

---

## P0：死代码修复

在分期实施前先修复 4 个已发现的死代码参数。

| 模块 | 参数 | 当前状态 | 处理方案 |
|------|------|---------|---------|
| `dimred.py` | `use_mde` | 声明但 run() 未读取 | **实现**：检测 `pynndescent` 可用时用 MDE 替代 UMAP 布局 |
| `annotation.py` | `method` | 声明但 run() 未读取，始终走 auto_marker | **实现**：`auto_marker`(现有逻辑) / `manual`(用户在 custom_markers 中提供完整 CellType→Cluster 映射) |
| `deg.py` | `reference` | dynamic_select 声明但 run() 未传入 rank_genes_groups | **实现**：将 reference 值传入 `sc.tl.rank_genes_groups(..., reference=ref)` |
| `trajectory.py` | `method` | 提取但只有 diffusion_map 代码路径 | **清理**：移除 schema 中未实现的 slingshot 选项，method 仅保留 diffusion_map |

---

## P1：QC + 标准化 + HVG

### 1.1 QC 模块增强（6 → 12 参数）

现有参数保留不变，新增：

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `hb_perc` | number | 0 | 基础 | 血红蛋白基因百分比上限（0=不过滤） |
| `doublet_method` | select | `scrublet` | 基础 | 双细胞检测方法：scrublet / none |
| `scrublet_threshold` | number | 0.25 | 高级 | Scrublet score 阈值 |
| `batch_adaptive_qc` | checkbox | False | 高级 | 按批次独立计算 MAD 阈值过滤 |
| `mad_multiplier` | number | 3.0 | 高级 | MAD 倍数（batch_adaptive_qc 启用时生效） |
| `save_counts_layer` | checkbox | True | 高级 | QC 过滤前保存原始 counts 到 layers["counts"] |

**实现要点**：
- 血红蛋白基因列表在 `io_utils.py` 中定义：`HB_GENES = {'HBA1', 'HBA2', 'HBB', 'HBD', 'HBG1', 'HBG2', 'HBE1', 'HBQ1', 'HBZ'}`
- `batch_adaptive_qc`：按 batch_key 分组，对 n_genes/nUMIs 使用 median ± k*MAD 计算上下界
- `save_counts_layer`：在过滤前执行 `adata.layers["counts"] = adata.X.copy()`

### 1.2 标准化模块（新：normalize.py，3 参数）

从 preprocess.py 拆出标准化逻辑：

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `method` | select | `log1p` | 基础 | 标准化方法：log1p(即 shiftlog) / pearson_residuals |
| `target_sum` | number | 10000 | 基础 | CPM 目标总数 |
| `clip_values` | checkbox | True | 高级 | Pearson 残差裁剪到 ±√n（仅 pearson_residuals 模式） |

**实现要点**：
- log1p 模式：调用 `ov.pp.preprocess(mode='shiftlog')`（现有逻辑）
- pearson_residuals 模式：调用 `scanpy.experimental.pp.normalize_pearson_residuals()` 或 omicverse 等价函数
- 标准化后保留 `adata.layers["normalized"]` 存储标准化数据

### 1.3 HVG 模块增强（新：hvg.py，3 → 12 参数）

从 preprocess.py 拆出 HVG 逻辑并大幅扩展：

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `n_top_genes` | number | 2000 | 基础 | HVG 数量（已有） |
| `batch_key` | text | `` | 基础 | 批次列名（已有） |
| `hvg_flavor` | select | `seurat_v3` | 基础 | HVG 选择方法：seurat_v3 / cell_ranger / seurat |
| `batch_hvg_strategy` | select | `intersection` | 高级 | 批次 HVG 合并策略：intersection / union |
| `exclude_mt_genes` | checkbox | False | 高级 | 排除线粒体基因（MT- 前缀） |
| `exclude_cc_genes` | checkbox | False | 高级 | 排除细胞周期基因 |
| `force_include_genes` | textarea | `` | 高级 | 强制包含的基因（换行或逗号分隔） |
| `cc_scoring` | checkbox | False | 高级 | 计算 S/G2M 期评分并存入 obs |
| `regress_cc` | checkbox | False | 高级 | 回归去除细胞周期效应（需先 cc_scoring） |

**实现要点**：
- `hvg_flavor`：传入 `sc.pp.highly_variable_genes(flavor=...)`
- `batch_hvg_strategy`：按 batch_key 分组选 HVG，然后 intersection（取交集）或 union（取并集）
- `exclude_mt_genes`：过滤 var_names 中以 "MT-" / "mt-" 开头的基因
- `exclude_cc_genes`：使用 scanpy 内置的细胞周期基因列表（`sc.datasets.cell_cycle_genes`）
- `force_include_genes`：在 HVG 结果中将指定基因的 `highly_variable` 设为 True
- `regress_cc`：调用 `sc.pp.regress_out(adata, ['S_score', 'G2M_score'])`

---

## P2：降维 + 批次校正

### 2.1 降维模块增强（2 → 10 参数）

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `n_comps` | number | 50 | 基础 | PCA 主成分数量（已有） |
| `auto_n_comps` | select | `none` | 高级 | 自动选 PC：none / elbow / kneedle |
| `umap_n_neighbors` | number | 15 | 基础 | UMAP 邻居数 |
| `umap_min_dist` | number | 0.5 | 基础 | UMAP 最小距离 |
| `umap_metric` | select | `euclidean` | 高级 | 距离度量：euclidean / cosine / correlation / manhattan |
| `umap_spread` | number | 1.0 | 高级 | UMAP spread |
| `enable_tsne` | checkbox | False | 高级 | 额外生成 t-SNE 嵌入 |
| `tsne_perplexity` | number | 30 | 高级 | t-SNE 困惑度 |
| `tsne_learning_rate` | number | 1000 | 高级 | t-SNE 学习率 |
| `use_mde` | checkbox | False | 高级 | 使用 MDE 替代 UMAP（P0 修复后） |

**实现要点**：
- `auto_n_comps=elbow`：PCA 后分析方差解释比例曲线，找拐点
- `auto_n_comps=kneedle`：使用 `kneed` 包自动检测肘部
- `umap_*`：传入 `sc.tl.umap(min_dist=..., spread=..., metric=...)`
- `enable_tsne`：调用 `sc.tl.tsne(perplexity=..., learning_rate=...)`
- `use_mde`：使用 `pynndescent` 的 MDE 算法生成 2D 嵌入

### 2.2 批次校正模块增强（4 → 12 参数）

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `method` | select | `harmony` | 基础 | 已有，扩展选项 |
| `batch_key` | text | `batch` | 基础 | 已有 |
| `n_pcs` | number | 50 | 基础 | 已有 |
| `max_epochs` | number | 200 | 高级 | 已有（仅 scVI/SysVI） |
| `harmony_theta` | number | 2.0 | 高级 | Harmony 多样性惩罚 |
| `harmony_lambda` | number | 1.0 | 高级 | Harmony 正则化 |
| `bbknn_neighbors_within_batch` | number | 3 | 高级 | BBKNN 每批次邻居数 |
| `scvi_n_latent` | number | 30 | 高级 | scVI 潜在维度 |
| `scvi_n_hidden` | number | 128 | 高级 | scVI 隐藏层大小 |
| `scvi_n_layers` | number | 1 | 高级 | scVI 层数 |
| `scvi_dropout_rate` | number | 0.1 | 高级 | scVI dropout |
| `evaluate_correction` | checkbox | False | 高级 | 计算批次校正评估指标 |

**新增方法**：
- `bbknn`：调用 `bbknn.bbknn(adata, batch_key=..., neighbors_within_batch=...)`
- `scanorama`：调用 `scanorama.integrate_scanpy()`
- scVI：扩展现有 SysVI 路径的参数暴露

**评估指标**（evaluate_correction=True 时）：
- ASW（平均轮廓宽度）：批次混合 + 细胞类型纯度
- LISI（局部逆辛普森指数）：计算 iLISI 和 cLISI
- 图连通性：校正后邻居图的连通性

---

## P3：聚类 + QC 重评估

### 3.1 聚类模块增强（2 → 9 参数）

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `resolutions` | text | `0.6,0.8,1.0` | 基础 | 已有 |
| `n_neighbors` | number | 15 | 基础 | 已有 |
| `clustering_method` | select | `leiden` | 基础 | 算法：leiden / louvain |
| `n_iterations` | number | 2 | 高级 | Leiden 迭代次数 |
| `distance_metric` | select | `euclidean` | 高级 | 邻居图距离度量 |
| `use_corrected` | checkbox | True | 高级 | 使用校正后表示（如有 batch_correct 输出） |
| `auto_select_resolution` | checkbox | False | 高级 | 自动选择最优分辨率 |
| `resolution_metric` | select | `silhouette` | 高级 | 评估指标：silhouette / calinski / davies_bouldin |
| `subcluster_key` | text | `` | 高级 | 指定父簇进行子聚类（留空=不子聚类） |

**实现要点**：
- `clustering_method=louvain`：调用 `sc.tl.louvain(resolution=...)`
- `auto_select_resolution`：对每个分辨率计算聚类质量指标，选最优
- `subcluster_key`：对指定簇的子集重跑 Leiden，结果存入新列

### 3.2 QC 重评估增强（3 → 6 参数）

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `cluster_key` | text | `leiden` | 基础 | 已有 |
| `doublet_threshold` | number | 0.3 | 基础 | 已有 |
| `mt_threshold` | number | 15.0 | 基础 | 已有 |
| `ribosomal_threshold` | number | 0 | 高级 | 核糖体比例阈值（0=不检查） |
| `min_cells_per_cluster` | number | 10 | 高级 | 最小细胞数阈值 |
| `auto_remove` | checkbox | False | 高级 | 自动移除低质量簇（输出过滤后 h5ad） |

---

## P4：注释 + DE

### 4.1 注释模块增强（5 → 14 参数）

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `method` | select | `auto_marker` | 基础 | P0 修复后：auto_marker / manual / celltypist |
| `cluster_key` | text | `leiden` | 基础 | 已有 |
| `resolution` | text | `0.8` | 基础 | 已有 |
| `marker_set` | select | `TME` | 基础 | 已有 |
| `custom_markers` | textarea | `` | 基础 | 已有 |
| `scoring_method` | select | `scanpy` | 高级 | 评分方法：scanpy / ucell |
| `assignment_strategy` | select | `best_score` | 高级 | 分配策略：best_score / threshold |
| `min_score_threshold` | number | 0.2 | 高级 | 最低分阈值（threshold 策略） |
| `celltypist_model` | select | `Immune_All_Low` | 高级 | CellTypist 模型（celltypist 模式） |
| `celltypist_threshold` | number | 0.5 | 高级 | CellTypist 概率阈值 |
| `celltypist_majority_voting` | checkbox | True | 高级 | CellTypist 多数投票 |
| `confidence_method` | select | `entropy` | 高级 | 置信度：entropy / score_margin |
| `mark_unknown` | checkbox | True | 高级 | 低置信度标 Unknown |
| `merge_similar_threshold` | number | 0 | 高级 | 相似簇合并（0=不合并） |

**新增方法**：
- `manual`（P0）：用户在 custom_markers 中提供完整映射，格式 `ClusterID:CellType` 每行一个
- `celltypist`：调用 `celltypist.annotate()`，支持模型选择、概率阈值、多数投票

### 4.2 DE 模块增强（7 → 13 参数）

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `groupby` | text | `` | 基础 | 已有 |
| `reference` | dynamic_select | `rest` | 基础 | P0 修复后生效 |
| `method` | select | `wilcoxon` | 基础 | 已有，扩展 t-test_overestim_var |
| `n_genes` | number | 20 | 基础 | 已有 |
| `show_dotplot` | checkbox | True | 基础 | 已有 |
| `plot_genes_umap` | text | `` | 基础 | 已有 |
| `custom_dotplot_genes` | textarea | `` | 基础 | 已有 |
| `pval_cutoff` | number | 0.05 | 高级 | padj 过滤阈值 |
| `logfc_cutoff` | number | 1.0 | 高级 | logFC 过滤阈值 |
| `min_pct` | number | 0.1 | 高级 | 最小表达比例 |
| `correction_method` | select | `benjamini_hochberg` | 高级 | 多重检验校正 |
| `volcano_top_n` | number | 10 | 高级 | 火山图标注基因数 |
| `volcano_genes` | textarea | `` | 高级 | 火山图自定义标注基因 |

---

## P5：轨迹 + 比例

### 5.1 轨迹模块增强（3 → 9 参数）

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `method` | select | `diffusion_map` | 基础 | P0 清理后仅 diffusion_map |
| `cluster_key` | text | `leiden` | 基础 | 已有 |
| `plot_genes` | textarea | `` | 基础 | 已有 |
| `enable_paga` | checkbox | False | 高级 | 启用 PAGA 轨迹图 |
| `paga_threshold` | number | 0.05 | 高级 | PAGA 连接阈值 |
| `n_diffcomps` | number | 15 | 高级 | 扩散图成分数量 |
| `start_cluster` | text | `` | 高级 | 起始簇（留空=自动） |
| `n_dcs` | number | 10 | 高级 | 伪时间计算用扩散成分 |
| `n_branchings` | number | 0 | 高级 | 分支点数量 |

**实现要点**：
- `enable_paga`：调用 `sc.tl.paga()` + `sc.pl.paga()`，叠加在 UMAP 上
- `start_cluster`：传入 `sc.tl.dpt(root=...)` 的 root cell 索引
- `n_branchings`：传入 `sc.tl.dpt(n_branchings=...)`

### 5.2 比例分析增强（3 → 6 参数）

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `groupby` | text | `celltype` | 基础 | 已有 |
| `batch_key` | text | `batch` | 基础 | 已有 |
| `compare_groups` | text | `` | 基础 | 已有 |
| `stat_test` | select | `chi_square` | 高级 | 统计检验：chi_square / fisher_exact / permutation |
| `n_permutations` | number | 1000 | 高级 | 置换检验次数 |
| `min_cells_per_group` | number | 10 | 高级 | 最小细胞数过滤 |

---

## P6：细胞通讯(新) + 可视化 + 导出

### 6.1 细胞通讯模块（新：cell_communication.py，8 参数）

| 参数 | 类型 | 默认值 | 分组 | 说明 |
|------|------|-------|------|------|
| `method` | select | `liana` | 基础 | 分析方法：liana |
| `resource` | select | `consensus` | 基础 | LR 数据库：consensus / cellcall / cellchatdb / omnipath |
| `organism` | select | `human` | 基础 | 物种：human / mouse |
| `cluster_key` | text | `celltype` | 基础 | 细胞类型列名 |
| `min_prop` | number | 0.1 | 高级 | 最小表达细胞比例 |
| `aggregate_method` | select | `cellchat` | 高级 | 聚合方法：cellchat / cellphone / geometric_mean |
| `top_n_interactions` | number | 20 | 高级 | 展示 Top N |
| `show_heatmap` | checkbox | True | 高级 | 通讯热图 |

**实现要点**：
- 调用 `liana.mt.rank_aggregate()` 获取配体-受体相互作用排名
- 生成交互气泡图、通讯热图、网络图
- 结果存入 adata.uns["liana_results"] 并导出 CSV

### 6.2 可视化全局参数

在 `analysis_select.html` 中新增"可视化设置"折叠组，参数存入各模块 params：

| 参数 | 类型 | 默认值 | 说明 |
|------|------|-------|------|
| `plot_theme` | select | `default` | 全局主题：default / dark / minimal |
| `color_palette` | select | `default` | 色板：default / tab20 / set3 / custom |
| `umap_point_size` | number | 5 | UMAP 点大小 |
| `umap_opacity` | number | 0.7 | UMAP 透明度 |
| `umap_legend_fontsize` | number | 10 | 图例字号 |

这些参数通过 `self.get_plotly_layout()` 统一应用到所有图表。

### 6.3 结果导出增强

在 `BaseAnalysis` 中新增 `export_results()` 方法，支持：

| 参数 | 类型 | 默认值 | 说明 |
|------|------|-------|------|
| `export_format` | select | `h5ad` | 输出格式：h5ad / csv / loom |
| `include_layers` | multiselect | counts | 导出层 |
| `include_obsm` | multiselect | X_pca, X_umap | 导出 obsm 键 |
| `compression` | select | `gzip` | 压缩：gzip / none |

---

## 前端组织

### 参数表单分组

每个模块的 `analysis_select.html` 表单分为两个折叠区域：
1. **基础设置**（默认展开）：最常用的 3-5 个参数
2. **高级设置**（默认折叠）：其余参数

### 模块切换

拆分后的 PIPELINE_ORDER：
```
qc → normalize → hvg → dimred → batch_correct → clustering → qc_reassess → annotation → deg → trajectory → proportion → cell_communication
```

### 预设模板

在 presets 系统中新增 SC 预设：
- **免疫分析**：QC 宽松 + log1p + seurat_v3 + harmony + leiden 0.8 + TME markers
- **肿瘤微环境**：QC 严格 + pearson + batch-aware HVG + bbknn + 高分辨率聚类
- **发育轨迹**：标准 QC + diffusion_map + PAGA + 伪时间

---

## 实施优先级与估算

| 阶段 | 涉及文件 | 新增参数 | 估算工作量 |
|------|---------|---------|-----------|
| P0 | dimred.py, annotation.py, deg.py, trajectory.py | 0（修复死代码） | 小 |
| P1 | qc.py, preprocess.py→normalize.py+hvg.py, routes/analysis.py, analysis_select.html | +17 | 中大 |
| P2 | dimred.py, batch_correct.py | +16 | 中 |
| P3 | clustering.py, qc_reassess.py | +10 | 中 |
| P4 | annotation.py, deg.py | +15 | 中大 |
| P5 | trajectory.py, proportion.py | +10 | 中 |
| P6 | cell_communication.py(新), visualization.py, base.py | +13 | 中 |

每个阶段独立可交付、可测试，完成后即可上线使用。
