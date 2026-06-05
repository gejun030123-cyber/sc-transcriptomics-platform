# Bulk RNA-seq 分析流程增强设计

**日期**: 2026-06-04
**状态**: Draft
**方案**: 方案 A（增量添加）

---

## 1. 背景与目标

当前平台的 Bulk RNA-seq 分析包含 5 个模块（QC、标准化、DEG、PCA/UMAP、热图），但缺少时序分析、基因集富集分析等主流平台标配功能。参数配置页面也没有说明文字，用户难以理解各参数的含义和建议范围。

**参考平台**: iDEP、Galaxy、DESeq2/edgeR 标准流程、clusterProfiler、OmicVerse

**目标**: 在不改动核心架构的前提下，新增 3 个模块 + 增强 4 个现有模块 + 添加参数提示系统。

---

## 2. 新增模块

### 2.1 `bulk_timecourse.py` — 时序差异分析

**功能**: 支持多时间点实验设计的差异基因检测和轨迹聚类。

**子功能**:

| 子功能 | 实现方式 | 说明 |
|--------|----------|------|
| 时序差异基因 | `ov.bulk.pyDEG().timecourse_deg(time_basis='spline')` | 自然样条基 + 调节 F 检验，q < 0.05 |
| 交互效应 | `timecourse_deg(group=..., time_basis='factor')` | 因子基（one-hot）交互 F 检验 |
| 轨迹聚类 | `ov.bulk.temporal_clusters()` | 模糊 c-means（Mfuzz），自动估计 fuzzifier m |
| 功能注释 | `ov.bulk.geneset_enrichment()` | 对每个聚类簇做 GO ORA |

**参数表单**:

| key | label | type | default | help |
|-----|-------|------|---------|------|
| `time_column` | 时间列名 | text | `minute` | obs 中表示时间点的列名。值应为数值（如 0, 15, 30, 60, 120, 180） |
| `group_column` | 分组列名（可选） | text | `` | 用于交互效应分析的分组列（如品系、处理条件）。留空则只做全时程分析 |
| `spline_df` | 样条自由度 | number | 3 | 控制时间曲线的平滑度。6 个时间点建议 3，4 个时间点建议 2。范围 2-5 |
| `n_clusters` | 轨迹聚类数 | number | 6 | 模糊 c-means 的聚类数。通常 4-8，取决于生物学复杂度 |
| `fuzzifier_m` | 模糊系数 | select | auto | auto 由算法估计，或手动设置 1.5-2.5。值越大聚类越模糊 |
| `fdr_threshold` | FDR 显著性阈值 | number | 0.05 | 时序差异基因的 q-value 截断值 |

**输出文件**:
- `timecourse_results.csv`: 全基因时序检验结果（F, pvalue, qvalue, sig）
- `timecourse_interaction.csv`: 交互效应结果（若指定 group_column）
- `timecourse_clusters.csv`: 基因聚类结果（cluster, membership）
- Plotly JSON: Q-Q 图、轨迹线图、基因×时间热图、聚类中心线图

**注册**: `modules/__init__.py` 中添加 `'bulk_timecourse': BulkTimecourseAnalysis`

---

### 2.2 `bulk_enrichment.py` — 基因集富集分析

**功能**: 对 DEG 结果做通路富集分析，支持 ORA 和 GSEA。

**实现**: 使用 `gseapy` 包（OmicVerse 底层也依赖 gseapy）。

**子功能**:

| 子功能 | 实现方式 | 说明 |
|--------|----------|------|
| ORA | `gseapy.enrichr()` 或 `ov.bulk.geneset_enrichment()` | 超几何检验，输入 DE 基因列表 |
| GSEA | `gseapy.prerank()` | 加权 Kolmogorov-Smirnov，输入全基因排序列表 |
| 数据库 | gseapy 内置库 | GO (BP/MF/CC)、KEGG、WikiPathways、Reactome |
| 物种 | Human / Mouse | 通过 organism 参数切换 |

**参数表单**:

| key | label | type | default | help |
|-----|-------|------|---------|------|
| `method` | 富集方法 | select | ORA | ORA：输入 DE 基因列表做超几何检验。GSEA：输入全基因按差异排序做秩检验 |
| `database` | 基因集数据库 | select | GO_BP | GO_BP：生物过程；GO_MF：分子功能；GO_CC：细胞组分；KEGG：代谢通路；WikiPathways：Wiki 通路；Reactome：反应组 |
| `organism` | 物种 | select | Human | Human：人类；Mouse：小鼠。基因 ID 需与所选物种匹配 |
| `pvalue_cutoff` | 显著性阈值 | number | 0.05 | 调整后 p-value 截断值，低于此值的通路被视为显著富集 |
| `top_n` | 展示通路数 | number | 20 | 可视化中显示的 Top N 显著通路 |
| `input_source` | 输入来源 | dynamic_select | 最近完成的 DEG 任务 | 从已完成的 DEG 任务中获取基因列表。ORA 模式取 sig!='normal' 的基因，GSEA 模式取全基因按 log2FC 排序 |

**输出文件**:
- `enrichment_results.csv`: 富集分析完整结果表
- Plotly JSON: 富集气泡图（x 轴 Gene Ratio，y 轴通路名，点大小=基因数，颜色=p-value）
- Plotly JSON: 条形图（Top N 通路 -log10(p-value)）
- 若 GSEA: running score 图

**注册**: `'bulk_enrichment': BulkEnrichmentAnalysis`

---

### 2.3 增强 `bulk_deg.py` — 更多统计方法

**新增方法**（在现有 t-test / Mann-Whitney 基础上）:

| 方法 | 实现 | 说明 |
|------|------|------|
| DESeq2 | `ov.bulk.pyDEG().deg_analysis(method='deseq2')` | 负二项分布模型，适用于小样本（n<10/组），RNA-seq 金标准 |

**新增参数**:

| key | label | type | default | help |
|-----|-------|------|---------|------|
| `method` | 统计方法 | select | t-test | t-test：参数检验，适合正态分布数据。Mann-Whitney：非参数检验，不假设正态。DESeq2：负二项分布模型，RNA-seq 金标准，适合小样本 |
| `fc_threshold` | Fold Change 阈值 | number | 2.0 | 基因表达变化倍数阈值。log2FC > log2(fc_threshold) 为上调，< -log2(fc_threshold) 为下调。常用值：1.5（宽松）、2.0（标准）、4.0（严格） |
| `pval_threshold` | 显著性阈值 | number | 0.05 | 调整后 p-value 截断值。0.05 为标准，0.01 为严格，0.1 为宽松 |
| `base_mean_filter` | 最低平均表达量 | number | 1 | 过滤低表达基因。BaseMean 低于此值的基因不参与分析和绘图。建议 1-10 |

---

## 3. 现有模块增强

### 3.1 `bulk_pca.py`

- 新增 **PCA 载荷图**：Top 10 基因在 PC1/PC2 上的载荷柱状图
- 新增 **肘部图**：方差解释比例累积曲线
- 新增 **降维方法切换**：参数 `dimred_method`（PCA / UMAP / t-SNE）
- 优化样本标注：n_obs > 50 时隐藏文本标签

### 3.2 `bulk_heatmap.py`

- 修复 **分组注释条**：当前代码有 `annotation_colors` 变量但未渲染到 Plotly 图中，需添加为 heatmap 的 `row_colors` 或独立 trace
- 新增 **基因聚类树状图**：在热图左侧/顶部显示 dendrogram
- 优化 **hover 信息**：显示基因名 + 样本名 + Z-score

### 3.3 `bulk_deg.py`（增强部分）

- **火山图基因标注**：自动标注 Top N 差异基因名称（Plotly annotation）
- **差异基因箱线图**：新参数 `plot_genes`，选取特定基因按分组绘制表达箱线图
- **BaseMean 过滤**：在 DEG 分析前过滤低表达基因

### 3.4 `bulk_qc.py`

- 新增 **样本相关性热图**：Pearson 相关矩阵热图
- 新增 **检测基因数 vs 文库大小散点图**（已有，优化 hover 标注样本名）

---

## 4. 参数提示系统

### 设计

在 `routes/analysis.py` 的 `PARAM_SCHEMAS` 中，为每个参数添加 `help` 字段：

```python
{
    'key': 'mito_perc',
    'label': '最大线粒体比例',
    'type': 'number',
    'default': 0.2,
    'step': 0.01,
    'help': '过滤线粒体基因比例高于此阈值的细胞/样本。人类样本建议 0.1-0.2，小鼠可放宽至 0.25。过高保留低质量细胞，过低丢失应激细胞。'
}
```

### 前端渲染

在 `templates/analysis_select.html` 的参数循环中，每个控件下方添加：

```html
{% if param.help %}
<small class="text-muted d-block mt-1" style="font-size:0.82rem; line-height:1.4;">
    {{ param.help }}
</small>
{% endif %}
```

### 覆盖范围

所有 14 个模块（9 单细胞 + 5 Bulk）的全部参数都需要补充 `help` 字段。

---

## 5. 路由变更

### `routes/analysis.py`

- `BULK_MODULE_LIST` 新增：
  - `{'name': 'bulk_timecourse', 'display': '时序分析', 'desc': '多时间点差异基因 + 轨迹聚类'}`
  - `{'name': 'bulk_enrichment', 'display': '通路富集', 'desc': 'GO/KEGG 通路富集分析 (ORA/GSEA)'}`
- `PARAM_SCHEMAS` 新增 `bulk_timecourse` 和 `bulk_enrichment` 的参数定义
- 所有现有参数的 `PARAM_SCHEMAS` 条目添加 `help` 字段

### `routes/api.py`

- 新增 `GET /api/enrichment-result/<task_id>`：返回富集分析 JSON（供前端 Plotly 渲染）
- 新增 `GET /api/obs-columns` 增强：在返回结果中新增 `time_candidates` 字段，列出 obs 中数值型且唯一值 < 20 的列名，供时序分析模块自动检测时间列

---

## 6. 模板变更

- `templates/analysis_select.html`：添加参数 `help` 渲染（见第 4 节）
- `templates/analysis_result.html`：无结构性变更，Plotly 图表已通过通用逻辑渲染

---

## 7. 依赖

| 包 | 用途 | 安装方式 |
|----|------|----------|
| `gseapy` | 基因集富集分析（ORA/GSEA） | `pip install gseapy` |
| `omicverse` | 时序分析、DESeq2 方法 | 已安装 |
| `pymfuzz` | 模糊 c-means 聚类 | `pip install pymfuzz`（omicverse 依赖） |

---

## 8. 不做的事情（YAGNI）

- 不引入 Pipeline 引擎（方案 C 的内容）
- 不重构 `BaseAnalysis` 基类
- 不添加用户登录系统
- 不添加结果对比功能
- 不支持自定义 GMT 文件上传（后续迭代）
- 不添加 t-SNE 以外的降维方法（仅在 bulk_pca 中添加 t-SNE 选项）

---

## 9. 参考来源

- OmicVerse 时序分析教程: https://omicverse.readthedocs.io/en/latest/Tutorials-bulk/t_timecourse.html
- OmicVerse DEG 教程: https://www.cnblogs.com/starlitnightly/p/18260566
- iDEP 平台: http://bioinformatics.sdstate.edu/idep/
- clusterProfiler: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
