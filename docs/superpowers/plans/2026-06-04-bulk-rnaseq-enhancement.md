# Bulk RNA-seq 分析流程增强实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 补齐 Bulk RNA-seq 分析能力：新增时序分析、通路富集模块，增强现有模块可视化，添加参数说明系统。

**Architecture:** 在现有 Flask + 模块注册架构上增量添加 2 个新模块（`bulk_timecourse`、`bulk_enrichment`），增强 4 个现有模块，为所有参数添加 `help` 字段并在前端渲染。时序分析基于 patsy 样条基 + statsmodels 线性模型实现，富集分析基于 OmicVerse 内置 `geneset_enrichment` / `geneset_enrichment_GSEA`。

**Tech Stack:** Python 3.13, Flask, scanpy, omicverse 2.1.2, plotly, patsy, statsmodels, scipy, fcmeans

---

## 前置依赖安装

在开始实现前，确认以下依赖已安装：

```bash
pip install gseapy fcmeans patsy statsmodels
```

---

## 文件结构

### 新建文件
- `modules/bulk_timecourse.py` — 时序差异基因 + 轨迹聚类模块
- `modules/bulk_enrichment.py` — 通路富集分析模块（ORA + GSEA）
- `tests/test_bulk_timecourse.py` — 时序分析模块测试
- `tests/test_bulk_enrichment.py` — 富集分析模块测试
- `tests/test_param_help.py` — 参数帮助系统测试

### 修改文件
- `routes/analysis.py` — 新增模块注册、参数定义（含 help）、模块列表
- `routes/api.py` — 新增 enrichment-result 接口、增强 obs-columns 接口
- `templates/analysis_select.html` — 参数 help 渲染
- `modules/__init__.py` — 注册新模块
- `modules/bulk_deg.py` — 增强 DEG 方法和可视化
- `modules/bulk_pca.py` — 增强 PCA 可视化
- `modules/bulk_heatmap.py` — 修复注释条、增强热图
- `modules/bulk_qc.py` — 新增样本相关性热图

---

## Task 1: 参数帮助系统

**Files:**
- Modify: `routes/analysis.py` — 为所有 `PARAM_SCHEMAS` 条目添加 `help` 字段
- Modify: `templates/analysis_select.html` — 渲染 `help` 文本

### Step 1: 为单细胞模块参数添加 help 字段

在 `routes/analysis.py` 的 `PARAM_SCHEMAS` 中，为每个参数对象添加 `help` 键。以下列出所有需要添加的 help 文本：

```python
# qc 模块
{'key': 'mito_perc', ..., 'help': '过滤线粒体基因比例高于此阈值的细胞。人类样本建议 0.1-0.2，小鼠可放宽至 0.25。过高保留低质量细胞，过低丢失应激细胞。'},
{'key': 'nUMIs', ..., 'help': '每个细胞的最小 UMI 总数。低于此值的细胞被视为碎片或死细胞。常用范围 500-1000。'},
{'key': 'detected_genes', ..., 'help': '每个细胞检测到的最小基因数。低于此值的细胞可能为低质量或空液滴。常用范围 200-500。'},
{'key': 'batch_key', ..., 'help': 'adata.obs 中标识实验批次的列名。用于分批次运行 Scrublet 双细胞检测。'},

# preprocess 模块
{'key': 'n_top_genes', ..., 'help': '选择的高变异基因数量。2000 为标准值，适合大多数分析。基因数过少会丢失生物学信号，过多会引入噪声。'},
{'key': 'target_sum', ..., 'help': '每个细胞标准化后的总计数目标。10000 为 scanpy 默认值，设为 None 则中位数标准化。'},

# dimred 模块
{'key': 'n_comps', ..., 'help': 'PCA 主成分数量。通常 30-50 即可捕获大部分方差。过多会引入噪声维度。可通过肘部图选择。'},
{'key': 'use_mde', ..., 'help': '使用 Minimum Dystortion Embedding 替代标准 UMAP，速度更快但结果略有差异。'},

# batch_correct 模块
{'key': 'method', ..., 'help': '批次校正方法。Harmony：速度快，推荐首选。ComBat：适用于已知批次的情况。SysVI：深度学习方法，适合复杂批次效应，需要 GPU。'},
{'key': 'batch_key', ..., 'help': 'adata.obs 中标识批次的列名。必须存在且包含至少 2 个不同的批次值。'},
{'key': 'n_pcs', ..., 'help': '用于批次校正的主成分数量。通常与降维步骤使用的 n_comps 一致。'},
{'key': 'max_epochs', ..., 'help': 'SysVI 方法的最大训练轮数。仅对 SysVI 方法有效。200 为默认值，复杂数据可能需要更多。'},

# clustering 模块
{'key': 'resolutions', ..., 'help': 'Leiden 聚类分辨率，多个值用逗号分隔。值越大聚类越细（簇越多）。0.4-0.6 适合粗分，0.8-1.0 标准，>1.0 细分。建议测试多个值。'},
{'key': 'n_neighbors', ..., 'help': '构建 KNN 图时的邻居数量。值越大聚类越平滑，越小越敏感。15 为默认值，小数据集可降至 10。'},

# annotation 模块
{'key': 'method', ..., 'help': '注释方法。auto_marker：使用内置 TME marker 基因自动打分。manual：手动指定 marker 基因。'},
{'key': 'cluster_key', ..., 'help': '用于分组的聚类列名。通常为 leiden 或 leiden_0.8 等。'},
{'key': 'resolution', ..., 'help': '对应的 Leiden 分辨率，用于定位正确的聚类列。'},

# deg 模块
{'key': 'groupby', ..., 'help': '差异分析的分组依据列名。如 celltype、leiden、condition 等。留空则自动使用 leiden。'},
{'key': 'reference', ..., 'help': '参考组。rest 表示以所有其他组为对照。选择特定组则以该组为对照。'},
{'key': 'method', ..., 'help': '统计检验方法。Wilcoxon：非参数检验，最常用。t-test：参数检验。logreg：逻辑回归。'},
{'key': 'n_genes', ..., 'help': '每个簇显示的 Top N 差异基因数。20 为标准值。'},

# trajectory 模块
{'key': 'method', ..., 'help': '轨迹推断方法。diffusion_map：基于扩散图的拟时序，适合连续过渡。slingshot：基于 MST 的轨迹，适合分支结构。'},
{'key': 'cluster_key', ..., 'help': '用于轨迹推断的聚类列名。'},

# proportion 模块
{'key': 'groupby', ..., 'help': '统计比例的细胞类型列名。通常为 celltype 或 leiden。'},
{'key': 'batch_key', ..., 'help': '用于比较的分组列名。如 batch、condition、treatment 等。'},
```

### Step 2: 为 Bulk 模块参数添加 help 字段

```python
# bulk_qc 模块
{'key': 'min_counts', ..., 'help': '最小文库 reads 数。低于此值的样本被过滤。人类/小鼠 RNA-seq 通常要求 ≥100000，小样本可降至 50000。'},
{'key': 'min_genes', ..., 'help': '每个样本检测到的最小基因数。低于此值的样本可能质量差。通常 5000-8000。'},
{'key': 'max_mt_pct', ..., 'help': '最大线粒体基因比例（%）。高于此值的样本可能降解严重。RNA-seq 通常 15-20%。'},

# bulk_normalize 模块
{'key': 'method', ..., 'help': '标准化方法。DESeq2：中位比率法，适用于差异分析前标准化，RNA-seq 金标准。CPM：每百万计数，简单但不考虑组成偏差。log2_quantile：分位数标准化，适合样本间可比性要求高的场景。'},

# bulk_deg 模块
{'key': 'groupby', ..., 'help': '分组列名。adata.obs 中用于区分实验组和对照组的列。如 condition、treatment、group。'},
{'key': 'group1', ..., 'help': '实验组名称。将与对照组比较计算差异基因。'},
{'key': 'group2', ..., 'help': '对照组名称。rest 表示以所有其他样本为对照。'},
{'key': 'method', ..., 'help': '统计方法。t-test：参数检验，适合正态分布数据，速度快。Mann-Whitney：非参数检验，不假设正态分布，更稳健。'},
{'key': 'fc_threshold', ..., 'help': 'Fold Change 阈值。log2FC > log2(fc) 为上调，< -log2(fc) 为下调。常用值：1.5（宽松）、2.0（标准）、4.0（严格）。'},
{'key': 'pval_threshold', ..., 'help': '调整后 p-value 显著性阈值。0.05 为标准，0.01 为严格，0.1 为宽松探索性分析。'},
{'key': 'top_n', ..., 'help': '结果中展示的 Top N 差异基因数。用于火山图标注和 Top 基因 CSV 导出。'},

# bulk_pca 模块
{'key': 'n_comps', ..., 'help': 'PCA 主成分数量。通常 5-10 即可。样本数少时自动降至 n_samples-1。'},
{'key': 'color_by', ..., 'help': '用于着色的 obs 列名。留空则不着色。如 condition、batch、celltype 等。'},
{'key': 'run_umap', ..., 'help': '是否同时运行 UMAP 降维。需要至少 10 个样本。'},

# bulk_heatmap 模块
{'key': 'heatmap_type', ..., 'help': '热图类型。top_var：显示最高变异的基因。deg：显示差异表达基因（需先运行 DEG 分析）。'},
{'key': 'top_n', ..., 'help': '热图中显示的基因数量。通常 30-100。过多会导致热图难以阅读。'},
{'key': 'groupby', ..., 'help': '样本分组列名，用于在热图旁添加分组注释条。留空则不添加。'},
```

### Step 3: 在前端渲染 help 文本

在 `templates/analysis_select.html` 中，参数表单循环内，在每个控件的 `</div>` 结束前添加：

```html
{% if param.help %}
<small class="text-muted d-block mt-1" style="font-size:0.82rem; line-height:1.4;">{{ param.help }}</small>
{% endif %}
```

在 `base.html` 的 `<style>` 中添加：

```css
.param-form small.text-muted { color: #6c757d !important; border-left: 2px solid #dee2e6; padding-left: 8px; margin-top: 4px; }
```

### Step 4: 验证

启动 Flask 应用，访问任一分析模块页面，确认参数下方显示灰色说明文字。

```bash
cd /data/GJ/platform && python app.py &
curl -s http://localhost:5000/projects/ 2>&1 | head -5
```

### Step 5: 提交

```bash
git add routes/analysis.py templates/analysis_select.html templates/base.html
git commit -m "feat: add parameter help text to all Bulk and SC analysis modules"
```

---

## Task 2: 增强 bulk_deg.py — DESeq2 方法 + 可视化增强

**Files:**
- Modify: `modules/bulk_deg.py`
- Modify: `routes/analysis.py` — 更新 bulk_deg 的 PARAM_SCHEMAS

### Step 1: 添加 DESeq2 方法支持

在 `modules/bulk_deg.py` 的 `run()` 方法中，在现有的 `method` 分支后添加 DESeq2 支持。在统计检验循环处修改：

```python
# 替换现有的逐基因检验循环
if method == 'deseq2':
    import omicverse as ov
    # 构造 DESeq2 所需的 count matrix
    count_df = pd.DataFrame(counts.T, index=gene_ids, columns=adata.obs.index)
    dds_obj = ov.bulk.pyDEG(count_df)
    dds_obj.drop_duplicates_index()
    dds_obj.normalize()
    result_df = dds_obj.deg_analysis(
        list(adata.obs.index[mask1]), list(adata.obs.index[mask2]),
        method='DEseq2', alpha=pval_threshold
    )
    log2fc = result_df['log2FC'].values if 'log2FC' in result_df.columns else result_df.iloc[:, 1].values
    pvalues = result_df['pvalue'].values if 'pvalue' in result_df.columns else result_df.iloc[:, 3].values
    gene_names = result_df.index.tolist()
    gene_ids = gene_names
    # Recalculate means
    mean1 = np.mean(data1, axis=0)
    mean2 = np.mean(data2, axis=0)
else:
    # 原有的 t-test / mann-whitney 逻辑
    ...
```

### Step 2: 添加 BaseMean 过滤参数

在 `PARAM_SCHEMAS` 的 `bulk_deg` 条目中添加：

```python
{'key': 'base_mean_filter', 'label': '最低平均表达量', 'type': 'number', 'default': 1, 'step': 0.5,
 'help': '过滤低表达基因。BaseMean 低于此值的基因不参与分析和绘图。建议 1-10。'}
```

在 `run()` 方法中，DEG 计算前添加过滤：

```python
base_mean_filter = float(self.params.get('base_mean_filter', 1))
raw_mean = counts.mean(axis=0)
expr_mask = raw_mean > base_mean_filter
counts = counts[:, expr_mask]
gene_ids = [gene_ids[i] for i in range(len(gene_ids)) if expr_mask[i]]
gene_names = [gene_names[i] for i in range(len(gene_names)) if expr_mask[i]]
```

### Step 3: 火山图添加基因标注

在火山图生成代码后，添加 Top N 基因的标注：

```python
# 在 fig_vol 的 traces 添加后
top_genes = deg_df[deg_df['regulation'] != 'NS'].head(top_n)
for _, row in top_genes.iterrows():
    fig_vol.add_annotation(
        x=row['log2FC'], y=-np.log10(max(row['padj'], 1e-300)),
        text=row['gene'], showarrow=True, arrowhead=2,
        font=dict(size=9, color='#333'), ax=20, ay=-30
    )
```

### Step 4: 添加差异基因箱线图

在 `PARAM_SCHEMAS` 的 `bulk_deg` 条目中添加：

```python
{'key': 'plot_genes', 'label': '展示基因（逗号分隔，可选）', 'type': 'text', 'default': '',
 'help': '指定要绘制箱线图的基因名，多个用逗号分隔。留空则不生成箱线图。'}
```

在 `run()` 方法的结果文件生成部分添加：

```python
plot_genes_str = self.params.get('plot_genes', '').strip()
if plot_genes_str:
    plot_gene_list = [g.strip() for g in plot_genes_str.split(',') if g.strip()]
    for pg in plot_gene_list:
        if pg in adata.var_names:
            fig_box = go.Figure()
            for grp_name, grp_mask in [(group1, mask1), (group2, mask2)]:
                vals = counts[grp_mask.values, list(adata.var_names).index(pg)]
                fig_box.add_trace(go.Box(y=vals, name=str(grp_name), boxpoints='all', jitter=0.3))
            fig_box.update_layout(title=f'{pg} 表达', yaxis_title='Expression',
                                 plot_bgcolor='white', width=400, height=350)
            fpath = os.path.join(plots_dir, f'bulk_deg_box_{pg}.json')
            with open(fpath, 'w') as f: f.write(fig_box.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'boxplot', 'label': f'{pg} Boxplot'})
```

### Step 5: 更新 PARAM_SCHEMAS 中 method 的 options

```python
{'key': 'method', 'label': '统计方法', 'type': 'select',
 'options': ['t-test', 'mann-whitney', 'deseq2'], 'default': 't-test',
 'help': '...'}
```

### Step 6: 验证

```bash
cd /data/GJ/platform
python -c "
from modules.bulk_deg import BulkDEGAnalysis
print('BulkDEGAnalysis loaded successfully')
print('Module name:', BulkDEGAnalysis.MODULE_NAME)
"
```

### Step 7: 提交

```bash
git add modules/bulk_deg.py routes/analysis.py
git commit -m "feat(bulk_deg): add DESeq2 method, BaseMean filter, volcano gene labels, boxplot"
```

---

## Task 3: 增强 bulk_pca.py

**Files:**
- Modify: `modules/bulk_pca.py`
- Modify: `routes/analysis.py` — 更新 bulk_pca 的 PARAM_SCHEMAS

### Step 1: 添加 PCA 载荷图

在 `run()` 方法的 PCA 方差图后添加：

```python
# PCA 载荷图
loadings = adata.varm['PCs'][:, :2] if 'PCs' in adata.varm else None
if loadings is not None:
    for pc_idx, pc_name in enumerate(['PC1', 'PC2']):
        top_idx = np.argsort(np.abs(loadings[:, pc_idx]))[::-1][:10]
        fig_load = go.Figure()
        fig_load.add_trace(go.Bar(
            x=[adata.var_names[i] for i in top_idx],
            y=loadings[top_idx, pc_idx],
            marker_color=['#e53935' if v > 0 else '#1a237e' for v in loadings[top_idx, pc_idx]]
        ))
        fig_load.update_layout(title=f'{pc_name} Top 10 载荷基因',
                              xaxis_title='Gene', yaxis_title='Loading',
                              plot_bgcolor='white', width=600, height=350,
                              xaxis_tickangle=45)
        fpath = os.path.join(plots_dir, f'bulk_pca_loadings_{pc_name.lower()}.json')
        with open(fpath, 'w') as f: f.write(fig_load.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': f'{pc_name} 载荷图'})
```

### Step 2: 添加肘部图

在方差图后添加累积方差曲线：

```python
fig_elbow = go.Figure()
cumvar = np.cumsum(pca_variance) * 100
fig_elbow.add_trace(go.Scatter(
    x=[f'PC{i+1}' for i in range(actual_comps)], y=cumvar,
    mode='lines+markers', marker=dict(size=6, color='#1a237e'),
    line=dict(width=2)
))
fig_elbow.add_hline(y=80, line_dash='dash', line_color='gray', annotation_text='80%')
fig_elbow.update_layout(title='PCA 方差累积曲线（肘部图）',
                       xaxis_title='主成分', yaxis_title='累积方差 (%)',
                       plot_bgcolor='white', width=600, height=350)
fpath = os.path.join(plots_dir, 'bulk_pca_elbow.json')
with open(fpath, 'w') as f: f.write(fig_elbow.to_json(engine="json"))
result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': '肘部图'})
```

### Step 3: 添加 t-SNE 选项

在 `PARAM_SCHEMAS` 的 `bulk_pca` 条目中添加：

```python
{'key': 'dimred_method', 'label': '降维方法', 'type': 'select', 'options': ['pca', 'umap', 'tsne'], 'default': 'pca',
 'help': '降维可视化方法。PCA：线性降维，保留全局结构。UMAP：非线性降维，保留局部结构。t-SNE：非线性降维，适合发现聚类。'}
```

在 `run()` 方法中，UMAP 部分改为根据参数选择方法：

```python
dimred_method = self.params.get('dimred_method', 'pca')
if dimred_method == 'tsne' and adata.n_obs >= 10:
    from sklearn.manifold import TSNE
    self.progress(75, "运行 t-SNE...")
    tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, adata.n_obs - 1))
    tsne_coords = tsne.fit_transform(adata.obsm['X_pca'])
    # ... 生成 t-SNE scatter plot
elif dimred_method == 'umap' and adata.n_obs >= 10:
    # 现有 UMAP 逻辑
    ...
```

### Step 4: 优化大样本标注

在所有 scatter 图中，当 `n_obs > 50` 时隐藏文本标注：

```python
text_labels = hover if adata.n_obs <= 50 else None
show_text = adata.n_obs <= 50
# 在 go.Scattergl 中: text=text_labels, textposition='top center' if show_text else None
```

### Step 5: 验证

```bash
python -c "
from modules.bulk_pca import BulkPCAAnalysis
print('BulkPCAAnalysis loaded, MODULE_NAME:', BulkPCAAnalysis.MODULE_NAME)
"
```

### Step 6: 提交

```bash
git add modules/bulk_pca.py routes/analysis.py
git commit -m "feat(bulk_pca): add loadings plot, elbow plot, t-SNE option, optimize large-sample labels"
```

---

## Task 4: 增强 bulk_heatmap.py 和 bulk_qc.py

**Files:**
- Modify: `modules/bulk_heatmap.py`
- Modify: `modules/bulk_qc.py`

### Step 1: 修复 bulk_heatmap.py 分组注释条

在热图生成部分，当 `annotation_colors` 不为 None 时，添加一个颜色条 trace：

```python
if annotation_colors:
    # 添加分组注释条作为独立 heatmap trace
    color_indices = [unique_groups.index(str(groups[i])) for i in range(len(groups))]
    fig_annotation = go.Figure()
    fig_annotation.add_trace(go.Heatmap(
        z=[[i] for i in color_indices],
        y=sample_ordered,
        x=['Group'],
        colorscale=[[i/(len(unique_groups)-1), palette[i % len(palette)]] for i in range(len(unique_groups))],
        showscale=False,
        text=[[unique_groups[i]] for i in color_indices],
        hovertemplate='%{y}: %{text}<extra></extra>'
    ))
    fig_annotation.update_layout(
        height=max(400, adata.n_obs * 25 + 150), width=100,
        margin=dict(l=0, r=0, t=30, b=40)
    )
    fpath_annot = os.path.join(plots_dir, 'bulk_heatmap_annotation.json')
    with open(fpath_annot, 'w') as f: f.write(fig_annotation.to_json(engine="json"))
    result_files.append({'file_path': fpath_annot, 'file_type': 'plotly_json', 'category': 'annotation', 'label': '分组注释条'})
```

### Step 2: 增强热图 hover 信息

修改热图 trace 的 hovertemplate：

```python
hovertemplate='样本: %{y}<br>基因: %{x}<br>Z-score: %{z:.2f}<extra></extra>'
```

### Step 3: bulk_qc.py 添加样本相关性热图

在 `run()` 方法的 QC 图生成部分添加：

```python
# 样本相关性热图
norm_for_corr = adata_filtered.copy()
sc.pp.normalize_total(norm_for_corr, target_sum=1e6)
sc.pp.log1p(norm_for_corr)
corr_data = norm_for_corr.X if not hasattr(norm_for_corr.X, 'toarray') else norm_for_corr.X.toarray()
corr_matrix = np.corrcoef(corr_data)
sample_labels_corr = norm_for_corr.obs.index.tolist()

fig_corr = go.Figure()
fig_corr.add_trace(go.Heatmap(
    z=corr_matrix.tolist(), x=sample_labels_corr, y=sample_labels_corr,
    colorscale='Blues', zmin=0, zmax=1,
    colorbar=dict(title='Pearson r'),
    hovertemplate='%{y} vs %{x}<br>r = %{z:.3f}<extra></extra>'
))
fig_corr.update_layout(
    title='样本相关性热图 (Pearson)',
    height=max(400, n_after * 30 + 100), width=max(400, n_after * 30 + 100),
    plot_bgcolor='white'
)
fpath = os.path.join(plots_dir, 'bulk_qc_corr.json')
with open(fpath, 'w') as f: f.write(fig_corr.to_json(engine="json"))
result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': '样本相关性热图'})
```

### Step 4: 验证

```bash
python -c "
from modules.bulk_heatmap import BulkHeatmapAnalysis
from modules.bulk_qc import BulkQCAnalysis
print('Both modules loaded successfully')
"
```

### Step 5: 提交

```bash
git add modules/bulk_heatmap.py modules/bulk_qc.py
git commit -m "feat: enhance heatmap annotation bars, add QC sample correlation heatmap"
```

---

## Task 5: 新增 bulk_enrichment.py — 通路富集分析模块

**Files:**
- Create: `modules/bulk_enrichment.py`
- Modify: `modules/__init__.py` — 注册模块
- Modify: `routes/analysis.py` — 添加模块定义和参数

### Step 1: 创建模块文件

创建 `modules/bulk_enrichment.py`：

```python
import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class BulkEnrichmentAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_enrichment"
    DISPLAY_NAME = "通路富集分析"
    DESCRIPTION = "GO/KEGG/WikiPathways 通路富集分析（ORA / GSEA）"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import omicverse as ov
        import plotly.graph_objects as go

        method = self.params.get('method', 'ORA')
        database = self.params.get('database', 'GO_BP')
        organism = self.params.get('organism', 'Human')
        pvalue_cutoff = float(self.params.get('pvalue_cutoff', 0.05))
        top_n = int(self.params.get('top_n', 20))
        input_source = self.params.get('input_source', '')

        self.progress(5, "加载基因列表...")

        # 从 DEG 结果 CSV 加载基因列表
        deg_genes = []
        gene_rnk = None

        if input_source and os.path.exists(input_source):
            deg_df = pd.read_csv(input_source)
            if 'regulation' in deg_df.columns:
                deg_genes = deg_df[deg_df['regulation'] != 'NS']['gene'].tolist()
                # GSEA 需要排序的基因列表
                if 'log2FC' in deg_df.columns:
                    gene_rnk = deg_df[['gene', 'log2FC']].dropna().sort_values('log2FC', ascending=False)
                    gene_rnk.columns = ['gene_name', 'rank']
            elif 'gene' in deg_df.columns:
                deg_genes = deg_df['gene'].tolist()

        if not deg_genes and method == 'ORA':
            raise ValueError("未找到差异基因列表。请先运行 DEG 分析，并在 input_source 中指定结果 CSV 文件路径。")

        self.progress(15, f"加载 {database} 基因集数据库...")

        # 映射数据库名称到文件
        db_map = {
            'GO_BP': 'GO_Biological_Process_2021',
            'GO_MF': 'GO_Molecular_Function_2021',
            'GO_CC': 'GO_Cellular_Component_2021',
            'KEGG': 'KEGG_2021_Human' if organism == 'Human' else 'KEGG_2019_Mouse',
            'WikiPathways': 'WikiPathways_2019_Human' if organism == 'Human' else 'WikiPathways_2019_Mouse',
            'Reactome': 'Reactome_2022',
        }

        db_filename = db_map.get(database, 'GO_Biological_Process_2021')
        organism_lower = organism.lower()

        # 下载基因集数据库（如果不存在）
        ov.utils.download_pathway_database()

        db_path = f'genesets/{db_filename}.txt'
        if not os.path.exists(db_path):
            # 尝试不同的路径格式
            alt_paths = [
                f'genesets/{db_filename}_{organism}.txt',
                f'genesets/{db_filename}_{organism_lower}.txt',
            ]
            for p in alt_paths:
                if os.path.exists(p):
                    db_path = p
                    break

        pathways_dict = ov.utils.geneset_prepare(db_path, organism=organism)

        self.progress(30, f"运行 {method} 富集分析...")

        result_files = []
        plots_dir = os.path.join(self.project_dir, 'plots')
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)

        if method == 'ORA':
            enr = ov.bulk.geneset_enrichment(
                gene_list=deg_genes,
                pathways_dict=pathways_dict,
                pvalue_type='adjust',
                pvalue_threshold=pvalue_cutoff,
                organism=organism_lower,
                outdir=os.path.join(self.project_dir, 'enrichr_tmp')
            )

            self.progress(60, "保存结果表...")
            csv_path = os.path.join(results_dir, 'enrichment_ora_results.csv')
            enr.to_csv(csv_path, index=False)
            result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': 'ORA 富集结果'})

            self.progress(70, "生成气泡图...")
            # 气泡图
            top_enr = enr.head(top_n)
            if len(top_enr) > 0:
                fig_bubble = go.Figure()
                fig_bubble.add_trace(go.Scatter(
                    x=top_enr['Fractions'] if 'Fractions' in top_enr.columns else top_enr.get('Odds Ratio', [1]*len(top_enr)),
                    y=top_enr['Term'] if 'Term' in top_enr.columns else top_enr.index.tolist(),
                    mode='markers',
                    marker=dict(
                        size=top_enr['Overlap'].apply(lambda x: int(str(x).split('/')[0]) * 3 + 5) if 'Overlap' in top_enr.columns else 10,
                        color=-np.log10(top_enr['P-value'].clip(lower=1e-300)),
                        colorscale='YlOrRd', showscale=True,
                        colorbar=dict(title='-log10(p)')
                    ),
                    text=top_enr['Term'] if 'Term' in top_enr.columns else top_enr.index.tolist(),
                    hovertemplate='%{text}<br>Fraction: %{x:.3f}<br>-log10(p): %{marker.color:.1f}<extra></extra>'
                ))
                fig_bubble.update_layout(
                    title=f'{database} ORA 富集分析 ({organism})',
                    xaxis_title='Gene Fraction',
                    yaxis=dict(autorange='reversed'),
                    plot_bgcolor='white', width=800, height=max(400, top_n * 25 + 100)
                )
                fpath = os.path.join(plots_dir, 'enrichment_ora_bubble.json')
                with open(fpath, 'w') as f: f.write(fig_bubble.to_json(engine="json"))
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'enrichment', 'label': 'ORA 气泡图'})

                # 条形图
                fig_bar = go.Figure()
                fig_bar.add_trace(go.Bar(
                    x=-np.log10(top_enr['P-value'].clip(lower=1e-300)),
                    y=top_enr['Term'] if 'Term' in top_enr.columns else top_enr.index.tolist(),
                    orientation='h',
                    marker_color='#e53935'
                ))
                fig_bar.update_layout(
                    title=f'Top {top_n} 富集通路',
                    xaxis_title='-log10(P-value)',
                    yaxis=dict(autorange='reversed'),
                    plot_bgcolor='white', width=700, height=max(400, top_n * 25 + 100)
                )
                fpath = os.path.join(plots_dir, 'enrichment_ora_bar.json')
                with open(fpath, 'w') as f: f.write(fig_bar.to_json(engine="json"))
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'enrichment', 'label': 'ORA 条形图'})

            n_sig = len(enr[enr['P-value'] < pvalue_cutoff]) if 'P-value' in enr.columns else len(enr)

        elif method == 'GSEA':
            if gene_rnk is None or len(gene_rnk) == 0:
                raise ValueError("GSEA 需要排序的基因列表。请确保 DEG 结果包含 log2FC 列。")

            pre_res = ov.bulk.geneset_enrichment_GSEA(
                gene_rnk=gene_rnk,
                pathways_dict=pathways_dict,
                processes=4,
                permutation_num=100,
                outdir=os.path.join(self.project_dir, 'enrichr_gsea_tmp')
            )

            self.progress(60, "处理 GSEA 结果...")
            enr = pre_res.res2d
            enr_sig = enr[enr['fdr'] < pvalue_cutoff].copy()

            csv_path = os.path.join(results_dir, 'enrichment_gsea_results.csv')
            enr.to_csv(csv_path, index=False)
            result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': 'GSEA 富集结果'})

            self.progress(70, "生成 GSEA 图表...")
            if len(enr_sig) > 0:
                top_gsea = enr_sig.head(top_n)
                fig_bar = go.Figure()
                fig_bar.add_trace(go.Bar(
                    x=top_gsea['nes'],
                    y=top_gsea.index.tolist(),
                    orientation='h',
                    marker_color=['#e53935' if v > 0 else '#1a237e' for v in top_gsea['nes']]
                ))
                fig_bar.update_layout(
                    title=f'Top {top_n} GSEA 通路 (NES)',
                    xaxis_title='Normalized Enrichment Score',
                    yaxis=dict(autorange='reversed'),
                    plot_bgcolor='white', width=700, height=max(400, top_n * 25 + 100)
                )
                fpath = os.path.join(plots_dir, 'enrichment_gsea_nes.json')
                with open(fpath, 'w') as f: f.write(fig_bar.to_json(engine="json"))
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'enrichment', 'label': 'GSEA NES 条形图'})

            n_sig = len(enr_sig)
        else:
            raise ValueError(f"不支持的富集方法: {method}")

        self.progress(90, "保存输出...")
        self.progress(100, "完成")

        return {
            'output_adata': None,
            'result_files': result_files,
            'summary': {
                'method': method,
                'database': database,
                'organism': organism,
                'n_input_genes': len(deg_genes) if method == 'ORA' else len(gene_rnk) if gene_rnk is not None else 0,
                'n_significant': n_sig,
                'pvalue_cutoff': pvalue_cutoff,
            }
        }
```

### Step 2: 注册模块

在 `modules/__init__.py` 中添加：

```python
from .bulk_enrichment import BulkEnrichmentAnalysis
```

在 `MODULE_REGISTRY` 中添加：

```python
'bulk_enrichment': BulkEnrichmentAnalysis,
```

在 `PIPELINE_ORDER` 末尾添加：

```python
'bulk_enrichment',
```

### Step 3: 添加路由配置

在 `routes/analysis.py` 的 `BULK_MODULE_LIST` 中添加：

```python
{'name': 'bulk_enrichment', 'display': '通路富集', 'desc': 'GO/KEGG/WikiPathways 通路富集分析（ORA / GSEA）'},
```

在 `PARAM_SCHEMAS` 中添加：

```python
'bulk_enrichment': [
    {'key': 'method', 'label': '富集方法', 'type': 'select', 'options': ['ORA', 'GSEA'], 'default': 'ORA',
     'help': 'ORA：超几何检验，输入 DE 基因列表，检验哪些通路过度代表。GSEA：秩检验，输入全基因按差异排序，检验通路在排序中的富集位置。'},
    {'key': 'database', 'label': '基因集数据库', 'type': 'select',
     'options': ['GO_BP', 'GO_MF', 'GO_CC', 'KEGG', 'WikiPathways', 'Reactome'], 'default': 'GO_BP',
     'help': 'GO_BP：生物过程。GO_MF：分子功能。GO_CC：细胞组分。KEGG：代谢和信号通路。WikiPathways：社区维护通路。Reactome：反应组数据库。'},
    {'key': 'organism', 'label': '物种', 'type': 'select', 'options': ['Human', 'Mouse'], 'default': 'Human',
     'help': 'Human：人类基因。Mouse：小鼠基因。基因 ID 需与所选物种匹配。'},
    {'key': 'pvalue_cutoff', 'label': '显著性阈值', 'type': 'number', 'default': 0.05, 'step': 0.01,
     'help': '调整后 p-value 截断值。0.05 为标准，0.01 为严格。'},
    {'key': 'top_n', 'label': '展示通路数', 'type': 'number', 'default': 20,
     'help': '可视化中显示的 Top N 显著通路数。'},
    {'key': 'input_source', 'label': 'DEG 结果文件', 'type': 'text', 'default': '',
     'help': '来自已完成的 DEG 分析的 CSV 结果文件路径。包含 gene 和 regulation/log2FC 列。'},
],
```

### Step 4: 验证

```bash
python -c "
from modules.bulk_enrichment import BulkEnrichmentAnalysis
print('BulkEnrichmentAnalysis loaded:', BulkEnrichmentAnalysis.MODULE_NAME)
"
```

### Step 5: 提交

```bash
git add modules/bulk_enrichment.py modules/__init__.py routes/analysis.py
git commit -m "feat: add bulk_enrichment module (ORA/GSEA with GO/KEGG/WikiPathways)"
```

---

## Task 6: 新增 bulk_timecourse.py — 时序分析模块

**Files:**
- Create: `modules/bulk_timecourse.py`
- Modify: `modules/__init__.py` — 注册模块
- Modify: `routes/analysis.py` — 添加模块定义和参数

### Step 1: 创建模块文件

创建 `modules/bulk_timecourse.py`。注意：OmicVerse 2.1.2 没有 `timecourse_deg`，需用 patsy + statsmodels 实现：

```python
import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class BulkTimecourseAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_timecourse"
    DISPLAY_NAME = "时序分析"
    DESCRIPTION = "多时间点差异基因检测 + 轨迹聚类"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def _build_spline_basis(self, time_vals, df=3):
        """构建自然样条基矩阵"""
        from patsy import dmatrix
        basis = dmatrix(f'bs(x, df={df}, degree=3) - 1',
                       {'x': time_vals}, return_type='dataframe')
        return basis.values

    def _run_temporal_f_test(self, counts, time_basis, n_spline_cols):
        """对每个基因做调节 F 检验（简化版）"""
        from scipy import stats
        import statsmodels.api as sm

        n_genes = counts.shape[1]
        n_obs = counts.shape[0]
        F_stats = np.zeros(n_genes)
        pvalues = np.zeros(n_genes)

        for i in range(n_genes):
            y = counts[:, i]
            try:
                # Full model: intercept + spline basis
                X_full = sm.add_constant(time_basis)
                model_full = sm.OLS(y, X_full).fit()
                # Reduced model: intercept only
                X_red = sm.add_constant(np.ones(n_obs))
                model_red = sm.OLS(y, X_red).fit()
                # F-test
                RSS_full = model_full.ssr
                RSS_red = model_red.ssr
                df_diff = n_spline_cols
                df_resid = n_obs - n_spline_cols - 1
                if RSS_full > 0 and df_resid > 0:
                    F = ((RSS_red - RSS_full) / df_diff) / (RSS_full / df_resid)
                    F_stats[i] = max(F, 0)
                    pvalues[i] = 1 - stats.f.cdf(F, df_diff, df_resid)
                else:
                    F_stats[i] = 0
                    pvalues[i] = 1.0
            except Exception:
                F_stats[i] = 0
                pvalues[i] = 1.0

        return F_stats, pvalues

    def run(self, input_path):
        import scanpy as sc
        import plotly.graph_objects as go
        from statsmodels.stats.multitest import multipletests

        self.progress(5, "加载数据...")
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(input_path)

        time_col = self.params.get('time_column', 'minute')
        group_col = self.params.get('group_column', '')
        spline_df = int(self.params.get('spline_df', 3))
        n_clusters = int(self.params.get('n_clusters', 6))
        fdr_threshold = float(self.params.get('fdr_threshold', 0.05))

        # 标准化
        counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
        counts = counts.astype(float)

        # 过滤低表达基因
        gene_mask = counts.mean(axis=0) > 1
        counts = counts[:, gene_mask]
        gene_names = [adata.var_names[i] for i in range(adata.n_vars) if gene_mask[i]]

        self.progress(10, "解析时间信息...")

        # 从 obs 获取时间向量
        if time_col in adata.obs.columns:
            time_vals = adata.obs[time_col].astype(float).values
        else:
            # 尝试从 obs 名解析
            raise ValueError(f"未找到时间列 '{time_col}'。请检查 obs 列名。")

        times_sorted = sorted(np.unique(time_vals))
        n_times = len(times_sorted)

        self.progress(20, f"构建样条基 (df={spline_df})...")
        time_basis = self._build_spline_basis(time_vals, df=spline_df)
        n_spline_cols = time_basis.shape[1]

        # 对数标准化（log2 CPM 风格）
        lib_sizes = counts.sum(axis=1, keepdims=True)
        cpm = counts / lib_sizes * 1e6
        lognorm = np.log2(cpm + 1)

        self.progress(30, "时序差异基因检验...")
        F_stats, pvalues = self._run_temporal_f_test(lognorm, time_basis, n_spline_cols)

        # BH FDR 校正
        try:
            _, qvalues, _, _ = multipletests(pvalues, method='fdr_bh')
        except Exception:
            qvalues = pvalues

        sig_mask = qvalues < fdr_threshold
        temporal_genes = [gene_names[i] for i in range(len(gene_names)) if sig_mask[i]]

        result_files = []
        plots_dir = os.path.join(self.project_dir, 'plots')
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)

        self.progress(50, "保存时序结果表...")
        tc_df = pd.DataFrame({
            'gene': gene_names,
            'F': np.round(F_stats, 4),
            'pvalue': pvalues,
            'qvalue': qvalues,
            'sig': ['temporal' if sig_mask[i] else 'normal' for i in range(len(gene_names))]
        })
        tc_df = tc_df.sort_values('qvalue')
        csv_path = os.path.join(results_dir, 'timecourse_results.csv')
        tc_df.to_csv(csv_path, index=False)
        result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': '时序差异基因结果'})

        # Q-Q 图
        self.progress(55, "生成 Q-Q 图...")
        from scipy import stats as sp_stats
        obs_F = np.sort(F_stats)
        quant = (np.arange(1, len(obs_F) + 1) - 0.5) / len(obs_F)
        df_red = n_obs - n_spline_cols - 1
        theo_F = sp_stats.f.ppf(quant, n_spline_cols, max(df_red, 1))

        fig_qq = go.Figure()
        fig_qq.add_trace(go.Scattergl(x=theo_F[:len(obs_F)], y=obs_F, mode='markers',
                                       marker=dict(size=3, color='#4878d0', opacity=0.5)))
        lim = np.nanpercentile(theo_F, 99.5)
        fig_qq.add_trace(go.Scatter(x=[0, lim], y=[0, lim], mode='lines',
                                    line=dict(color='gray', dash='dash', width=1), showlegend=False))
        fig_qq.update_layout(title='Q-Q: 时序调控检验',
                            xaxis_title='理论 F 分位数', yaxis_title='观测调节 F',
                            plot_bgcolor='white', width=500, height=500)
        fpath = os.path.join(plots_dir, 'timecourse_qq.json')
        with open(fpath, 'w') as f: f.write(fig_qq.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qq', 'label': 'Q-Q 图'})

        # 轨迹聚类
        if len(temporal_genes) >= n_clusters:
            self.progress(65, f"轨迹聚类 ({n_clusters} 簇)...")

            # 计算每个时间点的均值轨迹
            traj = np.zeros((len(temporal_genes), n_times))
            for ti, t in enumerate(times_sorted):
                mask_t = time_vals == t
                traj[:, ti] = lognorm[mask_t][:, [gene_names.index(g) for g in temporal_genes]].mean(axis=0)

            # Z-score 标准化
            traj_mean = traj.mean(axis=1, keepdims=True)
            traj_std = traj.std(axis=1, keepdims=True) + 1e-10
            traj_z = (traj - traj_mean) / traj_std

            # 模糊 c-means
            try:
                from fcmeans import FCM
                fcm = FCM(n_clusters=n_clusters, random_state=42)
                fcm.fit(traj_z)
                cluster_labels = fcm.predict(traj_z)
                membership = np.max(fcm.u, axis=1)
                centers = fcm.centers
            except ImportError:
                # Fallback: 使用 sklearn KMeans
                from sklearn.cluster import KMeans
                km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
                cluster_labels = km.fit_predict(traj_z)
                membership = np.ones(len(temporal_genes))
                centers = km.cluster_centers_

            self.progress(75, "保存聚类结果...")
            cluster_df = pd.DataFrame({
                'gene': temporal_genes,
                'cluster': cluster_labels + 1,
                'membership': np.round(membership, 3)
            })
            cluster_df = cluster_df.sort_values(['cluster', 'membership'], ascending=[True, False])
            csv_path = os.path.join(results_dir, 'timecourse_clusters.csv')
            cluster_df.to_csv(csv_path, index=False)
            result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': '轨迹聚类结果'})

            # 聚类中心线图
            self.progress(80, "生成轨迹图...")
            fig_centers = go.Figure()
            for ci in range(n_clusters):
                fig_centers.add_trace(go.Scatter(
                    x=times_sorted, y=centers[ci], mode='lines+markers',
                    name=f'C{ci+1}', marker=dict(size=5)
                ))
            fig_centers.add_hline(y=0, line_dash='dash', line_color='gray', line_width=0.7)
            fig_centers.update_layout(
                title='表达波聚类中心',
                xaxis_title='时间 (min)', yaxis_title='标准化表达',
                plot_bgcolor='white', width=700, height=400
            )
            fpath = os.path.join(plots_dir, 'timecourse_centers.json')
            with open(fpath, 'w') as f: f.write(fig_centers.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'trajectory', 'label': '聚类中心线图'})

            # 基因×时间热图（按聚类排序）
            self.progress(85, "生成热图...")
            order = np.argsort(cluster_labels * 1000 - membership)
            heat_ordered = traj_z[order]
            gene_ordered = [temporal_genes[i] for i in order]

            fig_heat = go.Figure()
            fig_heat.add_trace(go.Heatmap(
                z=heat_ordered.tolist(),
                x=[str(int(t)) for t in times_sorted],
                y=gene_ordered,
                colorscale='RdBu_r', zmid=0,
                colorbar=dict(title='Z-score'),
                hovertemplate='基因: %{y}<br>时间: %{x} min<br>Z-score: %{z:.2f}<extra></extra>'
            ))
            fig_heat.update_layout(
                title=f'时序基因热图 ({len(temporal_genes)} genes)',
                xaxis_title='时间 (min)', yaxis_title='',
                height=max(400, len(temporal_genes) * 3 + 150), width=500,
                plot_bgcolor='white'
            )
            fpath = os.path.join(plots_dir, 'timecourse_heatmap.json')
            with open(fpath, 'w') as f: f.write(fig_heat.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': '时序基因热图'})

        # 交互效应分析
        interaction_result = None
        if group_col and group_col in adata.obs.columns:
            self.progress(90, "交互效应分析...")
            groups = adata.obs[group_col].unique().tolist()
            if len(groups) >= 2:
                from patsy import dmatrix as dm
                # 构建因子基交互项
                group_dummies = pd.get_dummies(adata.obs[group_col]).values
                interaction_basis = np.hstack([time_basis, group_dummies[:, 1:]])
                # 添加交互项
                for gi in range(group_dummies.shape[1] - 1):
                    for ti in range(time_basis.shape[1]):
                        interaction_basis = np.hstack([interaction_basis,
                                                      (group_dummies[:, gi+1] * time_basis[:, ti]).reshape(-1, 1)])

                # 交互效应 F 检验
                n_interaction_cols = (group_dummies.shape[1] - 1) * time_basis.shape[1]
                F_int, p_int = self._run_interaction_f_test(lognorm, interaction_basis,
                                                           n_interaction_cols, time_basis.shape[1] + group_dummies.shape[1] - 1)
                try:
                    _, q_int, _, _ = multipletests(p_int, method='fdr_bh')
                except Exception:
                    q_int = p_int

                int_df = pd.DataFrame({
                    'gene': gene_names,
                    'F_interaction': np.round(F_int, 4),
                    'pvalue': p_int,
                    'qvalue': q_int,
                }).sort_values('pvalue')

                csv_path = os.path.join(results_dir, 'timecourse_interaction.csv')
                int_df.to_csv(csv_path, index=False)
                result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': '交互效应结果'})
                interaction_result = {
                    'n_sig_interaction': int(np.sum(q_int < fdr_threshold)),
                    'groups': groups,
                }

        self.progress(95, "保存 h5ad...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'timecourse_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "完成")

        summary = {
            'n_genes_total': len(gene_names),
            'n_temporal': len(temporal_genes),
            'pct_temporal': round(len(temporal_genes) / max(len(gene_names), 1) * 100, 1),
            'n_clusters': n_clusters,
            'spline_df': spline_df,
            'n_timepoints': n_times,
            'timepoints': [int(t) for t in times_sorted],
        }
        if interaction_result:
            summary['interaction'] = interaction_result

        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }

    def _run_interaction_f_test(self, counts, full_basis, n_interaction_cols, n_full_model_cols):
        """交互效应 F 检验"""
        from scipy import stats
        import statsmodels.api as sm

        n_genes = counts.shape[1]
        n_obs = counts.shape[0]
        F_stats = np.zeros(n_genes)
        pvalues = np.zeros(n_genes)

        # Reduced model: time + group (no interaction)
        n_reduced_cols = n_full_model_cols - n_interaction_cols

        for i in range(n_genes):
            y = counts[:, i]
            try:
                X_full = sm.add_constant(full_basis)
                model_full = sm.OLS(y, X_full).fit()
                X_red = sm.add_constant(full_basis[:, :n_reduced_cols])
                model_red = sm.OLS(y, X_red).fit()
                RSS_full = model_full.ssr
                RSS_red = model_red.ssr
                df_resid = n_obs - n_full_model_cols - 1
                if RSS_full > 0 and df_resid > 0:
                    F = ((RSS_red - RSS_full) / n_interaction_cols) / (RSS_full / df_resid)
                    F_stats[i] = max(F, 0)
                    pvalues[i] = 1 - stats.f.cdf(F, n_interaction_cols, df_resid)
            except Exception:
                F_stats[i] = 0
                pvalues[i] = 1.0

        return F_stats, pvalues
```

### Step 2: 注册模块

在 `modules/__init__.py` 中添加：

```python
from .bulk_timecourse import BulkTimecourseAnalysis
```

在 `MODULE_REGISTRY` 中添加：

```python
'bulk_timecourse': BulkTimecourseAnalysis,
```

### Step 3: 添加路由配置

在 `routes/analysis.py` 的 `BULK_MODULE_LIST` 中添加：

```python
{'name': 'bulk_timecourse', 'display': '时序分析', 'desc': '多时间点差异基因检测 + 轨迹聚类'},
```

在 `PARAM_SCHEMAS` 中添加：

```python
'bulk_timecourse': [
    {'key': 'time_column', 'label': '时间列名', 'type': 'text', 'default': 'minute',
     'help': 'obs 中表示时间点的列名。值应为数值（如 0, 15, 30, 60, 120, 180）。'},
    {'key': 'group_column', 'label': '分组列名（可选）', 'type': 'text', 'default': '',
     'help': '用于交互效应分析的分组列（如品系、处理条件）。留空则只做全时程分析。'},
    {'key': 'spline_df', 'label': '样条自由度', 'type': 'number', 'default': 3,
     'help': '控制时间曲线的平滑度。6 个时间点建议 3，4 个时间点建议 2。范围 2-5。'},
    {'key': 'n_clusters', 'label': '轨迹聚类数', 'type': 'number', 'default': 6,
     'help': '模糊 c-means 的聚类数。通常 4-8，取决于生物学复杂度。'},
    {'key': 'fdr_threshold', 'label': 'FDR 显著性阈值', 'type': 'number', 'default': 0.05, 'step': 0.01,
     'help': '时序差异基因的 q-value 截断值。0.05 为标准。'},
],
```

### Step 4: 验证

```bash
python -c "
from modules.bulk_timecourse import BulkTimecourseAnalysis
print('BulkTimecourseAnalysis loaded:', BulkTimecourseAnalysis.MODULE_NAME)
"
```

### Step 5: 提交

```bash
git add modules/bulk_timecourse.py modules/__init__.py routes/analysis.py
git commit -m "feat: add bulk_timecourse module (temporal DEG + fuzzy c-means trajectory clustering)"
```

---

## Task 7: API 增强 + 最终集成

**Files:**
- Modify: `routes/api.py` — 新增 enrichment-result 接口
- Modify: `routes/analysis.py` — 确保 enrichment 的 input_source 动态选择逻辑

### Step 1: 在 api.py 新增 enrichment 结果接口

在 `routes/api.py` 末尾添加：

```python
@api_bp.route('/enrichment-result/<task_id>')
def enrichment_result(task_id):
    """返回富集分析的 Plotly JSON 结果"""
    from models import ResultFile
    files = ResultFile.get_by_task(task_id)
    enrichment_files = [f for f in files if f.category == 'enrichment']
    result = []
    for f in enrichment_files:
        with open(f.file_path, 'r') as fh:
            data = json.load(fh)
        result.append({'id': f.id, 'label': f.label, 'data': data})
    return jsonify(result)
```

### Step 2: 在 analysis.py 中为 enrichment 添加 input_source 动态逻辑

在 `analyze` 路由的 GET 分支中，当 `module_name == 'bulk_enrichment'` 时，自动查找已完成的 DEG 任务的 CSV 结果路径：

```python
# 在 analyze 函数的 GET 分支，构建 uploaded_files 后
if module_name == 'bulk_enrichment':
    # 为 enrichment 模块提供 DEG 结果文件路径
    deg_tasks = [t for t in tasks if t.status == 'completed'
                 and t.module_name in ('bulk_deg', 'deg')]
    for t in deg_tasks:
        result_dir = os.path.join(Config.DATA_DIR, 'projects', pid, 'results')
        deg_csv = os.path.join(result_dir, 'bulk_deg_results.csv')
        if os.path.exists(deg_csv):
            # 通过 schema 的 default 值传递
            for param in schema:
                if param['key'] == 'input_source':
                    param['default'] = deg_csv
```

### Step 3: 验证完整流程

```bash
python -c "
from modules import MODULE_REGISTRY
print('Registered modules:', list(MODULE_REGISTRY.keys()))
print('Total:', len(MODULE_REGISTRY))
"
```

预期输出包含 `bulk_timecourse` 和 `bulk_enrichment`。

### Step 4: 启动应用验证

```bash
cd /data/GJ/platform
python app.py &
sleep 2
curl -s http://localhost:5000/ | grep -o '<title>[^<]*</title>'
# 预期: <title>控制台 - 生信分析平台</title>
curl -s http://localhost:5000/api/system/status
# 预期: JSON with system status
kill %1
```

### Step 5: 提交

```bash
git add routes/api.py routes/analysis.py
git commit -m "feat: add enrichment-result API and dynamic input_source for enrichment module"
```

---

## 实现顺序总结

| 顺序 | Task | 预计时间 | 依赖 |
|------|------|----------|------|
| 1 | 参数帮助系统 | 30min | 无 |
| 2 | 增强 bulk_deg.py | 45min | 无 |
| 3 | 增强 bulk_pca.py | 30min | 无 |
| 4 | 增强 bulk_heatmap.py + bulk_qc.py | 30min | 无 |
| 5 | 新增 bulk_enrichment.py | 45min | Task 2 (需要 DEG 结果) |
| 6 | 新增 bulk_timecourse.py | 60min | 无 |
| 7 | API 增强 + 最终集成 | 20min | Task 5, 6 |

**总计**: ~4 小时

Task 1-4 可并行实现，Task 5-7 需要按顺序。
