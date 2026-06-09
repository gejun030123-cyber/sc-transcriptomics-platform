# Bulk 热图增强设计文档 — 第一阶段

## Context

当前 `bulk_heatmap.py` 功能有限：仅支持 top_var/deg 两种基因选择、硬编码 z-score + Euclidean/Ward 聚类、单一注释条、不可配置的颜色方案。用户需要更灵活的热图配置，覆盖基因选择策略、数据变换、聚类方法、视觉样式四大类，并支持从表达式筛选结果导入基因。

## 架构：方案 B — 提取工具函数 + 增强模块

- 将可复用逻辑抽取到 `visualization.py`
- `bulk_heatmap.py` 调用工具函数，保持 run() 简洁
- 新参数添加到 `routes/analysis.py` 的 PARAM_SCHEMAS

---

## 一、基因选择策略扩展

### 新增参数

| key | type | default | 说明 |
|-----|------|---------|------|
| `var_metric` | select | `var` | 变异度量：var / mad / cv / range |
| `deg_source` | text | `''` | DEG 来源文件名（如 `results_0`），留空=自动选最新 |
| `deg_direction` | select | `both` | 方向过滤：both / up / down |
| `deg_sortby` | select | `padj` | 排序方式：padj / abs_logfc / logfc / cluster |
| `gene_import_source` | select | `manual` | 基因来源：manual / top_var / deg / expression_filter |
| `filter_expression` | textarea | `''` | 表达式筛选（如 `ALL:up`），仅 gene_import_source=expression_filter 时生效 |

### 实现逻辑

**var_metric 变异度量**：
```python
# visualization.py
def compute_gene_variability(data, metric='var'):
    if metric == 'var':     return np.var(data, axis=0)
    elif metric == 'mad':   return np.median(np.abs(data - np.median(data, axis=0)), axis=0)
    elif metric == 'cv':    return np.std(data, axis=0) / (np.mean(data, axis=0) + 1e-10)
    elif metric == 'range': return np.max(data, axis=0) - np.min(data, axis=0)
```

**deg_source / deg_direction / deg_sortby**：
- 扫描 `results/` 目录下的 `bulk_deg_results*.csv`
- 按 deg_direction 过滤 regulation 列
- 按 deg_sortby 排序后取 top_n

**gene_import_source = expression_filter**：
- 调用 `expression_parser.validate()` + `evaluate()` 解析表达式
- 从 DEG 矩阵中提取匹配基因
- 需要先加载所有 DEG 结果构建 padj/logfc 矩阵（复用 bulk_deg_integration 的逻辑）

---

## 二、数据变换与标准化

### 新增参数

| key | type | default | 说明 |
|-----|------|---------|------|
| `row_scaling` | select | `zscore` | 行标准化：none / zscore / center |
| `pseudocount` | number | `1` | log2 转换前加值（0.1 / 1 / 自定义） |
| `winsorize` | select | `none` | 极端值截断：none / 1pct / 5pct / custom |
| `winsorize_custom` | number | `''` | 自定义截断百分位（winsorize=custom 时生效） |
| `missing_value` | select | `ignore` | 缺失值处理：ignore / mean_fill / zero_fill |
| `clip_range` | text | `-3,3` | z-score 截断范围（逗号分隔的 min,max） |

### 实现逻辑

```python
# visualization.py
def transform_heatmap_data(data, row_scaling='zscore', pseudocount=1,
                           winsorize='none', clip_range=(-3, 3), missing_value='ignore'):
    """对热图数据进行标准化和变换。data: (samples, genes)"""
    # 1. 缺失值处理
    if missing_value == 'mean_fill':
        col_means = np.nanmean(data, axis=0)
        for j in range(data.shape[1]):
            mask = np.isnan(data[:, j])
            data[mask, j] = col_means[j]
    elif missing_value == 'zero_fill':
        data = np.nan_to_num(data, nan=0.0)

    # 2. pseudocount + log2（仅在原始计数时）
    # 由调用方决定是否需要 log 转换

    # 3. Winsorize
    if winsorize != 'none':
        pct = float(winsorize.replace('pct', '')) / 100
        from scipy.stats.mstats import winsorize as sp_winsorize
        data = sp_winsorize(data, limits=[pct, pct], axis=0)

    # 4. 行标准化
    if row_scaling == 'zscore':
        mean = data.mean(axis=0)
        std = data.std(axis=0) + 1e-10
        data = (data - mean) / std
    elif row_scaling == 'center':
        data = data - data.mean(axis=0)

    # 5. 截断
    if clip_range:
        data = np.clip(data, clip_range[0], clip_range[1])

    return data
```

---

## 三、聚类方法精细化

### 新增参数

| key | type | default | 说明 |
|-----|------|---------|------|
| `row_cluster` | select | `yes` | 行聚类（基因）：yes / no / preset |
| `col_cluster` | select | `yes` | 列聚类（样本）：yes / no / group_order |
| `row_method` | select | `ward` | 行聚类方法：ward / complete / average / single / mcquitty |
| `col_method` | select | `ward` | 列聚类方法（同上） |
| `row_metric` | select | `euclidean` | 行距离度量：euclidean / pearson / spearman / cosine |
| `col_metric` | select | `euclidean` | 列距离度量（同上） |
| `row_order_by` | select | `preset` | 不聚类时行排序：input / logfc / padj |
| `col_order_by` | select | `group` | 不聚类时列排序：input / group / batch |

### 实现逻辑

```python
# visualization.py
def cluster_heatmap(data, method='ward', metric='euclidean'):
    """层次聚类，返回排序后的行索引。"""
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist

    if metric == 'pearson':
        dist = pdist(data, metric=lambda u, v: 1 - np.corrcoef(u, v)[0, 1])
    elif metric == 'spearman':
        from scipy.stats import spearmanr
        dist = pdist(data, metric=lambda u, v: 1 - spearmanr(u, v)[0])
    elif metric == 'cosine':
        dist = pdist(data, metric='cosine')
    else:
        dist = pdist(data, metric='euclidean')

    link = linkage(dist, method=method)
    dendro = dendrogram(link, no_plot=True)
    return dendro['leaves']
```

**row_order_by / col_order_by**：
- `input`：保持原始顺序
- `logfc`：按 DEG 结果的 log2FC 排序
- `padj`：按 padj 排序
- `group`：按分组列排序（同组样本相邻）
- `batch`：按批次列排序

---

## 四、颜色与视觉样式

### 新增参数

| key | type | default | 说明 |
|-----|------|---------|------|
| `colorscale` | select | `RdBu_r` | 颜色方案：RdBu_r / RdYlBu / viridis / magma / Blues / custom |
| `reverse_color` | checkbox | `False` | 反转颜色 |
| `zmin` | text | `auto` | 最小值（auto=自动） |
| `zmax` | text | `auto` | 最大值（auto=自动） |
| `show_gene_labels` | select | `all` | 基因名显示：all / top_n / none |
| `show_sample_labels` | select | `all` | 样本名显示：all / group_only / none |
| `gene_font_size` | number | `8` | 基因名字体大小 |
| `sample_font_size` | number | `9` | 样本名字体大小 |
| `cell_border` | select | `none` | 单元格边框：none / thin / group_sep |
| `annotation_columns` | text | `''` | 多注释条列名（逗号分隔，如 `group,batch`） |

### 实现逻辑

colorscale 直接传给 `go.Heatmap(colorscale=...)`，Plotly 支持所有内置色图名。

`reverse_color`：`colorscale = colorscale + '_r'` 或反转自定义色图。

`zmin/zmax`：解析 `auto` → 不设置（Plotly 自动），否则转 float。

`show_gene_labels`：
- `all`：显示全部
- `top_n`：只显示前 N 个标签（通过 `ticktext` + `tickvals` 控制）
- `none`：`showticklabels=False`

`annotation_columns`：逗号分隔的列名列表，每列生成一个独立的注释条 figure。

---

## 五、visualization.py 新增工具函数清单

| 函数 | 作用 |
|------|------|
| `compute_gene_variability(data, metric)` | 计算基因变异度量 |
| `transform_heatmap_data(data, ...)` | 数据标准化/变换 |
| `cluster_heatmap(data, method, metric)` | 层次聚类返回排序索引 |
| `build_annotation_bar(obs, columns, palette)` | 多注释条生成 |
| `save_plotly_json(fig, plots_dir, filename, result_files, ...)` | 统一保存 |

---

## 六、PARAM_SCHEMAS 完整方案

```python
'bulk_heatmap': [
    # 基因选择
    {'key': 'gene_import_source', 'type': 'select',
     'options': ['manual', 'top_var', 'deg', 'expression_filter'],
     'default': 'top_var', ...},
    {'key': 'var_metric', 'type': 'select',
     'options': ['var', 'mad', 'cv', 'range'], 'default': 'var', ...},
    {'key': 'deg_source', 'type': 'text', 'default': '', ...},
    {'key': 'deg_direction', 'type': 'select',
     'options': ['both', 'up', 'down'], 'default': 'both', ...},
    {'key': 'deg_sortby', 'type': 'select',
     'options': ['padj', 'abs_logfc', 'logfc', 'cluster'], 'default': 'padj', ...},
    {'key': 'filter_expression', 'type': 'textarea', 'default': '', ...},
    {'key': 'top_n', 'type': 'number', 'default': 50, ...},
    {'key': 'custom_genes', 'type': 'textarea', 'default': '', ...},
    # 数据变换
    {'key': 'row_scaling', 'type': 'select',
     'options': ['zscore', 'center', 'none'], 'default': 'zscore', ...},
    {'key': 'pseudocount', 'type': 'number', 'default': 1, ...},
    {'key': 'winsorize', 'type': 'select',
     'options': ['none', '1pct', '5pct', 'custom'], 'default': 'none', ...},
    {'key': 'clip_range', 'type': 'text', 'default': '-3,3', ...},
    {'key': 'missing_value', 'type': 'select',
     'options': ['ignore', 'mean_fill', 'zero_fill'], 'default': 'ignore', ...},
    # 聚类
    {'key': 'row_cluster', 'type': 'select',
     'options': ['yes', 'no', 'preset'], 'default': 'yes', ...},
    {'key': 'col_cluster', 'type': 'select',
     'options': ['yes', 'no', 'group_order'], 'default': 'yes', ...},
    {'key': 'row_method', 'type': 'select',
     'options': ['ward', 'complete', 'average', 'single', 'mcquitty'], 'default': 'ward', ...},
    {'key': 'col_method', 'type': 'select',
     'options': ['ward', 'complete', 'average', 'single', 'mcquitty'], 'default': 'ward', ...},
    {'key': 'row_metric', 'type': 'select',
     'options': ['euclidean', 'pearson', 'spearman', 'cosine'], 'default': 'euclidean', ...},
    {'key': 'col_metric', 'type': 'select',
     'options': ['euclidean', 'pearson', 'spearman', 'cosine'], 'default': 'euclidean', ...},
    {'key': 'row_order_by', 'type': 'select',
     'options': ['input', 'logfc', 'padj'], 'default': 'input', ...},
    {'key': 'col_order_by', 'type': 'select',
     'options': ['input', 'group', 'batch'], 'default': 'input', ...},
    # 视觉样式
    {'key': 'colorscale', 'type': 'select',
     'options': ['RdBu_r', 'RdYlBu', 'viridis', 'magma', 'Blues', 'custom'],
     'default': 'RdBu_r', ...},
    {'key': 'reverse_color', 'type': 'checkbox', 'default': False, ...},
    {'key': 'zmin', 'type': 'text', 'default': 'auto', ...},
    {'key': 'zmax', 'type': 'text', 'default': 'auto', ...},
    {'key': 'show_gene_labels', 'type': 'select',
     'options': ['all', 'top_n', 'none'], 'default': 'all', ...},
    {'key': 'show_sample_labels', 'type': 'select',
     'options': ['all', 'group_only', 'none'], 'default': 'all', ...},
    {'key': 'gene_font_size', 'type': 'number', 'default': 8, ...},
    {'key': 'sample_font_size', 'type': 'number', 'default': 9, ...},
    {'key': 'annotation_columns', 'type': 'text', 'default': '', ...},
    {'key': 'groupby', 'type': 'text', 'default': '', ...},
]
```

---

## 七、改动文件清单

| 文件 | 改动 |
|------|------|
| `modules/visualization.py` | 新增 5 个工具函数 |
| `modules/bulk_heatmap.py` | 重构 run() 使用工具函数，添加新参数处理逻辑 |
| `routes/analysis.py` | 替换 bulk_heatmap PARAM_SCHEMAS（4→30+ 参数） |
| `tests/test_heatmap_helpers.py` | 新建，测试工具函数 |

## 八、向后兼容

- 旧参数 `heatmap_type` 保留作为 `gene_import_source` 的别名：`top_var` → `gene_import_source=top_var`，`deg` → `gene_import_source=deg`
- 旧参数 `custom_genes` 保留，等价于 `gene_import_source=manual`
- 旧预设文件中如含 `heatmap_type`，run() 中做映射兼容

## 九、验证方式

1. `python -m pytest tests/test_heatmap_helpers.py -v` — 工具函数测试通过
2. 启动 `python app.py`，进入 bulk_heatmap 配置页验证：
   - 基因选择：top_var / deg / expression_filter 各模式正常
   - 数据变换：z-score / center / none 切换生效
   - 聚类：不同方法/距离度量切换正常
   - 颜色：不同 colorscale 切换正常
3. 提交运行，检查结果页热图和注释条正确显示
