# Bulk RNA-seq QC 模块增强设计

**日期**: 2026-06-08
**状态**: 待实施
**方案**: 方案 A — 原地增强（在现有 bulk_qc.py 中增量添加）

---

## 背景与目标

当前 [modules/bulk_qc.py](modules/bulk_qc.py) 实现了基础质控：3 个阈值过滤（最小文库大小、最小基因数、最大线粒体比例）+ 4 张图（总览 2x2、相关性热图、PCA、过滤前后对比）。

用户数据为 4 组处理（ctrl, hmc3, moclel, rapa）× 3 重复的 Excel 计数矩阵，包含 count/FPKM + 基因注释（GeneSymbol, Description, GO, Pathway）。需要更全面的质控能力：

- 核糖体基因比例检测（RNA-seq QC 黄金搭档）
- 文库复杂度评估（Gini 系数）
- 基因层面过滤（低表达基因清除）
- 批次效应与样本关系可视化（组内/组间距离、PCA elbow）
- 离群样本自动检测（PCA 马氏距离）
- 过滤策略预设与决策日志

所有增强均在计数矩阵层面完成，不依赖上游比对文件（STAR/BAM）。

---

## 新增参数总览

| 参数 | type | default | 说明 |
|------|------|---------|------|
| `group_column` | dynamic_select | `` | 分组列名。填写后启用组内/组间距离、分组着色 PCA、violin 图；为空时自动从样本名推断前缀 |
| `max_ribo_pct` | number | 40.0 | 最大核糖体基因比例 (%)，高于此值的样本被过滤 |
| `min_gini` | number | 0 | 最小 Gini 系数（0 = 不过滤），低于此值提示低复杂度 |
| `min_sample_expr` | number | 0 | 基因在至少 N 个样本中 CPM > 1 才保留（0 = 不过滤） |
| `detect_outliers` | checkbox | true | 是否基于 PCA 马氏距离检测离群样本（仅告警，不自动剔除） |
| `filter_strategy` | select | `standard` | 过滤策略预设：`conservative`/`standard`/`custom` |

**filter_strategy 预设对应阈值：**

| 指标 | conservative | standard（默认） | custom |
|------|-------------|-----------------|--------|
| min_counts | 50,000 | 100,000 | 用户自填 |
| min_genes | 3,000 | 5,000 | 用户自填 |
| max_mt_pct | 30 | 20 | 用户自填 |
| max_ribo_pct | 60 | 40 | 用户自填 |

选择 `conservative` 或 `standard` 时自动覆盖对应数值参数；选择 `custom` 时用户可自由填写。

---

## 第 1 部分：新增样本级 QC 指标

### 核糖体基因比例

检测 `RPL`/`RPS` 前缀基因（核糖体蛋白基因）：

```python
ribo_mask = adata.var_names.str.startswith(('RPL', 'RPS'))
# 若有 gene_name 列，优先用 gene_name 匹配
if 'gene_name' in adata.var.columns:
    ribo_mask = adata.var['gene_name'].str.startswith(('RPL', 'RPS'))
sc.pp.calculate_qc_metrics(adata, qc_vars={'ribo': ribo_mask}, ...)
# → adata.obs['pct_counts_ribo']
```

与 MT 类似，`sc.pp.calculate_qc_metrics` 同时计算 `total_counts_ribo` 和 `pct_counts_ribo`。

### 文库复杂度 (Gini 系数)

对每个样本的非零 count 分布计算 Gini 系数：

```python
def _gini(values):
    """计算 Gini 系数（0=完全均匀，1=完全不均匀）"""
    sorted_vals = np.sort(values[values > 0])
    n = len(sorted_vals)
    if n == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * sorted_vals) / (n * np.sum(sorted_vals))) - (n + 1) / n
```

Gini > 0.8 提示文库复杂度低（PCR 过度扩增或低质量样本）。

### 基因检测新颖度 (Novelty)

```python
novelty = np.log10(n_genes) / np.log10(total_counts)
```

新颖度 < 0.8 可能提示高比例重复序列或污染。

### 汇总指标表

输出 `bulk_qc_sample_metrics.csv`：

| sample | lib_size | n_genes | mt_pct | ribo_pct | gini | novelty | group |
|--------|----------|---------|--------|----------|------|---------|-------|
| ctrl-1 | 2500000 | 15000 | 5.2 | 3.1 | 0.65 | 0.88 | ctrl |

---

## 第 2 部分：基因层面过滤

### 最低表达基因过滤

```python
if min_sample_expr > 0:
    # CPM 归一化
    lib_sizes = adata.X.sum(axis=1, keepdims=True)
    cpm = adata.X / lib_sizes * 1e6
    # 基因在至少 N 个样本中 CPM > 1
    expr_count_per_gene = (cpm > 1).sum(axis=0)
    gene_mask = expr_count_per_gene >= min_sample_expr
    # 记录被过滤基因
    removed_genes = adata.var_names[~gene_mask]
    adata = adata[:, gene_mask].copy()
```

输出 `bulk_qc_gene_filter.csv`：

| gene | expressed_in_n_samples | action |
|------|----------------------|--------|
| OR4F5 | 1 | removed |
| AL627309.5 | 0 | removed |

### 管家基因稳定性检查

自动检测常见管家基因（`GAPDH`, `ACTB`, `B2M`, `HPRT1`, `TBP`, `UBC`, `YWHAZ`, `SDHA`, `HMBS`, `RPLP0`），在 `var_names` 中查找匹配项。

对找到的管家基因，计算每个基因在样本间的 CV 系数（`std/mean`），生成热图：
- X 轴：样本
- Y 轴：管家基因
- 颜色：log10(CPM + 1)
- 标注每个基因的 CV 值

此检查不参与过滤，仅作为 QC 诊断输出。

---

## 第 3 部分：批次效应与样本关系分析

### 分组自动推断

若 `group_column` 为空，从样本名自动推断：

```python
def _infer_groups(sample_names):
    """从样本名推断分组：取分隔符前的前缀"""
    groups = []
    for name in sample_names:
        # 尝试 '-' 或 '_' 分割，取第一段
        for sep in ['-', '_']:
            if sep in str(name):
                groups.append(str(name).split(sep)[0])
                break
        else:
            groups.append(str(name))
    return groups
```

示例：`ctrl-1` → `ctrl`，`hmc3-2` → `hmc3`

### PCA 方差解释 Elbow 图

```python
from sklearn.decomposition import PCA
pca = PCA(n_components=min(10, n_samples - 1))
pca_result = pca.fit_transform(X_normalized)
# 柱状图：每个 PC 的方差解释比例
# 折线：累积方差解释比例
```

输出 `bulk_qc_pca_elbow.json`。

### 组内 vs 组间距离

```python
from scipy.spatial.distance import pdist, squareform
corr_matrix = np.corrcoef(X_normalized)  # Pearson
dist_matrix = 1 - corr_matrix  # 转为距离

# 按分组标签分组
for i, j in combinations(range(n_samples), 2):
    if groups[i] == groups[j]:
        intra_distances.append(dist_matrix[i, j])
    else:
        inter_distances.append(dist_matrix[i, j])
```

输出 `bulk_qc_group_distance.json`：箱线图展示组内距离 vs 组间距离，附 Welch t-test p-value。

### 增强相关性热图

在现有 Pearson 相关性热图基础上，Y 轴旁添加分组 annotation bar（色带），每个样本标记所属组。

### 离群样本检测

```python
from scipy.spatial.distance import mahalanobis
pca_coords = pca_result[:, :3]  # 前 3 个 PC
mean = pca_coords.mean(axis=0)
cov = np.cov(pca_coords.T)
cov_inv = np.linalg.pinv(cov)

mahal_distances = []
for i in range(n_samples):
    d = mahalanobis(pca_coords[i], mean, cov_inv)
    mahal_distances.append(d)

# MAD-based threshold
mad = np.median(np.abs(mahal_distances - np.median(mahal_distances)))
threshold = np.median(mahal_distances) + 3 * 1.4826 * mad
outliers = [samples[i] for i, d in enumerate(mahal_distances) if d > threshold]
```

离群样本在 PCA 图上标注红色，在 summary 中报告，**不自动剔除**。

---

## 第 4 部分：增强可视化

### QC 总览图扩展（3x2 子图）

从现有 2x2 扩展为 3x2：

| 位置 | 内容 | 状态 |
|------|------|------|
| (1,1) | 文库大小柱状图 | 保留 |
| (1,2) | 检测基因数柱状图 | 保留 |
| (2,1) | 线粒体比例柱状图 | 保留 |
| (2,2) | 核糖体比例柱状图 | **新增** |
| (3,1) | Gini 系数柱状图 | **新增** |
| (3,2) | 文库大小 vs 基因数散点（绿/红着色） | 保留 |

### QC 指标散点矩阵 (Pairs Plot)

`library_size × n_genes × mt_pct × ribo_pct` 四变量矩阵散点图，按分组着色。使用 Plotly `splom` (Scatterplot Matrix)。

```python
fig = go.Figure(data=go.Splom(
    dimensions=[
        dict(label='Library Size', values=obs['total_counts']),
        dict(label='N Genes', values=obs['n_genes_by_counts']),
        dict(label='MT%', values=obs['pct_counts_mt']),
        dict(label='Ribo%', values=obs['pct_counts_ribo']),
    ],
    marker=dict(color=group_colors, size=5),
))
```

### 各组 QC 指标小提琴图

若 `group_column` 已填写，为 mt_pct、ribo_pct、library_size、n_genes 各生成小提琴图（按组分面），用于检查组间系统性偏差。

### 管家基因稳定性热图

X 轴：样本，Y 轴：管家基因，颜色：log10(CPM + 1)，标注 CV 值。

### 输出文件汇总

| 文件名 | 类型 | 内容 | 条件 |
|--------|------|------|------|
| `bulk_qc_overview.json` | plotly_json | 3x2 QC 总览 | 始终生成 |
| `bulk_qc_corr.json` | plotly_json | 带分组条的相关性热图 | 始终生成 |
| `bulk_qc_pca.json` | plotly_json | PCA 散点（分组着色 + 离群标注） | 始终生成 |
| `bulk_qc_pca_elbow.json` | plotly_json | PCA 方差解释 elbow 图 | 始终生成 |
| `bulk_qc_group_distance.json` | plotly_json | 组内/组间距离箱线图 | group_column 非空 |
| `bulk_qc_pairs_plot.json` | plotly_json | QC 指标散点矩阵 | 始终生成 |
| `bulk_qc_violin_by_group.json` | plotly_json | 各组 QC 指标小提琴图 | group_column 非空 |
| `bulk_qc_filter.json` | plotly_json | 过滤前后对比柱状图 | 有样本被过滤时 |
| `bulk_qc_housekeeping.json` | plotly_json | 管家基因稳定性热图 | 找到管家基因时 |
| `bulk_qc_sample_metrics.csv` | csv | 全部样本 QC 指标表 | 始终生成 |
| `bulk_qc_filter_log.csv` | csv | 样本过滤日志（含 fail_reasons） | 始终生成 |
| `bulk_qc_gene_filter.csv` | csv | 基因过滤日志 | min_sample_expr > 0 |

---

## 第 5 部分：模块执行流程

```
输入 count matrix (Excel/CSV/H5AD)
  ↓
1. 读取数据
   - read_expression_matrix(input_path)
   - 自动推断分组（若 group_column 为空）
  ↓
2. 应用声明式过滤（apply_filters）
  ↓
3. 计算样本级 QC 指标
   - MT%, Ribo%, 文库大小, 基因数
   - Gini 系数, 新颖度
  ↓
4. 按 filter_strategy 加载阈值
   - conservative/standard 自动覆盖数值参数
   - custom 使用用户填写值
  ↓
5. 过滤样本
   - 应用所有阈值（min_counts, min_genes, max_mt_pct, max_ribo_pct）
   - 生成 bulk_qc_filter_log.csv
  ↓
6. 过滤基因
   - min_sample_expr 阈值
   - 生成 bulk_qc_gene_filter.csv
  ↓
7. 管家基因稳定性检查（不参与过滤）
  ↓
8. 生成所有图表
   - 3x2 QC 总览
   - 相关性热图（带分组条）
   - PCA + elbow 图
   - 散点矩阵 (pairs plot)
   - 各组小提琴图（条件）
   - 组内/组间距离（条件）
   - 管家基因热图（条件）
  ↓
9. 离群检测（PCA 马氏距离，仅告警）
  ↓
10. 保存 filtered adata + summary
```

---

## 第 6 部分：routes/analysis.py 参数更新

在 `PARAM_SCHEMAS['bulk_qc']` 末尾追加：

```python
        {'key': 'group_column', 'label': '分组列名（可选）', 'type': 'dynamic_select', 'default': '',
         'help': '样本分组列名。填写后启用组内/组间距离分析、分组着色 PCA 和小提琴图。留空则自动从样本名推断。'},
        {'key': 'max_ribo_pct', 'label': '最大核糖体基因比例 (%)', 'type': 'number', 'default': 40.0, 'step': 0.1,
         'help': '最大核糖体基因比例（%）。RPL/RPS 基因比例过高提示 rRNA 污染。PolyA 建库通常 < 5-10%，rRNA 去除建库可至 40-50%。'},
        {'key': 'min_gini', 'label': '最小文库复杂度 (Gini)', 'type': 'number', 'default': 0, 'step': 0.01,
         'help': '最小 Gini 系数（0 = 不过滤）。Gini > 0.8 提示文库复杂度低。'},
        {'key': 'min_sample_expr', 'label': '基因最低表达样本数', 'type': 'number', 'default': 0, 'step': 1,
         'help': '基因在至少 N 个样本中 CPM > 1 才保留。0 = 不过滤。建议设为样本总数的 10-20%。'},
        {'key': 'detect_outliers', 'label': '检测离群样本', 'type': 'checkbox', 'default': True,
         'help': '基于 PCA 马氏距离检测离群样本。仅告警，不自动剔除。'},
        {'key': 'filter_strategy', 'label': '过滤策略', 'type': 'select', 'default': 'standard',
         'options': ['conservative', 'standard', 'custom'],
         'help': 'conservative：宽松阈值（适合小样本）；standard：推荐阈值；custom：自定义所有阈值。'},
```

---

## 向后兼容

| 场景 | 行为 |
|------|------|
| 所有新参数留空/默认 | 行为与当前基本一致（多输出几张诊断图，不额外过滤） |
| `max_ribo_pct=40` | 典型数据 Ribo% < 10%，不影响 |
| `min_sample_expr=0` | 不过滤基因 |
| `group_column` 为空 | 自动推断分组，不生成分组特异性图 |
| `detect_outliers=true` | 仅在 summary 中报告，不剔除 |
| `filter_strategy=standard` | 使用当前默认阈值，无行为变化 |

---

## 错误处理

| 场景 | 处理方式 |
|------|----------|
| 无核糖体基因可检测 | Ribo% 全部记为 0，在 summary 中提示 |
| PCA 马氏距离计算失败（样本数 < 4） | 跳过离群检测，在 summary 中提示 |
| 自动推断分组失败（样本名无分隔符） | 所有样本归为一组，不生成分组相关图表 |
| 管家基因全部不在 var_names 中 | 跳过管家基因检查，不输出 housekeeping 图 |
| 基因过滤后基因数为 0 | 抛出 ValueError，提示 min_sample_expr 设置过高 |

---

## 改动文件汇总

| 文件 | 改动类型 | 说明 |
|------|----------|------|
| `modules/bulk_qc.py` | 大幅修改 | 新增 ~300 行：核糖体/Gini/novelty 指标、基因过滤、离群检测、分组分析、增强图表 |
| `routes/analysis.py` | 修改 | PARAM_SCHEMAS['bulk_qc'] 新增 6 个参数 |

---

## 不做的事情（YAGNI）

- 不支持上传 STAR/BAM 文件（仅计数矩阵层面）
- 不添加基因类型过滤（Description 列不够可靠判断 protein_coding vs lncRNA）
- 不添加动态阈值实时预览滑块（需要前端 JS 重构，超出本轮范围）
- 不添加 PDF/HTML 报告导出（需要额外依赖，后续迭代）
- 不自动剔除离群样本（仅告警）
- 不添加插入片段大小、5'/3' 覆盖偏差（需要上游比对数据）

---

## 实施状态

| 任务 | 内容 | 状态 |
|------|------|------|
| 1 | bulk_qc.py：样本级新指标（Ribo%、Gini、novelty） | 待实施 |
| 2 | bulk_qc.py：基因过滤 + 管家基因检查 | 待实施 |
| 3 | bulk_qc.py：分组推断 + 离群检测 + 组距离分析 | 待实施 |
| 4 | bulk_qc.py：增强图表（3x2 总览、pairs plot、elbow、violin、housekeeping） | 待实施 |
| 5 | routes/analysis.py：PARAM_SCHEMAS 更新 | 待实施 |

---

## 验证命令

```bash
# 验证模块可导入
python -c "from modules.bulk_qc import BulkQCAnalysis; print('OK')"

# 验证辅助函数
python -c "
from modules.bulk_qc import _gini, _infer_groups
import numpy as np
assert _gini(np.array([1, 1, 1])) == 0.0  # 完全均匀
assert _gini(np.array([0, 0, 100])) > 0.9  # 极不均匀
assert _infer_groups(['ctrl-1', 'ctrl-2', 'hmc3-1']) == ['ctrl', 'ctrl', 'hmc3']
print('Helper functions OK')
"

# 验证 PARAM_SCHEMAS 完整性
python -c "
from routes.analysis import PARAM_SCHEMAS
keys = [p['key'] for p in PARAM_SCHEMAS['bulk_qc']]
for k in ['max_ribo_pct', 'min_sample_expr', 'group_column', 'detect_outliers', 'filter_strategy']:
    assert k in keys, f'{k} missing'
print('PARAM_SCHEMAS OK')
"
```
