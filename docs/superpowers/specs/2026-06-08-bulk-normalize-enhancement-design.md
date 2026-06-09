# Bulk RNA-seq 标准化模块增强设计

**日期**: 2026-06-08
**状态**: 待实施
**方案**: 方案 A — 原地增强 + 修复下游重复标准化

---

## 背景与目标

当前 [modules/bulk_normalize.py](modules/bulk_normalize.py)（133 行）提供 3 种标准化方法（DESeq2/CPM/log2_quantile），1 个参数（method）。下游模块（bulk_heatmap、bulk_pca）各自独立执行标准化，导致用户串联 normalize→heatmap 时存在重复标准化风险。

本次增强目标：
- 新增 TMM、VST、rlog 三种标准化方法
- 新增低表达基因前置过滤
- 新增标准化前后诊断图（4 张）
- 修复下游模块重复标准化问题

---

## 新增参数总览

| 参数 | type | default | 说明 |
|------|------|---------|------|
| `method` | select | `deseq2` | 标准化方法（新增 `tmm`/`vst`/`rlog` 选项） |
| `min_expr_value` | number | 1 | 最小表达阈值（CPM），前置过滤用 |
| `min_expr_samples` | number | 3 | 基因在至少 N 个样本中达到阈值才保留（0=不过滤） |
| `max_zero_pct` | number | 0 | 最大零值比例 (%)，0=不过滤 |

---

## 第 1 部分：新增标准化方法

### 方法总览

| 方法 | key | 算法 | 输出 adata.X | 适用场景 |
|------|-----|------|-------------|---------|
| DESeq2 (现有) | `deseq2` | 中位比率法 size factor | log2(SF_norm + 1) | 差异分析 |
| CPM (现有) | `cpm` | 文库大小归一 | log2(CPM + 1) | 快速预览 |
| log2 分位数 (现有) | `log2_quantile` | log2 后分位数归一 | quantile_norm(log2) | 样本间可比 |
| **TMM (新增)** | `tmm` | edgeR trimmed mean of M-values | log2(TMM_norm + 1) | 差异分析（组成偏差大） |
| **VST (新增)** | `vst` | 方差稳定变换 | VST 值 | PCA/聚类/热图 |
| **rlog (新增)** | `rlog` | 正则化对数 | rlog 值 | 小样本 (<30) 高维分析 |

### TMM 实现

```python
def _tmm_normalize(counts):
    """TMM 标准化（edgeR 风格）。counts: (n_samples, n_genes) ndarray."""
    n_samples, n_genes = counts.shape
    lib_sizes = counts.sum(axis=1)

    # 选择参考样本（文库大小最接近中位数的样本）
    median_lib = np.median(lib_sizes)
    ref_idx = np.argmin(np.abs(lib_sizes - median_lib))

    factors = np.ones(n_samples)
    for i in range(n_samples):
        if i == ref_idx:
            continue
        # M = log2 ratio, A = average log expression
        mask = (counts[i] > 0) & (counts[ref_idx] > 0)
        Mi = np.log2((counts[i, mask] / lib_sizes[i]) / (counts[ref_idx, mask] / lib_sizes[ref_idx]))
        Ai = 0.5 * (np.log2(counts[i, mask] / lib_sizes[i]) + np.log2(counts[ref_idx, mask] / lib_sizes[ref_idx]))

        # 去除 top/bottom 30% M 和 top 5% A
        m_lo, m_hi = np.percentile(Mi, [30, 70])
        a_hi = np.percentile(Ai, 95)
        keep = (Mi >= m_lo) & (Mi <= m_hi) & (Ai <= a_hi)

        if keep.sum() > 0:
            # 加权均值（权重基于近似方差）
            weights = 1 / ((lib_sizes[i] - counts[i, mask]) / (lib_sizes[i] * counts[i, mask]) +
                           (lib_sizes[ref_idx] - counts[ref_idx, mask]) / (lib_sizes[ref_idx] * counts[ref_idx, mask]))
            factors[i] = 2 ** np.average(Mi[keep], weights=weights[keep])

    return factors
```

### VST 实现

```python
def _vst_transform(counts, size_factors):
    """近似 VST。counts: raw counts, size_factors: per-sample."""
    # 1. 归一化
    normed = counts / size_factors[:, None]
    # 2. 拟合 log(mean) ~ log(var) 关系
    gene_means = normed.mean(axis=0)
    gene_vars = normed.var(axis=0)
    valid = gene_means > 0
    log_mean = np.log2(gene_means[valid] + 1)
    log_var = np.log2(gene_vars[valid] + 1)
    # 线性回归 log_var = a + b * log_mean
    coeffs = np.polyfit(log_mean, log_var, 1)
    # 3. 积分变换：h(x) = integral_0^x 1/sqrt(f(t)) dt
    # 简化近似：VST ≈ log2( (normed + alpha) / (1 + alpha) ) * beta
    # 使用 log2(normed + 0.5) 作为快速近似
    vst = np.log2(normed + 0.5)
    return vst
```

### rlog 实现

与 VST 类似但对先验做正则化。对于样本数 < 30 的情况，正则化避免过度拟合：

```python
def _rlog_transform(counts, size_factors, prior_mean=0, prior_var=1):
    """近似 rlog。"""
    normed = counts / size_factors[:, None]
    gene_means = normed.mean(axis=0)
    # 正则化：收缩基因特异性 log 向全局先验
    log_normed = np.log2(normed + 0.5)
    gene_effects = log_normed.mean(axis=0) - prior_mean
    # 加权平均：样本越多，越信赖数据；样本越少，越信赖先验
    shrinkage = min(1.0, counts.shape[0] / 30.0)
    rlog = prior_mean + shrinkage * gene_effects[None, :] + (1 - shrinkage) * (log_normed - log_normed.mean(axis=0, keepdims=True))
    return rlog
```

### 统一输出约定

所有方法均写入：

```python
adata.layers['raw'] = 原始 counts 副本
adata.layers['normalized'] = 线性尺度标准化值（deseq2/tmm/cpm 输出）
adata.X = 最终分析用矩阵（log 转换或 VST/rlog）
adata.uns['normalization'] = {
    'method': method,
    'is_log_transformed': True,  # adata.X 是否已 log 转换
}
```

VST 和 rlog 方法不设置 `adata.layers['normalized']`（VST/rlog 直接输出到 adata.X）。

---

## 第 2 部分：低表达基因前置过滤

在标准化方法执行之前，先进行基因过滤。

### 过滤逻辑

```python
# 在标准化之前
if min_expr_samples > 0:
    lib_sizes = adata.X.sum(axis=1, keepdims=True)
    lib_sizes[lib_sizes == 0] = 1
    cpm = adata.X / lib_sizes * 1e6
    n_expr = np.array((cpm >= min_expr_value).sum(axis=0)).flatten()
    gene_mask = n_expr >= min_expr_samples
    if max_zero_pct > 0:
        zero_pct = np.array((adata.X == 0).sum(axis=0)).flatten() / adata.n_obs * 100
        gene_mask = gene_mask & (zero_pct <= max_zero_pct)
    adata = adata[:, gene_mask].copy()
```

### 阈值说明

| 参数 | 默认 | 说明 |
|------|------|------|
| min_expr_value | 1 | CPM 阈值。CPM > 1 意味着每百万 reads 中至少 1 条来自该基因 |
| min_expr_samples | 3 | 至少在 3 个样本中有效表达（典型最小组大小） |
| max_zero_pct | 0 | 默认不过滤。设为 50 表示允许最多 50% 样本中该基因为零 |

---

## 第 3 部分：标准化前后诊断图

| 图表 | 文件名 | 内容 | 条件 |
|------|--------|------|------|
| 文库大小对比 | `normalize_libsize.json` | 原始 vs 标准化后各列总和 | 始终 |
| Size Factor 分布 | `normalize_size_factors.json` | 各样本 size factor 柱状图 | deseq2/tmm |
| 箱线图前后对比 | `normalize_boxplot_compare.json` | 所有样本 log2 表达分布并排 | 始终 |
| PCA 前后对比 | `normalize_pca_compare.json` | 原始 vs 标准化后 PCA 左右并排 | 始终 |

### 箱线图对比实现

```python
fig = make_subplots(rows=1, cols=2, subplot_titles=['标准化前', '标准化后'])
# 左：原始 log2(counts + 1) 分布
for i in range(n_samples):
    fig.add_trace(go.Box(y=np.log2(raw_counts[i] + 1), name=sample_names[i], showlegend=False), row=1, col=1)
# 右：标准化后 adata.X 分布
for i in range(n_samples):
    fig.add_trace(go.Box(y=adata.X[i], name=sample_names[i], showlegend=False), row=1, col=2)
```

### PCA 对比实现

```python
fig = make_subplots(rows=1, cols=2, subplot_titles=['原始数据 PCA', '标准化后 PCA'])
# 左：原始数据
raw_norm = normalize_for_pca(raw_counts)  # CPM + log1p + scale
sc.pp.pca(raw_norm, n_comps=2)
# 右：标准化后
sc.pp.pca(adata_normalized, n_comps=2)
# 各自画散点...
```

---

## 第 4 部分：修复下游模块重复标准化

### 约定

normalize 模块在 `adata.uns['normalization']` 中写入标记：

```python
adata.uns['normalization'] = {
    'method': 'deseq2',
    'is_log_transformed': True,
}
```

### 下游模块检测逻辑

在 `bulk_heatmap.py` 和 `bulk_pca.py` 中，标准化代码块之前添加：

```python
if 'normalization' not in adata.uns:
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
```

### 需要修改的文件

| 文件 | 修改 |
|------|------|
| `modules/bulk_heatmap.py` | 添加 already_normalized 条件检测（~2 行） |
| `modules/bulk_pca.py` | 同上（~2 行） |

### 不需要修改的

- `bulk_deg.py`：用原始 counts 做差异分析
- `bulk_enrichment.py`：读基因列表，不涉及表达矩阵
- `bulk_timecourse.py`：独立读取原始数据
- `bulk_qc.py`：QC 从原始数据开始，不需修改

---

## 第 5 部分：Summary 与向后兼容

### Summary

```python
'summary': {
    'method': method,
    'n_samples': n_samples,
    'n_genes_before_filter': n_genes_original,
    'n_genes_after_filter': n_genes_filtered,
    'genes_filtered': n_genes_original - n_genes_filtered,
    'median_size_factor': float or None,
    'is_log_transformed': True/False,
}
```

### 向后兼容

| 场景 | 行为 |
|------|------|
| `method=deseq2`（默认） | 行为与当前一致 |
| 前置过滤参数全为 0/默认 | 不过滤，与当前一致 |
| 下游模块接收无标记旧数据 | `normalization` not in uns → 正常执行自己的标准化 |

### 参数选项更新

`method` 的 options 从 `['deseq2', 'cpm', 'log2_quantile']` 扩展为：
`['deseq2', 'tmm', 'cpm', 'vst', 'rlog', 'log2_quantile']`

---

## 改动文件汇总

| 文件 | 操作 | 说明 |
|------|------|------|
| `modules/bulk_normalize.py` | 大幅修改 | 新增 TMM/VST/rlog + 前置过滤 + 诊断图 + 输出标记 |
| `routes/analysis.py` | 修改 | PARAM_SCHEMAS 扩展（4 个参数） |
| `modules/bulk_heatmap.py` | 小改 | 添加 already_normalized 检测（~2 行） |
| `modules/bulk_pca.py` | 小改 | 添加 already_normalized 检测（~2 行） |

---

## 不做的事情（YAGNI）

- 不添加批次校正（ComBat-seq/RUVseq）→ 下一轮
- 不添加 TPM/FPKM（需要基因长度信息）
- 不添加 sctransform（单细胞方法）
- 不添加交互式实时预览（需前端 JS 重构）
- 不添加自动推荐方法（需下游分析意图信息）
- 不添加 RLE 图/mean-variance 图（后续迭代）

---

## 实施状态

| 任务 | 内容 | 状态 |
|------|------|------|
| 1 | bulk_normalize.py：TMM 方法实现 | 待实施 |
| 2 | bulk_normalize.py：VST/rlog 方法实现 | 待实施 |
| 3 | bulk_normalize.py：低表达基因前置过滤 | 待实施 |
| 4 | bulk_normalize.py：诊断图（4 张） | 待实施 |
| 5 | bulk_normalize.py：统一输出约定 + summary | 待实施 |
| 6 | routes/analysis.py：PARAM_SCHEMAS 更新 | 待实施 |
| 7 | bulk_heatmap.py + bulk_pca.py：重复标准化修复 | 待实施 |
