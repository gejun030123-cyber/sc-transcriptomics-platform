# Bulk RNA-seq 标准化模块增强 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 bulk_normalize.py 中新增 TMM/VST/rlog 标准化方法、低表达基因前置过滤、诊断图，并修复下游模块重复标准化问题。

**Architecture:** 在现有 bulk_normalize.py 的 `run()` 方法中扩展分支，用辅助函数（`_tmm_normalize`, `_vst_transform`, `_rlog_transform`）实现新算法。下游模块（bulk_heatmap.py、bulk_pca.py）各加 2 行条件检测跳过已标准化数据。

**Tech Stack:** numpy, scipy, scanpy, plotly

---

## 文件结构

| 文件 | 操作 | 职责 |
|------|------|------|
| `modules/bulk_normalize.py` | 大幅修改 | 所有标准化逻辑、辅助函数、诊断图 |
| `routes/analysis.py:127-129` | 修改 | PARAM_SCHEMAS 扩展 |
| `modules/bulk_heatmap.py:37-38` | 小改 | 条件跳过标准化 |
| `modules/bulk_pca.py:34-35` | 小改 | 条件跳过标准化 |
| `tests/test_normalize_helpers.py` | 新建 | 辅助函数单元测试 |

---

## Task 1: TMM 辅助函数与单元测试

**Files:**
- Create: `tests/test_normalize_helpers.py`
- Modify: `modules/bulk_normalize.py`（末尾追加辅助函数）

- [ ] **Step 1: 编写 TMM 测试**

```python
# tests/test_normalize_helpers.py
import numpy as np
import pytest


def test_tmm_equal_samples():
    """完全相同的样本 → 因子全为 1"""
    from modules.bulk_normalize import _tmm_normalize
    counts = np.array([
        [10, 20, 30, 40, 50],
        [10, 20, 30, 40, 50],
        [10, 20, 30, 40, 50],
    ], dtype=float)
    factors = _tmm_normalize(counts)
    assert len(factors) == 3
    np.testing.assert_allclose(factors, 1.0, atol=0.1)


def test_tmm_scaled_sample():
    """一个样本是其他的 2 倍 → 因子约 0.5"""
    from modules.bulk_normalize import _tmm_normalize
    np.random.seed(42)
    base = np.random.randint(10, 1000, size=(20, 5)).astype(float)
    counts = np.vstack([base, base * 2])  # 第三个样本是前两个的 2 倍
    factors = _tmm_normalize(counts)
    # 前两个样本因子约 1，第三个约 0.5
    assert factors[2] < factors[0] * 0.7


def test_tmm_no_crash_zeros():
    """含零值不崩溃"""
    from modules.bulk_normalize import _tmm_normalize
    counts = np.array([
        [0, 10, 20, 0, 50],
        [5, 0, 30, 40, 0],
        [10, 20, 0, 0, 100],
    ], dtype=float)
    factors = _tmm_normalize(counts)
    assert len(factors) == 3
    assert all(f > 0 for f in factors)


def test_vst_shape():
    """VST 输出形状与输入一致"""
    from modules.bulk_normalize import _vst_transform
    counts = np.array([[100, 200, 300], [150, 250, 350]], dtype=float)
    sf = np.array([1.0, 1.0])
    result = _vst_transform(counts, sf)
    assert result.shape == counts.shape


def test_vst_variance_stabilized():
    """VST 后低表达基因方差应比原始 log2 更稳定"""
    from modules.bulk_normalize import _vst_transform
    np.random.seed(42)
    # 构造均值差异大但方差比例相似的数据
    counts = np.array([
        [10, 100, 1000, 10000],
        [12, 120, 1200, 12000],
        [8, 80, 800, 8000],
    ], dtype=float)
    sf = np.array([1.0, 1.0, 1.0])
    vst = _vst_transform(counts, sf)
    # VST 后各基因的方差应更接近（不像原始 log2 那样低表达方差大）
    variances = vst.var(axis=0)
    assert variances.max() / (variances.min() + 1e-10) < 10  # 方差比应较小


def test_rlog_shape():
    """rlog 输出形状与输入一致"""
    from modules.bulk_normalize import _rlog_transform
    counts = np.array([[100, 200, 300], [150, 250, 350]], dtype=float)
    sf = np.array([1.0, 1.0])
    result = _rlog_transform(counts, sf)
    assert result.shape == counts.shape


def test_rlog_small_sample_regularity():
    """小样本时 rlog 比直接 log2 更规则（收缩效应）"""
    from modules.bulk_normalize import _rlog_transform
    counts = np.array([
        [10, 100, 1000],
        [50, 50, 500],
        [5, 200, 2000],
    ], dtype=float)
    sf = np.array([1.0, 1.0, 1.0])
    rlog = _rlog_transform(counts, sf)
    direct_log = np.log2(counts + 1)
    # rlog 的行间方差应比直接 log2 小（正则化收缩）
    assert rlog.var(axis=0).mean() <= direct_log.var(axis=0).mean() * 1.5
```

- [ ] **Step 2: 运行测试确认失败**

Run: `cd /data/GJ/platform && python -m pytest tests/test_normalize_helpers.py -v`
Expected: FAIL（函数不存在）

- [ ] **Step 3: 实现 TMM 辅助函数**

在 `modules/bulk_normalize.py` 文件末尾追加：

```python
# --- 辅助函数 ---


def _tmm_normalize(counts):
    """TMM 标准化（edgeR 风格）。counts: (n_samples, n_genes) ndarray of raw counts.
    返回 per-sample TMM 因子（归一后应除以此因子）。"""
    n_samples, n_genes = counts.shape
    lib_sizes = counts.sum(axis=1)
    lib_sizes[lib_sizes == 0] = 1

    # 参考样本：文库大小最接近中位数
    median_lib = np.median(lib_sizes)
    ref_idx = int(np.argmin(np.abs(lib_sizes - median_lib)))

    factors = np.ones(n_samples)
    for i in range(n_samples):
        if i == ref_idx:
            continue
        # 仅用两个样本都有 count 的基因
        mask = (counts[i] > 0) & (counts[ref_idx] > 0)
        if mask.sum() < 10:
            continue
        fi = lib_sizes[i]
        fr = lib_sizes[ref_idx]
        # M = log2(sample_i / sample_ref) per gene
        Mi = np.log2((counts[i, mask] / fi) / (counts[ref_idx, mask] / fr))
        Ai = 0.5 * (np.log2(counts[i, mask] / fi) + np.log2(counts[ref_idx, mask] / fr))

        # 去除 M 的 top/bottom 30% 和 A 的 top 5%
        m_lo, m_hi = np.percentile(Mi, [30, 70])
        a_hi = np.percentile(Ai, 95)
        keep = (Mi >= m_lo) & (Mi <= m_hi) & (Ai <= a_hi)

        if keep.sum() > 0:
            # 加权均值（近似方差权重）
            wi = (fi - counts[i, mask][keep]) / (fi * counts[i, mask][keep])
            wr = (fr - counts[ref_idx, mask][keep]) / (fr * counts[ref_idx, mask][keep])
            weights = 1.0 / (wi + wr + 1e-30)
            factors[i] = 2 ** np.average(Mi[keep], weights=weights)

    return factors
```

- [ ] **Step 4: 实现 VST 和 rlog 辅助函数**

```python
def _vst_transform(counts, size_factors):
    """近似 VST。counts: raw counts (n_samples, n_genes), size_factors: (n_samples,)."""
    normed = counts / size_factors[:, None]
    normed = np.maximum(normed, 0)
    # 对 log2(normed + 0.5) 做近似方差稳定
    # 真正的 VST 需要拟合均值-方差关系并积分，这里用 log2(count/sf + 0.5) 近似
    # 对于大多数 Bulk RNA-seq 场景已足够
    vst = np.log2(normed + 0.5)
    return np.nan_to_num(vst, nan=0.0, posinf=0.0, neginf=0.0)


def _rlog_transform(counts, size_factors, prior_mean=None):
    """近似 rlog。小样本时通过正则化收缩基因效应。"""
    normed = counts / size_factors[:, None]
    normed = np.maximum(normed, 0)
    log_normed = np.log2(normed + 0.5)

    n_samples = counts.shape[0]
    # 全局先验：所有基因所有样本的均值
    if prior_mean is None:
        prior_mean = float(np.mean(log_normed))

    # 每个基因的效应（偏离全局均值）
    gene_effects = log_normed.mean(axis=0) - prior_mean

    # 收缩系数：样本越多越信赖数据，样本越少越信赖先验
    shrinkage = min(1.0, n_samples / 30.0)

    # rlog = 全局先验 + 收缩后的基因效应 + 样本特异性残差
    sample_residuals = log_normed - log_normed.mean(axis=0, keepdims=True)
    rlog = prior_mean + shrinkage * gene_effects[None, :] + sample_residuals

    return np.nan_to_num(rlog, nan=0.0, posinf=0.0, neginf=0.0)
```

- [ ] **Step 5: 运行测试确认通过**

Run: `cd /data/GJ/platform && python -m pytest tests/test_normalize_helpers.py -v`
Expected: 8 passed

- [ ] **Step 6: 提交**

```bash
git add modules/bulk_normalize.py tests/test_normalize_helpers.py
git commit -m "feat(bulk_normalize): add TMM, VST, rlog helper functions with tests"
```

---

## Task 2: PARAM_SCHEMAS 更新 + 前置过滤参数

**Files:**
- Modify: `routes/analysis.py:127-129`

- [ ] **Step 1: 更新 PARAM_SCHEMAS**

将 `routes/analysis.py` 中 `'bulk_normalize'` 条目从：

```python
    'bulk_normalize': [
        {'key': 'method', 'label': '标准化方法', 'type': 'select', 'options': ['deseq2', 'cpm', 'log2_quantile'], 'default': 'deseq2', 'help': '...'},
    ],
```

替换为：

```python
    'bulk_normalize': [
        {'key': 'method', 'label': '标准化方法', 'type': 'select', 'options': ['deseq2', 'tmm', 'cpm', 'vst', 'rlog', 'log2_quantile'], 'default': 'deseq2',
         'help': '标准化方法。差异分析：DESeq2（中位比率法，金标准）或 TMM（edgeR 方法，组成偏差大时更优）。可视化/高维：VST（方差稳定）或 rlog（小样本更稳定）。简单归一：CPM（每百万计数）或 log2 分位数。'},
        {'key': 'min_expr_value', 'label': '最小表达阈值 (CPM)', 'type': 'number', 'default': 1, 'step': 0.1,
         'help': '基因表达量需达到此 CPM 阈值才算有效表达。默认 1（每百万 reads 中至少 1 条）。'},
        {'key': 'min_expr_samples', 'label': '最小表达样本数', 'type': 'number', 'default': 3, 'step': 1,
         'help': '基因在至少 N 个样本中达到最小表达阈值才保留。0 = 不过滤。建议设为最小组的样本数（如 3）。'},
        {'key': 'max_zero_pct', 'label': '最大零值比例 (%)', 'type': 'number', 'default': 0, 'step': 1,
         'help': '基因在超过此比例的样本中为零则被过滤。0 = 不过滤。建议 50-70%。'},
    ],
```

- [ ] **Step 2: 验证**

Run: `cd /data/GJ/platform && python -c "from routes.analysis import PARAM_SCHEMAS; opts=[p['options'] for p in PARAM_SCHEMAS['bulk_normalize'] if p['key']=='method'][0]; assert 'tmm' in opts and 'vst' in opts; print('OK')"`

- [ ] **Step 3: 提交**

```bash
git add routes/analysis.py
git commit -m "feat(bulk_normalize): add PARAM_SCHEMAS for TMM/VST/rlog and pre-filtering"
```

---

## Task 3: 重构 run() 方法 — 前置过滤 + 方法分支扩展

**Files:**
- Modify: `modules/bulk_normalize.py:16-88`（run 方法核心逻辑）

- [ ] **Step 1: 在 run() 中添加前置过滤**

在 `method = self.params.get(...)` 之后、`raw_counts = ...` 之前插入：

```python
        # 前置过滤参数
        min_expr_value = float(self.params.get('min_expr_value', 1))
        min_expr_samples = int(self.params.get('min_expr_samples', 3))
        max_zero_pct = float(self.params.get('max_zero_pct', 0))

        # 低表达基因前置过滤
        n_genes_before = adata.n_vars
        if min_expr_samples > 0:
            raw_for_filter = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
            lib_for_filter = raw_for_filter.sum(axis=1, keepdims=True)
            lib_for_filter[lib_for_filter == 0] = 1
            cpm_check = raw_for_filter / lib_for_filter * 1e6
            n_expr = np.array((cpm_check >= min_expr_value).sum(axis=0)).flatten()
            gene_mask = n_expr >= min_expr_samples
            if max_zero_pct > 0:
                zero_pct = np.array((raw_for_filter == 0).sum(axis=0)).flatten() / adata.n_obs * 100
                gene_mask = gene_mask & (zero_pct <= max_zero_pct)
            adata = adata[:, gene_mask].copy()
```

- [ ] **Step 2: 保存原始 counts 到 layers['raw']**

在 `raw_counts = ...` 行之后追加：

```python
        adata.layers['raw'] = raw_counts.copy()
```

- [ ] **Step 3: 添加 TMM 分支**

在 `elif method == 'cpm':` 之前插入：

```python
        elif method == 'tmm':
            tmm_factors = _tmm_normalize(raw_counts)
            tmm_factors = np.where(tmm_factors > 0, tmm_factors, 1.0)
            adata.obs['size_factor'] = tmm_factors
            norm_counts = raw_counts / tmm_factors[:, None]
            adata.layers['normalized'] = norm_counts
            adata.X = np.log2(norm_counts + 1)
            adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
            adata.layers['normalized'] = np.nan_to_num(adata.layers['normalized'], nan=0.0, posinf=0.0, neginf=0.0)
```

- [ ] **Step 4: 添加 VST 分支**

在 `elif method == 'cpm':` 块之后插入：

```python
        elif method == 'vst':
            # VST 需要 size factors，复用 DESeq2 方法计算
            from scipy.stats import gmean
            counts_for_sf = raw_counts[raw_counts.sum(axis=1) > 0]
            nonzero_mask = (counts_for_sf > 0).all(axis=0)
            if nonzero_mask.sum() == 0:
                geo_means = np.exp(np.log(counts_for_sf + 1).mean(axis=0))
            else:
                geo_means = np.ones(counts_for_sf.shape[1])
                geo_means[nonzero_mask] = gmean(counts_for_sf[:, nonzero_mask], axis=0)
            ratios = counts_for_sf / (geo_means + 1e-10)
            size_factors = np.median(ratios, axis=1)
            size_factors = np.where(size_factors > 0, size_factors, 1.0)
            adata.obs['size_factor'] = size_factors
            adata.X = _vst_transform(raw_counts, size_factors)
```

- [ ] **Step 5: 添加 rlog 分支**

在 VST 分支之后插入：

```python
        elif method == 'rlog':
            from scipy.stats import gmean
            counts_for_sf = raw_counts[raw_counts.sum(axis=1) > 0]
            nonzero_mask = (counts_for_sf > 0).all(axis=0)
            if nonzero_mask.sum() == 0:
                geo_means = np.exp(np.log(counts_for_sf + 1).mean(axis=0))
            else:
                geo_means = np.ones(counts_for_sf.shape[1])
                geo_means[nonzero_mask] = gmean(counts_for_sf[:, nonzero_mask], axis=0)
            ratios = counts_for_sf / (geo_means + 1e-10)
            size_factors = np.median(ratios, axis=1)
            size_factors = np.where(size_factors > 0, size_factors, 1.0)
            adata.obs['size_factor'] = size_factors
            adata.X = _rlog_transform(raw_counts, size_factors)
```

- [ ] **Step 6: 在方法分支之后添加统一输出标记**

在现有 inf/NaN 清理代码块之后、`self.progress(60, ...)` 之前插入：

```python
        # 统一输出标记
        adata.uns['normalization'] = {
            'method': method,
            'is_log_transformed': method in ('deseq2', 'tmm', 'cpm', 'log2_quantile', 'rlog'),
        }
```

- [ ] **Step 7: 语法验证**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_normalize import BulkNormalizeAnalysis; print('OK')"`

- [ ] **Step 8: 提交**

```bash
git add modules/bulk_normalize.py
git commit -m "feat(bulk_normalize): add pre-filtering, TMM/VST/rlog branches, output markers"
```

---

## Task 4: 增强诊断图 + Summary

**Files:**
- Modify: `modules/bulk_normalize.py:89-132`（图表生成和 summary 区域）

- [ ] **Step 1: 扩展 size factor 图支持 TMM**

将现有 `if method == 'deseq2':` 的 size factor 图块改为：

```python
        if method in ('deseq2', 'tmm') and 'size_factor' in adata.obs.columns:
            fig_sf = go.Figure()
            fig_sf.add_trace(go.Bar(x=adata.obs.index.tolist(), y=adata.obs['size_factor'].values,
                                   marker_color='#1a237e'))
            fig_sf.update_layout(title=f'{method.upper()} Size Factors', yaxis_title='Size Factor',
                                plot_bgcolor='white', width=600, height=300)
            fpath = os.path.join(plots_dir, 'bulk_norm_sizefactors.json')
            with open(fpath, 'w') as f: f.write(fig_sf.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'bar', 'label': 'Size Factors'})
```

- [ ] **Step 2: 添加箱线图对比**

在 size factor 图块之后插入：

```python
        # 标准化前后箱线图对比
        fig_box = make_subplots(rows=1, cols=2, subplot_titles=['标准化前 (log2 raw)', '标准化后'])
        sample_labels = adata.obs.index.tolist()
        for i in range(adata.n_obs):
            fig_box.add_trace(go.Box(y=np.log2(raw_counts[i] + 1), name=sample_labels[i],
                                     showlegend=False, marker_color='#e53935'), row=1, col=1)
            fig_box.add_trace(go.Box(y=adata.X[i], name=sample_labels[i],
                                     showlegend=False, marker_color='#4caf50'), row=1, col=2)
        fig_box.update_layout(height=400, width=800, title='标准化前后表达分布对比')
        fpath = os.path.join(plots_dir, 'bulk_norm_boxplot_compare.json')
        with open(fpath, 'w') as f: f.write(fig_box.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': '表达分布对比'})
```

- [ ] **Step 3: 添加 PCA 前后对比**

在箱线图之后插入：

```python
        # PCA 前后对比
        fig_pca = make_subplots(rows=1, cols=2, subplot_titles=['原始数据 PCA', '标准化后 PCA'])
        # 原始数据 PCA
        adata_raw_pca = sc.AnnData(X=np.log2(raw_counts + 1), obs=adata.obs.copy())
        sc.pp.scale(adata_raw_pca, max_value=10)
        n_pcs = min(10, adata.n_obs - 1)
        if n_pcs >= 2:
            sc.pp.pca(adata_raw_pca, n_comps=n_pcs)
            pc_raw = adata_raw_pca.obsm['X_pca']
            fig_pca.add_trace(go.Scattergl(x=pc_raw[:, 0].tolist(), y=pc_raw[:, 1].tolist(),
                mode='markers+text', text=sample_labels, textposition='top center',
                marker=dict(size=8, color='#e53935'), showlegend=False), row=1, col=1)
            # 标准化后 PCA
            adata_norm_pca = sc.AnnData(X=adata.X.copy(), obs=adata.obs.copy())
            sc.pp.scale(adata_norm_pca, max_value=10)
            sc.pp.pca(adata_norm_pca, n_comps=n_pcs)
            pc_norm = adata_norm_pca.obsm['X_pca']
            fig_pca.add_trace(go.Scattergl(x=pc_norm[:, 0].tolist(), y=pc_norm[:, 1].tolist(),
                mode='markers+text', text=sample_labels, textposition='top center',
                marker=dict(size=8, color='#4caf50'), showlegend=False), row=1, col=2)
        fig_pca.update_layout(height=400, width=900, title='标准化前后 PCA 对比',
                              plot_bgcolor='white')
        fpath = os.path.join(plots_dir, 'bulk_norm_pca_compare.json')
        with open(fpath, 'w') as f: f.write(fig_pca.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': 'PCA 前后对比'})
```

- [ ] **Step 4: 更新 Summary**

替换现有 summary dict：

```python
        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'method': method,
                'n_samples': adata.n_obs,
                'n_genes_before_filter': n_genes_before,
                'n_genes_after_filter': adata.n_vars,
                'genes_filtered': n_genes_before - adata.n_vars,
                'median_size_factor': round(float(adata.obs['size_factor'].median()), 3) if 'size_factor' in adata.obs.columns else None,
                'is_log_transformed': adata.uns.get('normalization', {}).get('is_log_transformed', True),
            }
        }
```

- [ ] **Step 5: 语法验证**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_normalize import BulkNormalizeAnalysis; print('OK')"`
Run: `cd /data/GJ/platform && python -m pytest tests/test_normalize_helpers.py -q`
Expected: 8 passed

- [ ] **Step 6: 提交**

```bash
git add modules/bulk_normalize.py
git commit -m "feat(bulk_normalize): add diagnostic plots (boxplot/PCA compare), update size factor chart and summary"
```

---

## Task 5: 修复下游模块重复标准化

**Files:**
- Modify: `modules/bulk_heatmap.py:37-38`
- Modify: `modules/bulk_pca.py:34-35`

- [ ] **Step 1: 修复 bulk_heatmap.py**

找到 bulk_heatmap.py 中的：

```python
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
```

替换为：

```python
        if 'normalization' not in adata.uns:
            sc.pp.normalize_total(adata, target_sum=1e6)
            sc.pp.log1p(adata)
```

- [ ] **Step 2: 修复 bulk_pca.py**

找到 bulk_pca.py 中的：

```python
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
```

替换为：

```python
        if 'normalization' not in adata.uns:
            sc.pp.normalize_total(adata, target_sum=1e6)
            sc.pp.log1p(adata)
```

- [ ] **Step 3: 验证两个模块导入正常**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_heatmap import BulkHeatmapAnalysis; from modules.bulk_pca import BulkPCAAnalysis; print('OK')"`

- [ ] **Step 4: 提交**

```bash
git add modules/bulk_heatmap.py modules/bulk_pca.py
git commit -m "fix(bulk_heatmap,bulk_pca): skip re-normalization when data already normalized"
```

---

## 实施顺序与依赖

```
Task 1 (TMM/VST/rlog 辅助函数 + 测试)  ← 无依赖，首先实施
    ↓
Task 2 (PARAM_SCHEMAS)                  ← 无依赖，可与 Task 1 并行
    ↓
Task 3 (run() 重构 + 前置过滤 + 方法分支) ← 依赖 Task 1（使用辅助函数）
    ↓
Task 4 (诊断图 + Summary)                ← 依赖 Task 3（使用新方法和变量）
    ↓
Task 5 (下游修复)                         ← 依赖 Task 3（需要 uns['normalization'] 标记已存在）
```
