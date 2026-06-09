# Bulk RNA-seq QC 模块增强 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 bulk_qc.py 中增量添加核糖体比例、Gini 系数、基因过滤、批次分析、离群检测等高级 QC 功能。

**Architecture:** 在现有 bulk_qc.py 的 `run()` 方法中按顺序插入新逻辑，用辅助函数（`_gini`, `_infer_groups`, `_detect_outliers_mahal`）保持模块化。PARAM_SCHEMAS 在 routes/analysis.py 中追加 6 个新参数。输出从 4 张图扩展为 ~12 个文件。

**Tech Stack:** scanpy, plotly, numpy, scipy, pandas, scikit-learn (PCA)

---

## 文件结构

| 文件 | 操作 | 职责 |
|------|------|------|
| `modules/bulk_qc.py` | 大幅修改 | 所有 QC 逻辑、辅助函数、图表生成 |
| `routes/analysis.py:116-120` | 修改 | PARAM_SCHEMAS['bulk_qc'] 追加 6 个参数 |
| `tests/test_bulk_qc_helpers.py` | 新建 | 辅助函数单元测试 |

---

## Task 1: 辅助函数与单元测试

**Files:**
- Create: `tests/test_bulk_qc_helpers.py`
- Modify: `modules/bulk_qc.py`（末尾追加辅助函数）

- [ ] **Step 1: 编写辅助函数测试**

```python
# tests/test_bulk_qc_helpers.py
import numpy as np
import pytest


def test_gini_uniform():
    """完全均匀分布 → Gini = 0"""
    from modules.bulk_qc import _gini
    vals = np.array([100.0, 100.0, 100.0, 100.0])
    assert abs(_gini(vals)) < 1e-6


def test_gini_concentrated():
    """一个值占主导 → Gini 接近 1"""
    from modules.bulk_qc import _gini
    vals = np.array([0.0, 0.0, 0.0, 1000.0])
    assert _gini(vals) > 0.9


def test_gini_all_zeros():
    """全零输入 → Gini = 0"""
    from modules.bulk_qc import _gini
    vals = np.array([0.0, 0.0, 0.0])
    assert _gini(vals) == 0.0


def test_gini_empty():
    """空数组 → Gini = 0"""
    from modules.bulk_qc import _gini
    vals = np.array([])
    assert _gini(vals) == 0.0


def test_infer_groups_dash():
    """样本名含 '-' → 取前缀"""
    from modules.bulk_qc import _infer_groups
    names = ['ctrl-1', 'ctrl-2', 'hmc3-1', 'hmc3-2']
    assert _infer_groups(names) == ['ctrl', 'ctrl', 'hmc3', 'hmc3']


def test_infer_groups_underscore():
    """样本名含 '_' → 取前缀"""
    from modules.bulk_qc import _infer_groups
    names = ['ctrl_1', 'rapa_3']
    assert _infer_groups(names) == ['ctrl', 'rapa']


def test_infer_groups_no_sep():
    """无分隔符 → 样本名本身"""
    from modules.bulk_qc import _infer_groups
    names = ['sampleA', 'sampleB']
    assert _infer_groups(names) == ['sampleA', 'sampleB']


def test_detect_outliers_no_outlier():
    """正常聚类数据 → 无离群点"""
    from modules.bulk_qc import _detect_outliers_mahal
    np.random.seed(42)
    coords = np.random.randn(20, 3)
    outliers = _detect_outliers_mahal(coords, [f's{i}' for i in range(20)])
    assert len(outliers) == 0


def test_detect_outliers_with_outlier():
    """含明显离群点 → 检测到"""
    from modules.bulk_qc import _detect_outliers_mahal
    np.random.seed(42)
    coords = np.random.randn(20, 3)
    coords = np.vstack([coords, [[50, 50, 50]]])  # 离群点
    samples = [f's{i}' for i in range(21)]
    outliers = _detect_outliers_mahal(coords, samples)
    assert 's20' in outliers


def test_detect_outliers_too_few():
    """样本数 < 4 → 返回空"""
    from modules.bulk_qc import _detect_outliers_mahal
    coords = np.array([[0, 0], [1, 1], [2, 2]])
    outliers = _detect_outliers_mahal(coords, ['a', 'b', 'c'])
    assert outliers == []
```

- [ ] **Step 2: 运行测试确认失败**

Run: `cd /data/GJ/platform && python -m pytest tests/test_bulk_qc_helpers.py -v`
Expected: 全部 FAIL（函数尚未定义）

- [ ] **Step 3: 实现辅助函数**

在 `modules/bulk_qc.py` 文件末尾追加：

```python
# --- 辅助函数 ---


def _gini(values):
    """计算 Gini 系数。values 为原始 count 数组（非负）。"""
    sorted_vals = np.sort(values[values > 0])
    n = len(sorted_vals)
    if n == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return (2.0 * np.sum(index * sorted_vals) / (n * np.sum(sorted_vals))) - (n + 1) / n


def _infer_groups(sample_names):
    """从样本名推断分组，取第一个分隔符前的前缀。"""
    groups = []
    for name in sample_names:
        name = str(name)
        assigned = False
        for sep in ['-', '_']:
            if sep in name:
                groups.append(name.split(sep)[0])
                assigned = True
                break
        if not assigned:
            groups.append(name)
    return groups


def _detect_outliers_mahal(pca_coords, sample_names):
    """基于 PCA 坐标的马氏距离检测离群样本，返回离群样本名列表。"""
    if pca_coords.shape[0] < 4 or pca_coords.shape[1] < 2:
        return []
    n_components = min(3, pca_coords.shape[1])
    coords = pca_coords[:, :n_components]
    mean = coords.mean(axis=0)
    cov = np.cov(coords.T)
    try:
        cov_inv = np.linalg.pinv(cov)
    except np.linalg.LinAlgError:
        return []
    distances = []
    for i in range(coords.shape[0]):
        diff = coords[i] - mean
        d = np.sqrt(diff @ cov_inv @ diff)
        distances.append(d)
    distances = np.array(distances)
    med = np.median(distances)
    mad = np.median(np.abs(distances - med))
    if mad < 1e-10:
        return []
    threshold = med + 3 * 1.4826 * mad
    return [sample_names[i] for i, d in enumerate(distances) if d > threshold]
```

- [ ] **Step 4: 运行测试确认通过**

Run: `cd /data/GJ/platform && python -m pytest tests/test_bulk_qc_helpers.py -v`
Expected: 10 passed

- [ ] **Step 5: 提交**

```bash
git add modules/bulk_qc.py tests/test_bulk_qc_helpers.py
git commit -m "feat(bulk_qc): add helper functions for Gini, group inference, outlier detection"
```

---

## Task 2: 更新 PARAM_SCHEMAS

**Files:**
- Modify: `routes/analysis.py:116-120`

- [ ] **Step 1: 更新 PARAM_SCHEMAS['bulk_qc']**

将 `routes/analysis.py` 中的 `'bulk_qc'` 条目从：

```python
    'bulk_qc': [
        {'key': 'min_counts', 'label': '最小文库 reads 数', 'type': 'number', 'default': 100000, 'help': '最小文库 reads 数。低于此值的样本被过滤。人类/小鼠 RNA-seq 通常要求 ≥100000，小样本可降至 50000。'},
        {'key': 'min_genes', 'label': '最小检测基因数', 'type': 'number', 'default': 5000, 'help': '每个样本检测到的最小基因数。低于此值的样本可能质量差。通常 5000-8000。'},
        {'key': 'max_mt_pct', 'label': '最大线粒体基因比例 (%)', 'type': 'number', 'default': 20.0, 'step': 0.1, 'help': '最大线粒体基因比例（%）。高于此值的样本可能降解严重。RNA-seq 通常 15-20%。'},
    ],
```

替换为：

```python
    'bulk_qc': [
        {'key': 'min_counts', 'label': '最小文库 reads 数', 'type': 'number', 'default': 100000, 'help': '最小文库 reads 数。低于此值的样本被过滤。人类/小鼠 RNA-seq 通常要求 ≥100000，小样本可降至 50000。'},
        {'key': 'min_genes', 'label': '最小检测基因数', 'type': 'number', 'default': 5000, 'help': '每个样本检测到的最小基因数。低于此值的样本可能质量差。通常 5000-8000。'},
        {'key': 'max_mt_pct', 'label': '最大线粒体基因比例 (%)', 'type': 'number', 'default': 20.0, 'step': 0.1, 'help': '最大线粒体基因比例（%）。高于此值的样本可能降解严重。RNA-seq 通常 15-20%。'},
        {'key': 'max_ribo_pct', 'label': '最大核糖体基因比例 (%)', 'type': 'number', 'default': 40.0, 'step': 0.1, 'help': '最大核糖体基因比例（%）。RPL/RPS 基因比例过高提示 rRNA 污染。PolyA 建库通常 < 5-10%，rRNA 去除建库可至 40-50%。'},
        {'key': 'min_gini', 'label': '最小文库复杂度 (Gini)', 'type': 'number', 'default': 0, 'step': 0.01, 'help': '最小 Gini 系数（0 = 不过滤）。Gini > 0.8 提示文库复杂度低（PCR 过度扩增）。'},
        {'key': 'min_sample_expr', 'label': '基因最低表达样本数', 'type': 'number', 'default': 0, 'step': 1, 'help': '基因在至少 N 个样本中 CPM > 1 才保留。0 = 不过滤。建议设为样本总数的 10-20%。'},
        {'key': 'group_column', 'label': '分组列名（可选）', 'type': 'text', 'default': '', 'help': '样本分组列名（adata.obs 中的列）。留空则自动从样本名推断（取第一个分隔符前的前缀）。填写后启用组内/组间距离分析和分组着色图。'},
        {'key': 'detect_outliers', 'label': '检测离群样本', 'type': 'checkbox', 'default': True, 'help': '基于 PCA 马氏距离检测离群样本。仅在 summary 中告警，不自动剔除。'},
        {'key': 'filter_strategy', 'label': '过滤策略', 'type': 'select', 'options': ['conservative', 'standard', 'custom'], 'default': 'standard', 'help': 'conservative：宽松阈值（适合小样本）；standard：推荐阈值；custom：自定义所有阈值。'},
    ],
```

- [ ] **Step 2: 验证导入正常**

Run: `cd /data/GJ/platform && python -c "from routes.analysis import PARAM_SCHEMAS; keys=[p['key'] for p in PARAM_SCHEMAS['bulk_qc']]; print(keys)"`
Expected: 包含 `max_ribo_pct`, `min_gini`, `min_sample_expr`, `group_column`, `detect_outliers`, `filter_strategy`

- [ ] **Step 3: 提交**

```bash
git add routes/analysis.py
git commit -m "feat(bulk_qc): add PARAM_SCHEMAS for ribo%, Gini, gene filtering, batch, outlier detection"
```

---

## Task 3: 扩展样本级 QC 指标（Ribo%、Gini、novelty）

**Files:**
- Modify: `modules/bulk_qc.py:29-50`（QC 指标计算与阈值解析区域）

- [ ] **Step 1: 在现有 MT 指标计算之后追加 Ribo%、Gini、novelty**

在 `modules/bulk_qc.py` 中，找到 `sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], ...)` 这行（约第 41 行），在其后追加：

```python
        # 核糖体基因检测
        if 'gene_name' in adata.var.columns:
            gene_names_for_ribo = adata.var['gene_name'].fillna('').astype(str)
        else:
            gene_names_for_ribo = adata.var_names.astype(str)
        adata.var['ribo'] = gene_names_for_ribo.str.startswith(('RPL', 'RPS'))
        sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
        ribo_pct = adata.obs['pct_counts_ribo'].values if 'pct_counts_ribo' in adata.obs.columns else np.zeros(n_before)
```

然后在 `mt_pct = ...` 行之后追加 Gini 和 novelty 计算：

```python
        # 文库复杂度 (Gini) 和新颖度
        raw_counts = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
        gini_values = np.array([_gini(raw_counts[i]) for i in range(n_before)])
        novelty_values = np.log10(n_genes_detected + 1) / np.log10(lib_sizes + 1)
```

- [ ] **Step 2: 扩展阈值解析，加入新参数**

在现有阈值解析块（`min_counts = int(self.params.get(...))` 那三行）之后追加：

```python
        max_ribo_pct = float(self.params.get('max_ribo_pct', 40.0))
        min_gini = float(self.params.get('min_gini', 0))
        min_sample_expr = int(self.params.get('min_sample_expr', 0))
        detect_outliers = self.params.get('detect_outliers', True)
        filter_strategy = self.params.get('filter_strategy', 'standard')

        # 过滤策略预设覆盖阈值
        if filter_strategy == 'conservative':
            min_counts = min(min_counts, 50000)
            min_genes = min(min_genes, 3000)
            max_mt_pct = max(max_mt_pct, 30.0)
            max_ribo_pct = max(max_ribo_pct, 60.0)
        # standard 使用当前默认值，不覆盖
```

注意：`conservative` 策略使用 `min()`/`max()` 而非直接赋值，这样用户手动填写的值如果比预设更严格会被保留。

- [ ] **Step 3: 扩展样本过滤 mask**

将现有 `mask = (lib_sizes >= min_counts) & (n_genes_detected >= min_genes) & (mt_pct <= max_mt_pct)` 替换为：

```python
        # 分组推断
        sample_names = adata.obs.index.tolist()
        group_col = self.params.get('group_column', '').strip()
        if group_col and group_col in adata.obs.columns:
            groups = adata.obs[group_col].astype(str).tolist()
        else:
            groups = _infer_groups(sample_names)
            adata.obs['_auto_group'] = groups
            group_col = '_auto_group'

        # 样本过滤
        mask = (lib_sizes >= min_counts) & (n_genes_detected >= min_genes) & (mt_pct <= max_mt_pct) & (ribo_pct <= max_ribo_pct)
        if min_gini > 0:
            mask = mask & (gini_values >= min_gini)

        # 过滤日志
        filter_log_rows = []
        for i in range(n_before):
            fail_reasons = []
            if lib_sizes[i] < min_counts:
                fail_reasons.append(f'lib_size<{min_counts}')
            if n_genes_detected[i] < min_genes:
                fail_reasons.append(f'n_genes<{min_genes}')
            if mt_pct[i] > max_mt_pct:
                fail_reasons.append(f'mt_pct>{max_mt_pct}')
            if ribo_pct[i] > max_ribo_pct:
                fail_reasons.append(f'ribo_pct>{max_ribo_pct}')
            if min_gini > 0 and gini_values[i] < min_gini:
                fail_reasons.append(f'gini<{min_gini}')
            filter_log_rows.append({
                'sample': sample_names[i],
                'group': groups[i],
                'passed': len(fail_reasons) == 0,
                'lib_size': int(lib_sizes[i]),
                'n_genes': int(n_genes_detected[i]),
                'mt_pct': round(float(mt_pct[i]), 2),
                'ribo_pct': round(float(ribo_pct[i]), 2),
                'gini': round(float(gini_values[i]), 4),
                'novelty': round(float(novelty_values[i]), 4),
                'fail_reasons': '; '.join(fail_reasons) if fail_reasons else '',
            })

        adata_filtered = adata[mask].copy()
        n_after = adata_filtered.n_obs
```

- [ ] **Step 4: 生成 sample metrics CSV 和 filter log CSV**

在 `self.progress(60, "生成质控图表...")` 之前追加：

```python
        # 输出样本指标表
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)

        metrics_df = pd.DataFrame(filter_log_rows)
        metrics_csv = os.path.join(results_dir, 'bulk_qc_sample_metrics.csv')
        metrics_df.to_csv(metrics_csv, index=False)
        result_files.append({'file_path': metrics_csv, 'file_type': 'csv', 'category': 'table', 'label': '样本 QC 指标'})

        # 过滤日志
        filter_log_csv = os.path.join(results_dir, 'bulk_qc_filter_log.csv')
        metrics_df.to_csv(filter_log_csv, index=False)
        result_files.append({'file_path': filter_log_csv, 'file_type': 'csv', 'category': 'table', 'label': '过滤日志'})
```

- [ ] **Step 5: 验证导入和语法**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_qc import BulkQCAnalysis; print('OK')"`
Expected: OK

- [ ] **Step 6: 提交**

```bash
git add modules/bulk_qc.py
git commit -m "feat(bulk_qc): add ribosomal %, Gini coefficient, novelty score, sample filter log"
```

---

## Task 4: 基因过滤与管家基因检查

**Files:**
- Modify: `modules/bulk_qc.py`（在样本过滤之后插入基因过滤逻辑）

- [ ] **Step 1: 在 `adata_filtered` 生成之后、图表生成之前插入基因过滤逻辑**

找到 `adata_filtered = adata[mask].copy()` 和 `n_after = adata_filtered.n_obs` 之后，在 `self.progress(60, "生成质控图表...")` 之前插入：

```python
        # 基因层面过滤
        genes_before_filter = adata_filtered.n_vars
        gene_filter_rows = []
        if min_sample_expr > 0:
            raw_filt = adata_filtered.X.toarray() if hasattr(adata_filtered.X, 'toarray') else np.asarray(adata_filtered.X)
            lib_sizes_filt = raw_filt.sum(axis=1, keepdims=True)
            lib_sizes_filt[lib_sizes_filt == 0] = 1  # 避免除零
            cpm = raw_filt / lib_sizes_filt * 1e6
            expr_count_per_gene = (cpm > 1).sum(axis=0)
            gene_mask = expr_count_per_gene >= min_sample_expr

            removed_gene_names = adata_filtered.var_names[~gene_mask].tolist()
            for g in removed_gene_names:
                idx = list(adata_filtered.var_names).index(g)
                gene_filter_rows.append({
                    'gene': g,
                    'expressed_in_n_samples': int(expr_count_per_gene[idx]),
                    'action': 'removed',
                })
            adata_filtered = adata_filtered[:, gene_mask].copy()

        # 管家基因稳定性检查
        housekeeping_genes = ['GAPDH', 'ACTB', 'B2M', 'HPRT1', 'TBP', 'UBC', 'YWHAZ', 'SDHA', 'HMBS', 'RPLP0']
        found_hk = [g for g in housekeeping_genes if g in adata_filtered.var_names]

        # 基因过滤日志
        if min_sample_expr > 0:
            gene_filter_csv = os.path.join(results_dir, 'bulk_qc_gene_filter.csv')
            pd.DataFrame(gene_filter_rows).to_csv(gene_filter_csv, index=False)
            result_files.append({'file_path': gene_filter_csv, 'file_type': 'csv', 'category': 'table', 'label': '基因过滤日志'})
```

- [ ] **Step 2: 语法验证**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_qc import BulkQCAnalysis; print('OK')"`
Expected: OK

- [ ] **Step 3: 提交**

```bash
git add modules/bulk_qc.py
git commit -m "feat(bulk_qc): add gene-level filtering and housekeeping gene stability check"
```

---

## Task 5: 生成增强图表（3x2 总览、Ribo bar、Gini bar、PCA elbow）

**Files:**
- Modify: `modules/bulk_qc.py:53-70`（QC 总览图生成区域）

- [ ] **Step 1: 替换 2x2 总览图为 3x2 布局**

将现有的 `fig = make_subplots(rows=2, cols=2, ...)` 到 `fig.update_layout(...)` 整块替换为：

```python
        fig = make_subplots(rows=3, cols=2,
            subplot_titles=['文库大小分布', '检测基因数',
                           '线粒体基因比例', '核糖体基因比例',
                           'Gini 系数', '文库大小 vs 检测基因数'],
            vertical_spacing=0.08)
        sample_idx = list(range(n_before))
        fig.add_trace(go.Bar(x=sample_idx, y=lib_sizes.tolist(), marker_color='#1a237e', name='文库大小'), row=1, col=1)
        fig.add_trace(go.Bar(x=sample_idx, y=n_genes_detected.tolist(), marker_color='#283593', name='基因数'), row=1, col=2)
        fig.add_trace(go.Bar(x=sample_idx, y=mt_pct.tolist(), marker_color='#e53935', name='MT%'), row=2, col=1)
        fig.add_trace(go.Bar(x=sample_idx, y=ribo_pct.tolist(), marker_color='#ff8f00', name='Ribo%'), row=2, col=2)
        fig.add_trace(go.Bar(x=sample_idx, y=gini_values.tolist(), marker_color='#6a1b9a', name='Gini'), row=3, col=1)
        colors = ['#4caf50' if m else '#e53935' for m in mask]
        fig.add_trace(go.Scattergl(x=lib_sizes.tolist(), y=n_genes_detected.tolist(), mode='markers',
            marker=dict(color=colors, size=6), name='样本'), row=3, col=2)
        fig.update_layout(height=900, width=800, showlegend=False, title='Bulk RNA-seq 质控总览')
```

- [ ] **Step 2: 在 PCA 图生成之后追加 PCA elbow 图**

找到 PCA 图保存代码之后，追加：

```python
        # PCA 方差解释 elbow 图
        pca_variance = adata_filtered.uns.get('pca_variance', None)
        if pca_variance is None and 'X_pca' in adata_filtered.obsm:
            # 从 PCA 坐标计算方差解释比例
            pca_var = np.var(adata_filtered.obsm['X_pca'], axis=0)
            pca_variance = pca_var / pca_var.sum()
        if pca_variance is not None:
            n_pcs = len(pca_variance)
            pc_labels = [f'PC{i+1}' for i in range(n_pcs)]
            cumulative = np.cumsum(pca_variance)
            fig_elbow = go.Figure()
            fig_elbow.add_trace(go.Bar(x=pc_labels, y=pca_variance.tolist(),
                marker_color='#1a237e', name='方差比例'))
            fig_elbow.add_trace(go.Scatter(x=pc_labels, y=cumulative.tolist(),
                mode='lines+markers', marker_color='#e53935', name='累积比例', yaxis='y2'))
            fig_elbow.update_layout(
                title='PCA 方差解释比例',
                xaxis_title='主成分', yaxis_title='方差解释比例',
                yaxis2=dict(title='累积比例', overlaying='y', side='right', range=[0, 1.05]),
                width=600, height=400, plot_bgcolor='white')
            fpath = os.path.join(plots_dir, 'bulk_qc_pca_elbow.json')
            with open(fpath, 'w') as f: f.write(fig_elbow.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': 'PCA 方差解释'})
```

- [ ] **Step 3: 语法验证**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_qc import BulkQCAnalysis; print('OK')"`
Expected: OK

- [ ] **Step 4: 提交**

```bash
git add modules/bulk_qc.py
git commit -m "feat(bulk_qc): expand overview to 3x2, add PCA elbow plot"
```

---

## Task 6: 相关性热图增强（分组 annotation bar）+ 离群检测

**Files:**
- Modify: `modules/bulk_qc.py`（相关性热图生成区域 + PCA 区域）

- [ ] **Step 1: 增强相关性热图，添加分组 annotation bar**

将现有的相关性热图代码（`fig_corr = go.Figure()` 到 `fig_corr.update_layout(...)`）替换为：

```python
        # 样本相关性热图（带分组 annotation bar）
        norm_for_corr = adata_filtered.copy()
        sc.pp.normalize_total(norm_for_corr, target_sum=1e6)
        sc.pp.log1p(norm_for_corr)
        corr_data = norm_for_corr.X if not hasattr(norm_for_corr.X, 'toarray') else norm_for_corr.X.toarray()
        corr_matrix = np.corrcoef(corr_data)
        sample_labels_corr = norm_for_corr.obs.index.tolist()
        filtered_groups = [groups[sample_names.index(s)] for s in sample_labels_corr]
        unique_groups = sorted(set(filtered_groups))
        group_color_map = {g: f'hsl({i*360//len(unique_groups)},70%,50%)' for i, g in enumerate(unique_groups)}
        group_colors = [group_color_map[g] for g in filtered_groups]

        fig_corr = go.Figure()
        fig_corr.add_trace(go.Heatmap(
            z=corr_matrix.tolist(), x=sample_labels_corr, y=sample_labels_corr,
            colorscale='Blues', zmin=0, zmax=1,
            colorbar=dict(title='Pearson r'),
            hovertemplate='%{y} vs %{x}<br>r = %{z:.3f}<extra></extra>'
        ))
        # 分组 annotation bar
        for g in unique_groups:
            fig_corr.add_trace(go.Bar(
                x=[1], y=[g], orientation='h',
                marker_color=group_color_map[g], showlegend=True, name=g,
                xaxis='x2'))
        fig_corr.update_layout(
            title='样本相关性热图 (Pearson)',
            height=max(400, n_after * 30 + 100), width=max(500, n_after * 30 + 200),
            plot_bgcolor='white',
        )
```

- [ ] **Step 2: 在 PCA 图生成代码中添加分组着色和离群标注**

将现有 PCA 图代码替换为增强版：

```python
        self.progress(80, "运行 PCA 离群检测...")
        sc.pp.normalize_total(adata_filtered, target_sum=1e6)
        sc.pp.log1p(adata_filtered)
        n_comps = min(10, n_after - 1)
        if n_comps >= 2:
            sc.pp.pca(adata_filtered, n_comps=n_comps)
            pc = adata_filtered.obsm['X_pca']
            filtered_sample_names = adata_filtered.obs.index.tolist()
            filtered_groups_pca = [groups[sample_names.index(s)] for s in filtered_sample_names]
            unique_groups_pca = sorted(set(filtered_groups_pca))
            group_color_map_pca = {g: f'hsl({i*360//max(1,len(unique_groups_pca))},70%,50%)' for i, g in enumerate(unique_groups_pca)}

            # 离群检测
            outlier_samples = []
            if detect_outliers:
                outlier_samples = _detect_outliers_mahal(pc, filtered_sample_names)

            fig_pca = go.Figure()
            for g in unique_groups_pca:
                grp_mask = [filtered_groups_pca[i] == g for i in range(n_after)]
                grp_idx = [i for i in range(n_after) if grp_mask[i]]
                fig_pca.add_trace(go.Scattergl(
                    x=pc[grp_idx, 0].tolist(), y=pc[grp_idx, 1].tolist(),
                    mode='markers+text',
                    text=[filtered_sample_names[i] for i in grp_idx],
                    textposition='top center',
                    marker=dict(size=8, color=group_color_map_pca[g]),
                    name=g))
            # 离群点高亮
            if outlier_samples:
                out_idx = [filtered_sample_names.index(s) for s in outlier_samples if s in filtered_sample_names]
                if out_idx:
                    fig_pca.add_trace(go.Scattergl(
                        x=pc[out_idx, 0].tolist(), y=pc[out_idx, 1].tolist(),
                        mode='markers',
                        marker=dict(size=14, color='red', symbol='x', line=dict(width=2, color='darkred')),
                        name='离群样本'))
            fig_pca.update_layout(title='质控后样本 PCA', xaxis_title='PC1', yaxis_title='PC2',
                                 plot_bgcolor='white', width=600, height=500)
            fpath = os.path.join(plots_dir, 'bulk_qc_pca.json')
            with open(fpath, 'w') as f: f.write(fig_pca.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': '样本 PCA'})
```

- [ ] **Step 3: 语法验证**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_qc import BulkQCAnalysis; print('OK')"`
Expected: OK

- [ ] **Step 4: 提交**

```bash
git add modules/bulk_qc.py
git commit -m "feat(bulk_qc): enhance correlation heatmap with group bar, PCA with grouping and outlier highlight"
```

---

## Task 7: 新增图表（组距离、pairs plot、violin、housekeeping）

**Files:**
- Modify: `modules/bulk_qc.py`（在 PCA 图之后追加新图表）

- [ ] **Step 1: 在 PCA elbow 图代码之后追加组内/组间距离图**

```python
        # 组内 vs 组间距离箱线图
        if len(set(groups)) > 1 and n_after > 3:
            norm_arr = norm_for_corr.X if not hasattr(norm_for_corr.X, 'toarray') else norm_for_corr.X.toarray()
            corr_mat_all = np.corrcoef(norm_arr)
            dist_mat = 1 - corr_mat_all
            intra_dists, inter_dists = [], []
            for i in range(n_after):
                for j in range(i + 1, n_after):
                    if filtered_groups[i] == filtered_groups[j]:
                        intra_dists.append(float(dist_mat[i, j]))
                    else:
                        inter_dists.append(float(dist_mat[i, j]))

            fig_dist = go.Figure()
            if intra_dists:
                fig_dist.add_trace(go.Box(y=intra_dists, name='组内距离', marker_color='#4caf50'))
            if inter_dists:
                fig_dist.add_trace(go.Box(y=inter_dists, name='组间距离', marker_color='#e53935'))
            # Welch t-test
            try:
                from scipy.stats import ttest_ind
                _, pval = ttest_ind(intra_dists, inter_dists, equal_var=False)
                title_suffix = f' (p={pval:.2e})'
            except Exception:
                title_suffix = ''
            fig_dist.update_layout(
                title=f'组内 vs 组间距离{title_suffix}',
                yaxis_title='1 - Pearson r', width=500, height=400, plot_bgcolor='white')
            fpath = os.path.join(plots_dir, 'bulk_qc_group_distance.json')
            with open(fpath, 'w') as f: f.write(fig_dist.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': '组内/组间距离'})
```

- [ ] **Step 2: 追加 pairs plot (散点矩阵)**

```python
        # QC 指标散点矩阵 (Pairs Plot)
        obs_filtered = adata_filtered.obs
        pairs_groups = [groups[sample_names.index(s)] for s in obs_filtered.index.tolist()]
        unique_pg = sorted(set(pairs_groups))
        pg_color_map = {g: f'hsl({i*360//max(1,len(unique_pg))},70%,50%)' for i, g in enumerate(unique_pg)}
        pg_colors = [pg_color_map[g] for g in pairs_groups]

        fig_pairs = go.Figure(data=go.Splom(
            dimensions=[
                dict(label='Library Size', values=obs_filtered['total_counts'].tolist()),
                dict(label='N Genes', values=obs_filtered['n_genes_by_counts'].tolist()),
                dict(label='MT%', values=obs_filtered['pct_counts_mt'].tolist() if 'pct_counts_mt' in obs_filtered.columns else [0]*n_after),
                dict(label='Ribo%', values=obs_filtered['pct_counts_ribo'].tolist() if 'pct_counts_ribo' in obs_filtered.columns else [0]*n_after),
            ],
            marker=dict(color=pg_colors, size=5, line=dict(width=0.5, color='white')),
            text=obs_filtered.index.tolist(),
            showupperhalf=False,
        ))
        fig_pairs.update_layout(title='QC 指标散点矩阵', width=700, height=700)
        fpath = os.path.join(plots_dir, 'bulk_qc_pairs_plot.json')
        with open(fpath, 'w') as f: f.write(fig_pairs.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': 'QC 指标散点矩阵'})
```

- [ ] **Step 3: 追加各组小提琴图**

```python
        # 各组 QC 指标小提琴图
        if len(unique_groups) > 1:
            import plotly.express as px
            violin_data = []
            for i, s in enumerate(obs_filtered.index.tolist()):
                g = groups[sample_names.index(s)]
                violin_data.append({
                    'group': g,
                    'MT%': float(obs_filtered.loc[s, 'pct_counts_mt']) if 'pct_counts_mt' in obs_filtered.columns else 0,
                    'Ribo%': float(obs_filtered.loc[s, 'pct_counts_ribo']) if 'pct_counts_ribo' in obs_filtered.columns else 0,
                    'Library Size': float(obs_filtered.loc[s, 'total_counts']),
                    'N Genes': float(obs_filtered.loc[s, 'n_genes_by_counts']),
                })
            violin_df = pd.DataFrame(violin_data)
            fig_violin = make_subplots(rows=2, cols=2,
                subplot_titles=['MT%', 'Ribo%', 'Library Size', 'N Genes'],
                shared_xaxes=True)
            metrics = [('MT%', 1, 1), ('Ribo%', 1, 2), ('Library Size', 2, 1), ('N Genes', 2, 2)]
            for metric, row, col in metrics:
                for g in unique_groups:
                    vals = violin_df[violin_df['group'] == g][metric].tolist()
                    fig_violin.add_trace(go.Violin(
                        y=vals, name=g, box_visible=True, meanline_visible=True,
                        legendgroup=g, showlegend=(row == 1 and col == 1)),
                        row=row, col=col)
            fig_violin.update_layout(title='各组 QC 指标分布', height=600, width=700)
            fpath = os.path.join(plots_dir, 'bulk_qc_violin_by_group.json')
            with open(fpath, 'w') as f: f.write(fig_violin.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': '各组 QC 指标分布'})
```

- [ ] **Step 4: 追加管家基因热图**

```python
        # 管家基因稳定性热图
        if found_hk:
            norm_hk = adata_filtered[:, [g for g in found_hk if g in adata_filtered.var_names]].copy()
            if norm_hk.n_vars > 0:
                sc.pp.normalize_total(norm_hk, target_sum=1e6)
                sc.pp.log1p(norm_hk)
                hk_data = norm_hk.X if not hasattr(norm_hk.X, 'toarray') else norm_hk.X.toarray()
                hk_genes = norm_hk.var_names.tolist()
                hk_samples = norm_hk.obs.index.tolist()
                # CV 系数
                hk_cv = {}
                for j, g in enumerate(hk_genes):
                    vals = hk_data[:, j]
                    cv = float(np.std(vals) / (np.mean(vals) + 1e-10))
                    hk_cv[g] = cv
                cv_labels = [f'{g} (CV={hk_cv[g]:.2f})' for g in hk_genes]
                fig_hk = go.Figure(data=go.Heatmap(
                    z=hk_data.T.tolist(), x=hk_samples, y=cv_labels,
                    colorscale='YlOrRd', colorbar=dict(title='log10(CPM+1)')))
                fig_hk.update_layout(title='管家基因表达稳定性', width=max(400, len(hk_samples)*40+200), height=max(200, len(hk_genes)*30+100))
                fpath = os.path.join(plots_dir, 'bulk_qc_housekeeping.json')
                with open(fpath, 'w') as f: f.write(fig_hk.to_json(engine="json"))
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': '管家基因稳定性'})
```

- [ ] **Step 5: 语法验证**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_qc import BulkQCAnalysis; print('OK')"`
Expected: OK

- [ ] **Step 6: 提交**

```bash
git add modules/bulk_qc.py
git commit -m "feat(bulk_qc): add group distance, pairs plot, violin by group, housekeeping heatmap"
```

---

## Task 8: 更新 Summary 返回值 + 最终集成

**Files:**
- Modify: `modules/bulk_qc.py`（summary 和 return 区域）

- [ ] **Step 1: 替换现有 summary 为增强版**

将现有的 `return { 'output_adata': ..., 'result_files': ..., 'summary': {...} }` 替换为：

```python
        # 移除临时列
        if '_auto_group' in adata_filtered.obs.columns:
            adata_filtered.obs.drop(columns=['_auto_group'], inplace=True)

        removed_samples = [r['sample'] for r in filter_log_rows if not r['passed']]
        removed_reasons = {r['sample']: r['fail_reasons'] for r in filter_log_rows if not r['passed']}

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'samples_before': n_before,
                'samples_after': n_after,
                'samples_removed': n_before - n_after,
                'removed_samples': removed_samples,
                'removed_reasons': removed_reasons,
                'genes_before': genes_before_filter,
                'genes_after': adata_filtered.n_vars,
                'genes_removed': genes_before_filter - adata_filtered.n_vars,
                'outlier_samples': outlier_samples if detect_outliers else [],
                'filter_strategy': filter_strategy,
                'median_lib_size': int(np.median(adata_filtered.obs['total_counts'])),
                'median_genes': int(np.median(adata_filtered.obs['n_genes_by_counts'])),
                'median_ribo_pct': round(float(np.median(adata_filtered.obs['pct_counts_ribo'])), 2) if 'pct_counts_ribo' in adata_filtered.obs.columns else 0,
                'median_gini': round(float(np.median([_gini(raw_counts[sample_names.index(s)]) for s in adata_filtered.obs.index.tolist()])), 4),
            }
        }
```

- [ ] **Step 2: 确保 `outlier_samples` 变量在 n_comps < 2 时已定义**

在 Task 6 的 PCA 代码块外、`adata_filtered` 生成之后添加默认值：

```python
        outlier_samples = []
```

这样当 PCA 不可用时变量仍然存在。

- [ ] **Step 3: 语法验证**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_qc import BulkQCAnalysis; print('OK')"`
Expected: OK

- [ ] **Step 4: 运行完整测试套件**

Run: `cd /data/GJ/platform && python -m pytest tests/test_bulk_qc_helpers.py -v`
Expected: 10 passed

- [ ] **Step 5: 提交**

```bash
git add modules/bulk_qc.py
git commit -m "feat(bulk_qc): update summary with all new metrics, finalize integration"
```

---

## 实施顺序与依赖

```
Task 1 (辅助函数 + 测试)  ← 无依赖，首先实施
    ↓
Task 2 (PARAM_SCHEMAS)    ← 无依赖，可与 Task 1 并行
    ↓
Task 3 (样本级指标)        ← 依赖 Task 1（使用 _gini, _infer_groups）
    ↓
Task 4 (基因过滤)          ← 依赖 Task 3（使用过滤后的 adata）
    ↓
Task 5 (3x2 总览 + elbow) ← 依赖 Task 3（使用 ribo_pct, gini_values）
    ↓
Task 6 (热图增强 + PCA)    ← 依赖 Task 3, 5（使用 groups, _detect_outliers_mahal）
    ↓
Task 7 (新图表)            ← 依赖 Task 6（使用 filtered_groups, norm_for_corr）
    ↓
Task 8 (Summary 集成)      ← 依赖所有前序 Task
```
