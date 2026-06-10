# Bulk 热图增强实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 增强 bulk_heatmap 模块，支持灵活的基因选择策略、数据变换、聚类方法、颜色样式配置，并集成表达式筛选器。

**Architecture:** 提取 5 个可复用工具函数到 `visualization.py`，重构 `bulk_heatmap.py` 调用它们。新增 ~27 个参数到 PARAM_SCHEMAS。保持向后兼容旧参数名。

**Tech Stack:** Python, numpy, scipy, plotly, pandas

---

## 文件清单

| 文件 | 操作 | 职责 |
|------|------|------|
| `modules/visualization.py` | 修改 | 新增 5 个热图工具函数 |
| `tests/test_heatmap_helpers.py` | 新建 | 工具函数单元测试 |
| `routes/analysis.py` | 修改 | 替换 bulk_heatmap PARAM_SCHEMAS |
| `modules/bulk_heatmap.py` | 修改 | 重构 run() 使用工具函数，添加新功能 |

---

## Task 1: 添加 compute_gene_variability 工具函数

**Files:**
- Modify: `modules/visualization.py`
- Test: `tests/test_heatmap_helpers.py`

- [ ] **Step 1: 编写失败测试**

```python
# tests/test_heatmap_helpers.py
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from modules.visualization import compute_gene_variability


class TestComputeGeneVariability:
    @pytest.fixture
    def sample_data(self):
        # (5 samples, 3 genes)
        return np.array([
            [1.0, 10.0, 0.1],
            [2.0, 20.0, 0.2],
            [3.0, 30.0, 0.3],
            [4.0, 40.0, 0.4],
            [5.0, 50.0, 0.5],
        ])

    def test_var(self, sample_data):
        result = compute_gene_variability(sample_data, metric='var')
        assert len(result) == 3
        np.testing.assert_almost_equal(result[0], np.var([1, 2, 3, 4, 5]))
        np.testing.assert_almost_equal(result[1], np.var([10, 20, 30, 40, 50]))

    def test_mad(self, sample_data):
        result = compute_gene_variability(sample_data, metric='mad')
        assert len(result) == 3
        assert all(v >= 0 for v in result)
        # MAD of uniform sequence should be 0
        uniform = np.array([[1, 5], [1, 5], [1, 5], [1, 5]])
        result_u = compute_gene_variability(uniform, metric='mad')
        assert result_u[0] == 0.0

    def test_cv(self, sample_data):
        result = compute_gene_variability(sample_data, metric='cv')
        assert len(result) == 3
        # CV should be dimensionless, all positive
        assert all(v >= 0 for v in result)
        # Gene with higher relative variation should have higher CV
        assert result[0] > 0

    def test_range(self, sample_data):
        result = compute_gene_variability(sample_data, metric='range')
        np.testing.assert_almost_equal(result[0], 4.0)  # 5 - 1
        np.testing.assert_almost_equal(result[1], 40.0)  # 50 - 10

    def test_unknown_metric_raises(self, sample_data):
        with pytest.raises(ValueError):
            compute_gene_variability(sample_data, metric='unknown')
```

- [ ] **Step 2: 运行测试确认失败**

```bash
cd /data/GJ/platform && python -m pytest tests/test_heatmap_helpers.py::TestComputeGeneVariability -v
```

Expected: `ImportError: cannot import name 'compute_gene_variability'`

- [ ] **Step 3: 实现函数**

在 `modules/visualization.py` 末尾添加：

```python
def compute_gene_variability(data, metric='var'):
    """计算每个基因的变异度量。data: (samples, genes)"""
    if metric == 'var':
        return np.var(data, axis=0)
    elif metric == 'mad':
        median = np.median(data, axis=0)
        return np.median(np.abs(data - median), axis=0)
    elif metric == 'cv':
        mean = np.mean(data, axis=0)
        std = np.std(data, axis=0)
        return std / (np.abs(mean) + 1e-10)
    elif metric == 'range':
        return np.max(data, axis=0) - np.min(data, axis=0)
    else:
        raise ValueError(f"未知变异度量: {metric}，支持: var/mad/cv/range")
```

- [ ] **Step 4: 运行测试确认通过**

```bash
cd /data/GJ/platform && python -m pytest tests/test_heatmap_helpers.py::TestComputeGeneVariability -v
```

Expected: `5 passed`

- [ ] **Step 5: 提交**

```bash
git add modules/visualization.py tests/test_heatmap_helpers.py
git commit -m "feat(viz): add compute_gene_variability helper (var/mad/cv/range)"
```

---

## Task 2: 添加 transform_heatmap_data 工具函数

**Files:**
- Modify: `modules/visualization.py`
- Modify: `tests/test_heatmap_helpers.py`

- [ ] **Step 1: 编写失败测试**

```python
# 在 tests/test_heatmap_helpers.py 中追加
from modules.visualization import transform_heatmap_data


class TestTransformHeatmapData:
    @pytest.fixture
    def sample_data(self):
        return np.array([
            [10.0, 100.0, 5.0],
            [20.0, 200.0, np.nan],
            [30.0, 300.0, 15.0],
            [40.0, 400.0, 20.0],
        ])

    def test_zscore(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='zscore')
        # 每列均值应接近 0
        np.testing.assert_almost_equal(result[:, 0].mean(), 0.0, decimal=10)
        # 每列标准差应接近 1
        assert abs(result[:, 0].std() - 1.0) < 0.01

    def test_center(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='center')
        # 仅减均值，不除标准差
        np.testing.assert_almost_equal(result[:, 0].mean(), 0.0, decimal=10)
        # 标准差不为 1
        assert result[:, 0].std() != 1.0

    def test_none_scaling(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='none')
        # 不做变换，但 NaN 应被处理
        assert not np.any(np.isnan(result))

    def test_winsorize(self):
        data = np.array([[1, 100], [2, 200], [3, 300], [1000, 400], [5, 500]])
        result = transform_heatmap_data(data.copy(), row_scaling='none', winsorize='5pct')
        # 极端值应被截断
        assert result[3, 0] < 1000

    def test_missing_mean_fill(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='none', missing_value='mean_fill')
        assert not np.any(np.isnan(result))
        # NaN 被填充为该列均值
        col2_mean = np.nanmean(sample_data[:, 2])
        assert abs(result[1, 2] - col2_mean) < 0.01

    def test_missing_zero_fill(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='none', missing_value='zero_fill')
        assert result[1, 2] == 0.0

    def test_clip_range(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='zscore', clip_range=(-1, 1))
        assert np.all(result >= -1.0)
        assert np.all(result <= 1.0)

    def test_no_clip(self, sample_data):
        result = transform_heatmap_data(sample_data.copy(), row_scaling='zscore', clip_range=None)
        # 不截断时可能有超过 3 的值
        assert np.any(np.abs(result) > 0)
```

- [ ] **Step 2: 运行测试确认失败**

```bash
cd /data/GJ/platform && python -m pytest tests/test_heatmap_helpers.py::TestTransformHeatmapData -v
```

Expected: `ImportError`

- [ ] **Step 3: 实现函数**

在 `modules/visualization.py` 的 `compute_gene_variability` 之后添加：

```python
def transform_heatmap_data(data, row_scaling='zscore', pseudocount=1,
                           winsorize='none', clip_range=(-3, 3), missing_value='ignore'):
    """对热图数据进行标准化和变换。data: (samples, genes)，原地修改并返回。"""
    # 1. 缺失值处理
    if missing_value == 'mean_fill':
        col_means = np.nanmean(data, axis=0)
        for j in range(data.shape[1]):
            mask = np.isnan(data[:, j])
            if mask.any():
                data[mask, j] = col_means[j]
    elif missing_value == 'zero_fill':
        data = np.nan_to_num(data, nan=0.0)

    # 2. Winsorize（截断极端值）
    if winsorize != 'none':
        if winsorize == 'custom':
            pct = 0.01  # 由调用方自行处理自定义百分位
        else:
            pct = float(winsorize.replace('pct', '')) / 100
        from scipy.stats.mstats import winsorize as sp_winsorize
        data = sp_winsorize(data, limits=[pct, pct], axis=0).data if hasattr(
            sp_winsorize(data, limits=[pct, pct], axis=0), 'data') else data

    # 3. 行标准化
    if row_scaling == 'zscore':
        mean = np.mean(data, axis=0)
        std = np.std(data, axis=0) + 1e-10
        data = (data - mean) / std
    elif row_scaling == 'center':
        data = data - np.mean(data, axis=0)
    # 'none': 不做变换

    # 4. 截断
    if clip_range is not None:
        data = np.clip(data, clip_range[0], clip_range[1])

    return data
```

- [ ] **Step 4: 运行测试确认通过**

```bash
cd /data/GJ/platform && python -m pytest tests/test_heatmap_helpers.py::TestTransformHeatmapData -v
```

Expected: `8 passed`

- [ ] **Step 5: 提交**

```bash
git add modules/visualization.py tests/test_heatmap_helpers.py
git commit -m "feat(viz): add transform_heatmap_data helper (zscore/center/winsorize/missing)"
```

---

## Task 3: 添加 cluster_heatmap 工具函数

**Files:**
- Modify: `modules/visualization.py`
- Modify: `tests/test_heatmap_helpers.py`

- [ ] **Step 1: 编写失败测试**

```python
# 在 tests/test_heatmap_helpers.py 中追加
from modules.visualization import cluster_heatmap


class TestClusterHeatmap:
    @pytest.fixture
    def cluster_data(self):
        # 6 samples, 4 genes — 两组明显分离
        return np.array([
            [1, 1, -1, -1],  # group A
            [1, 1, -1, -1],
            [1, 1, -1, -1],
            [-1, -1, 1, 1],  # group B
            [-1, -1, 1, 1],
            [-1, -1, 1, 1],
        ], dtype=float)

    def test_euclidean_ward(self, cluster_data):
        order = cluster_heatmap(cluster_data, method='ward', metric='euclidean')
        assert len(order) == 6
        assert set(order) == {0, 1, 2, 3, 4, 5}
        # group A (0,1,2) 应该相邻，group B (3,4,5) 应该相邻
        a_positions = [order.index(i) for i in [0, 1, 2]]
        b_positions = [order.index(i) for i in [3, 4, 5]]
        assert max(a_positions) - min(a_positions) <= 2
        assert max(b_positions) - min(b_positions) <= 2

    def test_pearson_metric(self, cluster_data):
        order = cluster_heatmap(cluster_data, method='complete', metric='pearson')
        assert len(order) == 6
        assert set(order) == {0, 1, 2, 3, 4, 5}

    def test_cosine_metric(self, cluster_data):
        order = cluster_heatmap(cluster_data, method='average', metric='cosine')
        assert len(order) == 6

    def test_single_sample(self):
        data = np.array([[1, 2, 3]])
        order = cluster_heatmap(data)
        assert order == [0]

    def test_two_samples(self):
        data = np.array([[1, 2], [3, 4]])
        order = cluster_heatmap(data)
        assert len(order) == 2
```

- [ ] **Step 2: 运行测试确认失败**

```bash
cd /data/GJ/platform && python -m pytest tests/test_heatmap_helpers.py::TestClusterHeatmap -v
```

Expected: `ImportError`

- [ ] **Step 3: 实现函数**

在 `modules/visualization.py` 的 `transform_heatmap_data` 之后添加：

```python
def cluster_heatmap(data, method='ward', metric='euclidean'):
    """层次聚类，返回排序后的行索引列表。data: (n_samples, n_features)"""
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist

    if data.shape[0] <= 1:
        return list(range(data.shape[0]))

    if metric == 'pearson':
        dist = pdist(data, metric=lambda u, v: 1 - np.corrcoef(u, v)[0, 1])
    elif metric == 'spearman':
        from scipy.stats import spearmanr
        dist = pdist(data, metric=lambda u, v: 1 - spearmanr(u, v).correlation)
    elif metric == 'cosine':
        dist = pdist(data, metric='cosine')
    else:
        dist = pdist(data, metric='euclidean')

    link = linkage(dist, method=method)
    dendro = dendrogram(link, no_plot=True)
    return dendro['leaves']
```

- [ ] **Step 4: 运行测试确认通过**

```bash
cd /data/GJ/platform && python -m pytest tests/test_heatmap_helpers.py::TestClusterHeatmap -v
```

Expected: `5 passed`

- [ ] **Step 5: 提交**

```bash
git add modules/visualization.py tests/test_heatmap_helpers.py
git commit -m "feat(viz): add cluster_heatmap helper (ward/complete/average + euclidean/pearson/cosine)"
```

---

## Task 4: 添加 build_annotation_bar 和 save_plotly_json 工具函数

**Files:**
- Modify: `modules/visualization.py`
- Modify: `tests/test_heatmap_helpers.py`

- [ ] **Step 1: 编写失败测试**

```python
# 在 tests/test_heatmap_helpers.py 中追加
import pandas as pd
from modules.visualization import build_annotation_bar, save_plotly_json
import tempfile, json


class TestBuildAnnotationBar:
    def test_single_column(self):
        obs = pd.DataFrame({'group': ['A', 'A', 'B', 'B']}, index=['s1', 's2', 's3', 's4'])
        fig_data, unique_groups = build_annotation_bar(obs, ['group'], sample_order=[0, 1, 2, 3])
        assert 'group' in fig_data
        assert len(unique_groups['group']) == 2

    def test_multiple_columns(self):
        obs = pd.DataFrame({
            'group': ['A', 'A', 'B', 'B'],
            'batch': ['b1', 'b2', 'b1', 'b2']
        }, index=['s1', 's2', 's3', 's4'])
        fig_data, unique_groups = build_annotation_bar(obs, ['group', 'batch'])
        assert 'group' in fig_data
        assert 'batch' in fig_data

    def test_empty_columns(self):
        obs = pd.DataFrame({'group': ['A', 'B']}, index=['s1', 's2'])
        fig_data, unique_groups = build_annotation_bar(obs, [])
        assert len(fig_data) == 0


class TestSavePlotlyJson:
    def test_save_creates_file(self):
        import plotly.graph_objects as go
        fig = go.Figure(go.Heatmap(z=[[1, 2], [3, 4]]))
        with tempfile.TemporaryDirectory() as tmpdir:
            result_files = []
            save_plotly_json(fig, tmpdir, 'test_heatmap.json', result_files,
                            category='heatmap', label='测试热图')
            assert os.path.exists(os.path.join(tmpdir, 'test_heatmap.json'))
            assert len(result_files) == 1
            assert result_files[0]['category'] == 'heatmap'
            assert result_files[0]['label'] == '测试热图'

    def test_save_valid_json(self):
        import plotly.graph_objects as go
        fig = go.Figure(go.Heatmap(z=[[1, 2], [3, 4]]))
        with tempfile.TemporaryDirectory() as tmpdir:
            result_files = []
            save_plotly_json(fig, tmpdir, 'test.json', result_files)
            with open(os.path.join(tmpdir, 'test.json')) as f:
                data = json.load(f)
            assert 'data' in data
```

- [ ] **Step 2: 运行测试确认失败**

```bash
cd /data/GJ/platform && python -m pytest tests/test_heatmap_helpers.py::TestBuildAnnotationBar tests/test_heatmap_helpers.py::TestSavePlotlyJson -v
```

Expected: `ImportError`

- [ ] **Step 3: 实现函数**

在 `modules/visualization.py` 的 `cluster_heatmap` 之后添加：

```python
DEFAULT_PALETTE = ['#1a237e', '#e53935', '#4caf50', '#ff9800', '#9c27b0',
                   '#00bcd4', '#795548', '#607d8b', '#f44336', '#3f51b5']


def build_annotation_bar(obs, columns, sample_order=None, palette=None):
    """
    为多个注释列生成颜色映射数据。

    Returns:
        fig_data: dict {col_name: {'colors': [...], 'groups': [...], 'unique': [...]}}
        unique_groups: dict {col_name: [unique_values]}
    """
    if palette is None:
        palette = DEFAULT_PALETTE
    if sample_order is None:
        sample_order = list(range(len(obs)))

    fig_data = {}
    unique_groups = {}
    for col in columns:
        if col not in obs.columns:
            continue
        values = [str(obs[col].iloc[i]) for i in sample_order]
        uniq = sorted(set(values))
        color_map = {g: palette[i % len(palette)] for i, g in enumerate(uniq)}
        fig_data[col] = {
            'colors': [color_map[v] for v in values],
            'groups': values,
            'unique': uniq,
            'color_map': color_map,
        }
        unique_groups[col] = uniq
    return fig_data, unique_groups


def save_plotly_json(fig, plots_dir, filename, result_files,
                     file_type='plotly_json', category='heatmap', label=''):
    """保存 Plotly 图表为 JSON 并追加到 result_files 列表。"""
    import os
    fpath = os.path.join(plots_dir, filename)
    with open(fpath, 'w') as f:
        f.write(fig.to_json(engine="json"))
    result_files.append({
        'file_path': fpath, 'file_type': file_type,
        'category': category, 'label': label,
    })
```

- [ ] **Step 4: 运行测试确认通过**

```bash
cd /data/GJ/platform && python -m pytest tests/test_heatmap_helpers.py -v
```

Expected: `20 passed`

- [ ] **Step 5: 提交**

```bash
git add modules/visualization.py tests/test_heatmap_helpers.py
git commit -m "feat(viz): add build_annotation_bar and save_plotly_json helpers"
```

---

## Task 5: 更新 PARAM_SCHEMAS

**Files:**
- Modify: `routes/analysis.py:172-177`

- [ ] **Step 1: 替换 bulk_heatmap 的 PARAM_SCHEMAS**

找到 `routes/analysis.py` 中 `'bulk_heatmap': [` 开始的列表（约第 172 行），替换整个列表为：

```python
    'bulk_heatmap': [
        # === 基因选择 ===
        {'key': 'gene_import_source', 'label': '基因来源', 'type': 'select',
         'options': ['top_var', 'deg', 'expression_filter', 'manual'], 'default': 'top_var',
         'help': '基因来源。top_var：最高变异基因。deg：差异表达基因。expression_filter：表达式筛选。manual：手动输入。'},
        {'key': 'var_metric', 'label': '变异度量', 'type': 'select',
         'options': ['var', 'mad', 'cv', 'range'], 'default': 'var',
         'help': 'top_var 模式下的变异度量。var=方差，mad=中位绝对偏差，cv=变异系数，range=极差。'},
        {'key': 'deg_direction', 'label': 'DEG 方向过滤', 'type': 'select',
         'options': ['both', 'up', 'down'], 'default': 'both',
         'help': '仅显示上调/下调/全部差异基因。'},
        {'key': 'deg_sortby', 'label': 'DEG 排序方式', 'type': 'select',
         'options': ['padj', 'abs_logfc', 'logfc'], 'default': 'padj',
         'help': '差异基因在热图中的排列顺序。'},
        {'key': 'filter_expression', 'label': '表达式筛选', 'type': 'textarea', 'default': '',
         'help': '集合逻辑表达式（如 ALL:up），仅 gene_import_source=expression_filter 时生效。'},
        {'key': 'top_n', 'label': '显示基因数', 'type': 'number', 'default': 50, 'step': 5,
         'help': '热图中显示的基因数量。通常 30-100。'},
        {'key': 'custom_genes', 'label': '自定义基因列表', 'type': 'textarea', 'default': '',
         'help': '手动输入基因名，逗号或换行分隔。gene_import_source=manual 时生效。'},
        # === 数据变换 ===
        {'key': 'row_scaling', 'label': '行标准化', 'type': 'select',
         'options': ['zscore', 'center', 'none'], 'default': 'zscore',
         'help': 'zscore=Z 分数标准化（推荐），center=仅中心化，none=不变换。'},
        {'key': 'pseudocount', 'label': '伪计数', 'type': 'number', 'default': 1, 'step': 0.1,
         'help': 'log2 转换前加值，避免 log(0)。'},
        {'key': 'winsorize', 'label': '极端值截断', 'type': 'select',
         'options': ['none', '1pct', '5pct'], 'default': 'none',
         'help': '对极端值进行百分位截断，避免 outliers 主导颜色。'},
        {'key': 'clip_range', 'label': '颜色截断范围', 'type': 'text', 'default': '-3,3',
         'help': '标准化后的值截断范围，逗号分隔（如 -3,3）。留空不截断。'},
        {'key': 'missing_value', 'label': '缺失值处理', 'type': 'select',
         'options': ['ignore', 'mean_fill', 'zero_fill'], 'default': 'ignore',
         'help': '缺失值处理方式。ignore=忽略，mean_fill=行均值填充，zero_fill=填 0。'},
        # === 聚类 ===
        {'key': 'row_cluster', 'label': '基因聚类', 'type': 'select',
         'options': ['yes', 'no'], 'default': 'yes',
         'help': '是否对基因进行层次聚类。'},
        {'key': 'col_cluster', 'label': '样本聚类', 'type': 'select',
         'options': ['yes', 'no', 'group_order'], 'default': 'yes',
         'help': '样本聚类。group_order=按分组固定顺序排列。'},
        {'key': 'row_method', 'label': '基因聚类方法', 'type': 'select',
         'options': ['ward', 'complete', 'average', 'single', 'mcquitty'], 'default': 'ward',
         'help': '层次聚类算法。ward=最小方差（推荐），complete=最长距离，average=平均距离。'},
        {'key': 'col_method', 'label': '样本聚类方法', 'type': 'select',
         'options': ['ward', 'complete', 'average', 'single', 'mcquitty'], 'default': 'ward',
         'help': '样本层次聚类算法。'},
        {'key': 'row_metric', 'label': '基因距离度量', 'type': 'select',
         'options': ['euclidean', 'pearson', 'spearman', 'cosine'], 'default': 'euclidean',
         'help': '基因间距离计算方式。euclidean=欧氏距离，pearson=相关距离。'},
        {'key': 'col_metric', 'label': '样本距离度量', 'type': 'select',
         'options': ['euclidean', 'pearson', 'spearman', 'cosine'], 'default': 'euclidean',
         'help': '样本间距离计算方式。'},
        # === 视觉样式 ===
        {'key': 'colorscale', 'label': '颜色方案', 'type': 'select',
         'options': ['RdBu_r', 'RdYlBu', 'viridis', 'magma', 'Blues', 'PiYG'], 'default': 'RdBu_r',
         'help': '热图颜色方案。RdBu_r=红蓝发散（推荐），viridis=感知均匀。'},
        {'key': 'reverse_color', 'label': '反转颜色', 'type': 'checkbox', 'default': False,
         'help': '反转颜色映射方向。'},
        {'key': 'zmin', 'label': '最小值', 'type': 'text', 'default': 'auto',
         'help': '颜色映射最小值。auto=自动，或输入数字（如 -3）。'},
        {'key': 'zmax', 'label': '最大值', 'type': 'text', 'default': 'auto',
         'help': '颜色映射最大值。auto=自动，或输入数字（如 3）。'},
        {'key': 'show_gene_labels', 'label': '基因名显示', 'type': 'select',
         'options': ['all', 'top20', 'none'], 'default': 'all',
         'help': '基因名标签显示方式。'},
        {'key': 'show_sample_labels', 'label': '样本名显示', 'type': 'select',
         'options': ['all', 'none'], 'default': 'all',
         'help': '样本名标签显示方式。'},
        {'key': 'gene_font_size', 'label': '基因名字体大小', 'type': 'number', 'default': 8, 'step': 1,
         'help': '基因名标签字体大小。'},
        {'key': 'sample_font_size', 'label': '样本名字体大小', 'type': 'number', 'default': 9, 'step': 1,
         'help': '样本名标签字体大小。'},
        {'key': 'annotation_columns', 'label': '注释条列名', 'type': 'text', 'default': '',
         'help': '额外注释条列名，逗号分隔（如 group,batch）。groupby 列自动包含。'},
        {'key': 'groupby', 'label': '样本分组列名', 'type': 'text', 'default': '',
         'help': '主分组列名，用于默认注释条。留空则不添加。'},
    ],
```

- [ ] **Step 2: 验证导入无报错**

```bash
cd /data/GJ/platform && python -c "from routes.analysis import PARAM_SCHEMAS; print(len(PARAM_SCHEMAS['bulk_heatmap']), 'params'); print('OK')"
```

Expected: `29 params` (or similar count) and `OK`

- [ ] **Step 3: 提交**

```bash
git add routes/analysis.py
git commit -m "feat(heatmap): expand PARAM_SCHEMAS with 27 new parameters for heatmap enhancement"
```

---

## Task 6: 重构 bulk_heatmap.py — 基因选择部分

**Files:**
- Modify: `modules/bulk_heatmap.py`

- [ ] **Step 1: 重写 run() 的基因选择逻辑**

替换 `bulk_heatmap.py` 中 `self.progress(40, "选择基因...")` 到 `title = ...` 之间的代码（约第 48-70 行）为：

```python
        self.progress(40, "选择基因...")

        gene_import_source = self.params.get('gene_import_source',
                                              'manual' if self.params.get('custom_genes', '').strip() else
                                              self.params.get('heatmap_type', 'top_var'))
        top_n = int(self.params.get('top_n', 50))
        custom_genes_str = self.params.get('custom_genes', '').strip()
        var_metric = self.params.get('var_metric', 'var')

        if gene_import_source == 'manual' or custom_genes_str:
            # 手动输入基因列表
            gene_list = [g.strip() for g in custom_genes_str.replace('\n', ',').split(',') if g.strip()]
            var_names_list = list(adata.var_names)
            top_idx = [var_names_list.index(g) for g in gene_list if g in var_names_list]
            not_found = [g for g in gene_list if g not in var_names_list]
            if not top_idx:
                raise ValueError(f"自定义基因列表中没有找到任何匹配基因。请检查基因名是否正确。")
            title = f'自定义基因热图 ({len(top_idx)} genes)'
            if not_found:
                title += f'，{len(not_found)} 个未找到'

        elif gene_import_source == 'deg':
            # 从 DEG 结果选择基因
            deg_direction = self.params.get('deg_direction', 'both')
            deg_sortby = self.params.get('deg_sortby', 'padj')
            results_dir = os.path.join(self.project_dir, 'results')
            deg_files = sorted([f for f in os.listdir(results_dir)
                                if f.startswith('bulk_deg_results') and f.endswith('.csv')
                                and 'merged' not in f and 'all_comparisons' not in f
                                and 'lrt' not in f and 'top_genes' not in f])
            if not deg_files:
                raise ValueError("未找到 DEG 结果文件，请先运行 bulk_deg")
            deg_df = pd.read_csv(os.path.join(results_dir, deg_files[0]))
            if deg_direction == 'up':
                deg_df = deg_df[deg_df['regulation'] == 'Up']
            elif deg_direction == 'down':
                deg_df = deg_df[deg_df['regulation'] == 'Down']
            else:
                deg_df = deg_df[deg_df['regulation'] != 'NS']
            if deg_sortby == 'padj':
                deg_df = deg_df.sort_values('padj')
            elif deg_sortby == 'abs_logfc':
                deg_df = deg_df.sort_values('log2FC', key=abs, ascending=False)
            elif deg_sortby == 'logfc':
                deg_df = deg_df.sort_values('log2FC', ascending=False)
            gene_list = deg_df['gene'].head(top_n).tolist()
            top_idx = [i for i, g in enumerate(adata.var_names) if g in set(gene_list)]
            title = f'Top {len(top_idx)} 差异基因热图 ({deg_direction})'

        elif gene_import_source == 'expression_filter':
            # 从表达式筛选结果导入基因
            from modules.expression_parser import validate, evaluate
            filter_expr = self.params.get('filter_expression', '').strip()
            if not filter_expr:
                raise ValueError("expression_filter 模式需要填写筛选表达式")
            results_dir = os.path.join(self.project_dir, 'results')
            deg_files = sorted([f for f in os.listdir(results_dir)
                                if f.startswith('bulk_deg_results') and f.endswith('.csv')
                                and 'merged' not in f and 'all_comparisons' not in f
                                and 'lrt' not in f and 'top_genes' not in f])
            if not deg_files:
                raise ValueError("未找到 DEG 结果文件，请先运行 bulk_deg")
            # 构建 DEG 矩阵
            comparisons = {}
            for f in deg_files:
                df = pd.read_csv(os.path.join(results_dir, f))
                if 'gene' in df.columns and 'log2FC' in df.columns:
                    name = f.replace('bulk_deg_', '').replace('.csv', '')
                    comparisons[name] = df
            comp_names = sorted(comparisons.keys())
            all_genes_union = set()
            for df in comparisons.values():
                all_genes_union.update(df['gene'].tolist())
            all_genes_union = sorted(all_genes_union)
            logfc_m = pd.DataFrame(0.0, index=all_genes_union, columns=comp_names)
            padj_m = pd.DataFrame(1.0, index=all_genes_union, columns=comp_names)
            fc_thresh = float(self.params.get('fc_threshold', 2.0))
            pv_thresh = float(self.params.get('pval_threshold', 0.05))
            import numpy as _np
            for name, df in comparisons.items():
                df_idx = df.set_index('gene')
                common = pd.Index(all_genes_union).intersection(df_idx.index)
                if len(common) > 0:
                    logfc_m.loc[common, name] = df_idx.loc[common, 'log2FC']
                    padj_m.loc[common, name] = df_idx.loc[common, 'padj']
            log2fc_t = _np.log2(fc_thresh)
            gene_sets = {}
            for c in comp_names:
                sig = (padj_m[c] < pv_thresh) & (abs(logfc_m[c]) >= log2fc_t)
                gene_sets[c] = set(logfc_m.index[sig])
            ast, err, _ = validate(filter_expr, comp_names)
            if err:
                raise ValueError(f"表达式错误: {err}")
            filtered = evaluate(ast, comp_names, gene_sets, padj_m, logfc_m, pv_thresh, log2fc_t)
            filtered = sorted(filtered)[:top_n]
            top_idx = [i for i, g in enumerate(adata.var_names) if g in set(filtered)]
            title = f'筛选基因热图 ({len(top_idx)} genes)'

        else:
            # top_var: 按变异度量选择
            from modules.visualization import compute_gene_variability
            gene_scores = compute_gene_variability(norm_data, metric=var_metric)
            top_idx = np.argsort(gene_scores)[::-1][:top_n]
            metric_names = {'var': '方差', 'mad': 'MAD', 'cv': '变异系数', 'range': '极差'}
            title = f'Top {top_n} 高变异基因热图 ({metric_names.get(var_metric, var_metric)})'
```

- [ ] **Step 2: 运行 Python 语法检查**

```bash
cd /data/GJ/platform && python -c "from modules.bulk_heatmap import BulkHeatmapAnalysis; print('OK')"
```

Expected: `OK`

- [ ] **Step 3: 提交**

```bash
git add modules/bulk_heatmap.py
git commit -m "feat(heatmap): refactor gene selection with deg/expression_filter/var_metric support"
```

---

## Task 7: 重构 bulk_heatmap.py — 数据变换 + 聚类 + 视觉样式

**Files:**
- Modify: `modules/bulk_heatmap.py`

- [ ] **Step 1: 替换数据变换和聚类逻辑**

替换 `heat_data = norm_data[:, top_idx]` 到聚类完成之间的代码（约第 72-100 行）为：

```python
        heat_data = norm_data[:, top_idx].copy()
        gene_labels = [adata.var_names[i] for i in top_idx]
        sample_labels = adata.obs.index.tolist()

        # 数据变换
        from modules.visualization import transform_heatmap_data, cluster_heatmap
        row_scaling = self.params.get('row_scaling', 'zscore')
        pseudocount = float(self.params.get('pseudocount', 1))
        winsorize = self.params.get('winsorize', 'none')
        missing_value = self.params.get('missing_value', 'ignore')
        clip_str = self.params.get('clip_range', '-3,3').strip()
        clip_range = None
        if clip_str:
            parts = [float(x.strip()) for x in clip_str.split(',') if x.strip()]
            if len(parts) == 2:
                clip_range = (parts[0], parts[1])

        # pseudocount: 对原始计数加偏移后 log2 转换（仅在数据未经 log 变换时生效）
        if pseudocount != 1 and 'normalization' not in adata.uns:
            heat_data = np.log2(heat_data + pseudocount)

        heat_z = transform_heatmap_data(
            heat_data, row_scaling=row_scaling,
            winsorize=winsorize, clip_range=clip_range, missing_value=missing_value)

        # 聚类
        row_cluster = self.params.get('row_cluster', 'yes')
        col_cluster = self.params.get('col_cluster', 'yes')
        row_method = self.params.get('row_method', 'ward')
        col_method = self.params.get('col_method', 'ward')
        row_metric = self.params.get('row_metric', 'euclidean')
        col_metric = self.params.get('col_metric', 'euclidean')

        if col_cluster == 'yes' and adata.n_obs > 2:
            sample_order = cluster_heatmap(heat_z, method=col_method, metric=col_metric)
        elif col_cluster == 'group_order' and groupby and groupby in adata.obs.columns:
            # 按分组排序
            groups = adata.obs[groupby].astype(str)
            group_order = sorted(groups.unique())
            sample_order = []
            for g in group_order:
                sample_order.extend([i for i in range(len(groups)) if groups.iloc[i] == g])
        else:
            sample_order = list(range(adata.n_obs))

        if row_cluster == 'yes' and len(top_idx) > 2:
            gene_order = cluster_heatmap(heat_z.T, method=row_method, metric=row_metric)
        else:
            gene_order = list(range(len(top_idx)))

        heat_ordered = heat_z[np.ix_(sample_order, gene_order)]
        sample_ordered = [sample_labels[i] for i in sample_order]
        gene_ordered = [gene_labels[i] for i in gene_order]
```

- [ ] **Step 2: 替换热图渲染逻辑**

替换 `self.progress(75, "生成热图...")` 到热图保存之间的代码（约第 102-144 行）为：

```python
        self.progress(75, "生成热图...")

        # 视觉样式参数
        colorscale = self.params.get('colorscale', 'RdBu_r')
        reverse_color = self.params.get('reverse_color', False)
        if reverse_color:
            colorscale = colorscale + '_r' if not colorscale.endswith('_r') else colorscale[:-2]

        zmin_str = self.params.get('zmin', 'auto').strip()
        zmax_str = self.params.get('zmax', 'auto').strip()
        zmin = float(zmin_str) if zmin_str and zmin_str != 'auto' else None
        zmax = float(zmax_str) if zmax_str and zmax_str != 'auto' else None

        show_gene_labels = self.params.get('show_gene_labels', 'all')
        show_sample_labels = self.params.get('show_sample_labels', 'all')
        gene_font_size = int(self.params.get('gene_font_size', 8))
        sample_font_size = int(self.params.get('sample_font_size', 9))

        # 构建热图 trace
        heatmap_kwargs = dict(
            z=heat_ordered.tolist(),
            x=gene_ordered,
            y=sample_ordered,
            colorscale=colorscale,
            colorbar=dict(title='Z-score' if row_scaling == 'zscore' else 'Value'),
            hovertemplate='样本: %{y}<br>基因: %{x}<br>值: %{z:.2f}<extra></extra>'
        )
        if row_scaling in ('zscore', 'center') and zmin is None and zmax is None:
            heatmap_kwargs['zmid'] = 0
        if zmin is not None:
            heatmap_kwargs['zmin'] = zmin
        if zmax is not None:
            heatmap_kwargs['zmax'] = zmax

        fig = go.Figure()
        fig.add_trace(go.Heatmap(**heatmap_kwargs))

        # 标签显示控制
        xaxis_kwargs = dict(tickangle=45, tickfont=dict(size=gene_font_size))
        yaxis_kwargs = dict(tickfont=dict(size=sample_font_size))
        if show_gene_labels == 'top20':
            show_n = min(20, len(gene_ordered))
            xaxis_kwargs['tickvals'] = list(range(show_n))
            xaxis_kwargs['ticktext'] = gene_ordered[:show_n]
        elif show_gene_labels == 'none':
            xaxis_kwargs['showticklabels'] = False
        if show_sample_labels == 'none':
            yaxis_kwargs['showticklabels'] = False

        fig.update_layout(
            title=title,
            xaxis=xaxis_kwargs,
            yaxis=yaxis_kwargs,
            height=max(400, len(sample_ordered) * 25 + 150),
            width=max(600, len(gene_ordered) * 12 + 200),
            plot_bgcolor='white'
        )

        from modules.visualization import save_plotly_json
        save_plotly_json(fig, plots_dir, 'bulk_heatmap.json', result_files,
                        category='heatmap', label=title)
```

- [ ] **Step 3: 替换注释条生成逻辑**

替换注释条生成代码（原 `if annotation_colors:` 块，约第 147-165 行）为：

```python
        # 注释条（支持多列）
        from modules.visualization import build_annotation_bar, DEFAULT_PALETTE
        annot_cols_str = self.params.get('annotation_columns', '').strip()
        annot_cols = [c.strip() for c in annot_cols_str.split(',') if c.strip()]
        if groupby and groupby in adata.obs.columns and groupby not in annot_cols:
            annot_cols.insert(0, groupby)

        if annot_cols:
            annot_data, _ = build_annotation_bar(adata.obs, annot_cols,
                                                 sample_order=sample_order, palette=DEFAULT_PALETTE)
            for col_name, col_info in annot_data.items():
                color_indices = [col_info['unique'].index(g) for g in col_info['groups']]
                fig_annot = go.Figure()
                cs = [[i / max(len(col_info['unique']) - 1, 1), col_info['color_map'][g]]
                      for i, g in enumerate(col_info['unique'])]
                if len(col_info['unique']) == 1:
                    cs = [[0, list(col_info['color_map'].values())[0]]]
                fig_annot.add_trace(go.Heatmap(
                    z=[[i] for i in color_indices],
                    y=sample_ordered, x=[col_name],
                    colorscale=cs, showscale=False,
                    text=[[col_info['groups'][i]] for i in range(len(col_info['groups']))],
                    hovertemplate='%{y}: %{text}<extra></extra>'
                ))
                fig_annot.update_layout(
                    height=max(400, len(sample_ordered) * 25 + 150), width=100,
                    margin=dict(l=0, r=0, t=30, b=40)
                )
                save_plotly_json(fig_annot, plots_dir, f'bulk_heatmap_annotation_{col_name}.json',
                                result_files, category='annotation', label=f'{col_name} 注释条')
```

- [ ] **Step 4: 验证完整模块可导入**

```bash
cd /data/GJ/platform && python -c "from modules.bulk_heatmap import BulkHeatmapAnalysis; print('DISPLAY_NAME:', BulkHeatmapAnalysis.DISPLAY_NAME); print('OK')"
```

Expected: `DISPLAY_NAME: Bulk 热图可视化` and `OK`

- [ ] **Step 5: 运行所有测试**

```bash
cd /data/GJ/platform && python -m pytest tests/ -v 2>&1 | tail -15
```

Expected: `all passed`

- [ ] **Step 6: 提交**

```bash
git add modules/bulk_heatmap.py
git commit -m "feat(heatmap): refactor data transform, clustering, colorscale, multi-annotation bars"
```

---

## Task 8: 更新 README 并最终验证

**Files:**
- Modify: `README.md`

- [ ] **Step 1: 更新 README 中热图模块描述**

将 `README.md` 中热图分析行更新为：

```
| 热图分析 | Top 变异/差异基因热图、多种变异度量(MAD/CV)、z-score/中心化、灵活聚类(pearson/cosine)、多样式配置、表达式筛选导入 | scipy |
```

- [ ] **Step 2: 运行全部测试**

```bash
cd /data/GJ/platform && python -m pytest tests/ -v 2>&1 | tail -10
```

Expected: `all passed`

- [ ] **Step 3: 提交并推送**

```bash
git add README.md
git commit -m "docs: update README for heatmap enhancement features"
git push origin feature/bulk-rnaseq-sc-enhancement
```
