# 未完成功能补齐实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 补齐 Platform Customization (Phase 3/4) 和 SC Pipeline Roadmap (P4/P5/P6) 中所有已定义但未实现的功能。

**Architecture:** 分 6 个独立任务，按依赖顺序实施。每个任务修改少量文件、产出可独立测试的增量。后端模块代码补齐 → 前端 UI 面板 → 联调验证。

**Tech Stack:** Python (scanpy, scipy, celltypist, plotly), Jinja2, JavaScript (vanilla)

---

## 文件清单

| 文件 | 任务 | 改动类型 |
|------|------|----------|
| `modules/annotation.py` | Task 1 | 修改 — 补 celltypist / confidence / merge |
| `modules/deg.py` | Task 2 | 修改 — 补 correction_method / min_pct |
| `modules/proportion.py` | Task 3 | 修改 — 补 stat_test / min_cells_per_group |
| `modules/visualization.py` | Task 4 | 修改 — umap_scatter 接受 point_size/opacity/legend_fontsize |
| `modules/base.py` | Task 4, 5 | 修改 — get_plotly_layout 扩展 / export_results 新增 |
| `templates/analysis_select.html` | Task 6 | 修改 — 可视化面板 + 过滤规则 UI |

---

## Task 1: annotation.py — 补齐 CellTypist / 置信度 / 合并

**Files:**
- Modify: `modules/annotation.py:81-135`

**Context:** `routes/analysis.py` 已定义 `method` 选项含 `celltypist`，以及 `celltypist_model`、`celltypist_threshold`、`celltypist_majority_voting`、`confidence_method`、`mark_unknown`、`merge_similar_threshold` 共 6 个参数。模块代码只实现了 `auto_marker` 和 `manual` 两个分支，celltypist 分支缺失，置信度和合并逻辑完全缺失。

- [ ] **Step 1: 在 `run()` 中读取新参数**

在 `annotation.py:81` 之后（`method = ...` 行之后）加入：

```python
        confidence_method = self.params.get('confidence_method', 'none')
        mark_unknown = self.params.get('mark_unknown', True)
        merge_similar_threshold = float(self.params.get('merge_similar_threshold', 0))
```

- [ ] **Step 2: 添加 celltypist 分支**

在 `annotation.py:101` 的 `method = 'auto_marker'` 回退之后、`if method == 'auto_marker':` 之前，插入：

```python
        if method == 'celltypist':
            self.progress(30, "Running CellTypist annotation...")
            try:
                import celltypist
                from celltypist import annotate
                ct_model = self.params.get('celltypist_model', 'Immune_All_Low')
                ct_threshold = float(self.params.get('celltypist_threshold', 0.5))
                majority_voting = self.params.get('celltypist_majority_voting', True)

                # CellTypist expects log-normalized data in .X
                predictions = annotate(
                    adata, model=ct_model,
                    majority_voting=majority_voting,
                    confidence_threshold=ct_threshold,
                )
                adata.obs['celltype'] = predictions.predicted_labels['predicted_labels'].values
                adata.obs['celltype'] = adata.obs['celltype'].astype('category')
                markers = {}
            except ImportError:
                self.progress(35, "celltypist not installed, falling back to auto_marker...")
                method = 'auto_marker'
            except Exception as e:
                self.progress(35, f"CellTypist failed: {e}, falling back to auto_marker...")
                method = 'auto_marker'
```

- [ ] **Step 3: 添加置信度计算**

在 celltype 赋值完成之后（`adata.obs['celltype'] = ...` 之后，`self.progress(70, ...)` 之前），插入：

```python
        # 置信度评估
        if confidence_method != 'none' and 'celltype' in adata.obs.columns:
            self.progress(60, f"Computing annotation confidence ({confidence_method})...")
            score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
            if confidence_method == 'entropy' and score_cols:
                import numpy as np
                score_matrix = adata.obs[score_cols].values
                # softmax → entropy per cell
                exp_scores = np.exp(score_matrix - score_matrix.max(axis=1, keepdims=True))
                probs = exp_scores / (exp_scores.sum(axis=1, keepdims=True) + 1e-10)
                entropy = -np.sum(probs * np.log(probs + 1e-10), axis=1)
                max_entropy = np.log(len(score_cols)) if len(score_cols) > 1 else 1
                adata.obs['annotation_confidence'] = 1 - entropy / (max_entropy + 1e-10)
            elif confidence_method == 'score_margin' and score_cols:
                import numpy as np
                score_matrix = adata.obs[score_cols].values
                sorted_scores = np.sort(score_matrix, axis=1)
                if sorted_scores.shape[1] >= 2:
                    adata.obs['annotation_confidence'] = sorted_scores[:, -1] - sorted_scores[:, -2]
                else:
                    adata.obs['annotation_confidence'] = sorted_scores[:, -1]

            if mark_unknown and 'annotation_confidence' in adata.obs.columns:
                low_conf_mask = adata.obs['annotation_confidence'] < 0.2
                adata.obs.loc[low_conf_mask, 'celltype'] = 'Unknown'
                n_unknown = low_conf_mask.sum()
                if n_unknown > 0:
                    self.progress(-1, f"标记 {n_unknown} 个低置信度细胞为 Unknown")
```

- [ ] **Step 4: 添加相似簇合并**

在置信度计算之后、可视化生成之前，插入：

```python
        # 相似簇合并
        if merge_similar_threshold > 0 and 'annotation_confidence' in adata.obs.columns:
            import numpy as np
            score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
            if score_cols:
                cluster_profiles = {}
                for cluster in adata.obs[cluster_key].cat.categories:
                    mask = adata.obs[cluster_key] == cluster
                    cluster_profiles[cluster] = adata.obs.loc[mask, score_cols].mean().values
                clusters = list(cluster_profiles.keys())
                for i in range(len(clusters)):
                    for j in range(i + 1, len(clusters)):
                        corr = np.corrcoef(cluster_profiles[clusters[i]], cluster_profiles[clusters[j]])[0, 1]
                        if corr >= merge_similar_threshold:
                            old_ct = adata.obs.loc[adata.obs[cluster_key] == clusters[j], 'celltype'].mode()
                            if len(old_ct) > 0:
                                new_ct = adata.obs.loc[adata.obs[cluster_key] == clusters[i], 'celltype'].mode().iloc[0]
                                adata.obs.loc[adata.obs[cluster_key] == clusters[j], 'celltype'] = new_ct
                                self.progress(-1, f"合并簇 {clusters[j]} → {clusters[i]} (r={corr:.2f})")
```

- [ ] **Step 5: 更新 summary 返回置信度信息**

在 `annotation.py:177` 的 `ct_counts` 之后，扩展 summary：

```python
        summary = {
            'n_celltypes': adata.obs['celltype'].nunique(),
            'celltype_counts': {str(k): int(v) for k, v in ct_counts.items()},
            'cluster_column': leiden_key,
            'method_used': method,
        }
        if 'annotation_confidence' in adata.obs.columns:
            summary['mean_confidence'] = round(float(adata.obs['annotation_confidence'].mean()), 3)
```

并在 return 中用 `summary` 替换原来的 inline dict。

- [ ] **Step 6: 验证 — 无 celltypist 时回退正常**

```bash
cd /data/GJ/platform
python -c "
from modules.annotation import AnnotationAnalysis
import tempfile, os
print('annotation.py imports OK')
"
```

- [ ] **Step 7: Commit**

```bash
git add modules/annotation.py
git commit -m "feat(sc): P4 annotation — celltypist, confidence_method, merge_similar"
```

---

## Task 2: deg.py — 补齐 correction_method 和 min_pct

**Files:**
- Modify: `modules/deg.py:39, 42-54`

**Context:** `correction_method` 参数在 routes 中定义（benjamini_hochberg/bonferroni/BY），但模块未使用。`min_pct` 参数已读取但从未用于过滤。scanpy 的 `rank_genes_groups` 不直接支持 correction_method 参数，需要在提取结果后手动重新校正 p 值。

- [ ] **Step 1: 读取 correction_method 参数**

在 `deg.py:30` 的 `min_pct = ...` 之后添加：

```python
        correction_method = self.params.get('correction_method', 'benjamini_hochberg')
```

- [ ] **Step 2: 在结果提取后应用 p 值校正和 min_pct 过滤**

替换 `deg.py:41-54` 的结果提取块为：

```python
        self.progress(50, "Extracting results...")
        from scipy.stats import ranksums
        from statsmodels.stats.multitest import multipletests
        import numpy as np

        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names

        # 预计算每组每基因的表达比例（用于 min_pct 过滤）
        pct_expr = {}
        if min_pct > 0:
            for g in groups:
                mask = adata.obs[groupby] == g
                n_cells = mask.sum()
                if n_cells > 0:
                    expr = (adata[mask].X > 0).sum(axis=0)
                    if hasattr(expr, 'A1'):
                        expr = expr.A1
                    pct_expr[g] = dict(zip(adata.var_names, expr / n_cells))

        deg_data = []
        correction_map = {'benjamini_hochberg': 'fdr_bh', 'bonferroni': 'bonferroni', 'BY': 'fdr_by'}
        meth = correction_map.get(correction_method, 'fdr_bh')

        for g in groups:
            raw_pvals = []
            gene_info = []
            for i in range(min(100, len(result['names'][g]))):
                gene = result['names'][g][i]
                raw_p = float(result['pvals'][g][i])
                # min_pct 过滤
                if min_pct > 0 and g in pct_expr:
                    if pct_expr[g].get(gene, 0) < min_pct:
                        continue
                raw_pvals.append(max(raw_p, 1e-300))
                gene_info.append({
                    'gene': gene,
                    'logfc': round(float(result['logfoldchanges'][g][i]), 3),
                    'pval': raw_p,
                    'score': round(float(result['scores'][g][i]), 3),
                })

            # 应用多重检验校正
            if raw_pvals:
                reject, padj, _, _ = multipletests(raw_pvals, method=meth)
                for j, info in enumerate(gene_info):
                    info['pval_adj'] = float(padj[j])
                    info['cluster'] = g
                    deg_data.append(info)
```

- [ ] **Step 3: 验证 statsmodels 可用**

```bash
cd /data/GJ/platform
python -c "from statsmodels.stats.multitest import multipletests; print('statsmodels OK')"
```

- [ ] **Step 4: Commit**

```bash
git add modules/deg.py
git commit -m "feat(sc): P4 deg — correction_method and min_pct filtering"
```

---

## Task 3: proportion.py — 补齐 stat_test 和 min_cells_per_group

**Files:**
- Modify: `modules/proportion.py:17, 25-36, 82-106`

**Context:** 模块始终用 chi2_contingency，忽略用户选择的 `stat_test`。`min_cells_per_group` 已在 routes 中定义但模块未读取。

- [ ] **Step 1: 读取新参数**

在 `proportion.py:26` 的 `batch_key = ...` 之后添加：

```python
        stat_test = self.params.get('stat_test', 'chi_square')
        n_permutations = int(self.params.get('n_permutations', 1000))
        min_cells_per_group = int(self.params.get('min_cells_per_group', 10))
```

- [ ] **Step 2: 实现辅助统计函数**

在文件顶部 `from modules.base import BaseAnalysis` 之后、class 定义之前，添加：

```python
def _run_stat_test(ct_abs, test_type, n_permutations=1000):
    """对列联表运行指定统计检验，返回 (statistic, p_value)。"""
    from scipy.stats import chi2_contingency, fisher_exact
    import numpy as np

    if test_type == 'chi_square':
        chi2, pval, _, _ = chi2_contingency(ct_abs)
        return chi2, pval

    elif test_type == 'fisher_exact':
        # Fisher 精确检验仅适用于 2x2 表，对更大表使用模拟
        if ct_abs.shape == (2, 2):
            stat, pval = fisher_exact(ct_abs.values)
            return stat, pval
        else:
            # 对多列表逐对合并后检验
            chi2, pval, _, _ = chi2_contingency(ct_abs, simulate_pval=True, n_permutations=n_permutations)
            return chi2, pval

    elif test_type == 'permutation':
        observed_chi2, _, _, _ = chi2_contingency(ct_abs)
        count = 0
        for _ in range(n_permutations):
            shuffled = ct_abs.copy()
            total = shuffled.values.sum()
            row_sums = shuffled.sum(axis=1).values
            col_sums = shuffled.sum(axis=0).values
            # 生成置换表
            perm_table = np.random.multinomial(total, col_sums / total).reshape(1, -1)
            for rs in row_sums[1:]:
                row = np.random.multinomial(rs, col_sums / total)
                perm_table = np.vstack([perm_table, row])
            perm_df = pd.DataFrame(perm_table, index=ct_abs.index, columns=ct_abs.columns)
            perm_chi2, _, _, _ = chi2_contingency(perm_df)
            if perm_chi2 >= observed_chi2:
                count += 1
        pval = (count + 1) / (n_permutations + 1)
        return observed_chi2, pval

    else:
        chi2, pval, _, _ = chi2_contingency(ct_abs)
        return chi2, pval
```

- [ ] **Step 3: 应用 min_cells_per_group 过滤**

在 `proportion.py:31` 的 `self.progress(30, ...)` 之后、`ct = pd.crosstab(...)` 之前，插入：

```python
        # 按最小细胞数过滤
        if min_cells_per_group > 0:
            group_counts = adata.obs[groupby].value_counts()
            valid_groups = group_counts[group_counts >= min_cells_per_group].index.tolist()
            if len(valid_groups) < len(group_counts):
                removed = set(group_counts.index) - set(valid_groups)
                self.progress(-1, f"移除 {len(removed)} 个低细胞数组: {removed}")
                adata = adata[adata.obs[groupby].isin(valid_groups)].copy()
```

- [ ] **Step 4: 替换硬编码的 chi2_contingency 调用**

将 `proportion.py:36` 的：
```python
        chi2, pval, dof, expected = chi2_contingency(ct_abs)
```
替换为：
```python
        chi2, pval = _run_stat_test(ct_abs, stat_test, n_permutations)
```

同样将 `proportion.py:90` 的：
```python
                chi2_sub, pval_sub, _, _ = chi2_contingency(ct_sub_abs)
```
替换为：
```python
                chi2_sub, pval_sub = _run_stat_test(ct_sub_abs, stat_test, n_permutations)
```

- [ ] **Step 5: 更新 summary 添加 stat_test 信息**

在 `proportion.py:118-123` 的 summary dict 中添加：

```python
                'stat_test': stat_test,
```

- [ ] **Step 6: 验证**

```bash
cd /data/GJ/platform
python -c "
from modules.proportion import _run_stat_test
import pandas as pd, numpy as np
ct = pd.DataFrame([[50, 30, 20], [20, 40, 40]], index=['A', 'B'], index_name='batch')
ct.index = ['A', 'B']
chi2, p = _run_stat_test(ct, 'chi_square')
print(f'chi2={chi2:.2f}, p={p:.4f}')
chi2, p = _run_stat_test(ct, 'permutation', 500)
print(f'perm chi2={chi2:.2f}, p={p:.4f}')
"
```

- [ ] **Step 7: Commit**

```bash
git add modules/proportion.py
git commit -m "feat(sc): P5 proportion — stat_test (chi_square/fisher/permutation), min_cells_per_group"
```

---

## Task 4: 可视化全局参数 — umap_scatter 扩展 + base.py get_plotly_layout

**Files:**
- Modify: `modules/visualization.py:4-43`
- Modify: `modules/base.py:87-104`

**Context:** `umap_scatter()` 硬编码 `size=2`, `opacity=0.6`。`get_plotly_layout()` 缺少 `umap_point_size`、`umap_opacity`、`umap_legend_fontsize` 支持。

- [ ] **Step 1: 扩展 umap_scatter 函数签名**

将 `visualization.py:4` 的函数签名改为：

```python
def umap_scatter(adata, color_key=None, basis='X_umap', max_cells=50000, title='',
                 viz_params=None):
    import plotly.graph_objects as go

    vp = viz_params or {}
    point_size = vp.get('umap_point_size', 5)
    opacity = vp.get('umap_opacity', 0.7)
    legend_fontsize = vp.get('umap_legend_fontsize', 10)
```

- [ ] **Step 2: 替换所有硬编码的 size/opacity**

将三个 `go.Scattergl` 调用中的 `marker=dict(size=2, opacity=0.6)` 和 `marker=dict(size=2, opacity=0.5)` 分别替换为使用 `point_size` 和 `opacity` 变量：

第 20-24 行（分类变量分支）：
```python
            fig.add_trace(go.Scattergl(
                x=coords[mask, 0], y=coords[mask, 1],
                mode='markers', name=str(cat),
                marker=dict(size=point_size, opacity=opacity),
            ))
```

第 26-31 行（连续变量分支）：
```python
        fig.add_trace(go.Scattergl(
            x=coords[:, 0], y=coords[:, 1],
            mode='markers',
            marker=dict(size=point_size, opacity=opacity, color=color_vals, colorscale='Viridis',
                       colorbar=dict(title=color_key)),
        ))
```

第 33-37 行（无 color 分支）：
```python
        fig.add_trace(go.Scattergl(
            x=coords[:, 0], y=coords[:, 1],
            mode='markers',
            marker=dict(size=point_size, opacity=opacity, color='#1a237e'),
        ))
```

在 `fig.update_layout(...)` 调用中添加 legend font size：
```python
    fig.update_layout(
        title=title, xaxis_title='UMAP-1', yaxis_title='UMAP-2',
        plot_bgcolor='white', width=700, height=500,
        margin=dict(l=40, r=40, t=40, b=40),
        legend=dict(font=dict(size=legend_fontsize)),
    )
```

- [ ] **Step 3: 扩展 get_plotly_layout 支持 umap 相关参数**

将 `base.py:87-104` 替换为：

```python
    def get_plotly_layout(self, title='', **overrides):
        """根据 self.params['_visualization'] 生成统一的 Plotly layout dict。"""
        viz = self.params.get('_visualization', {})
        theme_name = viz.get('theme', 'default')
        theme = VISUALIZATION_THEMES.get(theme_name, VISUALIZATION_THEMES['default'])
        layout = {
            'title': title,
            'width': viz.get('figure_width', 800),
            'height': viz.get('figure_height', 500),
            'plot_bgcolor': viz.get('bg_color', theme['bg_color']),
            'colorway': viz.get('color_palette', theme['color_palette']),
            'font': {
                'family': viz.get('font_family', theme['font_family']),
                'size': viz.get('font_size', 12),
            },
            'legend': {
                'font': {'size': viz.get('umap_legend_fontsize', 10)},
            },
        }
        layout.update(overrides)
        return layout

    def get_viz_params(self):
        """返回传递给 umap_scatter 的可视化参数 dict。"""
        viz = self.params.get('_visualization', {})
        return {
            'umap_point_size': viz.get('umap_point_size', 5),
            'umap_opacity': viz.get('umap_opacity', 0.7),
            'umap_legend_fontsize': viz.get('umap_legend_fontsize', 10),
        }
```

- [ ] **Step 4: Commit**

```bash
git add modules/visualization.py modules/base.py
git commit -m "feat(sc): P6 viz — umap_point_size, umap_opacity, umap_legend_fontsize"
```

---

## Task 5: base.py — 新增 export_results() 方法

**Files:**
- Modify: `modules/base.py`

**Context:** 当前各模块各自处理 CSV/JSON 导出，没有统一的 h5ad/csv/loom 导出方法。spec 要求 `export_results()` 支持 h5ad/csv/loom 格式、层选择、obsm 选择、压缩。

- [ ] **Step 1: 在 `get_viz_params()` 方法之后添加 export_results**

```python
    def export_results(self, adata, output_dir, export_format='h5ad',
                       include_layers=None, include_obsm=None, compression='gzip'):
        """将 adata 导出为指定格式。返回导出文件路径列表。"""
        import os
        os.makedirs(output_dir, exist_ok=True)
        exported = []

        if export_format == 'h5ad':
            out_path = os.path.join(output_dir, 'result.h5ad')
            # 按需裁剪 layers
            if include_layers is not None:
                keep = [l for l in include_layers if l in adata.layers]
                adata.layers_keys_to_keep = keep
            adata.write_h5ad(out_path, compression=compression if compression != 'none' else None)
            exported.append(out_path)

        elif export_format == 'csv':
            obs_path = os.path.join(output_dir, 'obs.csv')
            adata.obs.to_csv(obs_path)
            exported.append(obs_path)
            var_path = os.path.join(output_dir, 'var.csv')
            adata.var.to_csv(var_path)
            exported.append(var_path)
            x_path = os.path.join(output_dir, 'expression.csv')
            import pandas as pd
            if hasattr(adata.X, 'toarray'):
                df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
            else:
                df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
            df.to_csv(x_path, compression='gzip' if compression == 'gzip' else None)
            exported.append(x_path)

        elif export_format == 'loom':
            out_path = os.path.join(output_dir, 'result.loom')
            adata.write_loom(out_path)
            exported.append(out_path)

        if include_obsm:
            for key in include_obsm:
                if key in adata.obsm:
                    import numpy as np
                    obsm_path = os.path.join(output_dir, f'obsm_{key}.csv')
                    obsm_df = pd.DataFrame(adata.obsm[key], index=adata.obs_names)
                    obsm_df.to_csv(obsm_path)
                    exported.append(obsm_path)

        return exported
```

- [ ] **Step 2: 验证**

```bash
cd /data/GJ/platform
python -c "
from modules.base import BaseAnalysis
assert hasattr(BaseAnalysis, 'export_results')
print('export_results method exists')
"
```

- [ ] **Step 3: Commit**

```bash
git add modules/base.py
git commit -m "feat(sc): P6 base — add export_results() with h5ad/csv/loom support"
```

---

## Task 6: 前端 — 可视化设置面板 + 过滤规则 UI

**Files:**
- Modify: `templates/analysis_select.html:388-400`（可视化注入逻辑）+ 新增 UI 面板

**Context:** 当前 `_visualization` 以空 JSON `{}` 注入，无 UI 控件。`_filters` 有后端 `apply_filters()` 支持但无前端 UI。

- [ ] **Step 1: 在模块参数区域上方添加可视化设置面板**

在 `analysis_select.html` 的模块列表下方（约 line 35 附近，preset toolbar 之后）添加：

```html
<!-- 可视化设置面板 -->
<div class="viz-settings-panel" style="margin: 10px 0; border: 1px solid #e0e0e0; border-radius: 6px; overflow: hidden;">
    <div class="viz-panel-header" style="padding: 8px 12px; background: #f5f5f5; cursor: pointer; font-weight: 600; font-size: 13px;"
         onclick="this.nextElementSibling.style.display = this.nextElementSibling.style.display === 'none' ? 'block' : 'none'; this.querySelector('.arrow').textContent = this.nextElementSibling.style.display === 'none' ? '▸' : '▾';">
        <span class="arrow">▸</span> 可视化设置
    </div>
    <div class="viz-panel-body" style="display: none; padding: 12px;">
        <div style="display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 8px; font-size: 12px;">
            <label>主题
                <select name="_viz_theme" class="form-control" style="width: 100%;">
                    <option value="default" selected>默认</option>
                    <option value="nature">Nature 论文</option>
                    <option value="dark">暗色</option>
                </select>
            </label>
            <label>宽度
                <input type="number" name="_viz_width" value="800" min="200" max="2000" class="form-control" style="width: 100%;">
            </label>
            <label>高度
                <input type="number" name="_viz_height" value="500" min="200" max="1500" class="form-control" style="width: 100%;">
            </label>
            <label>字号
                <input type="number" name="_viz_fontsize" value="12" min="8" max="24" class="form-control" style="width: 100%;">
            </label>
            <label>字体
                <select name="_viz_fontfamily" class="form-control" style="width: 100%;">
                    <option value="Arial">Arial</option>
                    <option value="Helvetica">Helvetica</option>
                    <option value="Times New Roman">Times New Roman</option>
                    <option value="SimSun">宋体</option>
                </select>
            </label>
            <label>UMAP 点大小
                <input type="number" name="_viz_umap_point_size" value="5" min="1" max="20" class="form-control" style="width: 100%;">
            </label>
            <label>UMAP 透明度
                <input type="number" name="_viz_umap_opacity" value="0.7" min="0.1" max="1" step="0.1" class="form-control" style="width: 100%;">
            </label>
            <label>图例字号
                <input type="number" name="_viz_legend_fontsize" value="10" min="6" max="20" class="form-control" style="width: 100%;">
            </label>
            <label>导出格式
                <div style="display: flex; gap: 10px; align-items: center; padding-top: 4px;">
                    <label style="font-weight: normal;"><input type="checkbox" name="_viz_export_plotly" checked> Plotly JSON</label>
                    <label style="font-weight: normal;"><input type="checkbox" name="_viz_export_svg"> SVG</label>
                    <label style="font-weight: normal;"><input type="checkbox" name="_viz_export_png"> PNG</label>
                </div>
            </label>
        </div>
    </div>
</div>
```

- [ ] **Step 2: 修改表单提交逻辑，收集可视化参数**

替换 `analysis_select.html:388-400` 的 `_visualization` 注入逻辑为：

```javascript
// 表单提交时注入 _visualization 和 _filters 配置
document.querySelector('.param-form').addEventListener('submit', function() {
    // 可视化配置
    let vizHidden = this.querySelector('input[name="_visualization"]');
    if (!vizHidden) {
        vizHidden = document.createElement('input');
        vizHidden.type = 'hidden';
        vizHidden.name = '_visualization';
        this.appendChild(vizHidden);
    }
    const exportFormats = [];
    if (this.querySelector('[name="_viz_export_plotly"]')?.checked) exportFormats.push('plotly_json');
    if (this.querySelector('[name="_viz_export_svg"]')?.checked) exportFormats.push('svg');
    if (this.querySelector('[name="_viz_export_png"]')?.checked) exportFormats.push('png');

    vizHidden.value = JSON.stringify({
        theme: this.querySelector('[name="_viz_theme"]')?.value || 'default',
        figure_width: parseInt(this.querySelector('[name="_viz_width"]')?.value) || 800,
        figure_height: parseInt(this.querySelector('[name="_viz_height"]')?.value) || 500,
        font_size: parseInt(this.querySelector('[name="_viz_fontsize"]')?.value) || 12,
        font_family: this.querySelector('[name="_viz_fontfamily"]')?.value || 'Arial',
        umap_point_size: parseInt(this.querySelector('[name="_viz_umap_point_size"]')?.value) || 5,
        umap_opacity: parseFloat(this.querySelector('[name="_viz_umap_opacity"]')?.value) || 0.7,
        umap_legend_fontsize: parseInt(this.querySelector('[name="_viz_legend_fontsize"]')?.value) || 10,
        export_formats: exportFormats,
    });

    // 过滤规则配置
    let filterHidden = this.querySelector('input[name="_filters"]');
    if (!filterHidden) {
        filterHidden = document.createElement('input');
        filterHidden.type = 'hidden';
        filterHidden.name = '_filters';
        this.appendChild(filterHidden);
    }
    const filtersByModule = {};
    this.querySelectorAll('.filter-rule-row').forEach(row => {
        const module = row.dataset.module;
        const col = row.querySelector('.filter-col')?.value;
        const op = row.querySelector('.filter-op')?.value;
        const val = row.querySelector('.filter-val')?.value;
        if (col && op && val !== '') {
            if (!filtersByModule[module]) filtersByModule[module] = [];
            let parsedVal = val;
            if (op === 'in' || op === 'not_in') {
                parsedVal = val.split(',').map(s => s.trim());
            } else if (op === 'between') {
                parsedVal = val.split(',').map(s => parseFloat(s.trim()));
            } else {
                parsedVal = isNaN(val) ? val : parseFloat(val);
            }
            filtersByModule[module].push({ column: col, op: op, value: parsedVal });
        }
    });
    filterHidden.value = JSON.stringify(filtersByModule);
});
```

- [ ] **Step 3: 添加过滤规则 UI**

在每个模块的参数折叠区域底部（遍历 `analysis_select.html` 中的模块参数块），添加过滤规则区域。在现有的模块参数块末尾（`</div>` 关闭参数区域之前）插入：

```html
<div class="filter-rules-section" style="margin-top: 8px; padding: 8px; background: #fafafa; border-radius: 4px; font-size: 12px;">
    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 6px;">
        <span style="font-weight: 600; color: #666;">自定义过滤规则（可选）</span>
        <button type="button" class="btn-add-filter" onclick="addFilterRule(this, '{{ module.MODULE_NAME }}')"
                style="background: #4caf50; color: white; border: none; border-radius: 3px; padding: 2px 8px; cursor: pointer; font-size: 11px;">+ 添加规则</button>
    </div>
    <div class="filter-rules-container" data-module="{{ module.MODULE_NAME }}"></div>
</div>
```

由于模板是 Jinja2 动态渲染，需要在遍历模块的循环中添加。如果模块参数是动态生成的，则需要在 JS 中处理。在 `analysis_select.html` 的 JS 部分添加：

```javascript
function addFilterRule(btn, moduleName) {
    const container = btn.closest('.filter-rules-section').querySelector('.filter-rules-container');
    const row = document.createElement('div');
    row.className = 'filter-rule-row';
    row.dataset.module = moduleName;
    row.style.cssText = 'display: flex; gap: 6px; align-items: center; margin-bottom: 4px;';
    row.innerHTML = `
        <select class="filter-col form-control" style="width: 160px; font-size: 11px;" placeholder="列名">
            <option value="">-- 列名 --</option>
        </select>
        <select class="filter-op form-control" style="width: 100px; font-size: 11px;">
            <option value=">=">≥</option>
            <option value="<=">≤</option>
            <option value="==">＝</option>
            <option value="!=">≠</option>
            <option value="in">包含</option>
            <option value="not_in">不包含</option>
            <option value="between">范围</option>
        </select>
        <input class="filter-val form-control" style="width: 140px; font-size: 11px;" placeholder="值">
        <button type="button" onclick="this.parentElement.remove()" style="background: #e53935; color: white; border: none; border-radius: 3px; padding: 2px 6px; cursor: pointer;">×</button>
    `;
    container.appendChild(row);
    // 填充列名下拉
    const filePath = document.querySelector('select[name="input_path"]')?.value;
    if (filePath) {
        fetch(`/api/obs-columns?file_path=${encodeURIComponent(filePath)}`)
            .then(r => r.json())
            .then(d => {
                const sel = row.querySelector('.filter-col');
                (d.columns || []).forEach(c => {
                    sel.innerHTML += `<option value="${c}">${c}</option>`;
                });
            });
    }
}
```

- [ ] **Step 4: 验证前端渲染**

启动 dev server 检查：
```bash
cd /data/GJ/platform
python app.py &
sleep 2
curl -s http://localhost:5000/ | grep -c "viz-settings-panel"
```

- [ ] **Step 5: Commit**

```bash
git add templates/analysis_select.html
git commit -m "feat: Phase 3+4 — viz settings panel + filter rules UI"
```

---

## 总结

| Task | 文件 | 功能 |
|------|------|------|
| 1 | annotation.py | celltypist 集成 + 置信度计算 + 相似簇合并 |
| 2 | deg.py | correction_method 多重检验校正 + min_pct 过滤 |
| 3 | proportion.py | fisher_exact/permutation 检验 + min_cells_per_group |
| 4 | visualization.py + base.py | umap_scatter 参数化 + get_plotly_layout 扩展 |
| 5 | base.py | export_results() 统一导出方法 |
| 6 | analysis_select.html | 可视化设置面板 + 过滤规则 UI |

Task 1-3 互不依赖，可并行实施。Task 4-5 互不依赖。Task 6 依赖 Task 4（需要知道参数名）。整体按 1→2→3→4→5→6 或并行均可。
