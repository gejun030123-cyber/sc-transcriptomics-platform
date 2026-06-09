# Bulk DEG 差异分析增强 + 多组整合分析 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 增强 bulk_deg.py（新增 edgeR/limma、LRT、自动比较、高级参数）+ 新建 bulk_deg_integration.py（Upset 图、一致性评分、logFC 矩阵分析）。

**Architecture:** bulk_deg.py 的 `_run_single_comparison()` 函数扩展 method_map 和参数传递；`run()` 方法增加 auto_comparisons 自动比较逻辑。新建 bulk_deg_integration.py 作为独立模块读取 bulk_deg 输出。

**Tech Stack:** omicverse, inmoose, patsy, plotly, numpy, pandas, scipy

---

## 文件结构

| 文件 | 操作 | 职责 |
|------|------|------|
| `modules/bulk_deg.py` | 大幅修改 | 新方法 + LRT + 自动比较 + 高级参数 + 结果整合 |
| `modules/bulk_deg_integration.py` | **新建** | 多组差异整合模块 |
| `routes/analysis.py` | 修改 | 两个模块的 PARAM_SCHEMAS |
| `modules/__init__.py` | 修改 | 注册新模块 |

---

## Task 1: 扩展统计方法 + pyDEG 高级参数

**Files:**
- Modify: `modules/bulk_deg.py:39-66`（`_run_single_comparison` 函数签名和核心调用）

- [ ] **Step 1: 更新 `_run_single_comparison` 函数签名**

将函数签名从：
```python
def _run_single_comparison(adata, counts, group1_samples, group2_samples, group1, group2,
                           method, fc_threshold, pval_threshold, top_n, gene_id_to_name,
                           plots_dir, results_dir, suffix='', viz_params=None):
```

改为：
```python
def _run_single_comparison(adata, counts, group1_samples, group2_samples, group1, group2,
                           method, fc_threshold, pval_threshold, top_n, gene_id_to_name,
                           plots_dir, results_dir, suffix='', viz_params=None,
                           cooks_filter=True, independent_filter=True, padj_method='fdr_bh',
                           base_mean_filter=0, regulation_filter='both'):
```

- [ ] **Step 2: 更新 method_map**

将：
```python
    method_map = {'t-test': 'ttest', 'mann-whitney': 'wilcox', 'deseq2': 'DEseq2'}
```

改为：
```python
    method_map = {
        't-test': 'ttest', 'mann-whitney': 'wilcox',
        'deseq2': 'DEseq2', 'edger': 'edgepy', 'limma': 'limma'
    }
    ov_method = method_map.get(method, 'ttest')
    if ov_method in ('edgepy', 'limma'):
        try:
            import inmoose
        except ImportError:
            raise ImportError(f"方法 {method} 需要安装 inmoose: pip install inmoose patsy")
```

- [ ] **Step 3: 更新 deg_analysis 调用传入高级参数**

将：
```python
    result = dds.deg_analysis(group1_samples, group2_samples, method=ov_method)
```

改为：
```python
    result = dds.deg_analysis(
        group1_samples, group2_samples, method=ov_method,
        cooks_filter=cooks_filter, independent_filter=independent_filter,
        multipletests_method=padj_method
    )
```

- [ ] **Step 4: 接入 base_mean_filter**

在 `deg_df = deg_df.sort_values('padj')` 之后插入：
```python
    # 基础表达量过滤
    if base_mean_filter > 0 and 'mean_group1' in deg_df.columns:
        base_mean = (deg_df['mean_group1'] + deg_df['mean_group2']) / 2
        deg_df = deg_df[base_mean >= base_mean_filter].copy()
```

- [ ] **Step 5: 接入 regulation_filter**

在 `deg_df` 最终输出之前，按方向过滤：
```python
    if regulation_filter == 'up':
        deg_df = deg_df[deg_df['regulation'] == 'Up'].copy()
    elif regulation_filter == 'down':
        deg_df = deg_df[deg_df['regulation'] == 'Down'].copy()
```

在返回值之前插入（在 `n_up = ...` 之后）。

- [ ] **Step 6: 语法验证**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_deg import BulkDEGAnalysis; print('OK')"`

- [ ] **Step 7: 提交**

```bash
git add modules/bulk_deg.py
git commit -m "feat(bulk_deg): add edgeR/limma methods, expose pyDEG advanced params, wire base_mean_filter and regulation_filter"
```

---

## Task 2: LRT 多组检验

**Files:**
- Modify: `modules/bulk_deg.py`（`_run_single_comparison` 或新函数）

- [ ] **Step 1: 添加 LRT 支持**

在 `_run_single_comparison` 函数中，当 `method='edger'` 且需要 LRT 时，使用 `glmLRT`。

在 `dds.deg_analysis()` 调用之后，检测是否需要 LRT。实际实现：在 `run()` 方法中添加 LRT 分支，不通过 `_run_single_comparison`。

在 `run()` 方法的多比较循环之前，插入 LRT 分支：

```python
        # LRT 多组检验
        if test_type == 'lrt' and method == 'edger':
            self.progress(10, "运行 LRT 多组检验...")
            try:
                import inmoose
            except ImportError:
                raise ImportError("LRT 检验需要安装 inmoose: pip install inmoose patsy")
            import omicverse as ov
            count_df = pd.DataFrame(counts.T, index=adata.var_names.tolist(), columns=adata.obs.index.tolist())
            dds_lrt = ov.bulk.pyDEG(count_df)
            dds_lrt.drop_duplicates_index()
            dds_lrt.normalize()
            # 使用 edgeR 的 glmLRT
            from inmoose.edgepy import DGEList, glmLRT, estimateGLMCommonDisp, estimateGLMTagwiseDisp
            from patsy import dmatrix
            # 构建设计矩阵
            all_groups = adata.obs[groupby].astype(str).tolist()
            design_matrix = dmatrix(f'C({groupby}, Treatment(reference="{reference_group or unique_groups[0]}"))',
                                     data=adata.obs, return_type='dataframe')
            # DGEList + 估计离散度
            y = DGEList(counts=count_df.values.T, genes=count_df.index.tolist())
            y = estimateGLMCommonDisp(y, design_matrix)
            y = estimateGLMTagwiseDisp(y, design_matrix)
            # 拟合模型 + LRT
            from inmoose.edgepy import glmFit
            fit = glmFit(y, design_matrix)
            lrt = glmLRT(fit)
            lrt_results = lrt.table
            # 保存 LRT 结果
            lrt_csv = os.path.join(results_dir, 'bulk_deg_lrt_results.csv')
            lrt_results.to_csv(lrt_csv)
            result_files.append({'file_path': lrt_csv, 'file_type': 'csv', 'category': 'table', 'label': 'LRT 多组检验结果'})
```

**注意**：LRT 实现依赖 inmoose 的底层 API。如果 inmoose 版本不支持直接的 `glmFit`/`glmLRT`，回退到 pyDEG 的 edgepy 方法（它内部已使用 LRT）。

- [ ] **Step 2: 语法验证 + 提交**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_deg import BulkDEGAnalysis; print('OK')"`
Commit: `git add modules/bulk_deg.py && git commit -m "feat(bulk_deg): add LRT multi-group test via edgeR glmLRT"`

---

## Task 3: 自动生成配对比较 + 参考组

**Files:**
- Modify: `modules/bulk_deg.py:230-330`（run() 方法中的比较循环）

- [ ] **Step 1: 添加 auto_comparisons 和 reference_group 参数解析**

在 `run()` 方法的参数解析区域追加：

```python
        auto_comparisons = self.params.get('auto_comparisons', 'manual')
        reference_group = self.params.get('reference_group', '').strip()
        test_type = self.params.get('test_type', 'pairwise')
```

- [ ] **Step 2: 自动生成比较列表**

在现有 `comparison_pairs = _parse_comparisons(...)` 之后插入：

```python
        # 自动比较模式
        if auto_comparisons == 'all_pairwise':
            unique_groups = sorted(adata.obs[groupby].astype(str).unique().tolist())
            comparison_pairs = []
            for i in range(len(unique_groups)):
                for j in range(i + 1, len(unique_groups)):
                    if reference_group and unique_groups[j] == reference_group:
                        comparison_pairs.append((unique_groups[i], unique_groups[j]))
                    elif reference_group and unique_groups[i] == reference_group:
                        comparison_pairs.append((unique_groups[j], unique_groups[i]))
                    else:
                        comparison_pairs.append((unique_groups[i], unique_groups[j]))
        elif auto_comparisons == 'vs_reference':
            if not reference_group:
                reference_group = sorted(adata.obs[groupby].astype(str).unique().tolist())[0]
            unique_groups = sorted(adata.obs[groupby].astype(str).unique().tolist())
            comparison_pairs = [(g, reference_group) for g in unique_groups if g != reference_group]
```

- [ ] **Step 3: 将新参数传递给 _run_single_comparison**

在所有 `_run_single_comparison(...)` 调用中追加新参数：

```python
            deg_df, comp_files, n_up, n_down = _run_single_comparison(
                adata, counts, g1_samples, g2_samples, g1, g2,
                method, fc_threshold, pval_threshold, top_n, gene_id_to_name,
                plots_dir, results_dir, suffix=str(idx), viz_params=viz_params,
                cooks_filter=cooks_filter, independent_filter=independent_filter,
                padj_method=padj_method, base_mean_filter=base_mean_filter,
                regulation_filter=regulation_filter
            )
```

对所有调用点（多比较模式和单比较模式）做同样修改。

- [ ] **Step 4: 语法验证 + 提交**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_deg import BulkDEGAnalysis; print('OK')"`
Commit: `git add modules/bulk_deg.py && git commit -m "feat(bulk_deg): add auto pairwise comparisons and reference group"`

---

## Task 4: 多比较结果整合（合并表 + 摘要 + Volcano 并排）

**Files:**
- Modify: `modules/bulk_deg.py:307-330`（多比较结果整合区域）

- [ ] **Step 1: 生成合并结果 CSV**

在现有 `all_deg_dfs` 合并逻辑之后，生成 logFC/padj 矩阵格式的合并表：

```python
        # 多比较合并表（logFC + padj 矩阵格式）
        if len(all_deg_dfs) > 1:
            merged_parts = []
            for idx, (comp_name, deg_df) in enumerate(zip(comparison_names, all_deg_dfs)):
                part = deg_df[['gene', 'log2FC', 'padj', 'regulation']].copy()
                part.columns = ['gene', f'{comp_name}_log2FC', f'{comp_name}_padj', f'{comp_name}_regulation']
                if merged_parts:
                    merged_parts.append(part.drop(columns=['gene']))
                else:
                    merged_parts.append(part)
            merged_df = pd.concat(merged_parts, axis=1)
            merged_csv = os.path.join(results_dir, 'bulk_deg_merged_comparisons.csv')
            merged_df.to_csv(merged_csv, index=False)
            result_files.append({'file_path': merged_csv, 'file_type': 'csv', 'category': 'table', 'label': '多比较合并结果'})
```

- [ ] **Step 2: 计算共享差异基因**

```python
            # 共享差异基因统计
            up_sets = [set(deg_df[deg_df['regulation']=='Up']['gene']) for deg_df in all_deg_dfs]
            down_sets = [set(deg_df[deg_df['regulation']=='Down']['gene']) for deg_df in all_deg_dfs]
            shared_up = len(set.intersection(*up_sets)) if len(up_sets) >= 2 else 0
            shared_down = len(set.intersection(*down_sets)) if len(down_sets) >= 2 else 0
```

- [ ] **Step 3: 生成 Volcano 并排展示**

```python
            # Volcano 并排展示
            from plotly.subplots import make_subplots
            n_comp = len(comparison_names)
            fig_multi = make_subplots(rows=1, cols=n_comp, subplot_titles=comparison_names,
                                       horizontal_spacing=0.05)
            for idx, (comp_name, deg_df) in enumerate(zip(comparison_names, all_deg_dfs)):
                colors = ['#e53935' if r == 'Up' else '#1565c0' if r == 'Down' else '#9e9e9e'
                          for r in deg_df['regulation']]
                neg_log_p = -np.log10(deg_df['padj'].values + 1e-300)
                fig_multi.add_trace(go.Scattergl(
                    x=deg_df['log2FC'].tolist(), y=neg_log_p.tolist(),
                    mode='markers', marker=dict(color=colors, size=4),
                    text=deg_df['gene'].tolist(), showlegend=False),
                    row=1, col=idx+1)
            fig_multi.update_layout(height=400, width=300*n_comp, title='多组比较 Volcano 图')
            fpath = os.path.join(plots_dir, 'bulk_deg_volcano_multi.json')
            with open(fpath, 'w') as f: f.write(fig_multi.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'volcano', 'label': '多组比较 Volcano'})
```

- [ ] **Step 4: 更新 summary**

在多比较 summary 中追加 `per_comparison`、`shared_up`、`shared_down`。

- [ ] **Step 5: 语法验证 + 提交**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_deg import BulkDEGAnalysis; print('OK')"`
Commit: `git add modules/bulk_deg.py && git commit -m "feat(bulk_deg): add merged comparison table, shared DEG stats, multi-volcano plot"`

---

## Task 5: PARAM_SCHEMAS 更新

**Files:**
- Modify: `routes/analysis.py:137-151`

- [ ] **Step 1: 更新 bulk_deg PARAM_SCHEMAS**

在 `PARAM_SCHEMAS['bulk_deg']` 末尾追加新参数：

```python
        {'key': 'test_type', 'label': '检验类型', 'type': 'select', 'options': ['pairwise', 'lrt'], 'default': 'pairwise',
         'help': 'pairwise：两两比较。lrt：似然比检验（仅 edger），一次性检验所有组间差异。'},
        {'key': 'auto_comparisons', 'label': '自动生成比较', 'type': 'select', 'options': ['manual', 'all_pairwise', 'vs_reference'], 'default': 'manual',
         'help': 'manual：手动输入比较。all_pairwise：自动生成所有两两配对。vs_reference：所有组 vs 参考组。'},
        {'key': 'reference_group', 'label': '参考组（可选）', 'type': 'dynamic_select', 'depends_on': 'groupby', 'default': '',
         'help': '指定参考组。确保 logFC 方向一致（正值=该组>参考组）。'},
        {'key': 'cooks_filter', 'label': "Cook's 距离过滤", 'type': 'checkbox', 'default': True,
         'help': "剔除 Cook's 距离过大的异常高表达基因。"},
        {'key': 'independent_filter', 'label': '独立过滤', 'type': 'checkbox', 'default': True,
         'help': '自动去除低表达基因，提升检测效力。'},
        {'key': 'padj_method', 'label': 'p 值校正方法', 'type': 'select', 'options': ['fdr_bh', 'bonferroni', 'holm', 'fdr_by'], 'default': 'fdr_bh',
         'help': '多重检验校正方法。BH：最常用。Bonferroni：最严格。Holm：逐步校正。BY：依赖性校正。'},
        {'key': 'regulation_filter', 'label': '差异方向', 'type': 'select', 'options': ['both', 'up', 'down'], 'default': 'both',
         'help': 'both：双向。up：仅上调。down：仅下调。'},
```

同时将现有 `method` 的 options 从 `['t-test', 'mann-whitney', 'deseq2']` 扩展为 `['t-test', 'mann-whitney', 'deseq2', 'edger', 'limma']`。

- [ ] **Step 2: 验证 + 提交**

Run: `cd /data/GJ/platform && python -c "from routes.analysis import PARAM_SCHEMAS; print('OK')"`
Commit: `git add routes/analysis.py && git commit -m "feat(bulk_deg): add PARAM_SCHEMAS for edgeR/limma, LRT, auto comparisons, advanced params"`

---

## Task 6: 新建 bulk_deg_integration 模块

**Files:**
- Create: `modules/bulk_deg_integration.py`
- Modify: `modules/__init__.py`（注册新模块）
- Modify: `routes/analysis.py`（添加 PARAM_SCHEMAS）

- [ ] **Step 1: 创建模块骨架**

```python
# modules/bulk_deg_integration.py
import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class BulkDEGIntegrationAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_deg_integration"
    DISPLAY_NAME = "多组差异整合分析"
    DESCRIPTION = "多组比较结果整合：Upset 图、一致性评分、logFC 矩阵分析"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        import scanpy as sc

        self.progress(5, "加载差异分析结果...")
        # 从项目 results 目录扫描 bulk_deg 输出
        results_dir = os.path.join(self.project_dir, 'results')
        deg_files = [f for f in os.listdir(results_dir)
                     if f.startswith('bulk_deg_') and f.endswith('.csv') and 'merged' not in f and 'lrt' not in f
                     and 'skipped' not in f and 'gene_filter' not in f]

        min_comparisons = int(self.params.get('min_comparisons', 2))
        consistency_n = int(self.params.get('consistency_n', 50))
        fc_threshold = float(self.params.get('fc_threshold', 2.0))
        pval_threshold = float(self.params.get('pval_threshold', 0.05))

        # 解析每个比较结果
        self.progress(20, "解析比较结果...")
        comparisons = {}
        for f in deg_files:
            df = pd.read_csv(os.path.join(results_dir, f))
            # 提取比较名（从文件名）
            name = f.replace('bulk_deg_', '').replace('.csv', '')
            comparisons[name] = df

        if len(comparisons) < 2:
            return {'output_adata': input_path, 'result_files': [], 'summary': {'error': '需要至少 2 个比较结果'}}

        # 构建 logFC 矩阵和 padj 矩阵
        self.progress(30, "构建 logFC/padj 矩阵...")
        all_genes = set()
        for df in comparisons.values():
            all_genes.update(df['gene'].tolist())
        all_genes = sorted(all_genes)
        comp_names = sorted(comparisons.keys())

        logfc_matrix = pd.DataFrame(index=all_genes, columns=comp_names, dtype=float)
        padj_matrix = pd.DataFrame(index=all_genes, columns=comp_names, dtype=float)
        regulation_matrix = pd.DataFrame(index=all_genes, columns=comp_names, dtype=float)

        for name, df in comparisons.items():
            df_indexed = df.set_index('gene')
            for g in all_genes:
                if g in df_indexed.index:
                    logfc_matrix.loc[g, name] = df_indexed.loc[g, 'log2FC']
                    padj_matrix.loc[g, name] = df_indexed.loc[g, 'padj']
                    reg = df_indexed.loc[g, 'regulation']
                    regulation_matrix.loc[g, name] = 1 if reg == 'Up' else -1 if reg == 'Down' else 0

        logfc_matrix = logfc_matrix.fillna(0)
        padj_matrix = padj_matrix.fillna(1)
        regulation_matrix = regulation_matrix.fillna(0)

        # ... 继续实现各个功能 ...
```

- [ ] **Step 2: 实现一致性评分**

```python
        # 一致性评分
        self.progress(50, "计算一致性评分...")
        consistency_scores = []
        for g in all_genes:
            sig_count = sum(1 for c in comp_names if padj_matrix.loc[g, c] < pval_threshold
                           and abs(logfc_matrix.loc[g, c]) >= np.log2(fc_threshold))
            if sig_count >= min_comparisons:
                signs = regulation_matrix.loc[g, comp_names].values
                neg_log_p = -np.log10(padj_matrix.loc[g, comp_names].values + 1e-300)
                score = float(np.mean(signs * neg_log_p))
                consistency_scores.append({
                    'gene': g,
                    'consistency_score': round(score, 4),
                    'n_significant': sig_count,
                    'direction': 'Up' if score > 0 else 'Down' if score < 0 else 'Mixed',
                })
        consistency_df = pd.DataFrame(consistency_scores).sort_values('consistency_score', ascending=False)
```

- [ ] **Step 3: 实现 Upset 图**

```python
        # Upset 图（用 plotly bar chart + matrix）
        self.progress(60, "生成 Upset 图...")
        sig_sets = {}
        for c in comp_names:
            sig_genes = set(logfc_matrix.index[(padj_matrix[c] < pval_threshold) &
                                                (abs(logfc_matrix[c]) >= np.log2(fc_threshold))])
            sig_sets[c] = sig_genes

        # 计算所有交集组合
        from itertools import combinations
        intersection_data = []
        for r in range(1, len(comp_names) + 1):
            for combo in combinations(comp_names, r):
                intersection = sig_sets[combo[0]]
                for c in combo[1:]:
                    intersection = intersection & sig_sets[c]
                if intersection:
                    intersection_data.append({
                        'sets': combo,
                        'count': len(intersection),
                        'label': ' ∩ '.join(combo),
                    })
        intersection_data.sort(key=lambda x: x['count'], reverse=True)
```

- [ ] **Step 4: 实现 logFC 矩阵热图**

```python
        # logFC 矩阵热图
        self.progress(75, "生成 logFC 矩阵热图...")
        top_genes = consistency_df.head(consistency_n)['gene'].tolist() if len(consistency_df) > 0 else all_genes[:consistency_n]
        logfc_subset = logfc_matrix.loc[logfc_matrix.index.isin(top_genes)]

        if len(logfc_subset) > 0:
            # Clip to [-5, 5]
            logfc_clipped = logfc_subset.clip(-5, 5)
            fig_heat = go.Figure(data=go.Heatmap(
                z=logfc_clipped.values.tolist(),
                x=comp_names,
                y=logfc_clipped.index.tolist(),
                colorscale='RdBu_r', zmid=0,
                colorbar=dict(title='log2FC')
            ))
            fig_heat.update_layout(
                title='Top 差异基因 logFC 矩阵',
                height=max(400, len(logfc_subset) * 12 + 100),
                width=max(500, len(comp_names) * 80 + 200),
            )
```

- [ ] **Step 5: 实现方向一致性热图 + 比较相关性**

```python
        # 方向一致性热图
        reg_subset = regulation_matrix.loc[regulation_matrix.index.isin(top_genes)]
        fig_dir = go.Figure(data=go.Heatmap(
            z=reg_subset.values.tolist(), x=comp_names, y=reg_subset.index.tolist(),
            colorscale=[[0, '#1565c0'], [0.5, '#ffffff'], [1, '#e53935']],
            zmid=0, showscale=False))

        # 比较间 logFC 相关性
        corr_matrix = logfc_matrix.corr(method='pearson')
        fig_corr = go.Figure(data=go.Heatmap(
            z=corr_matrix.values.tolist(), x=comp_names, y=comp_names,
            colorscale='RdBu_r', zmid=0, text=np.round(corr_matrix.values, 2).tolist(),
            texttemplate='%{text}', colorbar=dict(title='Pearson r')))
```

- [ ] **Step 6: 保存所有输出 + 注册模块**

```python
        # 保存输出
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        # 保存各图表和 CSV...
        consistency_csv = os.path.join(results_dir, 'deg_integration_consistency.csv')
        consistency_df.to_csv(consistency_csv, index=False)
        result_files.append(...)

        # 注册模块
        # 在 modules/__init__.py 的 MODULE_REGISTRY 中追加：
        # 'bulk_deg_integration': BulkDEGIntegrationAnalysis
```

- [ ] **Step 7: 语法验证 + 提交**

Run: `cd /data/GJ/platform && python -c "from modules.bulk_deg_integration import BulkDEGIntegrationAnalysis; print('OK')"`
Commit: `git add modules/bulk_deg_integration.py modules/__init__.py routes/analysis.py && git commit -m "feat: add bulk_deg_integration module (upset, consistency score, logFC matrix)"`

---

## 实施顺序与依赖

```
Task 1 (edgeR/limma + 高级参数)    ← 无依赖
    ↓
Task 2 (LRT 多组检验)              ← 依赖 Task 1（需要 edger 方法）
    ↓
Task 3 (自动比较 + 参考组)          ← 依赖 Task 1
    ↓
Task 4 (多比较结果整合)              ← 依赖 Task 3（需要多个比较结果）
    ↓
Task 5 (PARAM_SCHEMAS)             ← 依赖 Task 1-4
    ↓
Task 6 (bulk_deg_integration)      ← 依赖 Task 4（读取 bulk_deg 输出）
```
