# Bulk RNA-seq 自定义比较功能增强实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 为 Bulk RNA-seq 管线的 4 个模块增加自定义比较能力：多组两两比较、自定义合并组、上调/下调分开富集、时序分组配对比较、自定义基因列表。

**Architecture:** 在现有模块上增量扩展参数，不引入新模块。每个模块的新增功能通过新参数控制，留空时行为与当前完全一致（向后兼容）。前端模板新增 textarea 类型支持。

**Tech Stack:** Python, Flask, OmicVerse, scanpy, Plotly, scipy, Jinja2, Bootstrap 5

**设计文档:** `docs/superpowers/specs/2026-06-06-bulk-custom-comparison-design.md`

---

## 文件结构

| 文件 | 操作 | 职责 |
|------|------|------|
| `modules/bulk_heatmap.py` | 修改 | 新增 `custom_genes` 参数支持 |
| `modules/bulk_enrichment.py` | 修改 | 新增 `split_direction` + `custom_genes` 参数支持 |
| `modules/bulk_deg.py` | 修改 | 新增 `comparisons` + `custom_groups` 参数支持 |
| `modules/bulk_timecourse.py` | 修改 | 新增 `pairwise_groups` 参数支持 |
| `routes/analysis.py` | 修改 | 4 个模块的 PARAM_SCHEMAS 新增参数 |
| `templates/analysis_select.html` | 修改 | textarea 参数类型渲染 |

---

## Task 1: bulk_heatmap — 自定义基因列表

**Files:**
- Modify: `modules/bulk_heatmap.py`
- Modify: `routes/analysis.py`

- [ ] **Step 1: 在 `modules/bulk_heatmap.py` 的 `run()` 方法中，在 `self.progress(40, "选择基因...")` 之后添加 custom_genes 分支**

在第 46 行 `self.progress(40, "选择基因...")` 之后、第 48 行 `if hm_type == 'top_var'` 之前插入：

```python
        custom_genes_str = self.params.get('custom_genes', '').strip()
        if custom_genes_str:
            # 自定义基因列表模式
            gene_list = [g.strip() for g in custom_genes_str.replace('\n', ',').split(',') if g.strip()]
            var_names_list = list(adata.var_names)
            top_idx = [var_names_list.index(g) for g in gene_list if g in var_names_list]
            not_found = [g for g in gene_list if g not in var_names_list]
            if not top_idx:
                raise ValueError(f"自定义基因列表中没有找到任何匹配基因。请检查基因名是否正确。")
            title = f'自定义基因热图 ({len(top_idx)} genes)'
            if not_found:
                title += f'，{len(not_found)} 个未找到'
        elif hm_type == 'top_var' or hm_type not in ('deg',):
```

同时将原来的 `if hm_type == 'top_var'` 改为 `elif hm_type == 'top_var'`（因为前面加了 `if custom_genes_str` 分支）。

- [ ] **Step 2: 在 `routes/analysis.py` 的 `PARAM_SCHEMAS['bulk_heatmap']` 末尾添加 custom_genes 参数**

在第 137 行 `{'key': 'groupby', ...}` 之后添加：

```python
        {'key': 'custom_genes', 'label': '自定义基因列表（可选）', 'type': 'textarea', 'default': '', 'help': '手动输入基因名，逗号或换行分隔。填写后忽略热图类型和基因数参数，直接用此列表绘制热图。'},
```

- [ ] **Step 3: 验证模块可导入**

```bash
cd /data/GJ/platform && python -c "from modules.bulk_heatmap import BulkHeatmapAnalysis; print('OK')"
```

Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add modules/bulk_heatmap.py routes/analysis.py
git commit -m "feat(bulk_heatmap): add custom_genes parameter for user-defined gene list"
```

---

## Task 2: bulk_enrichment — 上调/下调分开富集 + 自定义基因列表

**Files:**
- Modify: `modules/bulk_enrichment.py`
- Modify: `routes/analysis.py`

- [ ] **Step 1: 在 `modules/bulk_enrichment.py` 的 `run()` 方法中，在读取参数后添加 custom_genes 和 split_direction 逻辑**

在第 26 行 `input_source = self.params.get('input_source', '')` 之后添加：

```python
        split_direction = self.params.get('split_direction', False) in (True, 'true', 'on', '1')
        custom_genes_str = self.params.get('custom_genes', '').strip()
```

在第 44 行 `if not deg_genes and method == 'ORA':` 之前，替换现有的基因列表加载逻辑（第 34-43 行）为：

```python
        # 自定义基因列表优先
        if custom_genes_str:
            deg_genes = [g.strip() for g in custom_genes_str.replace('\n', ',').split(',') if g.strip()]
            gene_rnk = None
        elif input_source and os.path.exists(input_source):
```

保留原有的 `elif input_source` 块不变。

- [ ] **Step 2: 在 ORA 分析完成后添加 split_direction 分支**

在第 99 行 `enr = ov.bulk.geneset_enrichment(...)` 调用之后、保存 CSV 之前（第 101 行），添加 split_direction 逻辑。具体位置：在 `if method == 'ORA':` 块内，替换整个 ORA 块为：

```python
        if method == 'ORA':
            if split_direction and not custom_genes_str and 'regulation' in pd.read_csv(input_source).columns if input_source else False:
                # 拆分 Up/Down 分别做 ORA
                deg_full = pd.read_csv(input_source)
                directions = {'Up': deg_full[deg_full['regulation'] == 'Up']['gene'].tolist(),
                              'Down': deg_full[deg_full['regulation'] == 'Down']['gene'].tolist()}
                all_enr = []
                for direction, dir_genes in directions.items():
                    if not dir_genes:
                        continue
                    enr_dir = ov.bulk.geneset_enrichment(
                        gene_list=dir_genes, pathways_dict=pathways_dict,
                        pvalue_type='adjust', pvalue_threshold=pvalue_cutoff,
                        organism=organism_lower,
                        outdir=os.path.join(self.project_dir, f'enrichr_{direction.lower()}_tmp')
                    )
                    enr_dir['direction'] = direction
                    all_enr.append(enr_dir)

                    top_enr = enr_dir.head(top_n)
                    if len(top_enr) > 0:
                        x_col = 'Fractions' if 'Fractions' in top_enr.columns else None
                        term_col = 'Term' if 'Term' in top_enr.columns else None
                        y_vals = top_enr[term_col].tolist() if term_col else top_enr.index.tolist()
                        colors_dir = '#e53935' if direction == 'Up' else '#1a237e'

                        fig_bubble = go.Figure()
                        fig_bubble.add_trace(go.Scatter(
                            x=top_enr[x_col].values if x_col else list(range(len(top_enr))),
                            y=y_vals, mode='markers',
                            marker=dict(
                                size=top_enr['Overlap'].apply(lambda x: int(str(x).split('/')[0]) * 3 + 5).values if 'Overlap' in top_enr.columns else [10] * len(top_enr),
                                color=-np.log10(top_enr['P-value'].clip(lower=1e-300).values),
                                colorscale='YlOrRd', showscale=True, colorbar=dict(title='-log10(p)')
                            ), hovertemplate='%{y}<br>-log10(p): %{marker.color:.1f}<extra></extra>'
                        ))
                        fig_bubble.update_layout(
                            title=f'{database} ORA ({direction} genes, {organism})',
                            xaxis_title='Gene Fraction' if x_col else 'Index',
                            yaxis=dict(autorange='reversed'),
                            plot_bgcolor='white', width=800, height=max(400, top_n * 25 + 100)
                        )
                        fpath = os.path.join(plots_dir, f'enrichment_ora_{direction.lower()}_bubble.json')
                        with open(fpath, 'w') as f: f.write(fig_bubble.to_json(engine="json"))
                        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'enrichment', 'label': f'ORA {direction} 气泡图'})

                if all_enr:
                    combined = pd.concat(all_enr, ignore_index=True)
                    csv_path = os.path.join(results_dir, 'enrichment_split_results.csv')
                    combined.to_csv(csv_path, index=False)
                    result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': '上调/下调富集结果'})
                    n_sig = len(combined[combined['P-value'] < pvalue_cutoff]) if 'P-value' in combined.columns else 0
            else:
                # 标准 ORA（原有逻辑，保持不变）
                enr = ov.bulk.geneset_enrichment(
                    gene_list=deg_genes, pathways_dict=pathways_dict,
                    pvalue_type='adjust', pvalue_threshold=pvalue_cutoff,
                    organism=organism_lower,
                    outdir=os.path.join(self.project_dir, 'enrichr_tmp')
                )
                # ... 原有的保存 CSV、气泡图、条形图逻辑 ...
```

在 `else:` 分支内保留原有的全部 ORA 代码（从保存 CSV 到 `n_sig = len(enr[...])`）。

- [ ] **Step 3: 更新 summary 返回值**

在 `return` 的 `summary` 字典中添加：

```python
                'input_mode': 'custom_genes' if custom_genes_str else ('split_direction' if split_direction else 'standard'),
```

- [ ] **Step 4: 在 `routes/analysis.py` 的 `PARAM_SCHEMAS['bulk_enrichment']` 末尾添加新参数**

在第 152 行 `{'key': 'input_source', ...}` 之后添加：

```python
        {'key': 'split_direction', 'label': '分开分析上调/下调基因', 'type': 'checkbox', 'default': False, 'help': '开启后将 DEG 结果按 Up/Down 拆分，分别做 ORA 富集分析，生成独立的气泡图。'},
        {'key': 'custom_genes', 'label': '自定义基因列表（可选）', 'type': 'textarea', 'default': '', 'help': '手动输入基因名，逗号或换行分隔。填写后忽略 DEG 结果文件，直接用此列表做 ORA。'},
```

- [ ] **Step 5: 验证模块可导入**

```bash
cd /data/GJ/platform && python -c "from modules.bulk_enrichment import BulkEnrichmentAnalysis; print('OK')"
```

Expected: `OK`

- [ ] **Step 6: Commit**

```bash
git add modules/bulk_enrichment.py routes/analysis.py
git commit -m "feat(bulk_enrichment): add split_direction up/down enrichment and custom_genes input"
```

---

## Task 3: bulk_deg — 多组两两比较 + 自定义合并组

**Files:**
- Modify: `modules/bulk_deg.py`
- Modify: `routes/analysis.py`

- [ ] **Step 1: 在 `modules/bulk_deg.py` 的 `run()` 方法开头添加参数解析辅助函数**

在文件顶部（第 6 行 `from modules.base import BaseAnalysis` 之后）添加：

```python


def _parse_comparisons(comp_str):
    """解析 'A-vs-B;C-vs-D' 为 [('A','B'), ('C','D')]"""
    if not comp_str or not comp_str.strip():
        return []
    pairs = []
    for item in comp_str.replace('\n', ';').split(';'):
        item = item.strip()
        if '-vs-' in item:
            parts = item.split('-vs-')
            if len(parts) == 2 and parts[0].strip() and parts[1].strip():
                pairs.append((parts[0].strip(), parts[1].strip()))
    return pairs


def _parse_custom_groups(cg_str):
    """解析 'High=Treated_1h+Treated_3h\nLow=Ctrl' 为 {'High': ['Treated_1h','Treated_3h'], 'Low': ['Ctrl']}"""
    if not cg_str or not cg_str.strip():
        return {}
    mapping = {}
    for line in cg_str.strip().split('\n'):
        line = line.strip()
        if '=' not in line:
            continue
        name, expr = line.split('=', 1)
        name = name.strip()
        members = [m.strip() for m in expr.split('+') if m.strip()]
        if name and members:
            mapping[name] = members
    return mapping
```

- [ ] **Step 2: 在 `run()` 方法中，参数读取之后、构建计数矩阵之前添加 custom_groups 和 comparisons 逻辑**

在第 31 行 `top_n = int(self.params.get('top_n', 20))` 之后，第 33 行 `self.progress(15, "构建计数矩阵...")` 之前插入：

```python
        comparisons_str = self.params.get('comparisons', '').strip()
        custom_groups_str = self.params.get('custom_groups', '').strip()
        custom_groups = _parse_custom_groups(custom_groups_str)
        comparison_pairs = _parse_comparisons(comparisons_str)

        # 自定义合并组：在 obs 中创建临时列
        if custom_groups and groupby in adata.obs.columns:
            new_col = '_custom_group'
            adata.obs[new_col] = adata.obs[groupby].astype(str)
            for new_name, members in custom_groups.items():
                mask = adata.obs[groupby].astype(str).isin(members)
                adata.obs.loc[mask, new_col] = new_name
            groupby = new_col
```

- [ ] **Step 3: 将现有的单次 DEG 逻辑提取为 `_run_single_comparison` 辅助函数**

在 `_parse_custom_groups` 函数之后添加：

```python


def _run_single_comparison(adata, counts, group1_samples, group2_samples, group1, group2,
                           method, fc_threshold, pval_threshold, top_n, gene_names, gene_id_to_name,
                           plots_dir, results_dir, suffix=''):
    """执行单次 DEG 比较，返回 (deg_df, result_files, n_up, n_down)"""
    import omicverse as ov
    import plotly.graph_objects as go

    count_df = pd.DataFrame(counts.T, index=adata.var_names.tolist(), columns=adata.obs.index.tolist())
    dds = ov.bulk.pyDEG(count_df)
    dds.drop_duplicates_index()
    method_map = {'t-test': 'ttest', 'mann-whitney': 'wilcox', 'deseq2': 'DEseq2'}
    ov_method = method_map.get(method, 'ttest')
    dds.normalize()
    result = dds.deg_analysis(group1_samples, group2_samples, method=ov_method)

    gene_ids_list = result.index.tolist()
    if gene_id_to_name:
        gene_names_mapped = [gene_id_to_name.get(g, g) for g in gene_ids_list]
    else:
        gene_names_mapped = gene_ids_list
    log2fc = result['log2FC'].values
    pvalues = result['pvalue'].values
    padj = result['qvalue'].values
    n_genes = len(gene_names_mapped)

    log2fc_threshold = np.log2(fc_threshold)
    regulation = []
    for i in range(n_genes):
        if padj[i] < pval_threshold and log2fc[i] >= log2fc_threshold:
            regulation.append('Up')
        elif padj[i] < pval_threshold and log2fc[i] <= -log2fc_threshold:
            regulation.append('Down')
        else:
            regulation.append('NS')

    mask1_arr = np.array([s in group1_samples for s in adata.obs.index])
    mask2_arr = np.array([s in group2_samples for s in adata.obs.index])
    mean1 = counts[mask1_arr].mean(axis=0)
    mean2 = counts[mask2_arr].mean(axis=0)

    deg_df = pd.DataFrame({
        'gene': gene_names_mapped, 'log2FC': np.round(log2fc, 4),
        'pvalue': pvalues, 'padj': padj,
        'mean_group1': np.round(mean1[:n_genes], 2),
        'mean_group2': np.round(mean2[:n_genes], 2),
        'regulation': regulation
    }).sort_values('padj')

    n_up = sum(1 for r in regulation if r == 'Up')
    n_down = sum(1 for r in regulation if r == 'Down')
    result_files = []

    # 火山图
    neg_log_padj = -np.log10(padj + 1e-300)
    color_map = {'Up': '#e53935', 'Down': '#1a237e', 'NS': '#bdbdbd'}
    fig_vol = go.Figure()
    for reg in ['NS', 'Up', 'Down']:
        idx = [i for i in range(n_genes) if regulation[i] == reg]
        fig_vol.add_trace(go.Scattergl(
            x=log2fc[idx], y=neg_log_padj[idx], mode='markers',
            marker=dict(color=color_map[reg], size=5, opacity=0.7),
            name=f'{reg} ({len(idx)})',
            text=[gene_names_mapped[i] for i in idx],
            hovertemplate='%{text}<br>log2FC: %{x:.2f}<br>-log10(padj): %{y:.2f}'
        ))
    fig_vol.add_hline(y=-np.log10(pval_threshold), line_dash='dash', line_color='gray')
    fig_vol.add_vline(x=log2fc_threshold, line_dash='dash', line_color='gray')
    fig_vol.add_vline(x=-log2fc_threshold, line_dash='dash', line_color='gray')
    top_genes_vol = deg_df[deg_df['regulation'] != 'NS'].head(top_n)
    for _, row in top_genes_vol.iterrows():
        fig_vol.add_annotation(
            x=row['log2FC'], y=-np.log10(max(row['padj'], 1e-300)),
            text=row['gene'], showarrow=True, arrowhead=2,
            font=dict(size=9, color='#333'), ax=20, ay=-30
        )
    fig_vol.update_layout(
        title=f'火山图 ({group1} vs {group2})',
        xaxis_title='log2(Fold Change)', yaxis_title='-log10(padj)',
        plot_bgcolor='white', width=700, height=500
    )
    safe_suffix = suffix.replace(' ', '_').replace('/', '_') if suffix else ''
    vol_file = f'bulk_deg_volcano_{safe_suffix}.json' if safe_suffix else 'bulk_deg_volcano.json'
    fpath = os.path.join(plots_dir, vol_file)
    with open(fpath, 'w') as f: f.write(fig_vol.to_json(engine="json"))
    result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'volcano', 'label': f'火山图 ({suffix})' if suffix else '火山图'})

    # 独立 CSV
    csv_name = f'bulk_deg_{safe_suffix}.csv' if safe_suffix else 'bulk_deg_results.csv'
    deg_df.to_csv(os.path.join(results_dir, csv_name), index=False)
    result_files.append({'file_path': os.path.join(results_dir, csv_name), 'file_type': 'csv', 'category': 'table', 'label': f'DEG ({suffix})' if suffix else '差异表达基因列表'})

    return deg_df, result_files, n_up, n_down
```

- [ ] **Step 4: 重写 `run()` 方法中的 DEG 计算和图表生成部分**

将第 84 行 `self.progress(25, "OmicVerse 归一化 + 差异分析...")` 之后到第 243 行（`top_genes.to_csv` 之后）的代码替换为：

```python
        if comparison_pairs:
            # 多组比较模式
            all_results = []
            all_files = []
            total_up = 0
            total_down = 0
            for i, (g1, g2) in enumerate(comparison_pairs):
                self.progress(25 + int(55 * i / len(comparison_pairs)), f"比较 {i+1}/{len(comparison_pairs)}: {g1} vs {g2}")
                if groupby in adata.obs.columns:
                    s1 = list(adata.obs.index[adata.obs[groupby] == g1])
                    s2 = list(adata.obs.index[adata.obs[groupby] == g2])
                else:
                    raise ValueError(f"分组列 '{groupby}' 中未找到组 '{g1}' 或 '{g2}'")
                if not s1 or not s2:
                    continue
                suffix = f'{g1}_vs_{g2}'
                deg_df_i, files_i, n_up_i, n_down_i = _run_single_comparison(
                    adata, counts, s1, s2, g1, g2, method, fc_threshold, pval_threshold,
                    top_n, gene_names, gene_id_to_name, plots_dir, results_dir, suffix
                )
                deg_df_i['comparison'] = f'{g1} vs {g2}'
                all_results.append(deg_df_i)
                all_files.extend(files_i)
                total_up += n_up_i
                total_down += n_down_i

            result_files = all_files
            n_up = total_up
            n_down = total_down

            # 合并所有比较结果
            if all_results:
                combined = pd.concat(all_results, ignore_index=True)
                csv_path = os.path.join(results_dir, 'bulk_deg_all_comparisons.csv')
                combined.to_csv(csv_path, index=False)
                result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': '所有比较结果合并'})
                deg_df = combined
        else:
            # 单次比较模式（原有逻辑）
            if groupby in adata.obs.columns:
                groups = adata.obs[groupby].unique().tolist()
                if not group1 or group1 not in groups:
                    group1 = groups[0]
                if group2 == 'rest' or (not group2 or group2 not in groups):
                    group2_samples = [s for s in adata.obs.index if adata.obs.loc[s, groupby] != group1]
                    group2 = f'rest (n={len(group2_samples)})'
                else:
                    group2_samples = list(adata.obs.index[adata.obs[groupby] == group2])
                group1_samples = list(adata.obs.index[adata.obs[groupby] == group1])
            else:
                n = counts.shape[0]
                half = n // 2
                group1_samples = list(adata.obs.index[:half])
                group2_samples = list(adata.obs.index[half:])
                group1, group2 = "Group1", "Group2"

            self.progress(25, "OmicVerse 归一化 + 差异分析...")
            deg_df, result_files, n_up, n_down = _run_single_comparison(
                adata, counts, group1_samples, group2_samples, group1, group2,
                method, fc_threshold, pval_threshold, top_n, gene_names, gene_id_to_name,
                plots_dir, results_dir
            )
```

- [ ] **Step 5: 在 `routes/analysis.py` 的 `PARAM_SCHEMAS['bulk_deg']` 末尾添加新参数**

在第 128 行 `{'key': 'plot_genes', ...}` 之后添加：

```python
        {'key': 'comparisons', 'label': '多组比较（可选）', 'type': 'text', 'default': '', 'help': '多个比较用分号分隔，格式：A-vs-B;C-vs-D。填写后实验组/对照组参数被忽略。示例：DrugA-vs-Control;DrugB-vs-Control;DrugA-vs-DrugB'},
        {'key': 'custom_groups', 'label': '自定义合并组（可选）', 'type': 'textarea', 'default': '', 'help': '每行一个定义，格式：新组名=原组1+原组2。示例：High=Treated_1h+Treated_3h。定义后可在多组比较中使用新组名。'},
```

- [ ] **Step 6: 验证模块可导入**

```bash
cd /data/GJ/platform && python -c "from modules.bulk_deg import BulkDEGAnalysis; print('OK')"
```

Expected: `OK`

- [ ] **Step 7: Commit**

```bash
git add modules/bulk_deg.py routes/analysis.py
git commit -m "feat(bulk_deg): add multi-pairwise comparisons and custom merged groups"
```

---

## Task 4: bulk_timecourse — 分组配对比较

**Files:**
- Modify: `modules/bulk_timecourse.py`
- Modify: `routes/analysis.py`

- [ ] **Step 1: 在 `modules/bulk_timecourse.py` 的 `run()` 方法中添加 pairwise_groups 参数读取**

在第 109 行 `fdr_threshold = float(self.params.get('fdr_threshold', 0.05))` 之后添加：

```python
        pairwise_groups_str = self.params.get('pairwise_groups', '').strip()
```

- [ ] **Step 2: 在交互效应分析之前（第 350 行 `self.progress(85, "交互效应分析...")` 之前）添加 pairwise 比较逻辑**

在第 348 行 `result_files.append(...)` 之后、第 350 行 `self.progress(85, ...)` 之前插入：

```python
        # Pairwise group comparison at each timepoint
        if pairwise_groups_str and group_column and group_column in adata.obs.columns:
            if '-vs-' in pairwise_groups_str:
                parts = pairwise_groups_str.split('-vs-')
                pair_a, pair_b = parts[0].strip(), parts[1].strip()
            else:
                pair_a, pair_b = '', ''

            if pair_a and pair_b:
                self.progress(82, f"分组配对比较: {pair_a} vs {pair_b}...")
                from scipy import stats as sp_stats

                time_unique_sorted = np.sort(unique_times)
                pairwise_rows = []
                for t in time_unique_sorted:
                    mask_t = time_vals == t
                    mask_a = mask_t & (adata.obs[group_column].astype(str) == pair_a)
                    mask_b = mask_t & (adata.obs[group_column].astype(str) == pair_b)
                    if mask_a.sum() < 2 or mask_b.sum() < 2:
                        continue
                    expr_a = lognorm[mask_a, :]
                    expr_b = lognorm[mask_b, :]
                    for gi in range(n_genes):
                        tstat, pval = sp_stats.ttest_ind(expr_a[:, gi], expr_b[:, gi], equal_var=False)
                        mean_a = expr_a[:, gi].mean()
                        mean_b = expr_b[:, gi].mean()
                        log2fc_val = mean_a - mean_b  # already log2 CPM
                        pairwise_rows.append({
                            'gene': gene_names[gi], 'time': t,
                            'tstat': round(float(tstat), 4), 'pvalue': float(pval),
                            'log2FC': round(float(log2fc_val), 4),
                            'mean_a': round(float(mean_a), 4), 'mean_b': round(float(mean_b), 4),
                        })

                if pairwise_rows:
                    pw_df = pd.DataFrame(pairwise_rows)
                    # BH FDR per timepoint
                    for t in pw_df['time'].unique():
                        mask_t = pw_df['time'] == t
                        pvals_t = pw_df.loc[mask_t, 'pvalue'].values
                        try:
                            _, qvals_t, _, _ = multipletests(pvals_t, method='fdr_bh')
                        except Exception:
                            qvals_t = pvals_t
                        pw_df.loc[mask_t, 'qvalue'] = qvals_t

                    pw_df['sig'] = (pw_df['qvalue'] < fdr_threshold).astype(int)
                    pw_csv = os.path.join(results_dir, 'timecourse_pairwise_results.csv')
                    pw_df.to_csv(pw_csv, index=False)
                    result_files.append({'file_path': pw_csv, 'file_type': 'csv', 'category': 'table', 'label': f'{pair_a} vs {pair_b} 配对比较'})

                    # Heatmap: -log10(qvalue) per gene x time
                    pw_pivot = pw_df.pivot(index='gene', columns='time', values='qvalue').fillna(1.0)
                    pw_log10 = -np.log10(pw_pivot.clip(lower=1e-300))
                    fig_pw = go.Figure(data=go.Heatmap(
                        z=pw_log10.values.tolist(),
                        x=[str(t) for t in pw_pivot.columns.tolist()],
                        y=pw_pivot.index.tolist(),
                        colorscale='Reds', colorbar=dict(title='-log10(q)')
                    ))
                    fig_pw.update_layout(
                        title=f'{pair_a} vs {pair_b}: -log10(qvalue) per timepoint',
                        xaxis_title=time_column, yaxis_title='Gene',
                        width=700, height=max(400, len(pw_pivot) * 12 + 100),
                        yaxis=dict(tickfont=dict(size=7))
                    )
                    fpath = os.path.join(plots_dir, 'timecourse_pairwise_heatmap.json')
                    with open(fpath, 'w') as f: f.write(fig_pw.to_json(engine="json"))
                    result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': f'{pair_a} vs {pair_b} qvalue 热图'})

                    n_pw_sig = int(pw_df['sig'].sum())
                    n_pw_total = len(pw_df)
```

- [ ] **Step 3: 更新 summary 返回值**

在 `return` 的 `summary` 字典中添加：

```python
                'pairwise_comparison': pairwise_groups_str if pairwise_groups_str else None,
                'n_pairwise_sig': n_pw_sig if pairwise_groups_str and 'n_pw_sig' in dir() else 0,
```

- [ ] **Step 4: 在 `routes/analysis.py` 的 `PARAM_SCHEMAS['bulk_timecourse']` 末尾添加新参数**

在第 164 行 `{'key': 'fdr_threshold', ...}` 之后添加：

```python
        {'key': 'pairwise_groups', 'label': '分组配对比较（可选）', 'type': 'text', 'default': '', 'help': '格式：GroupA-vs-GroupB。在每个时间点对两组做 Welch t-test，生成时序差异热图。需同时填写分组列名。'},
```

- [ ] **Step 5: 验证模块可导入**

```bash
cd /data/GJ/platform && python -c "from modules.bulk_timecourse import BulkTimecourseAnalysis; print('OK')"
```

Expected: `OK`

- [ ] **Step 6: Commit**

```bash
git add modules/bulk_timecourse.py routes/analysis.py
git commit -m "feat(bulk_timecourse): add pairwise group comparison at each timepoint"
```

---

## Task 5: 模板 textarea 支持 + 前端验证

**Files:**
- Modify: `templates/analysis_select.html`

- [ ] **Step 1: 在 `templates/analysis_select.html` 第 93 行 `checkbox` 块之后、第 94 行 `number` 块之前添加 textarea 渲染**

在第 93 行 `{% if param.default %}checked{% endif %}>` 之后、第 94 行 `{% elif param.type == 'number' %}` 之前插入：

```html
                    {% elif param.type == 'textarea' %}
                    <textarea name="{{ param.key }}" class="form-control" id="param-{{ param.key }}"
                              rows="3" placeholder="{{ param.help }}">{{ param.default }}</textarea>
```

- [ ] **Step 2: 验证应用可启动**

```bash
cd /data/GJ/platform && python -c "from app import create_app; app = create_app(); print('app OK')"
```

Expected: `app OK`

- [ ] **Step 3: 验证 textarea 渲染**

```bash
cd /data/GJ/platform && python -c "
from app import create_app
from models import Project
from config import Config
import os
app = create_app()
with app.test_client() as c:
    p = Project(name='test')
    proj_dir = os.path.join(Config.DATA_DIR, 'projects', p.id)
    for sub in ['uploads', 'intermediate', 'results', 'plots']:
        os.makedirs(os.path.join(proj_dir, sub), exist_ok=True)
    p.save()
    resp = c.get(f'/projects/{p.id}/analyze/bulk_deg')
    html = resp.data.decode()
    assert 'textarea' in html, 'textarea not rendered'
    assert 'comparisons' in html, 'comparisons param missing'
    assert 'custom_groups' in html, 'custom_groups param missing'
    print('textarea rendering OK')
    import shutil
    shutil.rmtree(proj_dir)
    p.delete()
"
```

Expected: `textarea rendering OK`

- [ ] **Step 4: Commit**

```bash
git add templates/analysis_select.html
git commit -m "feat: add textarea parameter type support for multi-line input"
```

---

## Task 6: 端到端验证

- [ ] **Step 1: 验证所有模块可导入**

```bash
cd /data/GJ/platform && python -c "
from modules.bulk_heatmap import BulkHeatmapAnalysis
from modules.bulk_enrichment import BulkEnrichmentAnalysis
from modules.bulk_deg import BulkDEGAnalysis
from modules.bulk_timecourse import BulkTimecourseAnalysis
print('All 4 modules loaded successfully')
"
```

Expected: `All 4 modules loaded successfully`

- [ ] **Step 2: 验证 PARAM_SCHEMAS 包含所有新参数**

```bash
cd /data/GJ/platform && python -c "
from routes.analysis import PARAM_SCHEMAS
assert 'custom_genes' in [p['key'] for p in PARAM_SCHEMAS['bulk_heatmap']]
assert 'split_direction' in [p['key'] for p in PARAM_SCHEMAS['bulk_enrichment']]
assert 'custom_genes' in [p['key'] for p in PARAM_SCHEMAS['bulk_enrichment']]
assert 'comparisons' in [p['key'] for p in PARAM_SCHEMAS['bulk_deg']]
assert 'custom_groups' in [p['key'] for p in PARAM_SCHEMAS['bulk_deg']]
assert 'pairwise_groups' in [p['key'] for p in PARAM_SCHEMAS['bulk_timecourse']]
print('All new params present in PARAM_SCHEMAS')
"
```

Expected: `All new params present in PARAM_SCHEMAS`

- [ ] **Step 3: 验证 bulk_deg 辅助函数**

```bash
cd /data/GJ/platform && python -c "
from modules.bulk_deg import _parse_comparisons, _parse_custom_groups
assert _parse_comparisons('A-vs-B;C-vs-D') == [('A','B'), ('C','D')]
assert _parse_comparisons('') == []
assert _parse_custom_groups('High=Treated_1h+Treated_3h\nLow=Ctrl') == {'High': ['Treated_1h', 'Treated_3h'], 'Low': ['Ctrl']}
assert _parse_custom_groups('') == {}
print('Helper functions OK')
"
```

Expected: `Helper functions OK`

- [ ] **Step 4: 提交（如果有修复）**

如有任何修复，在此提交。否则跳过。

```bash
git add -A
git commit -m "fix: address issues found during end-to-end verification"
```
