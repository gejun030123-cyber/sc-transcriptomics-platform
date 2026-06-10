# Bulk 热图增强第二阶段实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development

**Goal:** 在第一阶段基础上补充 8 个高优先级功能：DEG 来源选择器、上调/下调分开、基因注释条、相关性热图参数化、筛选预设模板、clip 校验、log 开关、注释条配色。

**Architecture:** 在现有 bulk_heatmap.py 的 run() 方法中增量添加逻辑，新增参数到 PARAM_SCHEMAS。

---

## Task 1: 新增 PARAM_SCHEMAS + clip_range 校验 + log_transform

**Files:**
- Modify: `routes/analysis.py` — 在 bulk_heatmap 列表末尾添加 7 个参数
- Modify: `modules/bulk_heatmap.py` — clip_range 校验增强 + log_transform 逻辑

在 PARAM_SCHEMAS 的 bulk_heatmap 列表末尾（`groupby` 之后）追加：

```python
        # === 第二阶段新增 ===
        {'key': 'deg_comparison', 'label': 'DEG 比较选择', 'type': 'text', 'default': '',
         'help': '指定使用哪组 DEG 结果（如 results_0），留空取最新。gene_import_source=deg 时生效。'},
        {'key': 'up_down_separate', 'label': '上下调分开排列', 'type': 'checkbox', 'default': False,
         'help': '上调基因在上、下调在下，中间留间隙。'},
        {'key': 'gene_annotation_columns', 'label': '基因注释列名', 'type': 'text', 'default': '',
         'help': '从 adata.var 中读取基因属性作为注释条，逗号分隔（如 gene_type,pathway）。'},
        {'key': 'log_transform', 'label': 'Log2 转换', 'type': 'select',
         'options': ['auto', 'yes', 'no'], 'default': 'auto',
         'help': 'auto=根据数据是否已标准化自动判断，yes=强制 log2(x+pseudocount)，no=不转换。'},
        {'key': 'corr_method', 'label': '相关性方法', 'type': 'select',
         'options': ['pearson', 'spearman'], 'default': 'pearson',
         'help': '样本相关性热图的计算方法。'},
        {'key': 'corr_colorscale', 'label': '相关性热图配色', 'type': 'select',
         'options': ['Blues', 'RdBu_r', 'viridis', 'YlOrRd'], 'default': 'Blues',
         'help': '相关性热图的颜色方案。'},
        {'key': 'annotation_palette', 'label': '注释条自定义配色', 'type': 'text', 'default': '',
         'help': '格式：ctrl=#1565c0,hmc3=#e53935。留空使用默认配色。'},
```

clip_range 校验增强：在 bulk_heatmap.py 的 clip_range 解析处替换为：

```python
        clip_str = self.params.get('clip_range', '-3,3').strip()
        clip_range = None
        if clip_str:
            parts = [x.strip() for x in clip_str.split(',') if x.strip()]
            if len(parts) != 2:
                raise ValueError(f"clip_range 格式错误: '{clip_str}'，应为 'min,max'，如 '-3,3'")
            try:
                clip_range = (float(parts[0]), float(parts[1]))
            except ValueError:
                raise ValueError(f"clip_range 数值解析失败: '{clip_str}'")
```

log_transform 逻辑替换：将现有的 pseudocount 判断替换为：

```python
        log_transform = self.params.get('log_transform', 'auto')
        do_log = False
        if log_transform == 'yes':
            do_log = True
        elif log_transform == 'auto':
            do_log = 'normalization' not in adata.uns
        if do_log:
            heat_data = np.log2(heat_data + pseudocount)
```

---

## Task 2: DEG 来源选择器 + 上调/下调分开排列

**Files:**
- Modify: `modules/bulk_heatmap.py`

**deg_comparison 支持：** 在 deg 分支中，读取 `deg_comparison` 参数匹配具体文件：

```python
        elif gene_import_source == 'deg':
            deg_comparison = self.params.get('deg_comparison', '').strip()
            # ... existing deg_files scan ...
            if deg_comparison:
                matched = [f for f in deg_files if deg_comparison in f.replace('bulk_deg_', '').replace('.csv', '')]
                if matched:
                    deg_files = matched
            # ... rest of existing deg logic using deg_files[0] ...
```

**up_down_separate 支持：** 在聚类排序后、热图渲染前添加：

```python
        up_down_separate = self.params.get('up_down_separate', False)
        if isinstance(up_down_separate, str):
            up_down_separate = up_down_separate.lower() in ('true', '1', 'yes', 'on')

        if up_down_separate and gene_import_source == 'deg':
            # 获取每个基因的 regulation 信息
            deg_df_for_reg = pd.read_csv(os.path.join(self.project_dir, 'results', deg_files[0]))
            gene_reg_map = dict(zip(deg_df_for_reg['gene'], deg_df_for_reg['regulation']))
            up_genes = [g for g in gene_ordered if gene_reg_map.get(g) == 'Up']
            down_genes = [g for g in gene_ordered if gene_reg_map.get(g) == 'Down']
            other_genes = [g for g in gene_ordered if g not in up_genes and g not in down_genes]
            if up_genes or down_genes:
                gene_ordered = up_genes + down_genes + other_genes
                # 在上调和下调之间插入 NaN 空隙行
                gap_row = [float('nan')] * len(gene_ordered)
                up_data = heat_ordered[:, [gene_labels.index(g) for g in up_genes]] if up_genes else np.empty((heat_ordered.shape[0], 0))
                down_data = heat_ordered[:, [gene_labels.index(g) for g in down_genes]] if down_genes else np.empty((heat_ordered.shape[0], 0))
                other_data = heat_ordered[:, [gene_labels.index(g) for g in other_genes]] if other_genes else np.empty((heat_ordered.shape[0], 0))
                gap = np.full((heat_ordered.shape[0], 1), np.nan)
                parts = [p for p in [up_data, gap, down_data, other_data] if p.shape[1] > 0]
                if len(parts) > 1:
                    heat_ordered = np.hstack(parts)
                    gap_idx = up_data.shape[1] if up_genes else 0
                    gene_ordered = up_genes + ['---'] + down_genes + other_genes
```

---

## Task 3: 基因维度注释条 + 自定义注释条配色

**Files:**
- Modify: `modules/bulk_heatmap.py`

在现有样本注释条代码之后，添加基因注释条逻辑：

```python
        # 基因维度注释条
        gene_annot_cols_str = self.params.get('gene_annotation_columns', '').strip()
        gene_annot_cols = [c.strip() for c in gene_annot_cols_str.split(',') if c.strip()]
        if gene_annot_cols:
            for col_name in gene_annot_cols:
                if col_name not in adata.var.columns:
                    continue
                gene_values = [str(adata.var.loc[g, col_name]) if g in adata.var.index else 'NA'
                               for g in gene_ordered if g != '---']
                uniq = sorted(set(gene_values))
                color_map = {g: DEFAULT_PALETTE[i % len(DEFAULT_PALETTE)] for i, g in enumerate(uniq)}
                color_indices = [uniq.index(v) for v in gene_values]
                fig_ga = go.Figure()
                n_groups = len(uniq)
                cs = [[0, list(color_map.values())[0]]] if n_groups <= 1 else \
                     [[i / (n_groups - 1), color_map[g]] for i, g in enumerate(uniq)]
                fig_ga.add_trace(go.Heatmap(
                    z=[color_indices], x=gene_values, y=[col_name],
                    colorscale=cs, showscale=False,
                    text=[gene_values], hovertemplate='%{x}: %{text}<extra></extra>'
                ))
                fig_ga.update_layout(
                    height=60, width=max(600, len(gene_values) * 12 + 200),
                    margin=dict(l=0, r=0, t=5, b=0)
                )
                save_plotly_json(fig_ga, plots_dir, f'bulk_heatmap_gene_annot_{col_name}.json',
                                result_files, category='annotation', label=f'{col_name} 基因注释条')
```

自定义配色：解析 `annotation_palette` 参数替换 DEFAULT_PALETTE：

```python
        annotation_palette_str = self.params.get('annotation_palette', '').strip()
        custom_palette = None
        if annotation_palette_str:
            custom_palette = {}
            for pair in annotation_palette_str.split(','):
                if '=' in pair:
                    k, v = pair.split('=', 1)
                    custom_palette[k.strip()] = v.strip()
```

在 `build_annotation_bar` 调用后，若 `custom_palette` 非空，替换对应分组的颜色：

```python
            if custom_palette:
                for col_name, col_info in annot_data.items():
                    for group, color in custom_palette.items():
                        if group in col_info['color_map']:
                            col_info['color_map'][group] = color
                    col_info['colors'] = [col_info['color_map'][v] for v in col_info['groups']]
```

---

## Task 4: 样本相关性热图参数化 + 筛选表达式预设模板

**Files:**
- Modify: `modules/bulk_heatmap.py` — 相关性热图参数化
- Modify: `templates/analysis_select.html` — 预设模板按钮

相关性热图：替换硬编码的 Pearson + Blues 为参数驱动：

```python
        self.progress(85, "生成样本相关性热图...")
        corr_method = self.params.get('corr_method', 'pearson')
        corr_colorscale = self.params.get('corr_colorscale', 'Blues')

        if corr_method == 'spearman':
            from scipy.stats import spearmanr
            corr_matrix, _ = spearmanr(norm_data, axis=1)
            if corr_matrix.ndim == 0:
                corr_matrix = np.array([[1.0]])
        else:
            corr_matrix = np.corrcoef(norm_data)

        fig_corr = go.Figure()
        fig_corr.add_trace(go.Heatmap(
            z=corr_matrix.tolist(),
            x=sample_labels, y=sample_labels,
            colorscale=corr_colorscale,
            zmin=0 if corr_colorscale == 'Blues' else None,
            zmax=1 if corr_colorscale == 'Blues' else None,
            colorbar=dict(title=f'{corr_method.capitalize()} r'),
            hovertemplate='%{y} vs %{x}<br>r = %{z:.3f}<extra></extra>'
        ))
        fig_corr.update_layout(
            title=f'样本相关性热图 ({corr_method.capitalize()})',
            height=max(400, adata.n_obs * 30 + 100),
            width=max(400, adata.n_obs * 30 + 100),
            plot_bgcolor='white'
        )
        save_plotly_json(fig_corr, plots_dir, 'bulk_corr_heatmap.json', result_files,
                        category='heatmap', label='样本相关性热图')
```

筛选表达式预设模板：在 analysis_select.html 的 `filter-expr-toolbar` 中扩展按钮：

```html
<button type="button" class="btn btn-outline-secondary btn-sm filter-tpl" data-expr="ALL:up">ALL:up</button>
<button type="button" class="btn btn-outline-secondary btn-sm filter-tpl" data-expr="ALL:down">ALL:down</button>
<button type="button" class="btn btn-outline-secondary btn-sm filter-tpl" data-expr="ANY:up">ANY:up</button>
<button type="button" class="btn btn-outline-secondary btn-sm filter-tpl" data-expr="ANY:down">ANY:down</button>
```

---

## Task 5: 测试 + 提交

**Files:**
- Modify: `tests/test_heatmap_helpers.py` — 新增 clip_range 校验测试
- Commit + push
