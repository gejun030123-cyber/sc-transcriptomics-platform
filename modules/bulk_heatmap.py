import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class BulkHeatmapAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_heatmap"
    DISPLAY_NAME = "Bulk 热图可视化"
    DESCRIPTION = "Top 差异基因热图、样本相关性热图"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import plotly.graph_objects as go

        self.progress(5, "加载数据...")
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(input_path)

        top_n = int(self.params.get('top_n', 50))
        groupby = self.params.get('groupby', '')
        hm_type = self.params.get('heatmap_type', 'top_var')

        self.progress(20, "计算数据矩阵...")
        counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
        counts = counts.astype(float)
        # 清理 inf/NaN
        counts = np.nan_to_num(counts, nan=0.0, posinf=0.0, neginf=0.0)

        if 'normalization' not in adata.uns:
            sc.pp.normalize_total(adata, target_sum=1e6)
            sc.pp.log1p(adata)
        norm_data = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
        norm_data = norm_data.astype(float)
        norm_data = np.nan_to_num(norm_data, nan=0.0, posinf=0.0, neginf=0.0)

        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        self.progress(40, "选择基因...")

        # 向后兼容：从旧参数映射到新参数
        gene_import_source = self.params.get('gene_import_source', '')
        if not gene_import_source:
            if self.params.get('custom_genes', '').strip():
                gene_import_source = 'manual'
            else:
                gene_import_source = self.params.get('heatmap_type', 'top_var')

        top_n = int(self.params.get('top_n', 50))
        custom_genes_str = self.params.get('custom_genes', '').strip()
        var_metric = self.params.get('var_metric', 'var')

        if gene_import_source == 'manual' or (custom_genes_str and gene_import_source not in ('deg', 'expression_filter')):
            gene_list = [g.strip() for g in custom_genes_str.replace('\n', ',').split(',') if g.strip()]
            var_idx_map = {g: i for i, g in enumerate(adata.var_names)}
            top_idx = [var_idx_map[g] for g in gene_list if g in var_idx_map]
            not_found = [g for g in gene_list if g not in var_idx_map]
            if not top_idx:
                raise ValueError(f"自定义基因列表中没有找到任何匹配基因。请检查基因名是否正确。")
            title = f'自定义基因热图 ({len(top_idx)} genes)'
            if not_found:
                title += f'，{len(not_found)} 个未找到'

        elif gene_import_source == 'deg':
            deg_direction = self.params.get('deg_direction', 'both')
            deg_sortby = self.params.get('deg_sortby', 'padj')
            results_dir = os.path.join(self.project_dir, 'results')
            deg_files = sorted([f for f in os.listdir(results_dir)
                                if f.startswith('bulk_deg_results') and f.endswith('.csv')
                                and 'merged' not in f and 'all_comparisons' not in f
                                and 'lrt' not in f and 'top_genes' not in f])
            deg_comparison = self.params.get('deg_comparison', '').strip()
            if deg_comparison:
                matched = [f for f in deg_files if deg_comparison in f.replace('bulk_deg_', '').replace('.csv', '')]
                if matched:
                    deg_files = matched[:1]
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
            for name, df in comparisons.items():
                df_idx = df.set_index('gene')
                common = pd.Index(all_genes_union).intersection(df_idx.index)
                if len(common) > 0:
                    logfc_m.loc[common, name] = df_idx.loc[common, 'log2FC']
                    padj_m.loc[common, name] = df_idx.loc[common, 'padj']
            fc_thresh = float(self.params.get('fc_threshold', 2.0))
            pv_thresh = float(self.params.get('pval_threshold', 0.05))
            log2fc_t = np.log2(fc_thresh)
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

        heat_data = norm_data[:, top_idx].copy()
        gene_labels = [adata.var_names[i] for i in top_idx]
        sample_labels = adata.obs.index.tolist()

        # 数据变换
        from modules.visualization import transform_heatmap_data, cluster_heatmap
        row_scaling = self.params.get('row_scaling', 'zscore')
        pseudocount = float(self.params.get('pseudocount', 1))
        winsorize_param = self.params.get('winsorize', 'none')
        missing_value = self.params.get('missing_value', 'ignore')
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

        # log2 转换
        log_transform = self.params.get('log_transform', 'auto')
        do_log = False
        if log_transform == 'yes':
            do_log = True
        elif log_transform == 'auto':
            do_log = 'normalization' not in adata.uns
        if do_log:
            heat_data = np.log2(heat_data + pseudocount)

        heat_z = transform_heatmap_data(
            heat_data, row_scaling=row_scaling,
            winsorize=winsorize_param, clip_range=clip_range, missing_value=missing_value)

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

        # 上调/下调分开排列
        up_down_separate = self.params.get('up_down_separate', False)
        if isinstance(up_down_separate, str):
            up_down_separate = up_down_separate.lower() in ('true', '1', 'yes', 'on')

        if up_down_separate and gene_import_source == 'deg':
            results_dir_sep = os.path.join(self.project_dir, 'results')
            sep_deg_files = sorted([f for f in os.listdir(results_dir_sep)
                                    if f.startswith('bulk_deg_results') and f.endswith('.csv')
                                    and 'merged' not in f and 'all_comparisons' not in f
                                    and 'lrt' not in f and 'top_genes' not in f])
            deg_comparison_sep = self.params.get('deg_comparison', '').strip()
            if deg_comparison_sep:
                sep_matched = [f for f in sep_deg_files if deg_comparison_sep in f.replace('bulk_deg_', '').replace('.csv', '')]
                if sep_matched:
                    sep_deg_files = sep_matched[:1]
            if sep_deg_files:
                sep_deg_df = pd.read_csv(os.path.join(results_dir_sep, sep_deg_files[0]))
                gene_reg_map = dict(zip(sep_deg_df['gene'], sep_deg_df['regulation']))
                up_genes = [g for g in gene_ordered if gene_reg_map.get(g) == 'Up']
                down_genes = [g for g in gene_ordered if gene_reg_map.get(g) == 'Down']
                other_genes = [g for g in gene_ordered if g not in up_genes and g not in down_genes]
                if up_genes and down_genes:
                    # Build reordered data with gap between up and down
                    up_idx = [gene_labels.index(g) for g in up_genes]
                    down_idx = [gene_labels.index(g) for g in down_genes]
                    other_idx = [gene_labels.index(g) for g in other_genes]
                    gap_col = np.full((heat_ordered.shape[0], 1), np.nan)
                    parts = []
                    if up_idx:
                        parts.append(heat_ordered[:, up_idx])
                    parts.append(gap_col)
                    if down_idx:
                        parts.append(heat_ordered[:, down_idx])
                    if other_idx:
                        parts.append(heat_ordered[:, other_idx])
                    heat_ordered = np.hstack(parts)
                    gene_ordered = up_genes + ['---'] + down_genes + other_genes

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

        # 注释条（支持多列）
        from modules.visualization import build_annotation_bar, DEFAULT_PALETTE
        annot_cols_str = self.params.get('annotation_columns', '').strip()
        annot_cols = [c.strip() for c in annot_cols_str.split(',') if c.strip()]
        if groupby and groupby in adata.obs.columns and groupby not in annot_cols:
            annot_cols.insert(0, groupby)

        # 自定义注释条配色
        annotation_palette_str = self.params.get('annotation_palette', '').strip()
        custom_palette = {}
        if annotation_palette_str:
            for pair in annotation_palette_str.split(','):
                if '=' in pair:
                    k, v = pair.split('=', 1)
                    custom_palette[k.strip()] = v.strip()

        if annot_cols:
            annot_data, _ = build_annotation_bar(adata.obs, annot_cols,
                                                 sample_order=sample_order, palette=DEFAULT_PALETTE)
            # 应用自定义配色
            if custom_palette:
                for col_name_c, col_info_c in annot_data.items():
                    for group_c, color_c in custom_palette.items():
                        if group_c in col_info_c['color_map']:
                            col_info_c['color_map'][group_c] = color_c
                    col_info_c['colors'] = [col_info_c['color_map'][v] for v in col_info_c['groups']]
            for col_name, col_info in annot_data.items():
                color_indices = [col_info['unique'].index(g) for g in col_info['groups']]
                fig_annot = go.Figure()
                n_groups = len(col_info['unique'])
                if n_groups <= 1:
                    cs = [[0, list(col_info['color_map'].values())[0]]]
                else:
                    cs = [[i / (n_groups - 1), col_info['color_map'][g]]
                          for i, g in enumerate(col_info['unique'])]
                fig_annot.add_trace(go.Heatmap(
                    z=[[i] for i in color_indices],
                    y=sample_ordered, x=[col_name],
                    colorscale=cs, showscale=False,
                    text=[[col_info['groups'][j]] for j in range(len(col_info['groups']))],
                    hovertemplate='%{y}: %{text}<extra></extra>'
                ))
                fig_annot.update_layout(
                    height=max(400, len(sample_ordered) * 25 + 150), width=100,
                    margin=dict(l=0, r=0, t=30, b=40)
                )
                save_plotly_json(fig_annot, plots_dir, f'bulk_heatmap_annotation_{col_name}.json',
                                result_files, category='annotation', label=f'{col_name} 注释条')

        # 基因维度注释条
        gene_annot_cols_str = self.params.get('gene_annotation_columns', '').strip()
        gene_annot_cols = [c.strip() for c in gene_annot_cols_str.split(',') if c.strip()]
        if gene_annot_cols:
            for gcol in gene_annot_cols:
                if gcol not in adata.var.columns:
                    continue
                # Get gene values in the display order (excluding gap markers)
                display_genes = [g for g in gene_ordered if g != '---']
                gene_values = []
                for g in display_genes:
                    if g in adata.var.index:
                        gene_values.append(str(adata.var.loc[g, gcol]))
                    else:
                        gene_values.append('NA')
                uniq = sorted(set(gene_values))
                g_color_map = {gv: DEFAULT_PALETTE[i % len(DEFAULT_PALETTE)] for i, gv in enumerate(uniq)}
                color_indices = [uniq.index(v) for v in gene_values]
                fig_ga = go.Figure()
                n_g = len(uniq)
                g_cs = [[0, list(g_color_map.values())[0]]] if n_g <= 1 else \
                       [[i / (n_g - 1), g_color_map[gv]] for i, gv in enumerate(uniq)]
                fig_ga.add_trace(go.Heatmap(
                    z=[color_indices], x=gene_values, y=[gcol],
                    colorscale=g_cs, showscale=False,
                    text=[gene_values], hovertemplate='%{x}: %{text}<extra></extra>'
                ))
                fig_ga.update_layout(
                    height=60, width=max(600, len(gene_values) * 12 + 200),
                    margin=dict(l=0, r=0, t=5, b=0)
                )
                save_plotly_json(fig_ga, plots_dir, f'bulk_heatmap_gene_annot_{gcol}.json',
                                result_files, category='annotation', label=f'{gcol} 基因注释条')

        self.progress(85, "生成样本相关性热图...")
        corr_method = self.params.get('corr_method', 'pearson')
        corr_colorscale = self.params.get('corr_colorscale', 'Blues')

        if corr_method == 'spearman':
            from scipy.stats import spearmanr as sp_spearmanr
            corr_result = sp_spearmanr(norm_data, axis=1)
            corr_matrix = corr_result.correlation if hasattr(corr_result, 'correlation') else np.array([[1.0]])
            if np.ndim(corr_matrix) == 0:
                corr_matrix = np.array([[1.0]])
        else:
            corr_matrix = np.corrcoef(norm_data)

        # 颜色范围：Blues 用 [0,1]，diverging 色图不设限
        corr_zmin = 0 if corr_colorscale == 'Blues' else None
        corr_zmax = 1 if corr_colorscale == 'Blues' else None
        corr_zmid = 0 if corr_colorscale != 'Blues' else None

        corr_kwargs = dict(
            z=corr_matrix.tolist(),
            x=sample_labels, y=sample_labels,
            colorscale=corr_colorscale,
            colorbar=dict(title=f'{corr_method.capitalize()} r'),
            hovertemplate='%{y} vs %{x}<br>r = %{z:.3f}<extra></extra>'
        )
        if corr_zmin is not None:
            corr_kwargs['zmin'] = corr_zmin
        if corr_zmax is not None:
            corr_kwargs['zmax'] = corr_zmax
        if corr_zmid is not None:
            corr_kwargs['zmid'] = corr_zmid

        fig_corr = go.Figure()
        fig_corr.add_trace(go.Heatmap(**corr_kwargs))
        fig_corr.update_layout(
            title=f'样本相关性热图 ({corr_method.capitalize()})',
            height=max(400, adata.n_obs * 30 + 100),
            width=max(400, adata.n_obs * 30 + 100),
            plot_bgcolor='white'
        )
        save_plotly_json(fig_corr, plots_dir, 'bulk_corr_heatmap.json', result_files,
                        category='heatmap', label='样本相关性热图')

        self.progress(95, "保存结果...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_heatmap_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'heatmap_type': hm_type,
                'n_genes_shown': len(top_idx),
                'n_samples': adata.n_obs,
                'groupby': groupby or '无',
            }
        }
