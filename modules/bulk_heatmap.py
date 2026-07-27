import os
import json
import re
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


def _selection_values(value):
    """Parse comma/newline selections while keeping user-entered sample IDs exact."""
    return list(dict.fromkeys(
        item.strip() for item in str(value or '').replace('\n', ',').split(',') if item.strip()
    ))


def _comparison_groups(comparison_label):
    """Extract the two contrast groups from a human-readable ``A vs B`` label."""
    match = re.search(r'^\s*(.+?)\s+vs\s+(.+?)\s*$', str(comparison_label or ''), flags=re.IGNORECASE)
    return [match.group(1).strip(), match.group(2).strip()] if match else []


def _select_heatmap_display_samples(obs, sample_display_mode='all', groupby='',
                                    selected_groups='', selected_samples='',
                                    comparison_label=''):
    """Return the selected sample positions and an auditable display contract.

    This is deliberately a visualization filter: it never changes the DEG
    table or tests significance for a newly displayed group.
    """
    allowed_modes = {'all', 'deg_groups', 'selected_groups', 'selected_samples'}
    mode = str(sample_display_mode or 'all').strip()
    if mode not in allowed_modes:
        mode = 'all'
    sample_names = [str(name) for name in obs.index]
    if mode == 'all':
        return list(range(len(sample_names))), {
            'mode': 'all', 'n_samples': len(sample_names), 'selected_groups': [],
        }

    if mode == 'selected_samples':
        requested = _selection_values(selected_samples)
        if not requested:
            raise ValueError('请选择至少一个要展示的样本。')
        available = set(sample_names)
        missing = [name for name in requested if name not in available]
        if missing:
            raise ValueError(f"未找到所选样本：{', '.join(missing[:5])}")
        selected = [index for index, name in enumerate(sample_names) if name in set(requested)]
        if len(selected) < 2:
            raise ValueError('热图至少需要展示 2 个样本。')
        return selected, {
            'mode': mode, 'n_samples': len(selected), 'selected_samples': requested,
            'selected_groups': [],
        }

    if not groupby or groupby not in obs.columns:
        raise ValueError('按分组展示热图时，请填写有效的样本分组列名。')
    group_values = obs[groupby].astype(str).tolist()
    available_groups = set(group_values)
    if mode == 'deg_groups':
        requested = _comparison_groups(comparison_label)
        if not requested:
            raise ValueError('仅显示 DEG 两组时，请先选择带有“组1 vs 组2”名称的 DEG 比较结果。')
    else:
        requested = _selection_values(selected_groups)
        if not requested:
            raise ValueError('请选择至少一个要展示的分组。')
    missing = [name for name in requested if name not in available_groups]
    if missing:
        raise ValueError(f"分组列 '{groupby}' 中未找到：{', '.join(missing[:5])}")
    selected_set = set(requested)
    selected = [index for index, value in enumerate(group_values) if value in selected_set]
    if len(selected) < 2:
        raise ValueError('热图至少需要展示 2 个样本。')
    return selected, {
        'mode': mode, 'n_samples': len(selected), 'selected_groups': requested,
        'groupby': groupby,
    }


def _resolve_deg_result_path(results_dir, requested_comparison=''):
    """Resolve a user-selected DEG CSV while retaining legacy filename inputs."""
    candidates = sorted(
        [name for name in os.listdir(results_dir)
         if name.startswith('bulk_deg_results') and name.endswith('.csv')
         and 'merged' not in name and 'all_comparisons' not in name
         and 'lrt' not in name and 'top_genes' not in name],
        key=lambda name: os.path.getmtime(os.path.join(results_dir, name)), reverse=True,
    )
    if not candidates:
        raise ValueError('未找到 DEG 结果文件，请先运行 bulk_deg')
    requested = str(requested_comparison or '').strip()
    if requested:
        candidate_path = os.path.realpath(requested)
        results_real = os.path.realpath(results_dir)
        if (candidate_path.startswith(results_real + os.sep)
                and os.path.isfile(candidate_path)
                and os.path.basename(candidate_path) in candidates):
            return candidate_path
        exact = [name for name in candidates if requested == name]
        if exact:
            return os.path.join(results_dir, exact[0])
        matched = [name for name in candidates
                   if requested in name.replace('bulk_deg_', '').replace('.csv', '')]
        if matched:
            return os.path.join(results_dir, matched[0])
    return os.path.join(results_dir, candidates[0])


class BulkHeatmapAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_heatmap"
    DISPLAY_NAME = "Bulk 热图可视化"
    DESCRIPTION = "Top 差异基因热图、样本相关性热图"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        from modules.native_figures import (
            heatmap_figure, annotation_strip_figure,
            correlation_heatmap_figure, correlation_pairwise_table,
            summarize_correlation_pairs,
        )

        self.progress(5, "加载数据...")
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(input_path)
        from modules.io_utils import infer_expression_measurement

        top_n = int(self.params.get('top_n', 50))
        groupby = self.params.get('groupby', '')
        hm_type = self.params.get('heatmap_type', 'top_var')

        self.progress(20, "计算数据矩阵...")
        counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
        counts = counts.astype(float)
        # 清理 inf/NaN
        counts = np.nan_to_num(counts, nan=0.0, posinf=0.0, neginf=0.0)

        if 'normalization' not in adata.uns:
            input_measurement = infer_expression_measurement(adata, input_path)
            if input_measurement == 'raw_counts':
                sc.pp.normalize_total(adata, target_sum=1e6)
                sc.pp.log1p(adata)
                adata.uns['normalization'] = {'method': 'heatmap_auto_cpm_log1p', 'is_log_transformed': True}
            else:
                adata.X = np.log2(np.maximum(adata.X, 0) + 1)
                adata.uns['normalization'] = {'method': 'heatmap_auto_log2', 'is_log_transformed': True}
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
        selected_deg_path = ''
        selected_deg_label = str(self.params.get('deg_comparison_label', '') or '').strip()
        # This value is also used by the sample-display contract.  Initialize
        # it for top_var/manual/expression_filter modes where the DEG branch is
        # not entered.
        deg_comparison = str(self.params.get('deg_comparison', '') or '').strip()

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
            selected_deg_path = _resolve_deg_result_path(results_dir, deg_comparison)
            deg_df = pd.read_csv(selected_deg_path)
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
            contrast_suffix = f' · {selected_deg_label}' if selected_deg_label else ''
            title = f'Top {len(top_idx)} 差异基因热图 ({deg_direction}){contrast_suffix}'

        elif gene_import_source == 'expression_filter':
            from modules.expression_parser import validate, evaluate
            filter_expr = self.params.get('filter_expression', '').strip()
            if not filter_expr:
                raise ValueError("expression_filter 模式需要填写筛选表达式")
            results_dir = os.path.join(self.project_dir, 'results')
            deg_files = sorted([f for f in os.listdir(results_dir)
                                if f.startswith('bulk_deg_results') and f.endswith('.csv')
                                and 'merged' not in f and 'all_comparisons' not in f
                                and 'lrt' not in f and 'top_genes' not in f],
                               key=lambda f: os.path.getmtime(os.path.join(results_dir, f)),
                               reverse=True)
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

        all_sample_labels = adata.obs.index.tolist()
        sample_display_mode = self.params.get('sample_display_mode', 'all')
        sample_display_groups = self.params.get('sample_display_groups', '')
        sample_display_names = self.params.get('sample_display_names', '')
        display_indices, display_contract = _select_heatmap_display_samples(
            adata.obs,
            sample_display_mode=sample_display_mode,
            groupby=groupby,
            selected_groups=sample_display_groups,
            selected_samples=sample_display_names,
            comparison_label=selected_deg_label or deg_comparison,
        )
        display_obs = adata.obs.iloc[display_indices].copy()
        heat_data = norm_data[np.ix_(display_indices, top_idx)].copy()
        gene_labels = [adata.var_names[i] for i in top_idx]
        sample_labels = display_obs.index.tolist()

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

        if col_cluster == 'yes' and len(sample_labels) > 2:
            sample_order = cluster_heatmap(heat_z, method=col_method, metric=col_metric)
        elif col_cluster == 'group_order' and groupby and groupby in display_obs.columns:
            groups = display_obs[groupby].astype(str)
            group_order = sorted(groups.unique())
            sample_order = []
            for g in group_order:
                sample_order.extend([i for i in range(len(groups)) if groups.iloc[i] == g])
        else:
            # The display matrix may contain only selected groups/samples.
            sample_order = list(range(len(sample_labels)))

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
            if selected_deg_path:
                sep_deg_df = pd.read_csv(selected_deg_path)
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

        # 视觉样式参数（展示型热图统一走原生 Matplotlib；不再生成同名 Plotly JSON）
        zmin_str = self.params.get('zmin', 'auto').strip()
        zmax_str = self.params.get('zmax', 'auto').strip()
        zmin = float(zmin_str) if zmin_str and zmin_str != 'auto' else None
        zmax = float(zmax_str) if zmax_str and zmax_str != 'auto' else None

        show_gene_labels = self.params.get('show_gene_labels', 'all')
        show_sample_labels = self.params.get('show_sample_labels', 'auto')
        gene_font_size = int(self.params.get('gene_font_size', 8))
        sample_font_size = int(self.params.get('sample_font_size', 9))

        if row_scaling in ('zscore', 'center') and zmin is None and zmax is None:
            zmin, zmax = -3.0, 3.0
        if zmin is None:
            finite = heat_ordered[np.isfinite(heat_ordered)]
            zmin = float(np.nanmin(finite)) if finite.size else -1.0
        if zmax is None:
            finite = heat_ordered[np.isfinite(heat_ordered)]
            zmax = float(np.nanmax(finite)) if finite.size else 1.0
        if zmin == zmax:
            zmin, zmax = zmin - 1.0, zmax + 1.0

        fig = heatmap_figure(
            heat_ordered,
            x_labels=gene_ordered,
            y_labels=sample_ordered,
            title=title,
            x_label='Gene',
            y_label='Sample',
            colorbar_label='Z-score' if row_scaling == 'zscore' else 'Value',
            vmin=zmin,
            vmax=zmax,
        )
        ax = fig.axes[0]
        if show_gene_labels == 'top20':
            show_n = min(20, len(gene_ordered))
            ax.set_xticks(np.arange(show_n), [str(x) for x in gene_ordered[:show_n]])
        elif show_gene_labels == 'none':
            ax.set_xticks([])
        if show_sample_labels == 'none':
            ax.set_yticks([])
        elif show_sample_labels == 'auto':
            from modules.native_figures import apply_sample_tick_labels
            apply_sample_tick_labels(
                ax, sample_ordered, axis='y', max_labels=18,
                font_size=sample_font_size,
            )
        ax.tick_params(axis='x', labelsize=gene_font_size)
        if show_sample_labels == 'all':
            ax.tick_params(axis='y', labelsize=sample_font_size)
        result_files.extend(self.save_matplotlib_figure(
            fig, plots_dir, 'bulk_heatmap.png', 'heatmap', title,
            formats=('png', 'svg'), dpi=300,
        ))

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
            annot_data, _ = build_annotation_bar(display_obs, annot_cols,
                                                 sample_order=sample_order, palette=DEFAULT_PALETTE)
            # 应用自定义配色
            if custom_palette:
                for col_name_c, col_info_c in annot_data.items():
                    for group_c, color_c in custom_palette.items():
                        if group_c in col_info_c['color_map']:
                            col_info_c['color_map'][group_c] = color_c
                    col_info_c['colors'] = [col_info_c['color_map'][v] for v in col_info_c['groups']]
            for col_name, col_info in annot_data.items():
                fig_annot = annotation_strip_figure(
                    col_info['groups'], title=col_name,
                    color_map=col_info['color_map'], orientation='vertical',
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_annot, plots_dir,
                    f'bulk_heatmap_annotation_{col_name}.png',
                    'annotation', f'{col_name} 注释条', formats=('png', 'svg'), dpi=300,
                ))

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
                fig_ga = annotation_strip_figure(
                    gene_values, title=gcol, color_map=g_color_map,
                    orientation='horizontal',
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_ga, plots_dir,
                    f'bulk_heatmap_gene_annot_{gcol}.png',
                    'annotation', f'{gcol} 基因注释条', formats=('png', 'svg'), dpi=300,
                ))

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

        corr_group_labels = (
            adata.obs[groupby].astype(str).tolist()
            if groupby and groupby in adata.obs.columns else None
        )
        fig_corr, corr_metadata = correlation_heatmap_figure(
            corr_matrix, all_sample_labels,
            title=f'样本相关性热图 ({corr_method.capitalize()})',
            method=corr_method.capitalize(), group_labels=corr_group_labels,
            colorscale=corr_colorscale, cluster=True, mask_diagonal=True,
        )
        result_files.extend(self.save_matplotlib_figure(
            fig_corr, plots_dir, 'bulk_corr_heatmap.png', 'heatmap',
            '样本相关性热图', formats=('png', 'svg'), dpi=300,
        ))
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        corr_pairs = correlation_pairwise_table(
            corr_matrix, sample_labels, corr_group_labels, method=corr_method)
        corr_pairs_csv = os.path.join(results_dir, 'bulk_correlation_pairs.csv')
        corr_pairs.to_csv(corr_pairs_csv, index=False)
        corr_summary = summarize_correlation_pairs(corr_pairs, method=corr_method)
        corr_summary_csv = os.path.join(results_dir, 'bulk_correlation_summary.csv')
        corr_summary.to_csv(corr_summary_csv, index=False)
        result_files.extend([
            {'file_path': corr_pairs_csv, 'file_type': 'csv', 'category': 'table',
             'label': f'样本两两 {corr_method.capitalize()} 相关性'},
            {'file_path': corr_summary_csv, 'file_type': 'csv', 'category': 'table',
             'label': '样本相关性组内/组间摘要'},
        ])

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
                'n_samples': len(sample_labels),
                'n_samples_total': adata.n_obs,
                'sample_display': display_contract,
                'groupby': groupby or '无',
                'correlation_off_diagonal': corr_metadata['off_diagonal'],
                'correlation_display_range': [
                    corr_metadata['display_vmin'], corr_metadata['display_vmax'],
                ],
            }
        }
