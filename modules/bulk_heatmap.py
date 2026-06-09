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
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import pdist

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
            var_names_list = list(adata.var_names)
            top_idx = [var_names_list.index(g) for g in gene_list if g in var_names_list]
            not_found = [g for g in gene_list if g not in var_names_list]
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

        heat_data = norm_data[:, top_idx]
        gene_labels = [adata.var_names[i] for i in top_idx]
        sample_labels = adata.obs.index.tolist()

        z_mean = heat_data.mean(axis=0)
        z_std = heat_data.std(axis=0) + 1e-10
        heat_z = (heat_data - z_mean) / z_std
        heat_z = np.clip(heat_z, -3, 3)

        self.progress(60, "聚类样本...")
        if adata.n_obs > 2:
            sample_dist = pdist(heat_z, metric='euclidean')
            sample_link = linkage(sample_dist, method='ward')
            dendro = dendrogram(sample_link, no_plot=True)
            sample_order = dendro['leaves']
        else:
            sample_order = list(range(adata.n_obs))

        if len(top_idx) > 2:
            gene_dist = pdist(heat_z.T, metric='euclidean')
            gene_link = linkage(gene_dist, method='ward')
            gene_dendro = dendrogram(gene_link, no_plot=True)
            gene_order = gene_dendro['leaves']
        else:
            gene_order = list(range(len(top_idx)))

        heat_ordered = heat_z[np.ix_(sample_order, gene_order)]
        sample_ordered = [sample_labels[i] for i in sample_order]
        gene_ordered = [gene_labels[i] for i in gene_order]

        self.progress(75, "生成热图...")
        annotation_colors = None
        colorbar_trace = None
        if groupby and groupby in adata.obs.columns:
            groups = [adata.obs[groupby].iloc[i] for i in sample_order]
            unique_groups = sorted(set(str(g) for g in groups))
            palette = ['#1a237e', '#e53935', '#4caf50', '#ff9800', '#9c27b0',
                        '#00bcd4', '#795548', '#607d8b', '#f44336', '#3f51b5']
            group_color_map = {g: palette[i % len(palette)] for i, g in enumerate(unique_groups)}
            annotation_colors = [group_color_map[str(g)] for g in groups]

        fig = go.Figure()
        fig.add_trace(go.Heatmap(
            z=heat_ordered.tolist(),
            x=gene_ordered,
            y=sample_ordered,
            colorscale='RdBu_r',
            zmid=0,
            colorbar=dict(title='Z-score'),
            hovertemplate='样本: %{y}<br>基因: %{x}<br>Z-score: %{z:.2f}<extra></extra>'
        ))

        annotations = []
        if annotation_colors:
            for i, color in enumerate(annotation_colors):
                annotations.append(dict(
                    x=-0.5, y=i, xref='x', yref='y',
                    text='', showarrow=False,
                    xanchor='right'
                ))

        fig.update_layout(
            title=title,
            xaxis=dict(tickangle=45, tickfont=dict(size=8)),
            yaxis=dict(tickfont=dict(size=9)),
            height=max(400, adata.n_obs * 25 + 150),
            width=max(600, len(top_idx) * 12 + 200),
            plot_bgcolor='white'
        )

        fpath = os.path.join(plots_dir, 'bulk_heatmap.json')
        with open(fpath, 'w') as f: f.write(fig.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': title})

        # 分组注释条
        if annotation_colors:
            color_indices = [unique_groups.index(str(groups[i])) for i in range(len(groups))]
            fig_annot = go.Figure()
            fig_annot.add_trace(go.Heatmap(
                z=[[i] for i in color_indices],
                y=sample_ordered,
                x=['Group'],
                colorscale=[[i / max(len(unique_groups) - 1, 1), palette[i % len(palette)]] for i in range(len(unique_groups))],
                showscale=False,
                text=[[unique_groups[i]] for i in color_indices],
                hovertemplate='%{y}: %{text}<extra></extra>'
            ))
            fig_annot.update_layout(
                height=max(400, adata.n_obs * 25 + 150), width=100,
                margin=dict(l=0, r=0, t=30, b=40)
            )
            fpath_annot = os.path.join(plots_dir, 'bulk_heatmap_annotation.json')
            with open(fpath_annot, 'w') as f: f.write(fig_annot.to_json(engine="json"))
            result_files.append({'file_path': fpath_annot, 'file_type': 'plotly_json', 'category': 'annotation', 'label': '分组注释条'})

        self.progress(85, "生成样本相关性热图...")
        corr_matrix = np.corrcoef(norm_data)
        fig_corr = go.Figure()
        fig_corr.add_trace(go.Heatmap(
            z=corr_matrix.tolist(),
            x=sample_labels,
            y=sample_labels,
            colorscale='Blues',
            zmin=0, zmax=1,
            colorbar=dict(title='Pearson r'),
            hovertemplate='%{y} vs %{x}<br>r = %{z:.3f}<extra></extra>'
        ))
        fig_corr.update_layout(
            title='样本相关性热图 (Pearson)',
            height=max(400, adata.n_obs * 30 + 100),
            width=max(400, adata.n_obs * 30 + 100),
            plot_bgcolor='white'
        )
        fpath = os.path.join(plots_dir, 'bulk_corr_heatmap.json')
        with open(fpath, 'w') as f: f.write(fig_corr.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': '样本相关性热图'})

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
