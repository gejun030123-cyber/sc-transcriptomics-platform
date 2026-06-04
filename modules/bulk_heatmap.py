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
        if input_path.endswith('.csv') or input_path.endswith('.txt'):
            df = pd.read_csv(input_path, sep=None if input_path.endswith('.csv') else '\t', index_col=0)
            adata = sc.AnnData(X=df.values.T, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))
        else:
            adata = sc.read_h5ad(input_path)

        top_n = int(self.params.get('top_n', 50))
        groupby = self.params.get('groupby', '')
        hm_type = self.params.get('heatmap_type', 'top_var')

        self.progress(20, "计算数据矩阵...")
        counts = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
        counts = counts.astype(float)

        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        norm_data = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
        norm_data = norm_data.astype(float)

        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        self.progress(40, "选择基因...")

        if hm_type == 'top_var':
            gene_var = np.var(norm_data, axis=0)
            top_idx = np.argsort(gene_var)[::-1][:top_n]
            title = f'Top {top_n} 高变异基因热图'
        elif hm_type == 'deg' and os.path.exists(os.path.join(self.project_dir, 'results', 'bulk_deg_results.csv')):
            deg_df = pd.read_csv(os.path.join(self.project_dir, 'results', 'bulk_deg_results.csv'))
            deg_sig = deg_df[deg_df['regulation'] != 'NS'].head(top_n)
            gene_list = deg_sig['gene'].tolist()
            top_idx = [i for i, g in enumerate(adata.var_names) if g in gene_list][:top_n]
            title = f'Top {top_n} 差异基因热图'
        else:
            gene_var = np.var(norm_data, axis=0)
            top_idx = np.argsort(gene_var)[::-1][:top_n]
            title = f'Top {top_n} 高变异基因热图'

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
            hovertemplate='样本: %{y}<br>基因: %{x}<br>Z-score: %{z:.2f}'
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
        with open(fpath, 'w') as f: json.dump(json.loads(fig.to_json()), f)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': title})

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
            hovertemplate='%{y} vs %{x}<br>r = %{z:.3f}'
        ))
        fig_corr.update_layout(
            title='样本相关性热图 (Pearson)',
            height=max(400, adata.n_obs * 30 + 100),
            width=max(400, adata.n_obs * 30 + 100),
            plot_bgcolor='white'
        )
        fpath = os.path.join(plots_dir, 'bulk_corr_heatmap.json')
        with open(fpath, 'w') as f: json.dump(json.loads(fig_corr.to_json()), f)
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
