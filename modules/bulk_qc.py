import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class BulkQCAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_qc"
    DISPLAY_NAME = "Bulk RNA-seq 质控"
    DESCRIPTION = "计数矩阵质控：文库大小、基因检测、离群值过滤"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        from modules.visualization import scatter_plot, bar_plot
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        self.progress(5, "加载计数矩阵...")
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(input_path)

        # 应用自定义过滤规则
        adata = self.apply_filters(adata, 'bulk_qc')

        min_counts = int(self.params.get('min_counts', 100000))
        min_genes = int(self.params.get('min_genes', 5000))
        max_mt_pct = float(self.params.get('max_mt_pct', 20.0))

        self.progress(20, "计算质控指标...")
        # 优先用 gene_name 检测线粒体基因（Ensembl ID 不以 MT- 开头）
        if 'gene_name' in adata.var.columns:
            gene_names_for_mt = adata.var['gene_name'].fillna('').astype(str)
        else:
            gene_names_for_mt = adata.var_names.astype(str)
        adata.var['mt'] = gene_names_for_mt.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        n_before = adata.n_obs
        lib_sizes = adata.obs['total_counts'].values
        n_genes_detected = adata.obs['n_genes_by_counts'].values
        mt_pct = adata.obs['pct_counts_mt'].values if 'pct_counts_mt' in adata.obs.columns else np.zeros(n_before)

        self.progress(40, "过滤样本...")
        mask = (lib_sizes >= min_counts) & (n_genes_detected >= min_genes) & (mt_pct <= max_mt_pct)
        adata_filtered = adata[mask].copy()
        n_after = adata_filtered.n_obs

        self.progress(60, "生成质控图表...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        fig = make_subplots(rows=2, cols=2,
            subplot_titles=['文库大小分布', '检测基因数',
                           '线粒体基因比例', '文库大小 vs 检测基因数'])
        fig.add_trace(go.Bar(x=list(range(n_before)), y=lib_sizes, marker_color='#1a237e', name='文库大小'), row=1, col=1)
        fig.add_trace(go.Bar(x=list(range(n_before)), y=n_genes_detected, marker_color='#283593', name='基因数'), row=1, col=2)
        fig.add_trace(go.Bar(x=list(range(n_before)), y=mt_pct, marker_color='#e53935', name='MT%'), row=2, col=1)
        colors = ['#4caf50' if m else '#e53935' for m in mask]
        fig.add_trace(go.Scattergl(x=lib_sizes, y=n_genes_detected, mode='markers',
            marker=dict(color=colors, size=6), name='样本'), row=2, col=2)
        fig.update_layout(height=600, width=800, showlegend=False, title='Bulk RNA-seq 质控总览')
        fpath = os.path.join(plots_dir, 'bulk_qc_overview.json')
        with open(fpath, 'w') as f: f.write(fig.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': '质控总览'})

        # 样本相关性热图
        norm_for_corr = adata_filtered.copy()
        sc.pp.normalize_total(norm_for_corr, target_sum=1e6)
        sc.pp.log1p(norm_for_corr)
        corr_data = norm_for_corr.X if not hasattr(norm_for_corr.X, 'toarray') else norm_for_corr.X.toarray()
        corr_matrix = np.corrcoef(corr_data)
        sample_labels_corr = norm_for_corr.obs.index.tolist()

        fig_corr = go.Figure()
        fig_corr.add_trace(go.Heatmap(
            z=corr_matrix.tolist(), x=sample_labels_corr, y=sample_labels_corr,
            colorscale='Blues', zmin=0, zmax=1,
            colorbar=dict(title='Pearson r'),
            hovertemplate='%{y} vs %{x}<br>r = %{z:.3f}<extra></extra>'
        ))
        fig_corr.update_layout(
            title='样本相关性热图 (Pearson)',
            height=max(400, n_after * 30 + 100), width=max(400, n_after * 30 + 100),
            plot_bgcolor='white'
        )
        fpath = os.path.join(plots_dir, 'bulk_qc_corr.json')
        with open(fpath, 'w') as f: f.write(fig_corr.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'heatmap', 'label': '样本相关性热图'})

        if n_before > n_after:
            fig_r = go.Figure()
            fig_r.add_trace(go.Bar(x=['过滤前', '过滤后'], y=[n_before, n_after], marker_color=['#e53935', '#4caf50']))
            fig_r.update_layout(title='样本过滤结果', yaxis_title='样本数', width=400, height=300)
            fpath = os.path.join(plots_dir, 'bulk_qc_filter.json')
            with open(fpath, 'w') as f: f.write(fig_r.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'qc', 'label': '样本过滤结果'})

        self.progress(80, "运行 PCA 离群检测...")
        sc.pp.normalize_total(adata_filtered, target_sum=1e6)
        sc.pp.log1p(adata_filtered)
        sc.pp.pca(adata_filtered, n_comps=min(10, n_after - 1))
        fig_pca = go.Figure()
        pc = adata_filtered.obsm['X_pca']
        fig_pca.add_trace(go.Scattergl(x=pc[:, 0], y=pc[:, 1], mode='markers+text',
            text=adata_filtered.obs.index.tolist(), textposition='top center',
            marker=dict(size=8, color='#1a237e')))
        fig_pca.update_layout(title='质控后样本 PCA', xaxis_title='PC1', yaxis_title='PC2',
                             plot_bgcolor='white', width=600, height=500)
        fpath = os.path.join(plots_dir, 'bulk_qc_pca.json')
        with open(fpath, 'w') as f: f.write(fig_pca.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': '样本 PCA'})

        self.progress(90, "保存输出...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_qc_output.h5ad')
        adata_filtered.write_h5ad(output_path)

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'samples_before': n_before,
                'samples_after': n_after,
                'samples_removed': n_before - n_after,
                'genes_total': adata_filtered.n_vars,
                'median_lib_size': int(np.median(adata_filtered.obs['total_counts'])),
                'median_genes': int(np.median(adata_filtered.obs['n_genes_by_counts'])),
            }
        }
