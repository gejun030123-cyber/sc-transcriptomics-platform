import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class BulkPCAAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_pca"
    DISPLAY_NAME = "Bulk PCA / UMAP 降维"
    DESCRIPTION = "对 Bulk RNA-seq 数据进行 PCA 和 UMAP 降维可视化"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import plotly.graph_objects as go

        self.progress(5, "加载数据...")
        from modules.io_utils import read_expression_matrix
        adata = read_expression_matrix(input_path)

        n_comps = int(self.params.get('n_comps', 10))
        color_by = self.params.get('color_by', '')
        run_umap = self.params.get('run_umap', True)

        self.progress(20, "标准化数据...")
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value=10)

        self.progress(40, "运行 PCA...")
        actual_comps = min(n_comps, adata.n_obs - 1, adata.n_vars - 1)
        if actual_comps < 2:
            raise ValueError(f"样本数不足 ({adata.n_obs})，至少需要 3 个样本才能进行 PCA 分析。")
        sc.pp.pca(adata, n_comps=actual_comps)
        pca_variance = adata.uns['pca']['variance_ratio']

        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        color_values = None
        if color_by and color_by in adata.obs.columns:
            color_values = adata.obs[color_by].tolist()
        elif color_by == '' or color_by is None:
            color_values = None

        self.progress(55, "生成 PCA 图...")
        pc = adata.obsm['X_pca']
        hover = adata.obs.index.tolist()

        fig_pca = go.Figure()
        if color_values is not None and not all(str(v) == str(color_values[0]) for v in color_values):
            unique_vals = sorted(set(str(v) for v in color_values))
            palette = ['#1a237e', '#e53935', '#4caf50', '#ff9800', '#9c27b0',
                        '#00bcd4', '#795548', '#607d8b', '#f44336', '#3f51b5']
            for i, val in enumerate(unique_vals):
                idx = [j for j, v in enumerate(color_values) if str(v) == val]
                fig_pca.add_trace(go.Scattergl(
                    x=pc[idx, 0], y=pc[idx, 1], mode='markers+text',
                    text=[adata.obs.index[j] for j in idx],
                    textposition='top center', textfont=dict(size=8),
                    marker=dict(size=8, color=palette[i % len(palette)]),
                    name=str(val)
                ))
        else:
            fig_pca.add_trace(go.Scattergl(
                x=pc[:, 0], y=pc[:, 1], mode='markers+text',
                text=hover, textposition='top center', textfont=dict(size=8),
                marker=dict(size=8, color='#1a237e'), name='样本'
            ))
        fig_pca.update_layout(
            title=f'PCA 分析 (n={adata.n_obs})',
            xaxis_title=f'PC1 ({pca_variance[0]*100:.1f}% variance)',
            yaxis_title=f'PC2 ({pca_variance[1]*100:.1f}% variance)',
            plot_bgcolor='white', width=700, height=550
        )
        fpath = os.path.join(plots_dir, 'bulk_pca.json')
        with open(fpath, 'w') as f: f.write(fig_pca.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': 'PCA 分析'})

        self.progress(65, "生成方差解释图...")
        fig_var = go.Figure()
        fig_var.add_trace(go.Bar(
            x=[f'PC{i+1}' for i in range(actual_comps)],
            y=pca_variance * 100, marker_color='#1a237e'
        ))
        fig_var.update_layout(
            title='PCA 方差解释比例', xaxis_title='主成分', yaxis_title='方差解释比例 (%)',
            plot_bgcolor='white', width=600, height=350
        )
        fpath = os.path.join(plots_dir, 'bulk_pca_variance.json')
        with open(fpath, 'w') as f: f.write(fig_var.to_json(engine="json"))
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'pca', 'label': '方差解释比例'})

        if run_umap and adata.n_obs >= 10:
            self.progress(75, "运行 UMAP...")
            sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1))
            sc.tl.umap(adata)
            umap_coords = adata.obsm['X_umap']

            fig_umap = go.Figure()
            if color_values is not None and not all(str(v) == str(color_values[0]) for v in color_values):
                unique_vals = sorted(set(str(v) for v in color_values))
                palette = ['#1a237e', '#e53935', '#4caf50', '#ff9800', '#9c27b0',
                            '#00bcd4', '#795548', '#607d8b', '#f44336', '#3f51b5']
                for i, val in enumerate(unique_vals):
                    idx = [j for j, v in enumerate(color_values) if str(v) == val]
                    fig_umap.add_trace(go.Scattergl(
                        x=umap_coords[idx, 0], y=umap_coords[idx, 1], mode='markers+text',
                        text=[adata.obs.index[j] for j in idx],
                        textposition='top center', textfont=dict(size=8),
                        marker=dict(size=8, color=palette[i % len(palette)]),
                        name=str(val)
                    ))
            else:
                fig_umap.add_trace(go.Scattergl(
                    x=umap_coords[:, 0], y=umap_coords[:, 1], mode='markers+text',
                    text=hover, textposition='top center', textfont=dict(size=8),
                    marker=dict(size=8, color='#1a237e'), name='样本'
                ))
            fig_umap.update_layout(
                title='UMAP 分析', xaxis_title='UMAP1', yaxis_title='UMAP2',
                plot_bgcolor='white', width=700, height=550
            )
            fpath = os.path.join(plots_dir, 'bulk_umap.json')
            with open(fpath, 'w') as f: f.write(fig_umap.to_json(engine="json"))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': 'UMAP 分析'})

        self.progress(90, "保存结果...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_pca_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_samples': adata.n_obs,
                'n_genes': adata.n_vars,
                'n_components': actual_comps,
                'pc1_variance_pct': round(float(pca_variance[0] * 100), 2),
                'pc2_variance_pct': round(float(pca_variance[1] * 100), 2),
                'umap_computed': bool(run_umap and adata.n_obs >= 10),
            }
        }
