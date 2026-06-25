import os
import json
from modules.base import BaseAnalysis

class TrajectoryAnalysis(BaseAnalysis):
    MODULE_NAME = "trajectory"
    DISPLAY_NAME = "轨迹分析"
    DESCRIPTION = "基于扩散图的拟时序分析"
    INPUT_REQUIRES = ['X_umap']

    def validate_input(self, adata):
        if 'neighbors' not in adata.uns:
            return "邻居图未找到，请先运行降维分析（dimred）"
        return None

    def run(self, input_path):
        import scanpy as sc
        import numpy as np
        from modules.visualization import umap_scatter
        import plotly.graph_objects as go

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        cluster_key = self.params.get('cluster_key', 'leiden')
        enable_paga = self.params.get('enable_paga', False)
        paga_threshold = float(self.params.get('paga_threshold', 0.05))
        n_diffcomps = int(self.params.get('n_diffcomps', 15))
        start_cluster = self.params.get('start_cluster', '').strip()
        n_dcs = int(self.params.get('n_dcs', 10))
        n_branchings = int(self.params.get('n_branchings', 0))

        self.progress(20, "Computing diffusion map...")
        sc.tl.diffmap(adata, n_comps=n_diffcomps)

        # Set root cell for pseudotime
        if start_cluster and cluster_key in adata.obs.columns:
            root_mask = adata.obs[cluster_key].astype(str) == start_cluster
            if root_mask.sum() > 0:
                adata.uns['iroot'] = np.where(root_mask)[0][0]

        self.progress(40, "Computing diffusion pseudotime...")
        sc.tl.dpt(adata, n_branchings=n_branchings, n_dcs=n_dcs)

        # PAGA
        if enable_paga and cluster_key in adata.obs.columns:
            self.progress(50, "Computing PAGA...")
            sc.tl.paga(adata, groups=cluster_key)

        self.progress(60, "Generating trajectory plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        color_key = cluster_key if cluster_key in adata.obs.columns else 'dpt_pseudotime'
        fig_json = json.dumps(umap_scatter(adata, color_key, title=f'Trajectory by {color_key}'))
        fpath = os.path.join(plots_dir, 'trajectory_pseudotime.json')
        with open(fpath, 'w') as f: f.write(fig_json)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': f'Trajectory by {color_key}'})

        # PAGA plot
        if enable_paga and 'paga' in adata.uns:
            try:
                paga_connectivities = adata.uns['paga']['connectivities'].toarray() if hasattr(adata.uns['paga']['connectivities'], 'toarray') else adata.uns['paga']['connectivities']
                fig_paga = go.Figure()
                if 'X_umap' in adata.obsm:
                    cluster_means = {}
                    for cat in adata.obs[cluster_key].unique():
                        mask = adata.obs[cluster_key] == cat
                        cluster_means[cat] = adata.obsm['X_umap'][mask].mean(axis=0)
                    for i, ci in enumerate(adata.obs[cluster_key].cat.categories if hasattr(adata.obs[cluster_key], 'cat') else sorted(adata.obs[cluster_key].unique())):
                        for j, cj in enumerate(adata.obs[cluster_key].cat.categories if hasattr(adata.obs[cluster_key], 'cat') else sorted(adata.obs[cluster_key].unique())):
                            if i < j and paga_connectivities[i, j] > paga_threshold:
                                xi, yi = cluster_means.get(ci, [0, 0])
                                xj, yj = cluster_means.get(cj, [0, 0])
                                fig_paga.add_trace(go.Scattergl(
                                    x=[xi, xj], y=[yi, yj], mode='lines',
                                    line=dict(width=paga_connectivities[i, j] * 5, color='gray'),
                                    showlegend=False
                                ))
                fig_paga.update_layout(title='PAGA Trajectory', xaxis_title='UMAP1', yaxis_title='UMAP2',
                                      plot_bgcolor='white', width=600, height=500)
                result_files.append(self.save_plotly_json(fig_paga, plots_dir, 'trajectory_paga.json', 'paga', 'PAGA Trajectory'))
            except Exception:
                pass

        if 'X_diffmap' in adata.obsm:
            dc = adata.obsm['X_diffmap'][:, :2]
            fig = go.Figure()
            color_vals = adata.obs['dpt_pseudotime'].values if 'dpt_pseudotime' in adata.obs.columns else None
            fig.add_trace(go.Scattergl(x=dc[:, 0], y=dc[:, 1], mode='markers',
                                       marker=dict(size=2, color=color_vals, colorscale='Viridis', colorbar=dict(title='Pseudotime'))))
            fig.update_layout(title='Diffusion Map', xaxis_title='DC1', yaxis_title='DC2',
                             plot_bgcolor='white', width=600, height=500)
            result_files.append(self.save_plotly_json(fig, plots_dir, 'trajectory_diffmap.json', 'diffusion_map', 'Diffusion Map'))

        # Gene expression along pseudotime
        plot_genes_str = self.params.get('plot_genes', '').strip()
        if plot_genes_str and 'dpt_pseudotime' in adata.obs.columns:
            import numpy as np
            gene_list = [g.strip() for g in plot_genes_str.replace('\n', ',').split(',') if g.strip()]
            gene_list = [g for g in gene_list if g in adata.var_names][:10]

            if gene_list:
                self.progress(75, f"Plotting {len(gene_list)} genes along pseudotime...")
                pt = adata.obs['dpt_pseudotime'].values
                sort_idx = np.argsort(pt)
                pt_sorted = pt[sort_idx]
                window = max(adata.n_obs // 50, 10)

                fig = go.Figure()
                for gene in gene_list:
                    expr = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()
                    expr_sorted = expr[sort_idx]
                    # Rolling mean
                    kernel = np.ones(window) / window
                    smoothed = np.convolve(expr_sorted, kernel, mode='valid')
                    pt_smooth = pt_sorted[window // 2: window // 2 + len(smoothed)]
                    fig.add_trace(go.Scattergl(x=pt_smooth, y=smoothed, mode='lines', name=gene))

                fig.update_layout(
                    title='Gene Expression Along Pseudotime',
                    xaxis_title='Pseudotime', yaxis_title='Expression',
                    plot_bgcolor='white', width=700, height=400
                )
                result_files.append(self.save_plotly_json(fig, plots_dir, 'trajectory_gene_expression.json', 'gene_expression', '基因拟时序表达'))

        self.progress(85, "Saving output...")
        output_path = self.save_output(adata, 'trajectory')

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'cluster_column': cluster_key,
                'max_pseudotime': round(float(adata.obs['dpt_pseudotime'].max()), 3) if 'dpt_pseudotime' in adata.obs.columns else None,
                'n_cells': adata.n_obs,
            }
        }
