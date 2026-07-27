import os
import logging
from modules.base import BaseAnalysis
from modules.io_utils import resolve_obs_grouping

logger = logging.getLogger(__name__)

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
        import pandas as pd
        import matplotlib.pyplot as plt
        from modules.native_figures import umap_figure, grouped_box_figure, histogram_figure, scatter_figure, line_figure
        from modules.figure_style import NATURE_PALETTE, NATURE_TEXT, NATURE_GRID

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        requested_cluster_key = str(self.params.get('cluster_key', 'leiden') or '').strip()
        cluster_key, cluster_info = resolve_obs_grouping(
            adata, requested_cluster_key,
            fallbacks=['leiden', 'leiden_0.8', 'leiden_0.6', 'leiden_1.0'],
            max_categories=50, max_numeric_categories=20,
            require_multiple=False,
        )
        if cluster_key is None:
            cluster_key = ''
        else:
            if cluster_key != requested_cluster_key:
                self.progress(-1, f"轨迹分组列已改用 '{cluster_key}'：{cluster_info.get('requested_reason', '')}")
            if not isinstance(adata.obs[cluster_key].dtype, pd.CategoricalDtype):
                adata.obs[cluster_key] = adata.obs[cluster_key].astype(str).astype('category')
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
        if 'iroot' not in adata.uns:
            if 'X_diffmap' in adata.obsm and adata.obsm['X_diffmap'].shape[1] > 0:
                adata.uns['iroot'] = int(np.argmin(adata.obsm['X_diffmap'][:, 0]))
            else:
                adata.uns['iroot'] = 0

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
        fig_umap = umap_figure(adata, color_key, title=f'Trajectory by {color_key}')
        result_files.extend(self.save_matplotlib_figure(
            fig_umap, plots_dir, 'trajectory_pseudotime.png', 'umap',
            f'Trajectory by {color_key}', formats=('png', 'svg'), dpi=300,
        ))

        if self.params.get('show_pseudotime_distribution', True) and 'dpt_pseudotime' in adata.obs.columns:
            pseudotime = adata.obs['dpt_pseudotime'].astype(float)
            if cluster_key in adata.obs.columns:
                labels = adata.obs[cluster_key].astype(str)
                fig_dist = grouped_box_figure(
                    pseudotime.values[:, None], labels, ["DPT pseudotime"],
                    title=f'Pseudotime Distribution by {cluster_key}',
                    x_label=cluster_key, y_label='DPT pseudotime', rotation=35,
                )
            else:
                fig_dist = histogram_figure(
                    pseudotime.values, labels=['Cells'], bins=60,
                    title='Pseudotime Distribution', x_label='DPT pseudotime',
                    colors=[NATURE_PALETTE[0]],
                )
            result_files.extend(self.save_matplotlib_figure(
                fig_dist, plots_dir, 'trajectory_pseudotime_distribution.png',
                'violin', 'Pseudotime Distribution', formats=('png', 'svg'), dpi=300,
            ))

        # PAGA plot
        if enable_paga and 'paga' in adata.uns:
            try:
                paga_connectivities = adata.uns['paga']['connectivities'].toarray() if hasattr(adata.uns['paga']['connectivities'], 'toarray') else adata.uns['paga']['connectivities']
                fig_paga, ax_paga = plt.subplots(figsize=(7.0, 5.5), dpi=150)
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
                                ax_paga.plot([xi, xj], [yi, yj], color='#98A2B3',
                                             linewidth=max(0.7, paga_connectivities[i, j] * 5),
                                             alpha=0.75, zorder=1)
                    for cat, (x_value, y_value) in cluster_means.items():
                        ax_paga.scatter(x_value, y_value, s=70, color=NATURE_PALETTE[0],
                                        edgecolor='white', linewidth=0.8, zorder=2)
                        ax_paga.text(x_value, y_value, str(cat), ha='center', va='center',
                                     fontsize=8, color='white', zorder=3)
                ax_paga.set_title('PAGA Trajectory', loc='left', fontsize=10,
                                  fontweight='semibold', color=NATURE_TEXT)
                ax_paga.set_xlabel('UMAP1', fontsize=9)
                ax_paga.set_ylabel('UMAP2', fontsize=9)
                ax_paga.grid(False)
                result_files.extend(self.save_matplotlib_figure(
                    fig_paga, plots_dir, 'trajectory_paga.png', 'paga',
                    'PAGA Trajectory', formats=('png', 'svg'), dpi=300,
                ))
            except Exception as e:
                logger.warning("生成 PAGA 轨迹图失败: %s", e)

        if 'X_diffmap' in adata.obsm:
            dc = adata.obsm['X_diffmap'][:, :2]
            color_vals = adata.obs['dpt_pseudotime'].values if 'dpt_pseudotime' in adata.obs.columns else None
            fig = scatter_figure(
                dc[:, 0], dc[:, 1], title='Diffusion Map', x_label='DC1', y_label='DC2',
                colors=color_vals if color_vals is not None else NATURE_PALETTE[0],
                size=8, alpha=0.78,
            )
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, 'trajectory_diffmap.png', 'diffusion_map',
                'Diffusion Map', formats=('png', 'svg'), dpi=300,
            ))

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

                fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=150)
                for gene in gene_list:
                    expr = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()
                    expr_sorted = expr[sort_idx]
                    # Rolling mean
                    kernel = np.ones(window) / window
                    smoothed = np.convolve(expr_sorted, kernel, mode='valid')
                    pt_smooth = pt_sorted[window // 2: window // 2 + len(smoothed)]
                    ax.plot(pt_smooth, smoothed, linewidth=1.8,
                            color=NATURE_PALETTE[gene_list.index(gene) % len(NATURE_PALETTE)],
                            label=gene)
                ax.set_title('Gene Expression Along Pseudotime', loc='left', fontsize=10,
                             fontweight='semibold', color=NATURE_TEXT)
                ax.set_xlabel('Pseudotime', fontsize=9)
                ax.set_ylabel('Expression', fontsize=9)
                ax.grid(axis='y', color=NATURE_GRID, linewidth=0.5, alpha=0.7)
                ax.legend(frameon=False, fontsize=8)
                result_files.extend(self.save_matplotlib_figure(
                    fig, plots_dir, 'trajectory_gene_expression.png', 'gene_expression',
                    '基因拟时序表达', formats=('png', 'svg'), dpi=300,
                ))

        self.progress(85, "Saving output...")
        output_path = self.save_output(adata, 'trajectory')

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'cluster_column': cluster_key,
                'requested_cluster_key': requested_cluster_key,
                'max_pseudotime': round(float(adata.obs['dpt_pseudotime'].max()), 3) if 'dpt_pseudotime' in adata.obs.columns else None,
                'n_cells': adata.n_obs,
            }
        }
