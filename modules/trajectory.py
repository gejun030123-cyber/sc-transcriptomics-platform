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
        scope_key = str(self.params.get('scope_key', '') or '').strip()
        if 'neighbors' not in adata.uns and not scope_key:
            return "邻居图未找到，请先运行降维分析（dimred）"
        if 'X_umap' not in adata.obsm:
            return "缺少 UMAP 嵌入（obsm['X_umap']），请先运行降维分析（dimred）"
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
        adata = self.apply_scope(adata)
        # A scope creates a new cell graph.  Reusing the full-dataset
        # ``connectivities`` after slicing changes the transition kernel and
        # can even leave a graph whose shape no longer matches the selected
        # cells.  Rebuild neighbours from the upstream representation instead.
        scope_key = str(self.params.get('scope_key', '') or '').strip()
        graph_recomputed = False
        representation = None
        if scope_key:
            if adata.n_obs < 3:
                raise ValueError('轨迹分析的范围子集至少需要 3 个细胞才能重建邻居图。')
            for candidate in (
                'X_pca', 'X_harmony', 'X_scanorama', 'X_scvi', 'X_scVI',
                'X_sysvi', 'X_combat', 'X_mnn',
            ):
                if candidate in adata.obsm:
                    representation = candidate
                    break
            if representation is None:
                raise ValueError(
                    "轨迹分析指定了分析范围，但子集没有可用于重建邻居图的 PCA/校正表示；"
                    "请先运行 dimred（或批次校正）后重试。"
                )
            n_neighbors = int(self.params.get('n_neighbors', 15) or 15)
            n_neighbors = max(2, min(n_neighbors, adata.n_obs - 1))
            neighbor_kwargs = {
                'n_neighbors': n_neighbors,
                'use_rep': representation,
            }
            if representation == 'X_pca':
                n_available_pcs = int(adata.obsm[representation].shape[1])
                requested_pcs = int(self.params.get('n_pcs', n_available_pcs) or n_available_pcs)
                neighbor_kwargs['n_pcs'] = max(1, min(requested_pcs, n_available_pcs))
            self.progress(12, f"重建范围子集邻居图（{representation}）...")
            sc.pp.neighbors(adata, **neighbor_kwargs)
            graph_recomputed = True
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
        root_component_used = None
        if start_cluster and cluster_key in adata.obs.columns:
            root_mask = adata.obs[cluster_key].astype(str) == start_cluster
            if root_mask.sum() > 0:
                adata.uns['iroot'] = np.where(root_mask)[0][0]
            else:
                self.progress(-1, f"起始簇 '{start_cluster}' 在 '{cluster_key}' 中不存在，"
                                  "已回退为自动选择根细胞（DC1 最小）。")
        if 'iroot' not in adata.uns:
            if 'X_diffmap' in adata.obsm and adata.obsm['X_diffmap'].shape[1] > 0:
                # Scanpy reserves DC0 for the stationary component; DC1 is the
                # first informative diffusion coordinate for root selection.
                root_component = 1 if adata.obsm['X_diffmap'].shape[1] > 1 else 0
                adata.uns['iroot'] = int(np.argmin(adata.obsm['X_diffmap'][:, root_component]))
                root_component_used = root_component
            else:
                adata.uns['iroot'] = 0

        # n_dcs 超过扩散成分数时 scanpy dpt 会直接抛 ValueError；提前钳制并提示。
        n_dcs_used = int(n_dcs)
        if 'X_diffmap' in adata.obsm:
            n_available_dcs = int(adata.obsm['X_diffmap'].shape[1])
            if n_dcs_used > n_available_dcs:
                n_dcs_used = n_available_dcs
                self.progress(-1, f"n_dcs={n_dcs} 超过扩散成分数 {n_available_dcs}，"
                                  f"已钳制为 {n_dcs_used}。")

        self.progress(40, "Computing diffusion pseudotime...")
        sc.tl.dpt(adata, n_branchings=n_branchings, n_dcs=n_dcs_used)

        # PAGA
        if enable_paga and cluster_key in adata.obs.columns:
            self.progress(50, "Computing PAGA...")
            sc.tl.paga(adata, groups=cluster_key)

        self.progress(60, "Generating trajectory plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        color_key = cluster_key if cluster_key in adata.obs.columns else 'dpt_pseudotime'
        try:
            fig_umap = umap_figure(adata, color_key, title=f'Trajectory by {color_key}')
            result_files.extend(self.save_matplotlib_figure(
                fig_umap, plots_dir, 'trajectory_pseudotime.png', 'umap',
                f'Trajectory by {color_key}', formats=('png', 'svg'), dpi=300,
            ))
        except KeyError as exc:
            self.progress(-1, f'UMAP 轨迹图失败（缺少嵌入坐标）: {exc}')

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
                # 不连通图中 pseudotime 可能为 inf，先过滤再画，避免 inf 污染
                # 排序与横轴。
                finite_mask = np.isfinite(pt)
                if finite_mask.sum() < 2:
                    self.progress(-1, '有效（非 inf）pseudotime 细胞不足，跳过基因趋势图。')
                    gene_list = []
                else:
                    sort_idx = np.argsort(pt[finite_mask])
                    pt_sorted = pt[finite_mask][sort_idx]
                    # 平滑窗口钳制为不超过有效细胞数，并取奇数避免中心偏移。
                    window = max(adata.n_obs // 50, 10)
                    window = min(window, int(finite_mask.sum()))
                    if window % 2 == 0:
                        window += 1
                    window = min(window, int(finite_mask.sum()))

                    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=150)
                    counts_layer = (
                        adata.layers.get('counts')
                        if hasattr(adata, 'layers') and 'counts' in adata.layers else None
                    )
                    use_counts_scale = counts_layer is not None
                    expression_label = 'Expression (log1p counts)' if use_counts_scale else 'Expression'
                    if use_counts_scale:
                        library_size = np.asarray(
                            counts_layer.sum(axis=1), dtype=float,
                        ).reshape(-1)
                        counts_scale = np.divide(
                            10000.0, library_size,
                            out=np.zeros_like(library_size, dtype=float),
                            where=library_size > 0,
                        )
                    for gene in gene_list:
                        if use_counts_scale and gene in adata.var_names:
                            position = adata.var_names.get_loc(gene)
                            column = counts_layer[:, position]
                            raw = (
                                column.toarray().ravel()
                                if hasattr(column, 'toarray')
                                else np.asarray(column).ravel()
                            )
                            expr = np.log1p(raw * counts_scale)
                        else:
                            expr = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()
                        expr_sorted = expr[finite_mask][sort_idx]
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
        max_pseudotime = None
        if 'dpt_pseudotime' in adata.obs.columns:
            finite_pt = adata.obs['dpt_pseudotime'].astype(float)
            finite_pt = finite_pt[np.isfinite(finite_pt)]
            if len(finite_pt) > 0:
                max_pseudotime = round(float(finite_pt.max()), 3)
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'cluster_column': cluster_key,
                'requested_cluster_key': requested_cluster_key,
                'max_pseudotime': max_pseudotime,
                'n_cells': adata.n_obs,
                'scope_key': str(self.params.get('scope_key', '') or '').strip() or None,
                'scope_values': ([v.strip() for v in str(self.params.get('scope_values', '') or '').split(',') if v.strip()] or None),
                'neighbor_graph_recomputed_for_scope': bool(graph_recomputed),
                'neighbor_representation': representation,
                'root_diffusion_component': root_component_used,
                'n_dcs_used': int(n_dcs_used),
            }
        }
