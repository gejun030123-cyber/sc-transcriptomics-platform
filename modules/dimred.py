from modules.base import BaseAnalysis

class DimredAnalysis(BaseAnalysis):
    MODULE_NAME = "dimred"
    DISPLAY_NAME = "降维分析"
    DESCRIPTION = "PCA、UMAP/t-SNE 降维，支持自动选 PC"
    INPUT_REQUIRES = []

    def run(self, input_path):
        import os
        import scanpy as sc
        import omicverse as ov
        import numpy as np
        import matplotlib.pyplot as plt
        from modules.native_figures import scatter_figure, bar_figure, line_figure
        from modules.figure_style import NATURE_PALETTE, nature_continuous_cmap, NATURE_TEXT

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        n_comps = int(self.params.get('n_comps', 50))
        use_mde = self.params.get('use_mde', False)
        auto_n_comps = self.params.get('auto_n_comps', 'none')

        # UMAP parameters
        umap_n_neighbors = int(self.params.get('umap_n_neighbors', 15))
        umap_min_dist = float(self.params.get('umap_min_dist', 0.5))
        umap_metric = self.params.get('umap_metric', 'euclidean')
        umap_spread = float(self.params.get('umap_spread', 1.0))

        # t-SNE parameters
        enable_tsne = self.params.get('enable_tsne', False)
        tsne_perplexity = float(self.params.get('tsne_perplexity', 30))
        tsne_learning_rate = float(self.params.get('tsne_learning_rate', 1000))

        self.progress(20, "Scaling data...")
        ov.pp.scale(adata, max_value=10)

        self.progress(35, f"Running PCA ({n_comps} components)...")
        try:
            sc.pp.pca(adata, n_comps=n_comps, layer='scaled')
        except TypeError:
            if 'scaled' in adata.layers:
                adata.X = adata.layers['scaled']
            sc.pp.pca(adata, n_comps=n_comps)

        # Auto-select number of PCs
        if auto_n_comps != 'none' and 'pca' in adata.uns:
            variance_ratio = adata.uns['pca']['variance_ratio']
            if auto_n_comps == 'elbow':
                # Elbow method: find point of maximum curvature
                cumvar = np.cumsum(variance_ratio)
                diffs = np.diff(cumvar)
                diffs2 = np.diff(diffs)
                n_comps_auto = int(np.argmax(diffs2) + 2)  # +2 for diff offsets
                n_comps_auto = max(5, min(n_comps_auto, n_comps))
                self.progress(40, f"Auto-selected {n_comps_auto} PCs (elbow)")
                n_comps = n_comps_auto
            elif auto_n_comps == 'kneedle':
                try:
                    from kneed import KneeLocator
                    x = range(1, len(variance_ratio) + 1)
                    kl = KneeLocator(x, variance_ratio, curve='convex', direction='decreasing')
                    if kl.knee:
                        n_comps_auto = max(5, min(int(kl.knee), n_comps))
                        self.progress(40, f"Auto-selected {n_comps_auto} PCs (kneedle)")
                        n_comps = n_comps_auto
                except ImportError:
                    self.progress(40, "kneed not installed, using default n_comps")

        self.progress(50, "Computing neighbors...")
        sc.pp.neighbors(adata, n_pcs=n_comps, n_neighbors=umap_n_neighbors, metric=umap_metric)

        self.progress(65, "Computing UMAP...")
        embedding_method = 'umap'
        if use_mde:
            try:
                import pymde
                mde = pymde.preserve_neighbors(adata.obsm['X_pca'], embedding_dim=2, device='cpu')
                embedding = mde.embed(verbose=False)
                adata.obsm['X_mde'] = embedding.cpu().numpy() if hasattr(embedding, 'cpu') else embedding.numpy()
                adata.obsm['X_umap'] = adata.obsm['X_mde']
                embedding_method = 'mde'
            except ImportError:
                self.progress(66, "pymde not installed, falling back to UMAP...")
                sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread)
        else:
            sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread)

        # Optional t-SNE
        if enable_tsne:
            self.progress(72, "Computing t-SNE...")
            sc.tl.tsne(adata, perplexity=tsne_perplexity, learning_rate=int(tsne_learning_rate), n_pcs=n_comps)

        self.progress(80, "Generating embedding plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        for color_key in ['batch', 'leiden', 'n_genes_by_counts']:
            if color_key in adata.obs.columns:
                try:
                    fig_static = self.build_publication_umap(
                        adata, color_key, title=f'UMAP colored by {color_key}'
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_static, plots_dir, f'dimred_umap_{color_key}.png',
                        'umap', f'UMAP by {color_key}'
                    ))
                    import matplotlib.pyplot as plt
                    plt.close(fig_static)
                except Exception as exc:
                    self.progress(-1, f'UMAP 静态导出失败：{exc}')

        # PCA scatter for early detection of outliers and batch/sample structure
        if self.params.get('show_pca_scatter', True) and 'X_pca' in adata.obsm:
            pca = adata.obsm['X_pca']
            color_key = next((k for k in ['batch', 'phase', 'n_genes_by_counts', 'total_counts'] if k in adata.obs.columns), None)
            color_values = None
            if color_key:
                values = adata.obs[color_key]
                try:
                    numeric_values = values.astype(float).values
                    color_values = numeric_values
                except (TypeError, ValueError):
                    labels = values.astype(str)
                    categories = sorted(labels.unique(), key=lambda x: (len(x), x))
                    color_values = np.asarray([
                        NATURE_PALETTE[categories.index(value) % len(NATURE_PALETTE)]
                        for value in labels
                    ], dtype=object)
            else:
                color_values = NATURE_PALETTE[0]
            fig_pca = scatter_figure(
                pca[:, 0], pca[:, 1], title='PCA Scatter', x_label='PC1', y_label='PC2',
                colors=color_values, size=10, alpha=0.75,
            )
            if color_key and np.issubdtype(np.asarray(color_values).dtype, np.number):
                artist = fig_pca.axes[0].collections[0]
                artist.set_cmap(nature_continuous_cmap())
                fig_pca.colorbar(artist, ax=fig_pca.axes[0], fraction=0.035,
                                 pad=0.025, label=color_key)
            elif color_key and isinstance(color_values, np.ndarray) and color_values.dtype == object:
                for i, cat in enumerate(sorted(adata.obs[color_key].astype(str).unique(),
                                               key=lambda x: (len(x), x))):
                    fig_pca.axes[0].scatter([], [], color=NATURE_PALETTE[i % len(NATURE_PALETTE)], label=cat)
                fig_pca.axes[0].legend(frameon=False, fontsize=8)
            result_files.extend(self.save_matplotlib_figure(
                fig_pca, plots_dir, 'dimred_pca_scatter.png', 'pca',
                'PCA Scatter', formats=('png', 'svg'), dpi=300,
            ))

        # t-SNE plot
        if enable_tsne and 'X_tsne' in adata.obsm:
            for color_key in ['batch', 'leiden']:
                if color_key in adata.obs.columns:
                    tsne = adata.obsm['X_tsne']
                    color_vals = adata.obs[color_key].astype(str).values if color_key in adata.obs.columns else None
                    categories = sorted(set(color_vals))
                    colors = np.asarray([
                        NATURE_PALETTE[categories.index(value) % len(NATURE_PALETTE)]
                        for value in color_vals
                    ], dtype=object)
                    fig = scatter_figure(
                        tsne[:, 0], tsne[:, 1], title=f't-SNE by {color_key}',
                        x_label='tSNE1', y_label='tSNE2', colors=colors,
                        size=9, alpha=0.68,
                    )
                    for i, category in enumerate(categories):
                        fig.axes[0].scatter([], [], color=NATURE_PALETTE[i % len(NATURE_PALETTE)], label=category)
                    fig.axes[0].legend(frameon=False, fontsize=8)
                    result_files.extend(self.save_matplotlib_figure(
                        fig, plots_dir, f'dimred_tsne_{color_key}.png', 'tsne',
                        f't-SNE by {color_key}', formats=('png', 'svg'), dpi=300,
                    ))

        # Variance ratio plot
        if 'pca' in adata.uns:
            vr = adata.uns['pca']['variance_ratio'][:min(50, n_comps)]
            fig, axes = plt.subplots(figsize=(7.5, 5.0), dpi=150)
            x_pc = np.arange(len(vr))
            axes.bar(x_pc, vr, color=NATURE_PALETTE[0], alpha=0.88, label='Individual')
            ax2 = axes.twinx()
            ax2.plot(x_pc, np.cumsum(vr), color=NATURE_PALETTE[3], linewidth=1.8,
                     marker='o', markersize=3.5, label='Cumulative')
            tick_step = max(1, int(np.ceil(len(vr) / 10)))
            tick_idx = x_pc[::tick_step]
            axes.set_xticks(tick_idx, [f'PC{i + 1}' for i in tick_idx], rotation=45)
            axes.set_xlabel('PC', fontsize=9)
            axes.set_ylabel('Variance Ratio', fontsize=9)
            ax2.set_ylabel('Cumulative', fontsize=9)
            axes.set_title('PCA Variance Ratio', loc='left', fontsize=10,
                           fontweight='semibold', color=NATURE_TEXT)
            handles, labels = axes.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            axes.legend(handles + h2, labels + l2, frameon=False, fontsize=8)
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, 'dimred_pca_variance.png', 'pca',
                'PCA Variance Ratio', formats=('png', 'svg'), dpi=300,
            ))

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'dimred')

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_cells': adata.n_obs,
                'n_pcs': n_comps,
                'pca_variance_ratio_top5': round(float(adata.uns['pca']['variance_ratio'][:5].sum()), 3) if 'pca' in adata.uns else None,
                'embedding_method': embedding_method,
                'tsne_enabled': enable_tsne,
            }
        }
