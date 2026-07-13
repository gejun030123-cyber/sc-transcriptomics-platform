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
        from modules.visualization import umap_scatter
        import json
        import numpy as np

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
                fig_json = json.dumps(umap_scatter(adata, color_key, title=f'UMAP colored by {color_key}'))
                fpath = os.path.join(plots_dir, f'dimred_umap_{color_key}.json')
                with open(fpath, 'w') as f: f.write(fig_json)
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': f'UMAP by {color_key}'})
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
                    self.progress(-1, f'静态 UMAP 导出失败（不影响交互图）：{exc}')

        # PCA scatter for early detection of outliers and batch/sample structure
        if self.params.get('show_pca_scatter', True) and 'X_pca' in adata.obsm:
            import plotly.graph_objects as go
            pca = adata.obsm['X_pca']
            color_key = next((k for k in ['batch', 'phase', 'n_genes_by_counts', 'total_counts'] if k in adata.obs.columns), None)
            fig_pca = go.Figure()
            if color_key:
                values = adata.obs[color_key]
                try:
                    numeric_values = values.astype(float).values
                    fig_pca.add_trace(go.Scattergl(
                        x=pca[:, 0],
                        y=pca[:, 1],
                        mode='markers',
                        marker=dict(
                            size=4,
                            color=numeric_values,
                            colorscale='Viridis',
                            opacity=0.75,
                            colorbar=dict(title=color_key),
                        ),
                        text=adata.obs_names.tolist(),
                        hovertemplate='%{text}<br>' + color_key + ': %{marker.color:.3f}<extra></extra>',
                    ))
                except (TypeError, ValueError):
                    labels = values.astype(str)
                    for cat in sorted(labels.unique(), key=lambda x: (len(x), x)):
                        mask = labels == cat
                        fig_pca.add_trace(go.Scattergl(
                            x=pca[mask.values, 0],
                            y=pca[mask.values, 1],
                            mode='markers',
                            marker=dict(size=4, opacity=0.7),
                            name=str(cat),
                            text=adata.obs_names[mask.values].tolist(),
                            hovertemplate='%{text}<br>' + color_key + ': ' + str(cat) + '<extra></extra>',
                        ))
            else:
                fig_pca.add_trace(go.Scattergl(
                    x=pca[:, 0],
                    y=pca[:, 1],
                    mode='markers',
                    marker=dict(size=4, color='#3949ab', opacity=0.7),
                    text=adata.obs_names.tolist(),
                    hovertemplate='%{text}<extra></extra>',
                ))
            fig_pca.update_layout(
                title='PCA Scatter',
                xaxis_title='PC1',
                yaxis_title='PC2',
                plot_bgcolor='white',
                width=700,
                height=520,
            )
            result_files.append(self.save_plotly_json(
                fig_pca, plots_dir, 'dimred_pca_scatter.json',
                'pca', 'PCA Scatter'
            ))

        # t-SNE plot
        if enable_tsne and 'X_tsne' in adata.obsm:
            for color_key in ['batch', 'leiden']:
                if color_key in adata.obs.columns:
                    import plotly.graph_objects as go
                    tsne = adata.obsm['X_tsne']
                    color_vals = adata.obs[color_key].astype(str).values if color_key in adata.obs.columns else None
                    fig = go.Figure()
                    fig.add_trace(go.Scattergl(x=tsne[:, 0], y=tsne[:, 1], mode='markers',
                                               marker=dict(size=3, opacity=0.6), text=color_vals))
                    fig.update_layout(title=f't-SNE by {color_key}', xaxis_title='tSNE1', yaxis_title='tSNE2',
                                     plot_bgcolor='white', width=600, height=500)
                    result_files.append(self.save_plotly_json(fig, plots_dir, f'dimred_tsne_{color_key}.json', 'tsne', f't-SNE by {color_key}'))

        # Variance ratio plot
        if 'pca' in adata.uns:
            import plotly.graph_objects as go
            vr = adata.uns['pca']['variance_ratio'][:min(50, n_comps)]
            fig = go.Figure()
            fig.add_trace(go.Bar(y=vr, name='Individual'))
            fig.add_trace(go.Scatter(y=np.cumsum(vr), mode='lines', name='Cumulative'))
            fig.update_layout(title='PCA Variance Ratio', xaxis_title='PC', yaxis_title='Variance Ratio',
                             plot_bgcolor='white', width=600, height=400)
            result_files.append(self.save_plotly_json(fig, plots_dir, 'dimred_pca_variance.json', 'pca', 'PCA Variance Ratio'))

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
