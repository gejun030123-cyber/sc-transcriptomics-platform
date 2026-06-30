from modules.base import BaseAnalysis
import logging

logger = logging.getLogger(__name__)

class ClusteringAnalysis(BaseAnalysis):
    MODULE_NAME = "clustering"
    DISPLAY_NAME = "聚类分析"
    DESCRIPTION = "KNN 图 + 多分辨率 Leiden 聚类"
    INPUT_REQUIRES = ['X_umap']

    def validate_input(self, adata):
        if 'X_umap' not in adata.obsm:
            return "UMAP not found. Run dimensionality reduction first."
        return None

    def run(self, input_path):
        import os
        import scanpy as sc
        from modules.visualization import umap_scatter
        import json
        import numpy as np

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        resolutions = [float(r.strip()) for r in str(self.params.get('resolutions', '0.6,0.8,1.0')).split(',')]
        n_neighbors = int(self.params.get('n_neighbors', 15))
        clustering_method = self.params.get('clustering_method', 'leiden')
        n_iterations = int(self.params.get('n_iterations', 2))
        distance_metric = self.params.get('distance_metric', 'euclidean')
        use_corrected = self.params.get('use_corrected', True)
        auto_select = self.params.get('auto_select_resolution', False)
        resolution_metric = self.params.get('resolution_metric', 'silhouette')

        # Determine representation to use
        use_rep = 'X_pca'
        if use_corrected:
            for key in ['X_pca_harmony', 'X_pca_combat', 'X_scanorama', 'X_sysvi', 'X_scVI', 'X_bbknn']:
                if key in adata.obsm:
                    use_rep = key
                    break

        self.progress(20, f"Computing KNN graph (n_neighbors={n_neighbors})...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep, metric=distance_metric)

        self.progress(40, f"Running {clustering_method} clustering at resolutions: {resolutions}...")
        for i, res in enumerate(resolutions):
            if clustering_method == 'louvain':
                sc.tl.louvain(adata, resolution=res, key_added=f'leiden_{res}')
            else:
                sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}', flavor="igraph", n_iterations=n_iterations)
            pct = 40 + int((i + 1) / len(resolutions) * 30)
            self.progress(pct, f"{clustering_method} resolution {res} done")

        if 'leiden' not in adata.obs.columns:
            adata.obs['leiden'] = adata.obs[f'leiden_{resolutions[0]}'].copy()

        self.progress(75, "Generating cluster UMAP plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        for res in resolutions:
            key = f'leiden_{res}'
            if key in adata.obs.columns:
                n_clusters = adata.obs[key].nunique()
                fig_json = json.dumps(umap_scatter(adata, key, title=f'Leiden (res={res}, {n_clusters} clusters)'))
                fpath = os.path.join(plots_dir, f'cluster_umap_{res}.json')
                with open(fpath, 'w') as f: f.write(fig_json)
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': f'Clusters (res={res})'})

        # 多分辨率 UMAP 比较图
        if 'X_umap' in adata.obsm:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
            n_res = len(resolutions)
            fig_multi = make_subplots(rows=1, cols=n_res, subplot_titles=[f'res={r}' for r in resolutions])
            for ci, res in enumerate(resolutions, 1):
                key = f'leiden_{res}'
                if key in adata.obs.columns:
                    coords = adata.obsm['X_umap'][:, :2]
                    cats = adata.obs[key].values
                    for cat in sorted(set(cats)):
                        mask = cats == cat
                        fig_multi.add_trace(go.Scattergl(
                            x=coords[mask, 0], y=coords[mask, 1], mode='markers',
                            marker=dict(size=2, opacity=0.6), name=str(cat), showlegend=(ci==1)
                        ), row=1, col=ci)
            fig_multi.update_layout(height=400, width=350*n_res, title='多分辨率聚类比较')
            for i in range(1, n_res+1):
                fig_multi.update_xaxes(title_text='UMAP-1', row=1, col=i)
                fig_multi.update_yaxes(title_text='UMAP-2', row=1, col=i)
            fpath = os.path.join(plots_dir, 'cluster_multi_res_umap.json')
            with open(fpath, 'w') as f: f.write(fig_multi.to_json())
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': '多分辨率聚类比较'})

        # 聚类 marker dotplot
        try:
            from modules.annotation import DEFAULT_TME_MARKERS
            import plotly.graph_objects as go
            first_key = f'leiden_{resolutions[0]}'
            if first_key in adata.obs.columns:
                dotplot_genes = []
                for genes in DEFAULT_TME_MARKERS.values():
                    dotplot_genes.extend([g for g in genes[:2] if g in adata.var_names])
                dotplot_genes = list(dict.fromkeys(dotplot_genes))[:20]
                if dotplot_genes:
                    import scanpy as sc2
                    sc2.tl.dendrogram(adata, groupby=first_key)
                    fig_dot = sc2.pl.dotplot(adata, var_names=dotplot_genes, groupby=first_key, return_fig=True)
                    import io, base64
                    buf = io.BytesIO()
                    fig_dot.savefig(buf, format='png', dpi=100, bbox_inches='tight')
                    import matplotlib.pyplot as plt
                    plt.close('all')
                    buf.seek(0)
                    img_b64 = base64.b64encode(buf.read()).decode()
                    fpath = os.path.join(plots_dir, 'cluster_dotplot.json')
                    with open(fpath, 'w') as f:
                        json.dump({'data': [{'type': 'image', 'source': f'data:image/png;base64,{img_b64}', 'xref': 'paper', 'yref': 'paper', 'x': 0, 'y': 1, 'sizex': 1, 'sizey': 1, 'sizing': 'stretch'}], 'layout': {'width': 800, 'height': 500, 'title': f'Marker Dotplot ({first_key})'}}, f)
                    result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'dotplot', 'label': f'Marker Dotplot'})
        except Exception as e:
            logger.warning("生成 marker dotplot 失败（注释模块可能不可用）: %s", e)

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'clustering')

        # Auto-select best resolution
        best_res = resolutions[0]
        if auto_select and len(resolutions) > 1:
            try:
                from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
                best_score = -float('inf')
                for res in resolutions:
                    key = f'leiden_{res}'
                    if key not in adata.obs.columns:
                        continue
                    labels = adata.obs[key].astype('category').cat.codes.values
                    rep_data = adata.obsm[use_rep]
                    if resolution_metric == 'silhouette':
                        score = silhouette_score(rep_data, labels)
                    elif resolution_metric == 'calinski':
                        score = calinski_harabasz_score(rep_data, labels)
                    elif resolution_metric == 'davies_bouldin':
                        score = -davies_bouldin_score(rep_data, labels)  # negate: lower is better
                    else:
                        score = silhouette_score(rep_data, labels)
                    if score > best_score:
                        best_score = score
                        best_res = res
            except Exception as e:
                logger.warning("自动选择分辨率失败: %s", e)

        self.progress(100, "Done")
        summary = {f'n_clusters_{res}': int(adata.obs[f'leiden_{res}'].nunique()) for res in resolutions if f'leiden_{res}' in adata.obs.columns}
        summary['resolutions'] = resolutions
        summary['best_resolution'] = best_res
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }
