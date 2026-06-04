import os
from modules.base import BaseAnalysis

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
        import scanpy as sc
        from modules.visualization import umap_scatter
        import json

        self.progress(5, "Loading data...")
        adata = sc.read_h5ad(input_path)
        resolutions = [float(r.strip()) for r in str(self.params.get('resolutions', '0.6,0.8,1.0')).split(',')]
        n_neighbors = int(self.params.get('n_neighbors', 15))

        self.progress(20, f"Computing KNN graph (n_neighbors={n_neighbors})...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors)

        self.progress(40, f"Running Leiden clustering at resolutions: {resolutions}...")
        for i, res in enumerate(resolutions):
            sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
            pct = 40 + int((i + 1) / len(resolutions) * 30)
            self.progress(pct, f"Leiden resolution {res} done")

        if 'leiden' not in adata.obs.columns:
            adata.obs['leiden'] = adata.obs[f'leiden_{resolutions[0]}'].copy()

        self.progress(75, "Generating cluster UMAP plots...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        for res in resolutions:
            key = f'leiden_{res}'
            if key in adata.obs.columns:
                n_clusters = adata.obs[key].nunique()
                fig_json = json.dumps(umap_scatter(adata, key, title=f'Leiden (res={res}, {n_clusters} clusters)'))
                fpath = os.path.join(plots_dir, f'cluster_umap_{res}.json')
                with open(fpath, 'w') as f: f.write(fig_json)
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': f'Clusters (res={res})'})

        self.progress(90, "Saving output...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'clustering_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "Done")
        summary = {f'n_clusters_{res}': int(adata.obs[f'leiden_{res}'].nunique()) for res in resolutions if f'leiden_{res}' in adata.obs.columns}
        summary['resolutions'] = resolutions
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }
