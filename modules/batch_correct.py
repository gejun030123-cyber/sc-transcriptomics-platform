import os
from modules.base import BaseAnalysis

class BatchCorrectAnalysis(BaseAnalysis):
    MODULE_NAME = "batch_correct"
    DISPLAY_NAME = "Batch Correction"
    DESCRIPTION = "Harmony, ComBat, or SysVI batch correction"
    INPUT_REQUIRES = ['X_pca']

    def validate_input(self, adata):
        if 'X_pca' not in adata.obsm:
            return "PCA not found. Run dimensionality reduction first."
        return None

    def run(self, input_path):
        import scanpy as sc
        import omicverse as ov
        from modules.visualization import umap_scatter
        import json

        self.progress(5, "Loading data...")
        adata = sc.read_h5ad(input_path)
        method = self.params.get('method', 'harmony')
        batch_key = self.params.get('batch_key', 'batch')
        n_pcs = int(self.params.get('n_pcs', 50))

        if method == 'harmony':
            self.progress(20, "Running Harmony batch correction...")
            ov.single.batch_correction(adata, batch_key=batch_key, methods=['harmony'], n_comps=n_pcs)
            adata.obsm['X_pca_harmony'] = adata.obsm['X_pca_harmony'].copy()
            corrected_key = 'X_pca_harmony'

        elif method == 'combat':
            self.progress(20, "Running ComBat batch correction...")
            ov.single.batch_correction(adata, batch_key=batch_key, methods=['combat'], n_comps=n_pcs)
            corrected_key = 'X_pca_combat'

        elif method == 'sysvi':
            max_epochs = int(self.params.get('max_epochs', 200))
            self.progress(20, f"Running SysVI batch correction ({max_epochs} epochs)...")
            import os as _os
            _os.environ['CUDA_VISIBLE_DEVICES'] = '0'
            from scvi.external import SysVI
            adata.layers['counts'] = adata.X.copy()
            SysVI.setup_anndata(adata, layer='counts', batch_key=batch_key)
            model = SysVI(adata)
            model.train(max_epochs=max_epochs, accelerator='gpu')
            latent = model.get_latent_representation()
            adata.obsm['X_sysvi'] = latent
            corrected_key = 'X_sysvi'
        else:
            return {'output_adata': input_path, 'result_files': [], 'summary': {'error': f'Unknown method: {method}'}}

        self.progress(60, "Computing UMAP on corrected embeddings...")
        sc.pp.neighbors(adata, use_rep=corrected_key, n_pcs=n_pcs)
        sc.tl.umap(adata)

        self.progress(75, "Generating comparison plots...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        fig_json = json.dumps(umap_scatter(adata, 'batch', title=f'UMAP after {method} correction (by batch)'))
        fpath = os.path.join(plots_dir, f'batch_umap_{method}.json')
        with open(fpath, 'w') as f: f.write(fig_json)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': f'UMAP after {method} (batch)'})

        self.progress(90, "Saving output...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, f'batch_correct_{method}_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'method': method,
                'n_cells': adata.n_obs,
                'embedding_key': corrected_key,
            }
        }
