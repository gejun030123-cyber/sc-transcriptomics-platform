import os
from modules.base import BaseAnalysis

class BatchCorrectAnalysis(BaseAnalysis):
    MODULE_NAME = "batch_correct"
    DISPLAY_NAME = "批次校正"
    DESCRIPTION = "Harmony、ComBat、BBKNN、Scanorama、SysVI 批次效应校正"
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
        import numpy as np

        self.progress(5, "Loading data...")
        adata = sc.read_h5ad(input_path)
        from modules.io_utils import remap_var_names
        adata = remap_var_names(adata)
        method = self.params.get('method', 'harmony')
        batch_key = self.params.get('batch_key', 'batch')
        n_pcs = int(self.params.get('n_pcs', 50))

        if method == 'harmony':
            self.progress(20, "Running Harmony batch correction...")
            theta = float(self.params.get('harmony_theta', 2.0))
            lam = float(self.params.get('harmony_lambda', 1.0))
            ov.single.batch_correction(adata, batch_key=batch_key, methods=['harmony'], n_comps=n_pcs)
            adata.obsm['X_pca_harmony'] = adata.obsm['X_pca_harmony'].copy()
            corrected_key = 'X_pca_harmony'

        elif method == 'combat':
            self.progress(20, "Running ComBat batch correction...")
            ov.single.batch_correction(adata, batch_key=batch_key, methods=['combat'], n_comps=n_pcs)
            corrected_key = 'X_pca_combat'

        elif method == 'bbknn':
            self.progress(20, "Running BBKNN batch correction...")
            neighbors_within_batch = int(self.params.get('bbknn_neighbors_within_batch', 3))
            import bbknn
            bbknn.bbknn(adata, batch_key=batch_key, neighbors_within_batch=neighbors_within_batch)
            corrected_key = 'X_pca'  # BBKNN modifies the neighbor graph, not PCA

        elif method == 'scanorama':
            self.progress(20, "Running Scanorama batch correction...")
            import scanorama
            # Split by batch, integrate
            batches = []
            batch_names = adata.obs[batch_key].unique()
            for b in batch_names:
                batches.append(adata[adata.obs[batch_key] == b])
            corrected = scanorama.integrate_scanpy(batches, dimred=n_pcs)
            # Reassemble
            import scipy.sparse as sp
            all_corrected = np.vstack([c for c in corrected])
            adata.obsm['X_scanorama'] = all_corrected
            corrected_key = 'X_scanorama'

        elif method == 'sysvi':
            max_epochs = int(self.params.get('max_epochs', 200))
            scvi_n_latent = int(self.params.get('scvi_n_latent', 30))
            scvi_n_hidden = int(self.params.get('scvi_n_hidden', 128))
            scvi_n_layers = int(self.params.get('scvi_n_layers', 1))
            scvi_dropout_rate = float(self.params.get('scvi_dropout_rate', 0.1))

            self.progress(20, f"Running SysVI batch correction ({max_epochs} epochs)...")
            import os as _os
            _os.environ['CUDA_VISIBLE_DEVICES'] = '0'
            from scvi.external import SysVI
            adata.layers['counts'] = adata.X.copy()
            SysVI.setup_anndata(adata, layer='counts', batch_key=batch_key)
            model = SysVI(adata, n_latent=scvi_n_latent, n_hidden=scvi_n_hidden,
                          n_layers=scvi_n_layers, dropout_rate=scvi_dropout_rate)
            model.train(max_epochs=max_epochs, accelerator='gpu')
            latent = model.get_latent_representation()
            adata.obsm['X_sysvi'] = latent
            corrected_key = 'X_sysvi'

        elif method == 'scvi':
            scvi_n_latent = int(self.params.get('scvi_n_latent', 30))
            scvi_n_hidden = int(self.params.get('scvi_n_hidden', 128))
            scvi_n_layers = int(self.params.get('scvi_n_layers', 1))
            scvi_dropout_rate = float(self.params.get('scvi_dropout_rate', 0.1))
            scvi_learning_rate = float(self.params.get('scvi_learning_rate', 1e-3))
            max_epochs = int(self.params.get('max_epochs', 200))

            self.progress(20, f"Running scVI batch correction ({max_epochs} epochs)...")
            import scvi
            adata.layers['counts'] = adata.X.copy()
            scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key=batch_key)
            model = scvi.model.SCVI(adata, n_latent=scvi_n_latent, n_hidden=scvi_n_hidden,
                                     n_layers=scvi_n_layers, dropout_rate=scvi_dropout_rate)
            model.train(max_epochs=max_epochs, lr=scvi_learning_rate)
            latent = model.get_latent_representation()
            adata.obsm['X_scVI'] = latent
            corrected_key = 'X_scVI'

        else:
            return {'output_adata': input_path, 'result_files': [], 'summary': {'error': f'Unknown method: {method}'}}

        self.progress(60, "Computing UMAP on corrected embeddings...")
        sc.pp.neighbors(adata, use_rep=corrected_key, n_pcs=n_pcs)
        sc.tl.umap(adata)

        self.progress(70, "Generating comparison plots...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        result_files = []

        fig_json = json.dumps(umap_scatter(adata, 'batch', title=f'UMAP after {method} correction (by batch)'))
        fpath = os.path.join(plots_dir, f'batch_umap_{method}.json')
        with open(fpath, 'w') as f: f.write(fig_json)
        result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': f'UMAP after {method} (batch)'})

        # Before/after comparison (if original UMAP exists)
        if 'X_umap' in adata.obsm:
            fig_before = json.dumps(umap_scatter(adata, 'batch', title='UMAP before correction (by batch)'))
            fpath_before = os.path.join(plots_dir, 'batch_umap_before.json')
            with open(fpath_before, 'w') as f: f.write(fig_before)
            result_files.append({'file_path': fpath_before, 'file_type': 'plotly_json', 'category': 'umap', 'label': 'UMAP before correction (batch)'})

        # Batch correction evaluation
        evaluate = self.params.get('evaluate_correction', False)
        eval_metrics = {}
        if evaluate and batch_key in adata.obs.columns:
            self.progress(80, "Evaluating batch correction...")
            try:
                from sklearn.metrics import silhouette_score
                # ASW for batch mixing (lower = better mixing)
                batch_labels = adata.obs[batch_key].astype('category').cat.codes.values
                asw_batch = silhouette_score(adata.obsm[corrected_key], batch_labels)
                eval_metrics['asw_batch'] = round(float(asw_batch), 4)

                # ASW for cell type purity (if available)
                for ct_key in ['celltype', 'leiden']:
                    if ct_key in adata.obs.columns:
                        ct_labels = adata.obs[ct_key].astype('category').cat.codes.values
                        if len(set(ct_labels)) > 1:
                            asw_ct = silhouette_score(adata.obsm[corrected_key], ct_labels)
                            eval_metrics[f'asw_{ct_key}'] = round(float(asw_ct), 4)
                            break

                # Graph connectivity
                import scipy.sparse as sp
                from scipy.sparse.csgraph import connected_components
                n_components, _ = connected_components(adata.obsp['connectivities'], directed=False)
                eval_metrics['graph_components'] = int(n_components)
            except Exception as e:
                eval_metrics['error'] = str(e)

        self.progress(90, "Saving output...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, f'batch_correct_{method}_output.h5ad')
        adata.write_h5ad(output_path)

        summary = {
            'method': method,
            'n_cells': adata.n_obs,
            'embedding_key': corrected_key,
        }
        if eval_metrics:
            summary['eval_metrics'] = eval_metrics

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }
