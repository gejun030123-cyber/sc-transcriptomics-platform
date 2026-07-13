from modules.base import BaseAnalysis


class BatchCorrectAnalysis(BaseAnalysis):
    MODULE_NAME = "batch_correct"
    DISPLAY_NAME = "批次校正"
    DESCRIPTION = "Harmony、ComBat、BBKNN、Scanorama、SysVI、scVI 批次效应整合与评价"
    INPUT_REQUIRES = ['X_pca']

    def validate_input(self, adata):
        if 'X_pca' not in adata.obsm:
            return "PCA not found. Run dimensionality reduction first."
        return None

    @staticmethod
    def _as_bool(value, default=False):
        if value is None:
            return default
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            return value.lower() in {'1', 'true', 'yes', 'on'}
        return bool(value)

    @staticmethod
    def _batch_entropy(counts):
        import numpy as np

        total = float(counts.sum())
        if total <= 0:
            return 0.0
        probs = counts[counts > 0] / total
        if len(probs) <= 1:
            return 0.0
        return float(-(probs * np.log(probs)).sum() / np.log(len(counts)))

    def _hvg_adata(self, adata):
        if 'highly_variable' in adata.var.columns and adata.var['highly_variable'].sum() > 0:
            return adata[:, adata.var['highly_variable'].astype(bool)].copy()
        return adata.copy()

    def _run_harmony(self, adata, batch_key, n_pcs):
        import harmonypy as hm
        import numpy as np

        theta = float(self.params.get('harmony_theta', 2.0))
        lamb = float(self.params.get('harmony_lambda', 1.0))
        pcs = np.asarray(adata.obsm['X_pca'][:, :n_pcs])
        try:
            result = hm.run_harmony(
                pcs,
                adata.obs,
                vars_use=[batch_key],
                theta=theta,
                lamb=lamb,
                max_iter_harmony=int(self.params.get('harmony_max_iter', 20)),
            )
        except TypeError:
            result = hm.run_harmony(
                pcs,
                adata.obs,
                vars_use=[batch_key],
                theta=theta,
                max_iter_harmony=int(self.params.get('harmony_max_iter', 20)),
            )
        embedding = result.Z_corr
        if embedding.shape[0] == pcs.shape[1]:
            embedding = embedding.T
        adata.obsm['X_pca_harmony'] = embedding
        return 'X_pca_harmony', {'theta': theta, 'lambda': lamb}

    def _run_combat(self, adata, batch_key, n_pcs):
        import scanpy as sc

        combat_adata = self._hvg_adata(adata)
        sc.pp.combat(combat_adata, key=batch_key)
        sc.pp.scale(combat_adata, max_value=10)
        actual_n_pcs = max(2, min(n_pcs, combat_adata.n_obs - 1, combat_adata.n_vars - 1))
        sc.tl.pca(combat_adata, n_comps=actual_n_pcs, svd_solver='arpack')
        adata.obsm['X_pca_combat'] = combat_adata.obsm['X_pca'][:, :actual_n_pcs].copy()
        return 'X_pca_combat', {'n_hvg_for_combat': int(combat_adata.n_vars)}

    def _run_bbknn(self, adata, batch_key):
        import bbknn

        requested_neighbors = max(1, int(self.params.get('bbknn_neighbors_within_batch', 3)))
        batch_counts = adata.obs[batch_key].astype(str).value_counts()
        min_batch_size = int(batch_counts.min())
        smallest_batches = [str(k) for k, v in batch_counts.items() if int(v) == min_batch_size]
        neighbors_within_batch = min(requested_neighbors, min_batch_size)
        method_info = {
            'requested_neighbors_within_batch': requested_neighbors,
            'neighbors_within_batch': neighbors_within_batch,
            'min_batch_size': min_batch_size,
            'smallest_batches': smallest_batches,
            'batch_counts': {str(k): int(v) for k, v in batch_counts.to_dict().items()},
            'graph_method': 'bbknn',
        }
        if neighbors_within_batch < requested_neighbors:
            method_info['warning'] = (
                'BBKNN neighbors_within_batch was automatically reduced from '
                f'{requested_neighbors} to {neighbors_within_batch} because the smallest '
                f'batch has {min_batch_size} cell(s).'
            )
        bbknn.bbknn(
            adata,
            batch_key=batch_key,
            neighbors_within_batch=neighbors_within_batch,
            use_rep='X_pca',
        )
        return 'X_pca', method_info

    def _run_scanorama(self, adata, batch_key, n_pcs):
        import numpy as np
        import scanorama

        work = self._hvg_adata(adata)
        batch_values = list(work.obs[batch_key].astype(str).unique())
        batches = [work[work.obs[batch_key].astype(str) == b].copy() for b in batch_values]
        scanorama.integrate_scanpy(batches, dimred=n_pcs)

        embedding = np.zeros((work.n_obs, n_pcs), dtype=np.float32)
        for batch_adata in batches:
            if 'X_scanorama' not in batch_adata.obsm:
                raise RuntimeError("Scanorama did not produce obsm['X_scanorama']")
            corrected = batch_adata.obsm['X_scanorama']
            positions = work.obs_names.get_indexer(batch_adata.obs_names)
            embedding[positions, :corrected.shape[1]] = corrected

        target = np.zeros((adata.n_obs, embedding.shape[1]), dtype=np.float32)
        target_positions = adata.obs_names.get_indexer(work.obs_names)
        target[target_positions, :] = embedding
        adata.obsm['X_scanorama'] = target
        return 'X_scanorama', {'n_hvg_for_scanorama': int(work.n_vars)}

    def _run_sysvi(self, adata, batch_key):
        import scvi
        from scvi.external import SysVI

        max_epochs = int(self.params.get('max_epochs', self.params.get('sysvi_epochs', 60)))
        scvi_n_latent = int(self.params.get('scvi_n_latent', 30))
        scvi_n_hidden = int(self.params.get('scvi_n_hidden', 128))
        scvi_n_layers = int(self.params.get('scvi_n_layers', 1))
        scvi_dropout_rate = float(self.params.get('scvi_dropout_rate', 0.1))
        sysvi_cycle_weight = float(self.params.get('sysvi_cycle_weight', 5.0))
        sysvi_kl_weight = float(self.params.get('sysvi_kl_weight', 1.0))
        sysvi_prior = self.params.get('sysvi_prior', 'vamp')
        sysvi_n_prior_components = int(self.params.get('sysvi_n_prior_components', 5))

        work = self._hvg_adata(adata)
        scvi.settings.seed = int(self.params.get('sysvi_seed', 0))
        SysVI.setup_anndata(work, batch_key=batch_key)
        model = SysVI(
            work,
            prior=sysvi_prior,
            n_prior_components=sysvi_n_prior_components,
            n_latent=scvi_n_latent,
            n_hidden=scvi_n_hidden,
            n_layers=scvi_n_layers,
            dropout_rate=scvi_dropout_rate,
            embed_categorical_covariates=self._as_bool(
                self.params.get('sysvi_embed_categorical_covariates', False)
            ),
        )
        train_kwargs = {
            'max_epochs': max_epochs,
            'check_val_every_n_epoch': 1,
            'plan_kwargs': {
                'kl_weight': sysvi_kl_weight,
                'z_distance_cycle_weight': sysvi_cycle_weight,
            },
        }
        try:
            train_kwargs.update({'accelerator': 'auto', 'devices': 'auto'})
            model.train(**train_kwargs)
        except TypeError:
            train_kwargs.pop('accelerator', None)
            train_kwargs.pop('devices', None)
            model.train(**train_kwargs)

        latent = model.get_latent_representation(adata=work)
        target = latent
        if work.n_obs != adata.n_obs or list(work.obs_names) != list(adata.obs_names):
            import numpy as np
            target = np.zeros((adata.n_obs, latent.shape[1]), dtype=latent.dtype)
            target_positions = adata.obs_names.get_indexer(work.obs_names)
            target[target_positions, :] = latent
        adata.obsm['X_sysvi'] = target
        return 'X_sysvi', {
            'max_epochs': max_epochs,
            'n_hvg_for_sysvi': int(work.n_vars),
            'cycle_weight': sysvi_cycle_weight,
            'kl_weight': sysvi_kl_weight,
            'prior': sysvi_prior,
        }

    def _run_scvi(self, adata, batch_key):
        import scvi

        scvi_n_latent = int(self.params.get('scvi_n_latent', 30))
        scvi_n_hidden = int(self.params.get('scvi_n_hidden', 128))
        scvi_n_layers = int(self.params.get('scvi_n_layers', 1))
        scvi_dropout_rate = float(self.params.get('scvi_dropout_rate', 0.1))
        scvi_learning_rate = float(self.params.get('scvi_learning_rate', 1e-3))
        max_epochs = int(self.params.get('max_epochs', 80))

        work = self._hvg_adata(adata)
        layer = 'counts' if 'counts' in work.layers else None
        scvi.model.SCVI.setup_anndata(work, layer=layer, batch_key=batch_key)
        model = scvi.model.SCVI(
            work,
            n_latent=scvi_n_latent,
            n_hidden=scvi_n_hidden,
            n_layers=scvi_n_layers,
            dropout_rate=scvi_dropout_rate,
        )
        try:
            model.train(max_epochs=max_epochs, lr=scvi_learning_rate, accelerator='auto', devices='auto')
        except TypeError:
            model.train(max_epochs=max_epochs, lr=scvi_learning_rate)
        latent = model.get_latent_representation()

        target = latent
        if work.n_obs != adata.n_obs or list(work.obs_names) != list(adata.obs_names):
            import numpy as np
            target = np.zeros((adata.n_obs, latent.shape[1]), dtype=latent.dtype)
            target_positions = adata.obs_names.get_indexer(work.obs_names)
            target[target_positions, :] = latent
        adata.obsm['X_scVI'] = target
        return 'X_scVI', {'max_epochs': max_epochs, 'n_hvg_for_scvi': int(work.n_vars), 'layer': layer}

    def _compute_evaluation(self, adata, corrected_key, batch_key):
        import numpy as np
        import pandas as pd
        import scanpy as sc
        from sklearn.metrics import silhouette_score
        from scipy.sparse.csgraph import connected_components

        metrics = {}
        sample_size = int(self.params.get('evaluation_sample_size', 10000))
        sample_size = max(1000, min(sample_size, adata.n_obs))
        batch_labels = adata.obs[batch_key].astype('category')
        metrics['n_batches'] = int(batch_labels.nunique())
        metrics['batch_counts'] = {
            str(k): int(v) for k, v in adata.obs[batch_key].astype(str).value_counts().to_dict().items()
        }

        rep = adata.obsm[corrected_key]
        if batch_labels.nunique() > 1:
            asw = silhouette_score(rep, batch_labels.cat.codes.values, sample_size=sample_size, random_state=7)
            metrics['asw_batch'] = round(float(asw), 4)
            metrics['abs_asw_batch'] = round(abs(float(asw)), 4)

        for label_key in [
            self.params.get('bio_label_key', ''),
            'celltype',
            'reference_celltype',
            'leiden',
        ]:
            if not label_key or label_key not in adata.obs.columns:
                continue
            labels = adata.obs[label_key].astype('category')
            if labels.nunique() <= 1:
                continue
            try:
                asw_label = silhouette_score(rep, labels.cat.codes.values, sample_size=sample_size, random_state=7)
                metrics[f'asw_{label_key}'] = round(float(asw_label), 4)
                break
            except Exception as exc:
                metrics[f'asw_{label_key}_error'] = str(exc)
                break

        if 'connectivities' in adata.obsp:
            conn = adata.obsp['connectivities']
            n_components, component_labels = connected_components(conn, directed=False)
            metrics['graph_components'] = int(n_components)
            if len(component_labels):
                largest = np.bincount(component_labels).max()
                metrics['largest_graph_component_fraction'] = round(float(largest / adata.n_obs), 4)

            try:
                csr = conn.tocsr()
                batch_codes = batch_labels.cat.codes.to_numpy()
                neighbor_entropies = []
                same_batch_fracs = []
                n_batches = max(1, batch_labels.nunique())
                for i in range(csr.shape[0]):
                    start, end = csr.indptr[i], csr.indptr[i + 1]
                    if start == end:
                        continue
                    indices = csr.indices[start:end]
                    weights = csr.data[start:end]
                    total_weight = weights.sum()
                    if total_weight <= 0:
                        continue
                    counts = np.zeros(n_batches, dtype=float)
                    np.add.at(counts, batch_codes[indices], weights)
                    same_batch_fracs.append(float(counts[batch_codes[i]] / total_weight))
                    neighbor_entropies.append(self._batch_entropy(counts))
                if neighbor_entropies:
                    metrics['mean_neighbor_batch_entropy'] = round(float(np.mean(neighbor_entropies)), 4)
                    metrics['mean_neighbor_same_batch_fraction'] = round(float(np.mean(same_batch_fracs)), 4)
            except Exception as exc:
                metrics['neighbor_mixing_error'] = str(exc)

        cluster_key = self.params.get('evaluation_cluster_key', '').strip()
        if not cluster_key:
            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else ''
        if not cluster_key or cluster_key not in adata.obs.columns:
            cluster_key = 'batch_eval_leiden'
            if cluster_key not in adata.obs.columns:
                eval_resolution = float(self.params.get('evaluation_resolution', 0.8))
                sc.tl.leiden(
                    adata,
                    resolution=eval_resolution,
                    key_added=cluster_key,
                    flavor='igraph',
                    n_iterations=2,
                )
                metrics['evaluation_resolution'] = eval_resolution

        table = pd.crosstab(adata.obs[cluster_key].astype(str), adata.obs[batch_key].astype(str))
        if not table.empty:
            cluster_sizes = table.sum(axis=1)
            fractions = table.div(cluster_sizes, axis=0)
            entropies = table.apply(lambda row: self._batch_entropy(row.values.astype(float)), axis=1)
            max_fracs = fractions.max(axis=1)
            weights = cluster_sizes / cluster_sizes.sum()
            metrics['evaluation_cluster_key'] = cluster_key
            metrics['n_eval_clusters'] = int(table.shape[0])
            metrics['weighted_batch_entropy'] = round(float((entropies * weights).sum()), 4)
            metrics['weighted_max_batch_fraction'] = round(float((max_fracs * weights).sum()), 4)
            metrics['batch_dominated_cluster_count_90pct'] = int((max_fracs >= 0.9).sum())
            metrics['single_batch_cluster_count'] = int((entropies == 0).sum())
            metrics['min_cluster_size'] = int(cluster_sizes.min())
            metrics['median_cluster_size'] = round(float(cluster_sizes.median()), 2)
            metrics['small_clusters_lt_50'] = int((cluster_sizes < 50).sum())

        return metrics

    def _make_eval_plots(self, adata, metrics, batch_key, method, plots_dir):
        import plotly.graph_objects as go

        result_files = []
        scalar_keys = [
            'abs_asw_batch',
            'weighted_batch_entropy',
            'weighted_max_batch_fraction',
            'mean_neighbor_batch_entropy',
            'mean_neighbor_same_batch_fraction',
            'largest_graph_component_fraction',
        ]
        x, y = [], []
        for key in scalar_keys:
            value = metrics.get(key)
            if isinstance(value, (int, float)):
                x.append(key)
                y.append(value)
        if x:
            fig = go.Figure(go.Bar(x=x, y=y, marker_color='#2563eb'))
            fig.update_layout(
                title=f'Batch Integration Metrics ({method})',
                xaxis_title='Metric',
                yaxis_title='Value',
                plot_bgcolor='white',
                width=850,
                height=420,
            )
            result_files.append(self.save_plotly_json(
                fig, plots_dir, f'batch_eval_metrics_{method}.json', 'metrics',
                f'Batch integration metrics ({method})'
            ))

        cluster_key = metrics.get('evaluation_cluster_key')
        if cluster_key and cluster_key in adata.obs.columns:
            import pandas as pd
            table = pd.crosstab(
                adata.obs[cluster_key].astype(str),
                adata.obs[batch_key].astype(str),
                normalize='index',
            )
            fig = go.Figure()
            for batch in table.columns:
                fig.add_trace(go.Bar(
                    x=table.index.tolist(),
                    y=table[batch].values,
                    name=str(batch),
                    hovertemplate='Cluster: %{x}<br>' + batch_key + ': ' + str(batch)
                    + '<br>Fraction: %{y:.2%}<extra></extra>',
                ))
            fig.update_layout(
                title=f'Evaluation Cluster Batch Composition ({method})',
                xaxis_title=cluster_key,
                yaxis_title='Fraction',
                barmode='stack',
                plot_bgcolor='white',
                width=max(750, 50 * max(1, len(table.index))),
                height=460,
            )
            result_files.append(self.save_plotly_json(
                fig, plots_dir, f'batch_eval_cluster_composition_{method}.json', 'bar',
                f'Evaluation cluster batch composition ({method})'
            ))

        table_rows = []
        for key, value in metrics.items():
            if isinstance(value, dict):
                value = ', '.join([f'{k}: {v}' for k, v in value.items()])
            table_rows.append([key, str(value)])
        if table_rows:
            fig = go.Figure(data=[go.Table(
                header=dict(values=['Metric', 'Value'], fill_color='#dbeafe', align='left'),
                cells=dict(values=list(zip(*table_rows)), fill_color='white', align='left'),
            )])
            fig.update_layout(title=f'Batch Integration Metric Table ({method})', width=850, height=520)
            result_files.append(self.save_plotly_json(
                fig, plots_dir, f'batch_eval_metric_table_{method}.json', 'table',
                f'Batch integration metric table ({method})'
            ))

        return result_files

    def run(self, input_path):
        import json
        import os
        import scanpy as sc
        from modules.visualization import umap_scatter

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        method = self.params.get('method', 'harmony')
        batch_key = self.params.get('batch_key', 'batch')
        if batch_key not in adata.obs.columns:
            raise ValueError(f"batch_key '{batch_key}' 不在 adata.obs 中，可用列: {list(adata.obs.columns)}")
        if adata.obs[batch_key].nunique() < 2:
            raise ValueError(f"batch_key '{batch_key}' 只有一个取值，无法评估或执行批次整合")

        n_pcs = int(self.params.get('n_pcs', 50))
        if 'X_pca' in adata.obsm:
            n_pcs = max(2, min(n_pcs, adata.obsm['X_pca'].shape[1]))
        umap_before = adata.obsm['X_umap'].copy() if 'X_umap' in adata.obsm else None
        method_info = {}
        graph_already_built = False

        if method == 'harmony':
            self.progress(20, "Running Harmony batch correction...")
            corrected_key, method_info = self._run_harmony(adata, batch_key, n_pcs)
        elif method == 'combat':
            self.progress(20, "Running Scanpy ComBat PCA correction...")
            corrected_key, method_info = self._run_combat(adata, batch_key, n_pcs)
        elif method == 'bbknn':
            self.progress(20, "Running BBKNN graph integration...")
            corrected_key, method_info = self._run_bbknn(adata, batch_key)
            graph_already_built = True
        elif method == 'scanorama':
            self.progress(20, "Running Scanorama integration...")
            corrected_key, method_info = self._run_scanorama(adata, batch_key, n_pcs)
        elif method == 'sysvi':
            self.progress(20, "Running SysVI integration...")
            corrected_key, method_info = self._run_sysvi(adata, batch_key)
        elif method == 'scvi':
            self.progress(20, "Running scVI integration...")
            corrected_key, method_info = self._run_scvi(adata, batch_key)
        else:
            return {'output_adata': input_path, 'result_files': [], 'summary': {'error': f'Unknown method: {method}'}}

        self.progress(60, "Computing UMAP on integrated representation...")
        if graph_already_built:
            sc.tl.umap(adata)
        else:
            n_pcs_kwargs = {}
            if corrected_key.startswith('X_pca'):
                n_pcs_kwargs['n_pcs'] = min(n_pcs, adata.obsm[corrected_key].shape[1])
            sc.pp.neighbors(adata, use_rep=corrected_key, **n_pcs_kwargs)
            sc.tl.umap(adata)

        self.progress(70, "Generating comparison plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        fig_json = json.dumps(umap_scatter(
            adata, batch_key, title=f'UMAP after {method} integration (by {batch_key})'
        ))
        fpath = os.path.join(plots_dir, f'batch_umap_{method}.json')
        with open(fpath, 'w') as handle:
            handle.write(fig_json)
        result_files.append({
            'file_path': fpath,
            'file_type': 'plotly_json',
            'category': 'umap',
            'label': f'UMAP after {method} (batch)',
        })

        if umap_before is not None:
            corrected_umap = adata.obsm['X_umap'].copy()
            adata.obsm['X_umap'] = umap_before
            fig_before = json.dumps(umap_scatter(
                adata, batch_key, title=f'UMAP before integration (by {batch_key})'
            ))
            adata.obsm['X_umap'] = corrected_umap
            fpath_before = os.path.join(plots_dir, 'batch_umap_before.json')
            with open(fpath_before, 'w') as handle:
                handle.write(fig_before)
            result_files.append({
                'file_path': fpath_before,
                'file_type': 'plotly_json',
                'category': 'umap',
                'label': 'UMAP before integration (batch)',
            })

        evaluate = self._as_bool(self.params.get('evaluate_correction', False))
        eval_metrics = {}
        if evaluate:
            self.progress(80, "Evaluating batch integration...")
            try:
                eval_metrics = self._compute_evaluation(adata, corrected_key, batch_key)
                result_files.extend(self._make_eval_plots(adata, eval_metrics, batch_key, method, plots_dir))
            except Exception as exc:
                eval_metrics = {'error': str(exc)}

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'batch_correct')

        summary = {
            'method': method,
            'n_cells': int(adata.n_obs),
            'n_genes': int(adata.n_vars),
            'embedding_key': corrected_key,
            'method_info': method_info,
        }
        if eval_metrics:
            summary['eval_metrics'] = eval_metrics

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }
