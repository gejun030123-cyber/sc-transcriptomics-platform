from modules.base import BaseAnalysis
from modules.io_utils import resolve_obs_grouping, obs_grouping_info


class BatchCorrectAnalysis(BaseAnalysis):
    MODULE_NAME = "batch_correct"
    DISPLAY_NAME = "批次校正"
    DESCRIPTION = "Harmony、ComBat、BBKNN、Scanorama、SysVI、scVI 批次效应整合与评价"
    INPUT_REQUIRES = ['X_pca']

    def validate_input(self, adata):
        if 'X_pca' not in adata.obsm:
            return "PCA not found. Run dimensionality reduction first."
        batch_key = str(self.params.get('batch_key', 'batch') or '').strip()
        # Keep the historical behavior for a missing batch column (run() emits
        # the actionable error), but reject an existing continuous/high-cardinality
        # column before an integration method can treat each cell as a batch.
        if batch_key in adata.obs.columns:
            info = obs_grouping_info(
                adata, batch_key, max_categories=50,
                max_numeric_categories=20, require_multiple=True,
            )
            if not info['valid']:
                return f"batch_key '{batch_key}' 不是有效的分类批次列：{info['reason']}"
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

    @staticmethod
    def _numeric_metric(value):
        """Return a finite Python float for a summary value, or ``None``."""
        import math

        try:
            value = float(value)
        except (TypeError, ValueError):
            return None
        return value if math.isfinite(value) else None

    @staticmethod
    def _evaluation_sample_indices(labels, requested_size, seed=7):
        """Pick a deterministic, batch-stratified metric sample.

        A plain random sample can omit a rare batch, making a before/after ASW
        comparison both unstable and occasionally invalid.  This helper first
        includes one observation from every non-missing batch, then fills the
        remaining budget at random.  The same returned indices are used before
        and after correction.
        """
        import numpy as np
        import pandas as pd

        categorical = pd.Series(labels).astype('category')
        n_obs = int(len(categorical))
        try:
            requested = max(1, int(requested_size))
        except (TypeError, ValueError):
            requested = 10000
        if n_obs == 0:
            return np.array([], dtype=int), {
                'requested_size': requested,
                'used_size': 0,
                'strategy': 'batch_stratified',
                'sampled_all_cells': True,
                'n_batches_in_sample': 0,
            }

        codes = categorical.cat.codes.to_numpy()
        valid_codes = np.unique(codes[codes >= 0])
        # At least one cell per observed batch is required to make the batch
        # comparison interpretable.  This is still bounded by n_obs.
        target_size = min(n_obs, max(min(requested, n_obs), int(len(valid_codes))))
        rng = np.random.RandomState(seed)
        selected = []
        for code in valid_codes:
            candidates = np.flatnonzero(codes == code)
            if len(candidates):
                selected.append(int(rng.choice(candidates)))

        if len(selected) < target_size:
            selected_set = set(selected)
            remaining = np.array([idx for idx in range(n_obs) if idx not in selected_set], dtype=int)
            n_extra = min(target_size - len(selected), len(remaining))
            if n_extra:
                selected.extend(int(idx) for idx in rng.choice(remaining, size=n_extra, replace=False))

        indices = np.asarray(sorted(selected), dtype=int)
        metadata = {
            'requested_size': requested,
            'used_size': int(len(indices)),
            'strategy': 'batch_stratified',
            'sampled_all_cells': bool(len(indices) == n_obs),
            'n_batches_in_sample': int(len(np.unique(codes[indices][codes[indices] >= 0]))) if len(indices) else 0,
            'n_missing_batch_labels': int((codes < 0).sum()),
        }
        return indices, metadata

    @staticmethod
    def _silhouette_on_indices(representation, labels, indices):
        """Calculate ASW safely on an explicit, shared set of observations."""
        import numpy as np
        import pandas as pd
        from sklearn.metrics import silhouette_score

        try:
            representation = np.asarray(representation, dtype=float)
        except (TypeError, ValueError):
            return None, '表示矩阵不是可用于距离计算的数值矩阵。'
        if representation.ndim != 2 or representation.shape[0] == 0:
            return None, '表示矩阵为空或不是二维矩阵。'
        indices = np.asarray(indices, dtype=int)
        if len(indices) == 0:
            return None, '评价抽样为空。'

        categorical = pd.Series(labels).astype('category')
        codes = categorical.cat.codes.to_numpy()[indices]
        subset = representation[indices]
        valid = (codes >= 0) & np.isfinite(subset).all(axis=1)
        subset = subset[valid]
        codes = codes[valid]
        n_labels = int(len(np.unique(codes)))
        if len(subset) < 3 or n_labels < 2 or n_labels >= len(subset):
            return None, '抽样后每个类别需要足够的非缺失观测，当前无法稳定计算 ASW。'
        try:
            return round(float(silhouette_score(subset, codes)), 4), None
        except Exception as exc:
            return None, str(exc)

    def _bio_label_key(self, adata):
        """Find one categorical biological label suitable for retention checks."""
        candidates = [
            self.params.get('bio_label_key', ''),
            'celltype',
            'reference_celltype',
            'annotation',
            'leiden',
        ]
        seen = set()
        for raw_key in candidates:
            key = str(raw_key or '').strip()
            if not key or key in seen or key not in adata.obs.columns:
                continue
            seen.add(key)
            info = obs_grouping_info(
                adata, key, max_categories=100,
                max_numeric_categories=20, require_multiple=True,
            )
            if info.get('valid'):
                return key
        return ''

    def _append_neighbor_mixing_metrics(self, metrics, representation, batch_codes,
                                        sample_indices, connectivities=None):
        """Append bounded kNN/graph mixing metrics for one representation.

        When a graph is available (notably after BBKNN), only the sampled rows
        are inspected, while all of each row's graph neighbours remain in the
        calculation.  Otherwise a bounded kNN graph is built on the same sample.
        """
        import numpy as np

        sample_indices = np.asarray(sample_indices, dtype=int)
        valid_indices = sample_indices[batch_codes[sample_indices] >= 0] if len(sample_indices) else sample_indices
        if len(valid_indices) < 2:
            metrics['neighbor_mixing_warning'] = '可用于邻居混合评价的非缺失批次细胞不足 2 个。'
            return

        n_batches = max(1, int(len(np.unique(batch_codes[batch_codes >= 0]))))

        def _summarize(neighbor_rows, source, neighbor_k=None):
            entropies, same_batch_fracs = [], []
            for cell_index, neighbor_indices, weights in neighbor_rows:
                neighbor_indices = np.asarray(neighbor_indices, dtype=int)
                weights = np.asarray(weights, dtype=float)
                in_bounds = (neighbor_indices >= 0) & (neighbor_indices < len(batch_codes))
                neighbor_indices = neighbor_indices[in_bounds]
                weights = weights[in_bounds]
                valid_batch = batch_codes[neighbor_indices] >= 0
                neighbor_indices = neighbor_indices[valid_batch]
                weights = weights[valid_batch]
                if not len(neighbor_indices):
                    continue
                if weights.size != neighbor_indices.size:
                    weights = np.ones(len(neighbor_indices), dtype=float)
                total_weight = float(weights.sum())
                if total_weight <= 0:
                    continue
                counts = np.zeros(n_batches, dtype=float)
                np.add.at(counts, batch_codes[neighbor_indices], weights)
                own_batch = batch_codes[cell_index]
                if own_batch < 0 or own_batch >= len(counts):
                    continue
                same_batch_fracs.append(float(counts[own_batch] / total_weight))
                entropies.append(self._batch_entropy(counts))

            if not entropies:
                metrics['neighbor_mixing_warning'] = '评价图中没有可用的加权邻居。'
                return
            metrics['neighbor_mixing_source'] = source
            if neighbor_k is not None:
                metrics['evaluation_neighbor_k'] = int(neighbor_k)
            metrics['n_neighbor_evaluated_cells'] = int(len(entropies))
            metrics['mean_neighbor_batch_entropy'] = round(float(np.mean(entropies)), 4)
            metrics['mean_neighbor_same_batch_fraction'] = round(float(np.mean(same_batch_fracs)), 4)

        if connectivities is not None:
            try:
                csr = connectivities.tocsr()
                if csr.shape[0] != len(batch_codes) or csr.shape[1] != len(batch_codes):
                    raise ValueError('connectivities 与 obs 行数不一致')
                rows = []
                for cell_index in valid_indices:
                    start, end = csr.indptr[cell_index], csr.indptr[cell_index + 1]
                    neighbor_indices = csr.indices[start:end]
                    weights = csr.data[start:end]
                    # A self-loop does not carry mixing evidence.
                    keep = neighbor_indices != cell_index
                    rows.append((cell_index, neighbor_indices[keep], weights[keep]))
                _summarize(rows, 'connectivities')
                return
            except Exception as exc:
                metrics['neighbor_graph_warning'] = f'无法使用 connectivities，改用表示空间 kNN：{exc}'

        try:
            from sklearn.neighbors import NearestNeighbors

            subset = np.asarray(representation, dtype=float)[valid_indices]
            finite = np.isfinite(subset).all(axis=1)
            subset = subset[finite]
            source_indices = valid_indices[finite]
            if len(source_indices) < 2:
                metrics['neighbor_mixing_warning'] = '表示空间中有效细胞不足 2 个。'
                return
            neighbor_k = min(15, len(source_indices) - 1)
            model = NearestNeighbors(n_neighbors=neighbor_k + 1, metric='euclidean')
            neighbor_locals = model.fit(subset).kneighbors(subset, return_distance=False)
            rows = []
            for local_index, local_neighbors in enumerate(neighbor_locals):
                local_neighbors = np.asarray(local_neighbors, dtype=int)
                local_neighbors = local_neighbors[local_neighbors != local_index][:neighbor_k]
                rows.append((
                    int(source_indices[local_index]),
                    source_indices[local_neighbors],
                    np.ones(len(local_neighbors), dtype=float),
                ))
            _summarize(rows, 'representation_knn', neighbor_k=neighbor_k)
        except Exception as exc:
            metrics['neighbor_mixing_error'] = str(exc)

    def _append_cluster_batch_metrics(self, adata, metrics, batch_key):
        """Append design-level cluster/batch composition checks after correction."""
        import pandas as pd

        requested_cluster_key = str(self.params.get('evaluation_cluster_key', '') or '').strip()
        cluster_key, cluster_info = resolve_obs_grouping(
            adata,
            requested_cluster_key,
            fallbacks=['leiden', 'leiden_0.8', 'leiden_0.6', 'leiden_1.0'],
            max_categories=50,
            max_numeric_categories=20,
            require_multiple=False,
        )
        if cluster_key is None:
            cluster_key = 'batch_eval_leiden'
            if cluster_key not in adata.obs.columns:
                try:
                    import scanpy as sc

                    eval_resolution = float(self.params.get('evaluation_resolution', 0.8))
                    sc.tl.leiden(
                        adata,
                        resolution=eval_resolution,
                        key_added=cluster_key,
                        flavor='igraph',
                        n_iterations=2,
                    )
                    metrics['evaluation_resolution'] = eval_resolution
                except Exception as exc:
                    metrics['cluster_evaluation_warning'] = (
                        f'无法生成评价分群，已保留表示与邻居混合指标：{exc}'
                    )
                    return

        metrics['requested_evaluation_cluster_key'] = requested_cluster_key
        metrics['evaluation_cluster_key'] = cluster_key
        if requested_cluster_key and not cluster_info.get('requested_valid', False):
            metrics['evaluation_cluster_key_warning'] = cluster_info.get(
                'requested_reason', '请求列不是有效的分类聚类列'
            )

        table = pd.crosstab(adata.obs[cluster_key].astype(str), adata.obs[batch_key].astype(str))
        if not table.empty:
            cluster_sizes = table.sum(axis=1)
            fractions = table.div(cluster_sizes, axis=0)
            entropies = table.apply(lambda row: self._batch_entropy(row.values.astype(float)), axis=1)
            max_fracs = fractions.max(axis=1)
            weights = cluster_sizes / cluster_sizes.sum()
            metrics['n_eval_clusters'] = int(table.shape[0])
            metrics['weighted_batch_entropy'] = round(float((entropies * weights).sum()), 4)
            metrics['weighted_max_batch_fraction'] = round(float((max_fracs * weights).sum()), 4)
            metrics['batch_dominated_cluster_count_90pct'] = int((max_fracs >= 0.9).sum())
            metrics['single_batch_cluster_count'] = int((entropies == 0).sum())
            metrics['min_cluster_size'] = int(cluster_sizes.min())
            metrics['median_cluster_size'] = round(float(cluster_sizes.median()), 2)
            metrics['small_clusters_lt_50'] = int((cluster_sizes < 50).sum())

    def _compute_evaluation(self, adata, representation, batch_key, sample_indices=None,
                            sampling=None, connectivities=None,
                            include_cluster_metrics=False):
        """Evaluate a representation on a shared, bounded sample.

        ``representation`` intentionally accepts either an obsm key or a matrix
        so callers can preserve the original PCA before an integration method
        changes the working AnnData object.
        """
        import numpy as np

        if isinstance(representation, str):
            representation = adata.obsm[representation]
        representation = np.asarray(representation)
        batch_labels = adata.obs[batch_key].astype('category')
        if sample_indices is None:
            sample_indices, sampling = self._evaluation_sample_indices(
                batch_labels, self.params.get('evaluation_sample_size', 10000)
            )
        else:
            sample_indices = np.asarray(sample_indices, dtype=int)
            sampling = dict(sampling or {})
            sampling.setdefault('requested_size', int(len(sample_indices)))
            sampling.setdefault('used_size', int(len(sample_indices)))
            sampling.setdefault('strategy', 'provided_shared_indices')
            sampling.setdefault('sampled_all_cells', bool(len(sample_indices) == len(batch_labels)))

        batch_codes = batch_labels.cat.codes.to_numpy()
        metrics = {
            'n_batches': int(batch_labels.nunique()),
            'batch_counts': {
                str(k): int(v) for k, v in adata.obs[batch_key].astype(str).value_counts().to_dict().items()
            },
            'evaluation_sample_size_requested': int(sampling.get('requested_size', len(sample_indices))),
            'n_evaluation_cells': int(len(sample_indices)),
            'evaluation_sample_strategy': str(sampling.get('strategy', 'batch_stratified')),
            'sampled_all_cells': bool(sampling.get('sampled_all_cells', False)),
        }
        if sampling.get('n_missing_batch_labels'):
            metrics['n_missing_batch_labels'] = int(sampling['n_missing_batch_labels'])

        asw_batch, batch_error = self._silhouette_on_indices(
            representation, batch_labels, sample_indices
        )
        if asw_batch is None:
            metrics['asw_batch_warning'] = batch_error
        else:
            metrics['asw_batch'] = asw_batch
            metrics['abs_asw_batch'] = round(abs(asw_batch), 4)

        bio_key = self._bio_label_key(adata)
        if not bio_key:
            metrics['bio_label_warning'] = '未找到可用的生物标签列，无法自动评价生物结构保留。'
        else:
            metrics['bio_label_key'] = bio_key
            asw_bio, bio_error = self._silhouette_on_indices(
                representation, adata.obs[bio_key], sample_indices
            )
            if asw_bio is None:
                metrics['asw_bio_warning'] = bio_error
            else:
                metrics['asw_bio'] = asw_bio
                # Retain the historical key shape for consumers that look up a
                # concrete label name, while providing a stable generic key for
                # the pre/post delta table.
                metrics[f'asw_{bio_key}'] = asw_bio

        self._append_neighbor_mixing_metrics(
            metrics, representation, batch_codes, sample_indices,
            connectivities=connectivities,
        )
        if include_cluster_metrics:
            self._append_cluster_batch_metrics(adata, metrics, batch_key)
        return metrics

    @classmethod
    def _build_evaluation_comparison(cls, before, after, sampling, method):
        """Create JSON-friendly before/after deltas with explicit directions."""
        metric_directions = {
            'asw_batch': 'toward_zero',
            'abs_asw_batch': 'lower',
            'mean_neighbor_batch_entropy': 'higher',
            'mean_neighbor_same_batch_fraction': 'lower',
            'asw_bio': 'higher',
        }
        delta = {}
        for key in metric_directions:
            before_value = cls._numeric_metric(before.get(key))
            after_value = cls._numeric_metric(after.get(key))
            if before_value is not None and after_value is not None:
                delta[key] = round(after_value - before_value, 4)

        scope = 'graph_only' if str(method).lower() == 'bbknn' else 'embedding_and_graph'
        warnings = []
        if scope == 'graph_only':
            warnings.append('BBKNN 主要改变邻居图而非 PCA 表示；应优先根据邻居混合指标和前后 UMAP 复核。')
        if not delta:
            warnings.append('没有足够的成对数值指标用于前后差值比较。')
        return {
            'before': before,
            'after': after,
            'delta': delta,
            'metric_directions': metric_directions,
            'sampling': dict(sampling or {}),
            'comparison_scope': scope,
            'warnings': warnings,
        }

    @classmethod
    def _comparison_rows(cls, comparison):
        """Flatten comparable metrics for the downloadable review table."""
        labels = {
            'asw_batch': 'Batch ASW',
            'abs_asw_batch': '|Batch ASW|',
            'mean_neighbor_batch_entropy': 'Mean neighbour batch entropy',
            'mean_neighbor_same_batch_fraction': 'Mean same-batch neighbour fraction',
            'asw_bio': 'Biological-label ASW',
        }
        before = comparison.get('before') or {}
        after = comparison.get('after') or {}
        directions = comparison.get('metric_directions') or {}
        rows = []
        for key, direction in directions.items():
            before_value = cls._numeric_metric(before.get(key))
            after_value = cls._numeric_metric(after.get(key))
            if before_value is None or after_value is None:
                continue
            delta = round(after_value - before_value, 4)
            rows.append({
                'metric': key,
                'label': labels.get(key, key),
                'before': before_value,
                'after': after_value,
                'delta_after_minus_before': delta,
                'desired_direction': direction,
            })
        return rows

    def _make_eval_plots(self, adata, metrics, batch_key, method, plots_dir, comparison=None):
        import os
        import re
        import pandas as pd
        import matplotlib.pyplot as plt
        from modules.native_figures import bar_figure, grouped_bar_figure
        from modules.figure_style import NATURE_TEXT

        result_files = []
        comparison_rows = self._comparison_rows(comparison) if comparison else []
        if comparison_rows:
            comparison_scope = str(comparison.get('comparison_scope', 'embedding_and_graph'))
            sampling = comparison.get('sampling') or {}
            for row in comparison_rows:
                row['comparison_scope'] = comparison_scope
                row['n_evaluation_cells'] = sampling.get('used_size')
                row['sampling_strategy'] = sampling.get('strategy')

            fig = grouped_bar_figure(
                [row['label'] for row in comparison_rows],
                [
                    ('Before correction', [row['before'] for row in comparison_rows]),
                    ('After correction', [row['after'] for row in comparison_rows]),
                ],
                title=f'Batch Integration: Before vs After ({method})',
                x_label='Metric', y_label='Value', rotation=35,
            )
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, f'batch_evaluation_pre_post_{method}.png', 'metrics',
                f'Batch integration pre/post metrics ({method})', formats=('png', 'svg'), dpi=300,
            ))

            results_dir = os.path.join(self.project_dir, 'results')
            os.makedirs(results_dir, exist_ok=True)
            safe_method = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(method))
            comparison_path = os.path.join(
                results_dir, f'batch_evaluation_pre_post_{safe_method}.csv'
            )
            pd.DataFrame(comparison_rows).to_csv(comparison_path, index=False)
            result_files.append({
                'file_path': comparison_path,
                'file_type': 'csv',
                'category': 'table',
                'label': f'Batch integration pre/post metrics ({method})',
            })

        scalar_keys = [
            'abs_asw_batch',
            'asw_bio',
            'weighted_batch_entropy',
            'weighted_max_batch_fraction',
            'mean_neighbor_batch_entropy',
            'mean_neighbor_same_batch_fraction',
        ]
        x, y = [], []
        for key in scalar_keys:
            value = metrics.get(key)
            if isinstance(value, (int, float)):
                x.append(key)
                y.append(value)
        if x:
            fig = bar_figure(
                x, y, title=f'Batch Integration Metrics ({method})',
                x_label='Metric', y_label='Value', rotation=45,
            )
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, f'batch_eval_metrics_{method}.png', 'metrics',
                f'Batch integration metrics ({method})', formats=('png', 'svg'), dpi=300,
            ))

        cluster_key = metrics.get('evaluation_cluster_key')
        if cluster_key and cluster_key in adata.obs.columns:
            import pandas as pd
            table = pd.crosstab(
                adata.obs[cluster_key].astype(str),
                adata.obs[batch_key].astype(str),
                normalize='index',
            )
            fig = grouped_bar_figure(
                table.index.tolist(), [(str(batch), table[batch].values) for batch in table.columns],
                title=f'Evaluation Cluster Batch Composition ({method})',
                x_label=cluster_key, y_label='Fraction', rotation=35, stacked=True,
            )
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, f'batch_eval_cluster_composition_{method}.png', 'bar',
                f'Evaluation cluster batch composition ({method})', formats=('png', 'svg'), dpi=300,
            ))

        table_rows = []
        for key, value in metrics.items():
            if isinstance(value, dict):
                value = ', '.join([f'{k}: {v}' for k, v in value.items()])
            table_rows.append([key, str(value)])
        if table_rows:
            fig, ax = plt.subplots(figsize=(9.0, max(4.8, 0.28 * len(table_rows) + 1.8)), dpi=150)
            ax.axis('off')
            table_artist = ax.table(cellText=table_rows, colLabels=['Metric', 'Value'],
                                    loc='center', cellLoc='left', colLoc='left')
            table_artist.auto_set_font_size(False)
            table_artist.set_fontsize(8)
            table_artist.scale(1, 1.35)
            for (row, col), cell in table_artist.get_celld().items():
                cell.set_edgecolor('#D0D5DD')
                cell.set_linewidth(0.45)
                if row == 0:
                    cell.set_facecolor('#E8F0F8')
                    cell.set_text_props(weight='semibold', color=NATURE_TEXT)
                else:
                    cell.set_facecolor('white')
            ax.set_title(f'Batch Integration Metric Table ({method})', loc='left', pad=12,
                         fontsize=10, fontweight='semibold', color=NATURE_TEXT)
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, f'batch_eval_metric_table_{method}.png', 'table',
                f'Batch integration metric table ({method})', formats=('png', 'svg'), dpi=300,
            ))

        return result_files

    def run(self, input_path):
        import os
        import scanpy as sc
        from modules.native_figures import umap_figure

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        method = self.params.get('method', 'harmony')
        requested_batch_key = str(self.params.get('batch_key', 'batch') or '').strip()
        batch_key, batch_info = resolve_obs_grouping(
            adata, requested_batch_key, max_categories=50,
            max_numeric_categories=20, require_multiple=True,
        )
        if batch_key is None:
            if requested_batch_key not in adata.obs.columns:
                raise ValueError(f"batch_key '{requested_batch_key}' 不在 adata.obs 中，可用列: {list(adata.obs.columns)}")
            raise ValueError(f"batch_key '{requested_batch_key}' 不是有效的分类批次列：{batch_info.get('reason', '')}")
        if not batch_info.get('requested_valid', False):
            raise ValueError(f"batch_key '{requested_batch_key}' 不是有效的分类批次列：{batch_info.get('requested_reason', '')}")
        if batch_key not in adata.obs.columns:
            raise ValueError(f"batch_key '{batch_key}' 不在 adata.obs 中，可用列: {list(adata.obs.columns)}")
        if not hasattr(adata.obs[batch_key].dtype, 'categories'):
            adata.obs[batch_key] = adata.obs[batch_key].astype(str).astype('category')
        if adata.obs[batch_key].nunique() < 2:
            raise ValueError(f"batch_key '{batch_key}' 只有一个取值，无法评估或执行批次整合")

        n_pcs = int(self.params.get('n_pcs', 50))
        if 'X_pca' in adata.obsm:
            n_pcs = max(2, min(n_pcs, adata.obsm['X_pca'].shape[1]))
        umap_before = adata.obsm['X_umap'].copy() if 'X_umap' in adata.obsm else None
        method_info = {}
        graph_already_built = False
        # Evaluate the unmodified representation before an integration method
        # can replace a graph or embedding.  The selected cells are shared with
        # the post-correction evaluation, so deltas are directly comparable.
        evaluate = self._as_bool(self.params.get('evaluate_correction', True))
        pre_eval_metrics = {}
        evaluation_sampling = None
        evaluation_indices = None
        if evaluate:
            self.progress(15, "Profiling batch structure before integration...")
            evaluation_indices, evaluation_sampling = self._evaluation_sample_indices(
                adata.obs[batch_key], self.params.get('evaluation_sample_size', 10000)
            )
            try:
                pre_connectivities = (
                    adata.obsp['connectivities'] if 'connectivities' in adata.obsp else None
                )
                pre_eval_metrics = self._compute_evaluation(
                    adata,
                    adata.obsm['X_pca'],
                    batch_key,
                    sample_indices=evaluation_indices,
                    sampling=evaluation_sampling,
                    connectivities=pre_connectivities,
                    include_cluster_metrics=False,
                )
            except Exception as exc:
                pre_eval_metrics = {'error': str(exc)}

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

        fig_after = umap_figure(
            adata, batch_key,
            title=f'UMAP after {method} integration (by {batch_key})',
        )
        result_files.extend(self.save_matplotlib_figure(
            fig_after, plots_dir, f'batch_umap_{method}.png', 'umap',
            f'UMAP after {method} (batch)', formats=('png', 'svg'), dpi=300,
        ))

        if umap_before is not None:
            corrected_umap = adata.obsm['X_umap'].copy()
            adata.obsm['X_umap'] = umap_before
            fig_before = umap_figure(
                adata, batch_key, title=f'UMAP before integration (by {batch_key})'
            )
            adata.obsm['X_umap'] = corrected_umap
            result_files.extend(self.save_matplotlib_figure(
                fig_before, plots_dir, 'batch_umap_before.png', 'umap',
                'UMAP before integration (batch)', formats=('png', 'svg'), dpi=300,
            ))

        eval_metrics = {}
        evaluation_comparison = None
        if evaluate:
            self.progress(80, "Evaluating batch integration...")
            try:
                post_connectivities = (
                    adata.obsp['connectivities'] if 'connectivities' in adata.obsp else None
                )
                eval_metrics = self._compute_evaluation(
                    adata,
                    adata.obsm[corrected_key],
                    batch_key,
                    sample_indices=evaluation_indices,
                    sampling=evaluation_sampling,
                    connectivities=post_connectivities,
                    include_cluster_metrics=True,
                )
                evaluation_comparison = self._build_evaluation_comparison(
                    pre_eval_metrics, eval_metrics, evaluation_sampling, method
                )
                result_files.extend(self._make_eval_plots(
                    adata, eval_metrics, batch_key, method, plots_dir,
                    comparison=evaluation_comparison,
                ))
            except Exception as exc:
                eval_metrics = {'error': str(exc)}
                if pre_eval_metrics:
                    evaluation_comparison = self._build_evaluation_comparison(
                        pre_eval_metrics, eval_metrics, evaluation_sampling, method
                    )

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'batch_correct')

        summary = {
            'method': method,
            'n_cells': int(adata.n_obs),
            'n_genes': int(adata.n_vars),
            'embedding_key': corrected_key,
            'requested_batch_key': requested_batch_key,
            'batch_key': batch_key,
            'method_info': method_info,
        }
        if eval_metrics:
            summary['eval_metrics'] = eval_metrics
        if evaluation_comparison:
            summary['evaluation_comparison'] = evaluation_comparison

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }
