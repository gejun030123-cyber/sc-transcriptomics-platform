from modules.base import BaseAnalysis
import copy
import logging
import math
import re

from modules.pc_strategy import resolution_filename_token, resolve_analysis_n_pcs

logger = logging.getLogger(__name__)


def parse_resolutions(value, default='0.6,0.8,1.0'):
    """Parse and validate clustering resolutions from UI or API input."""
    import numpy as np

    raw_value = default if value is None else value
    tokens = list(raw_value) if isinstance(raw_value, (list, tuple)) else str(raw_value).split(',')
    if not tokens or any(not str(token).strip() for token in tokens):
        raise ValueError('resolutions 不能为空，且必须是正数列表')

    parsed = []
    for token in tokens:
        try:
            resolution = float(str(token).strip())
        except (TypeError, ValueError) as exc:
            raise ValueError(f'无效的聚类分辨率：{token!r}') from exc
        if not np.isfinite(resolution) or resolution <= 0:
            raise ValueError(f'聚类分辨率必须是有限正数：{token!r}')
        if not any(math.isclose(resolution, existing, rel_tol=0.0, abs_tol=1e-12)
                   for existing in parsed):
            parsed.append(resolution)
    if not parsed:
        raise ValueError('至少需要一个聚类分辨率')
    return parsed


def cluster_key_for_method(method, resolution):
    """Return a method-specific per-resolution AnnData observation key."""
    prefix = 'louvain' if str(method).lower() == 'louvain' else 'leiden'
    return f'{prefix}_{resolution}'


def limited_representation(rep_data, n_pcs):
    """Limit a representation to the dimensions used for KNN/metrics."""
    import numpy as np

    array = np.asarray(rep_data)
    if array.ndim != 2:
        raise ValueError('聚类表示必须是二维矩阵')
    if n_pcs is None or int(n_pcs) <= 0:
        return array
    return array[:, :min(int(n_pcs), array.shape[1])]


def metric_subset(rep_data, labels, max_cells=10000, random_state=0):
    """Deterministically subsample representation/labels for resolution metrics."""
    import numpy as np

    array = np.asarray(rep_data)
    labels = np.asarray(labels)
    max_cells = int(max_cells or 0)
    # 非正数视为“未设置”，钳制到默认值而不是全量评分：
    # silhouette/calinski 全量是 O(n²)，10 万细胞会耗尽内存。
    if max_cells <= 0:
        max_cells = 10000
    if len(labels) <= max_cells:
        return array, labels

    rng = np.random.default_rng(int(random_state))
    unique = np.unique(labels)
    mandatory = []
    for label in unique:
        indices = np.flatnonzero(labels == label)
        mandatory.append(int(rng.choice(indices)))
    mandatory = np.asarray(mandatory, dtype=int)
    if len(mandatory) >= max_cells:
        selected = np.sort(mandatory[:max_cells])
    else:
        remaining = np.setdiff1d(np.arange(len(labels)), mandatory, assume_unique=False)
        extra = rng.choice(remaining, size=max_cells - len(mandatory), replace=False)
        selected = np.sort(np.concatenate([mandatory, extra]))
    return array[selected], labels[selected]

class ClusteringAnalysis(BaseAnalysis):
    MODULE_NAME = "clustering"
    DISPLAY_NAME = "聚类分析"
    DESCRIPTION = "KNN 图 + 多分辨率 Leiden 聚类"
    # These are representations produced by the batch-correction module.
    # Keep the order deterministic so Harmony is preferred when it exists,
    # while still supporting the other correction methods exposed by the UI.
    CORRECTED_REPRESENTATIONS = (
        'X_pca_harmony', 'X_pca_combat', 'X_scanorama', 'X_sysvi', 'X_scVI',
    )
    # UMAP is rebuilt from the representation selected below.  Requiring an
    # old X_umap here would allow a pre-Harmony embedding to remain attached
    # to a corrected clustering result.
    INPUT_REQUIRES = []

    @staticmethod
    def _as_bool(value, default=False):
        if value is None:
            return default
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            return value.strip().lower() in {'1', 'true', 'yes', 'on'}
        return bool(value)

    def _resolve_representation(self, adata, use_corrected):
        """Resolve the embedding used for both neighbors and clustering.

        ``use_corrected`` is an explicit request, not a preference hint.  If
        no corrected embedding is available, fail loudly instead of silently
        returning to ``X_pca`` and labelling the result as post-integration.
        """
        if use_corrected:
            for key in self.CORRECTED_REPRESENTATIONS:
                if key in adata.obsm:
                    return key, False
            bbknn_info = (adata.uns or {}).get('bbknn_corrected_graph')
            if (
                isinstance(bbknn_info, dict)
                and bbknn_info.get('method') == 'bbknn'
                and ('connectivities' in adata.obsp or 'neighbors' in adata.uns)
            ):
                # BBKNN 的整合结果保存在邻居图而不是 obsm；复用该图。
                return 'X_pca', True
            available = ', '.join(sorted(str(key) for key in adata.obsm.keys())) or '无'
            raise ValueError(
                'use_corrected=true 但输入数据没有校正后表示；'
                '不会静默回退到 X_pca。请先运行批次校正，或明确关闭 use_corrected。'
                f' 当前 obsm: {available}'
            )
        if 'X_pca' not in adata.obsm:
            raise ValueError(
                'No valid clustering representation: use_corrected=false '
                '需要输入数据包含 X_pca；请先运行降维。'
            )
        return 'X_pca', False

    def validate_input(self, adata):
        # A missing parameter means the programmatic API requested the
        # historical PCA-only path.  The web schema supplies True explicitly,
        # in which case the strict corrected-representation contract applies.
        use_corrected = self._as_bool(self.params.get('use_corrected', False), False)
        try:
            self._resolve_representation(adata, use_corrected)
        except ValueError as exc:
            return str(exc)
        return None

    def _build_knn(self, adata, use_rep, n_neighbors, n_pcs, distance_metric):
        """Build the graph on a bounded representation and return used dims."""
        import numpy as np
        representation = np.asarray(adata.obsm[use_rep])
        used_n_pcs = min(int(n_pcs), representation.shape[1])
        if used_n_pcs < int(n_pcs):
            self.progress(
                -1,
                f'请求 {n_pcs} 个聚类维度超过 {use_rep} 的实际维度，已使用 {used_n_pcs}',
            )

        import scanpy as sc
        # 始终显式传入 use_rep=X_pca：scanpy 在 n_vars <= 50 且 use_rep=None
        # 时会退回到 adata.X（tools/_utils.py _get_pca_or_small_x），导致 KNN
        # 建在表达矩阵而不是 PCA 坐标上，与模块契约和记录元数据不符。
        sc.pp.neighbors(
            adata, n_neighbors=n_neighbors, use_rep=use_rep,
            n_pcs=used_n_pcs, metric=distance_metric,
        )
        return used_n_pcs

    def _run_clusterings(self, adata, resolutions, clustering_method, n_iterations):
        import scanpy as sc

        for index, resolution in enumerate(resolutions):
            key = cluster_key_for_method(clustering_method, resolution)
            if clustering_method == 'louvain':
                sc.tl.louvain(adata, resolution=resolution, key_added=key)
            else:
                sc.tl.leiden(
                    adata, resolution=resolution, key_added=key,
                    flavor='igraph', n_iterations=n_iterations,
                )
            pct = 40 + int((index + 1) / len(resolutions) * 30)
            self.progress(pct, f"{clustering_method} resolution {resolution} done")

    @staticmethod
    def _as_dense_array(matrix):
        """Convert only the selected expression slice to a dense float array."""
        import numpy as np

        if hasattr(matrix, 'toarray'):
            matrix = matrix.toarray()
        return np.asarray(matrix, dtype=float)

    @staticmethod
    def _safe_cluster_filename_token(value):
        """Keep data-derived cluster labels inside the managed plots folder."""
        token = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(value or '')).strip('._')
        return token[:48] or 'cluster'

    @staticmethod
    def _sample_indices(indices, maximum, seed):
        """Return a deterministic, bounded cell sample for static violins."""
        import numpy as np

        indices = np.asarray(indices, dtype=int)
        maximum = max(1, int(maximum))
        if len(indices) <= maximum:
            return indices
        rng = np.random.default_rng(int(seed))
        return np.sort(rng.choice(indices, size=maximum, replace=False))

    @staticmethod
    def _marker_genes_by_cluster(adata, marker_selection, cluster_groups, top_n):
        """Keep selected markers ordered by their owning cluster and present in var."""
        available = {str(gene) for gene in adata.var_names}
        marker_map = (marker_selection or {}).get('cluster_markers', {}) or {}
        selected = {}
        for group in cluster_groups:
            genes = []
            for gene in marker_map.get(str(group), [])[:max(1, int(top_n))]:
                gene = str(gene)
                if gene in available and gene not in genes:
                    genes.append(gene)
            if genes:
                selected[str(group)] = genes
        return selected

    @classmethod
    def _marker_expression_source(cls, adata):
        """Prefer recoverable counts for display, without changing the analysis matrix."""
        import numpy as np

        if 'counts' in adata.layers:
            matrix = adata.layers['counts']
            library_size = np.asarray(matrix.sum(axis=1), dtype=float).reshape(-1)
            return {
                'matrix': matrix,
                'gene_index': {str(gene): index for index, gene in enumerate(adata.var_names)},
                'normalize_counts': True,
                'library_size': library_size,
                'label': 'Library-size normalized log1p expression',
                'source': 'layers["counts"] → library-size normalized log1p',
            }
        if adata.raw is not None:
            return {
                'matrix': adata.raw.X,
                'gene_index': {str(gene): index for index, gene in enumerate(adata.raw.var_names)},
                'normalize_counts': False,
                'library_size': None,
                'label': 'Expression (adata.raw display scale)',
                'source': 'adata.raw.X',
            }
        return {
            'matrix': adata.X,
            'gene_index': {str(gene): index for index, gene in enumerate(adata.var_names)},
            'normalize_counts': False,
            'library_size': None,
            'label': 'Expression (adata.X display scale)',
            'source': 'adata.X',
        }

    @classmethod
    def _marker_expression_values(cls, source, indices, genes):
        """Extract a small marker slice, normalizing UMI counts only for display."""
        import numpy as np

        gene_indices = [source['gene_index'].get(str(gene)) for gene in genes]
        if any(index is None for index in gene_indices):
            missing = [str(gene) for gene, index in zip(genes, gene_indices) if index is None]
            raise KeyError(f'Marker gene not found in display matrix: {missing}')
        values = cls._as_dense_array(source['matrix'][indices][:, gene_indices])
        if values.ndim == 1:
            values = values.reshape(-1, len(gene_indices))
        if source['normalize_counts']:
            library_size = np.asarray(source['library_size'], dtype=float)[indices]
            values = np.log1p(np.divide(
                values, library_size[:, None], out=np.zeros_like(values),
                where=library_size[:, None] > 0,
            ) * 10_000.0)
        # An already scaled X/raw matrix may contain negative values.  A
        # violin of those values is still useful as a relative display, but
        # clipping avoids presenting a negative count-like expression scale.
        return np.nan_to_num(np.clip(values, 0.0, None), nan=0.0, posinf=0.0, neginf=0.0)

    def _render_marker_heatmap(self, adata, cluster_key, marker_selection,
                               plots_dir, result_files):
        """Render per-cluster marker means as a gene-wise z-score heatmap."""
        import numpy as np
        from modules.native_figures import heatmap_figure

        if not self._as_bool(self.params.get('show_marker_heatmap', True), True):
            return False, 0
        labels = adata.obs[cluster_key].astype(str)
        groups = sorted(labels.unique().tolist(), key=lambda item: (len(item), item))
        try:
            top_n = max(1, int(self.params.get('marker_heatmap_top_n', 3)))
        except (TypeError, ValueError):
            top_n = 3
        genes_by_cluster = self._marker_genes_by_cluster(
            adata, marker_selection, groups, top_n,
        )
        genes = [
            gene for group in groups for gene in genes_by_cluster.get(str(group), [])
        ]
        genes = list(dict.fromkeys(genes))
        if not genes:
            return False, 0

        expression = self._as_dense_array(adata[:, genes].X)
        means = np.vstack([
            expression[labels.to_numpy() == str(group)].mean(axis=0)
            for group in groups
        ])
        z_scores = (means - means.mean(axis=0, keepdims=True)) / (
            means.std(axis=0, keepdims=True) + 1e-10
        )
        figure = heatmap_figure(
            np.clip(z_scores.T, -3, 3),
            x_labels=groups,
            y_labels=genes,
            title='Cluster Marker Heatmap',
            x_label='Cluster',
            y_label='Cluster-specific marker gene',
            colorbar_label='Mean expression z-score',
        )
        result_files.extend(self.save_matplotlib_figure(
            figure, plots_dir, 'cluster_marker_heatmap.png',
            'heatmap', 'Cluster Marker Heatmap', formats=('png', 'svg'), dpi=300,
        ))
        return True, len(genes)

    def _render_cluster_marker_violins(self, adata, cluster_key, marker_selection,
                                       plots_dir, result_files):
        """Generate one target-cluster-versus-rest marker violin panel per cluster."""
        import numpy as np
        from modules.native_figures import marker_violin_figure

        if not self._as_bool(self.params.get('show_marker_violin', True), True):
            return [], None
        labels = adata.obs[cluster_key].astype(str)
        groups = sorted(labels.unique().tolist(), key=lambda item: (len(item), item))
        if len(groups) < 2:
            return [], None
        try:
            top_n = max(1, int(self.params.get('marker_violin_top_n', 3)))
        except (TypeError, ValueError):
            top_n = 3
        try:
            max_cells = max(50, int(self.params.get('marker_violin_max_cells_per_group', 2000)))
        except (TypeError, ValueError):
            max_cells = 2000
        genes_by_cluster = self._marker_genes_by_cluster(
            adata, marker_selection, groups, top_n,
        )
        if not genes_by_cluster:
            return [], None
        source = self._marker_expression_source(adata)
        label_values = labels.to_numpy()
        completed = []
        for group_index, group in enumerate(groups):
            genes = genes_by_cluster.get(str(group), [])
            if not genes:
                continue
            target_indices = self._sample_indices(
                np.flatnonzero(label_values == str(group)), max_cells,
                seed=group_index,
            )
            rest_indices = self._sample_indices(
                np.flatnonzero(label_values != str(group)), max_cells,
                seed=10_000 + group_index,
            )
            if not len(target_indices) or not len(rest_indices):
                continue
            try:
                target_values = self._marker_expression_values(source, target_indices, genes)
                rest_values = self._marker_expression_values(source, rest_indices, genes)
                figure = marker_violin_figure(
                    [
                        [target_values[:, index] for index in range(len(genes))],
                        [rest_values[:, index] for index in range(len(genes))],
                    ],
                    [f'Cluster {group}', 'Other clusters'],
                    genes,
                    title=f'Cluster {group} Differential Marker Expression',
                    y_label=source['label'],
                )
                result_files.extend(self.save_matplotlib_figure(
                    figure,
                    plots_dir,
                    f'cluster_marker_violin_{group_index + 1:02d}_{self._safe_cluster_filename_token(group)}.png',
                    'violin',
                    f'Cluster {group} Differential Marker Violin (vs other clusters)',
                    formats=('png', 'svg'), dpi=300,
                ))
                completed.append(str(group))
            except Exception as exc:
                logger.warning('生成 cluster %s 差异 marker violin 失败: %s', group, exc)
                self.progress(-1, f'Cluster {group} 差异 marker violin 生成失败：{exc}')
        return completed, source['source']

    def _run_marker_preview(self, adata, cluster_key, plots_dir, result_files):
        """Render dotplot, heatmap, and cluster-versus-rest marker violins."""
        visualization_defaults = {
            'dotplot': False, 'heatmap': False, 'heatmap_n_genes': 0,
            'violin_clusters': [], 'violin_expression_source': None,
        }
        marker_selection = {
            'cluster_markers': {}, 'genes': [], 'warnings': [],
            'n_data_driven': 0, 'n_classic_anchor': 0, 'n_relaxed': 0,
            'parameters': {},
            'visualizations': visualization_defaults.copy(),
        }
        if not self._as_bool(self.params.get('compute_marker_preview', True), True):
            marker_selection['warnings'].append(
                'Marker preview disabled for proportions/expression-only batch run.'
            )
            return marker_selection

        try:
            from modules.annotation import DEFAULT_UNIVERSAL_MARKERS, select_cluster_marker_genes
            if cluster_key not in adata.obs.columns:
                return marker_selection
            marker_selection = select_cluster_marker_genes(
                adata,
                cluster_key,
                classic_markers=DEFAULT_UNIVERSAL_MARKERS,
                method=self.params.get('marker_selection_method', 'wilcoxon'),
                n_rank_genes=int(self.params.get('marker_rank_genes', 200)),
                min_markers=int(self.params.get('marker_min_per_cluster', 2)),
                max_markers=int(self.params.get('marker_max_per_cluster', 5)),
                padj_cutoff=float(self.params.get('marker_padj_cutoff', 0.05)),
                min_pct=float(self.params.get('marker_min_pct', 0.10)),
                min_delta_pct=float(self.params.get('marker_min_delta_pct', 0.05)),
            )
        except Exception as exc:
            logger.warning("选择 cluster marker 失败（注释模块可能不可用）: %s", exc)
            marker_selection['warnings'].append(f'cluster marker selection failed: {exc}')
            return marker_selection

        marker_selection.setdefault('visualizations', visualization_defaults.copy())
        if marker_selection.get('warnings'):
            self.progress(-1, '；'.join(marker_selection['warnings'][:2]))
        dotplot_genes = marker_selection.get('genes', [])[:40]
        if dotplot_genes:
            try:
                import matplotlib.pyplot as plt
                import scanpy as sc
                sc.tl.dendrogram(adata, groupby=cluster_key)
                fig_dot = sc.pl.dotplot(
                    adata, var_names=dotplot_genes, groupby=cluster_key,
                    return_fig=True,
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_dot, plots_dir, 'cluster_marker_dotplot.png', 'dotplot',
                    'Cluster-specific Marker Dotplot'
                ))
                plt.close(getattr(fig_dot, 'fig', fig_dot))
                marker_selection['visualizations']['dotplot'] = True
            except Exception as exc:
                logger.warning('生成 marker dotplot 失败: %s', exc)
                self.progress(-1, f'Marker dotplot 生成失败：{exc}')

        try:
            heatmap_drawn, heatmap_n_genes = self._render_marker_heatmap(
                adata, cluster_key, marker_selection, plots_dir, result_files,
            )
            marker_selection['visualizations']['heatmap'] = bool(heatmap_drawn)
            marker_selection['visualizations']['heatmap_n_genes'] = int(heatmap_n_genes)
        except Exception as exc:
            logger.warning('生成 cluster marker heatmap 失败: %s', exc)
            self.progress(-1, f'Cluster marker 热图生成失败：{exc}')

        try:
            violin_clusters, violin_source = self._render_cluster_marker_violins(
                adata, cluster_key, marker_selection, plots_dir, result_files,
            )
            marker_selection['visualizations']['violin_clusters'] = violin_clusters
            marker_selection['visualizations']['violin_expression_source'] = violin_source
        except Exception as exc:
            logger.warning('生成 cluster marker violin 失败: %s', exc)
            self.progress(-1, f'Cluster marker violin 生成失败：{exc}')
        return marker_selection

    def _select_best_resolution(self, adata, resolutions, clustering_method,
                                use_rep, used_n_pcs, auto_select,
                                resolution_metric, metric_subsample):
        """Select a resolution with bounded, reproducible quality metrics."""
        import numpy as np

        best_res = resolutions[0]
        resolution_scores = {}
        if not auto_select or len(resolutions) <= 1:
            return best_res, resolution_scores
        try:
            from sklearn.metrics import (
                calinski_harabasz_score,
                davies_bouldin_score,
                silhouette_score,
            )
            best_score = -float('inf')
            metric_rep = limited_representation(adata.obsm[use_rep], used_n_pcs)
            for resolution in resolutions:
                key = cluster_key_for_method(clustering_method, resolution)
                if key not in adata.obs.columns:
                    continue
                labels = adata.obs[key].astype('category').cat.codes.values
                rep_data, metric_labels = metric_subset(
                    metric_rep, labels, max_cells=metric_subsample,
                )
                n_labels = len(np.unique(metric_labels))
                if n_labels < 2 or len(metric_labels) <= n_labels:
                    continue
                if resolution_metric == 'calinski':
                    score = calinski_harabasz_score(rep_data, metric_labels)
                elif resolution_metric == 'davies_bouldin':
                    score = -davies_bouldin_score(rep_data, metric_labels)
                else:
                    score = silhouette_score(rep_data, metric_labels)
                resolution_scores[str(resolution)] = float(score)
                if score > best_score:
                    best_score = score
                    best_res = resolution
        except Exception as exc:
            logger.warning("自动选择分辨率失败: %s", exc)
        return best_res, resolution_scores

    def _plot_resolution_sankey(self, adata, resolutions, clustering_method,
                                plots_dir, result_files):
        """Render proportional, adjacent-resolution cluster flows.

        Cluster identifiers only identify a partition within one resolution;
        they are not persistent identities across resolutions.  The alluvial
        bands therefore encode the actual cells shared by adjacent partitions,
        instead of implying correspondence from matching cluster IDs.
        """
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        from matplotlib.patches import PathPatch, Rectangle
        from matplotlib.path import Path
        from modules.figure_style import NATURE_MUTED, NATURE_PALETTE, NATURE_TEXT, stable_category_colors

        if not self.params.get('show_resolution_sankey', True) or len(resolutions) <= 1:
            return
        sankey_res = sorted([
            resolution for resolution in resolutions
            if cluster_key_for_method(clustering_method, resolution) in adata.obs.columns
        ])
        if len(sankey_res) <= 1:
            return

        def _categories_and_counts(resolution):
            key = cluster_key_for_method(clustering_method, resolution)
            labels = adata.obs[key].astype(str)
            categories = sorted(labels.unique(), key=lambda item: (len(item), item))
            counts = labels.value_counts()
            return categories, {category: int(counts.get(category, 0)) for category in categories}

        categories_by_resolution = {}
        counts_by_resolution = {}
        for resolution in sankey_res:
            categories, counts = _categories_and_counts(resolution)
            categories_by_resolution[resolution] = categories
            counts_by_resolution[resolution] = counts

        flows = []
        for left, right in zip(sankey_res[:-1], sankey_res[1:]):
            left_key = cluster_key_for_method(clustering_method, left)
            right_key = cluster_key_for_method(clustering_method, right)
            flow = pd.crosstab(
                adata.obs[left_key].astype(str), adata.obs[right_key].astype(str),
            )
            for left_category in flow.index:
                for right_category in flow.columns:
                    count = int(flow.loc[left_category, right_category])
                    if count > 0:
                        flows.append((left, right, str(left_category), str(right_category), count))
        if not flows:
            return

        # Each column spans the same vertical scale.  A small fixed gap helps
        # distinguish adjacent clusters without distorting their cell counts.
        total_cells = max(int(adata.n_obs), 1)
        max_clusters = max(len(categories) for categories in categories_by_resolution.values())
        fig_sankey, ax_sankey = plt.subplots(
            figsize=(max(8.0, 2.8 * len(sankey_res)),
                     max(5.6, min(18.0, 2.8 + 0.38 * max_clusters))), dpi=150,
        )
        x_positions = np.linspace(0.10, 0.90, len(sankey_res))
        node_width = min(0.065, 0.18 / max(len(sankey_res), 1))
        node_bounds = {}
        for res_index, resolution in enumerate(sankey_res):
            categories = categories_by_resolution[resolution]
            gap = min(0.014, 0.12 / max(len(categories) - 1, 1)) if len(categories) > 1 else 0.0
            available_height = 0.86 - gap * max(len(categories) - 1, 0)
            cursor = 0.93
            for category in categories:
                height = available_height * counts_by_resolution[resolution][category] / total_cells
                node_bounds[(resolution, category)] = (cursor - height, cursor)
                cursor -= height + gap
            ax_sankey.text(
                x_positions[res_index], 0.99, f'Resolution {resolution}',
                ha='center', va='bottom', fontsize=9,
                color=NATURE_TEXT, fontweight='semibold',
            )

        resolution_index = {resolution: index for index, resolution in enumerate(sankey_res)}
        node_color_keys = [
            f'{resolution}:{category}'
            for resolution in sankey_res
            for category in categories_by_resolution[resolution]
        ]
        node_colors = stable_category_colors(node_color_keys, palette=NATURE_PALETTE)

        # Allocate the source and target portions of each band independently.
        # This conserves the full height of every node even when the number of
        # clusters (and therefore column gap) differs between resolutions.
        source_segments = {}
        target_segments = {}
        for left, right in zip(sankey_res[:-1], sankey_res[1:]):
            pair_flows = {
                (left_category, right_category): count
                for flow_left, flow_right, left_category, right_category, count in flows
                if flow_left == left and flow_right == right
            }
            for left_category in categories_by_resolution[left]:
                low, high = node_bounds[(left, left_category)]
                cursor = high
                for right_category in categories_by_resolution[right]:
                    count = pair_flows.get((left_category, right_category), 0)
                    if not count:
                        continue
                    height = (high - low) * count / counts_by_resolution[left][left_category]
                    source_segments[(left, right, left_category, right_category)] = (cursor - height, cursor)
                    cursor -= height
            for right_category in categories_by_resolution[right]:
                low, high = node_bounds[(right, right_category)]
                cursor = high
                for left_category in categories_by_resolution[left]:
                    count = pair_flows.get((left_category, right_category), 0)
                    if not count:
                        continue
                    height = (high - low) * count / counts_by_resolution[right][right_category]
                    target_segments[(left, right, left_category, right_category)] = (cursor - height, cursor)
                    cursor -= height

        # Draw true bands before nodes.  Cubic edges make each transition easy
        # to follow while the band width remains proportional to cell count.
        for left, right, left_category, right_category, _count in flows:
            source_low, source_high = source_segments[(left, right, left_category, right_category)]
            target_low, target_high = target_segments[(left, right, left_category, right_category)]
            left_x = x_positions[resolution_index[left]] + node_width / 2
            right_x = x_positions[resolution_index[right]] - node_width / 2
            bend = (right_x - left_x) * 0.42
            path = Path(
                [
                    (left_x, source_low),
                    (left_x + bend, source_low),
                    (right_x - bend, target_low),
                    (right_x, target_low),
                    (right_x, target_high),
                    (right_x - bend, target_high),
                    (left_x + bend, source_high),
                    (left_x, source_high),
                    (left_x, source_low),
                ],
                [
                    Path.MOVETO,
                    Path.CURVE4, Path.CURVE4, Path.CURVE4,
                    Path.LINETO,
                    Path.CURVE4, Path.CURVE4, Path.CURVE4,
                    Path.CLOSEPOLY,
                ],
            )
            color = node_colors[f'{left}:{left_category}']
            ax_sankey.add_patch(PathPatch(path, facecolor=color, edgecolor='none', alpha=0.38, zorder=1))

        for resolution in sankey_res:
            for category in categories_by_resolution[resolution]:
                low, high = node_bounds[(resolution, category)]
                height = high - low
                x_value = x_positions[resolution_index[resolution]]
                color = node_colors[f'{resolution}:{category}']
                ax_sankey.add_patch(Rectangle(
                    (x_value - node_width / 2, low), node_width, height,
                    facecolor=color, edgecolor='white', linewidth=0.7, zorder=3,
                ))
                if height >= 0.042:
                    label = f'{category}\n{counts_by_resolution[resolution][category]:,}'
                    fontsize = 6.8
                elif height >= 0.020:
                    label = category
                    fontsize = 5.8
                else:
                    continue
                ax_sankey.text(
                    x_value, (low + high) / 2, label,
                    ha='center', va='center', fontsize=fontsize, color='white',
                    linespacing=0.92, zorder=4,
                )
        ax_sankey.set_xlim(0, 1)
        ax_sankey.set_ylim(0, 1.04)
        ax_sankey.axis('off')
        ax_sankey.set_title(
            'Cluster Resolution Flow', loc='left', pad=12,
            fontsize=10, fontweight='semibold', color=NATURE_TEXT,
        )
        fig_sankey.text(
            0.5, 0.018,
            'Band width = shared cells between adjacent resolutions; cluster IDs are resolution-specific.',
            ha='center', va='bottom', fontsize=7.5, color=NATURE_MUTED,
        )
        layout_rect = (0, 0.08, 1, 0.93)
        fig_sankey.tight_layout(pad=1.1, rect=layout_rect)
        fig_sankey._native_layout_rect = layout_rect
        result_files.extend(self.save_matplotlib_figure(
            fig_sankey, plots_dir, 'cluster_resolution_sankey.png',
            'sankey', '分辨率分群流向图', formats=('png', 'svg'), dpi=300,
        ))

    def run(self, input_path):
        import os
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import scanpy as sc
        from modules.figure_style import NATURE_PALETTE
        from modules.native_figures import (
            bar_figure,
            grouped_bar_figure,
            umap_figure,
            umap_panel_figure,
        )

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        resolutions = parse_resolutions(self.params.get('resolutions', '0.6,0.8,1.0'))
        n_neighbors = int(self.params.get('n_neighbors', 15))
        clustering_method = str(self.params.get('clustering_method', 'leiden') or 'leiden').lower()
        if clustering_method not in {'leiden', 'louvain'}:
            raise ValueError(f'不支持的聚类算法：{clustering_method}')
        n_iterations = int(self.params.get('n_iterations', 2))
        distance_metric = self.params.get('distance_metric', 'euclidean')
        # The web schema sends its default (True) explicitly.  Keep the
        # programmatic API's omitted parameter compatible with the historical
        # PCA-only path, while an explicit True remains strict below.
        use_corrected = self._as_bool(self.params.get('use_corrected', False), False)
        requested_n_pcs = int(self.params.get('n_pcs', 25))
        if requested_n_pcs <= 0:
            raise ValueError('n_pcs 必须是正整数')
        auto_select = self._as_bool(self.params.get('auto_select_resolution', False), False)
        resolution_metric = self.params.get('resolution_metric', 'silhouette')
        metric_subsample = int(self.params.get('metric_subsample', 10000))
        batch_key = str(self.params.get('batch_key', 'batch') or '').strip()
        primary_resolution = self.params.get('primary_resolution', None)
        if primary_resolution == '':
            primary_resolution = None

        # Resolve once and use this same representation for neighbors, UMAP,
        # resolution metrics, and the persisted summary.  In particular,
        # ``use_corrected=true`` must never silently fall back to X_pca.
        use_rep, reuse_graph = self._resolve_representation(adata, use_corrected)
        representation = np.asarray(adata.obsm[use_rep])
        if representation.ndim != 2 or representation.shape[1] == 0:
            raise ValueError(f'聚类表示 {use_rep} 不是有效的二维矩阵')

        self.progress(20, f"Computing KNN graph (n_neighbors={n_neighbors})...")
        used_n_pcs, pc_diagnostics = resolve_analysis_n_pcs(
            adata, requested_n_pcs, representation_key=use_rep,
        )
        if used_n_pcs < requested_n_pcs or pc_diagnostics.get('upstream_selected_n_pcs'):
            self.progress(
                -1,
                f"聚类实际使用 {used_n_pcs} PCs（请求 {requested_n_pcs}；"
                f"上游选择 {pc_diagnostics.get('upstream_selected_n_pcs', '—')}，"
                f"表示维度 {pc_diagnostics['available_n_pcs']}）",
            )
        if reuse_graph:
            # BBKNN 校正结果保存在邻居图（obsp connectivities）中，直接复用，
            # 不再重建 KNN，否则 BBKNN 的批次整合会被丢弃。
            neighbors_params = (adata.uns.get('neighbors', {}) or {}).get('params', {}) or {}
            if 'connectivities' not in adata.obsp or not neighbors_params:
                raise RuntimeError(
                    '输入标记为 BBKNN 校正图，但缺少 connectivities/neighbors；'
                    '请重新运行批次校正。'
                )
            self.progress(20, '复用 BBKNN 校正后的邻居图，不再重建 KNN。')
        else:
            self._build_knn(
                adata, use_rep, n_neighbors, used_n_pcs, distance_metric,
            )
            neighbors_params = (adata.uns.get('neighbors', {}) or {}).get('params', {}) or {}
            recorded_neighbors_rep = neighbors_params.get('use_rep')
            if use_rep != 'X_pca' and recorded_neighbors_rep != use_rep:
                raise RuntimeError(
                    f'邻居图表示不一致：请求 {use_rep}，Scanpy 记录为 {recorded_neighbors_rep!r}'
                )
        # Scanpy records ``None`` when X_pca is selected by its automatic
        # default; expose the resolved representation rather than that
        # implementation detail in the platform summary.
        neighbors_use_rep = use_rep
        adata.uns['clustering_neighbors'] = {
            'requested_use_corrected': bool(use_corrected),
            'resolved_use_rep': use_rep,
            'neighbors_use_rep': neighbors_use_rep,
            'n_neighbors': int(n_neighbors),
            'n_pcs': int(used_n_pcs),
            'distance_metric': str(distance_metric),
        }

        # The clustering graph may now be built on Harmony/ComBat/other
        # corrected coordinates.  Recompute UMAP from that graph instead of
        # plotting an X_umap inherited from the pre-correction input.
        self.progress(30, f"Computing UMAP from {use_rep} neighbors...")
        sc.tl.umap(adata, random_state=0)
        adata.uns['clustering_umap'] = {
            'source_representation': use_rep,
            'n_pcs': int(used_n_pcs),
            'recomputed_after_neighbors': True,
            'random_state': 0,
        }

        self.progress(40, f"Running {clustering_method} clustering at resolutions: {resolutions}...")
        self._run_clusterings(adata, resolutions, clustering_method, n_iterations)

        self.progress(75, "Generating cluster UMAP plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        for res in resolutions:
            key = cluster_key_for_method(clustering_method, res)
            if key in adata.obs.columns:
                n_clusters = adata.obs[key].nunique()
                try:
                    fig_static = self.build_publication_umap(
                        adata, key,
                        title=f'{clustering_method.title()} (res={res}, {n_clusters} clusters)'
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_static,
                        plots_dir,
                        f'cluster_umap_res_{resolution_filename_token(res)}.png',
                        'umap', f'Clusters (res={res})'
                    ))
                    import matplotlib.pyplot as plt
                    plt.close(fig_static)
                except Exception as exc:
                    self.progress(-1, f'聚类 UMAP 导出失败：{exc}')

        # 多分辨率 UMAP 比较图
        if 'X_umap' in adata.obsm:
            res_keys = [cluster_key_for_method(clustering_method, res) for res in resolutions
                        if cluster_key_for_method(clustering_method, res) in adata.obs.columns]
            fig_multi = umap_panel_figure(
                adata, res_keys,
                titles=[f'res={res}' for res in resolutions
                        if cluster_key_for_method(clustering_method, res) in adata.obs.columns],
                point_size=4.5, opacity=0.68,
                ncols=min(3, max(1, len(res_keys))),
            )
            if fig_multi is not None:
                result_files.extend(self.save_matplotlib_figure(
                    fig_multi, plots_dir, 'cluster_multi_res_umap.png',
                    'umap', '多分辨率聚类比较', formats=('png', 'svg'), dpi=300,
                    preserve_aspect=True,
                ))

        self._plot_resolution_sankey(
            adata, resolutions, clustering_method, plots_dir, result_files,
        )

        # 聚类 marker dotplot：按每个 cluster 独立排名差异 marker，而不是
        # 从固定面板截取前几个基因。保留 Universal marker 作为少量人工
        # 验证锚点，且把最终基因和阈值写入 summary。
        first_key = cluster_key_for_method(clustering_method, resolutions[0])
        marker_selection = self._run_marker_preview(
            adata, first_key, plots_dir, result_files,
        )

        marker_selection_payload = {
            'cluster_key': first_key,
            'cluster_markers': marker_selection.get('cluster_markers', {}),
            'n_data_driven': int(marker_selection.get('n_data_driven', 0)),
            'n_classic_anchor': int(marker_selection.get('n_classic_anchor', 0)),
            'n_relaxed': int(marker_selection.get('n_relaxed', 0)),
            'warnings': marker_selection.get('warnings', []),
            'parameters': marker_selection.get('parameters', {}),
            'visualizations': marker_selection.get('visualizations', {}),
        }
        adata.uns['marker_selection'] = copy.deepcopy(marker_selection_payload)

        # Select the primary resolution.  ``best_resolution`` is only a valid
        # label when an actual metric comparison was performed.
        auto_best_res, resolution_scores = self._select_best_resolution(
            adata, resolutions, clustering_method, use_rep, used_n_pcs,
            auto_select, resolution_metric, metric_subsample,
        )
        auto_selection_valid = bool(auto_select and resolution_scores)
        selected_res = auto_best_res if auto_selection_valid else resolutions[0]
        resolution_selection_mode = 'auto_metric' if auto_selection_valid else 'default_first'

        if primary_resolution is not None:
            try:
                requested_res = float(primary_resolution)
                matched_res = next((res for res in resolutions if math.isclose(
                    requested_res, res, rel_tol=0.0, abs_tol=1e-8,
                )), None)
                if matched_res is not None and cluster_key_for_method(clustering_method, matched_res) in adata.obs.columns:
                    selected_res = matched_res
                    resolution_selection_mode = 'explicit_primary'
                else:
                    logger.warning("指定主分辨率 %s 不在已计算分辨率中，使用 %s", primary_resolution, selected_res)
            except (TypeError, ValueError):
                logger.warning("无法解析主分辨率 %s，使用 %s", primary_resolution, selected_res)

        primary_key = cluster_key_for_method(clustering_method, selected_res)
        if primary_key in adata.obs.columns:
            adata.obs['leiden'] = adata.obs[primary_key].copy()
            if clustering_method == 'louvain':
                adata.obs['louvain'] = adata.obs[primary_key].copy()
        elif 'leiden' not in adata.obs.columns:
            fallback_key = cluster_key_for_method(clustering_method, resolutions[0])
            adata.obs['leiden'] = adata.obs[fallback_key].copy()

        if self.params.get('show_cluster_size_bar', True) and 'leiden' in adata.obs.columns:
            cluster_labels = adata.obs['leiden'].astype(str)
            cluster_counts = cluster_labels.value_counts().sort_index(key=lambda idx: idx.map(lambda x: (len(x), x)))
            size_colors = [
                NATURE_PALETTE[3] if int(value) < 50 else NATURE_PALETTE[0]
                for value in cluster_counts.values
            ]
            fig_size = bar_figure(
                cluster_counts.index.tolist(), cluster_counts.values,
                title=f'Cluster Size Distribution (primary res={selected_res})',
                x_label='Cluster', y_label='Cell count',
                colors=size_colors,
                annotations=[f'{int(value):,}' + (' (<50 review)' if int(value) < 50 else '')
                             for value in cluster_counts.values],
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_size, plots_dir, 'cluster_size_bar.png',
                'bar', 'Cluster Size Distribution', formats=('png', 'svg'), dpi=300,
            ))

        batch_values = adata.obs[batch_key] if batch_key in adata.obs.columns else None
        cluster_donor_composition = []
        from modules.io_utils import obs_grouping_info
        batch_info = obs_grouping_info(
            adata, batch_key, max_categories=50,
            max_numeric_categories=20, require_multiple=False,
        ) if batch_values is not None else {'valid': False, 'reason': '列不存在'}
        batch_is_invalid = batch_values is not None and not batch_info['valid']
        if (
            self.params.get('show_cluster_batch_composition', True)
            and 'leiden' in adata.obs.columns
            and batch_values is not None
            and not batch_is_invalid
        ):
            batch_table = pd.crosstab(
                adata.obs['leiden'].astype(str),
                batch_values.astype(str),
                normalize='index',
            )
            fig_batch = grouped_bar_figure(
                batch_table.index.tolist(),
                [(str(batch), batch_table[batch].values) for batch in batch_table.columns],
                title=f'Batch Composition by Cluster ({batch_key})',
                x_label='Cluster', y_label='Fraction',
                stacked=True, max_series=12,
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_batch, plots_dir, 'cluster_batch_composition.png',
                'bar', 'Cluster Batch Composition', formats=('png', 'svg'), dpi=300,
            ))
            from modules.sc_figure_diagnostics import cluster_composition_table
            composition_frame = cluster_composition_table(adata.obs, 'leiden', batch_key)
            if not composition_frame.empty:
                cluster_donor_composition = composition_frame.to_dict(orient='records')
                composition_path = os.path.join(
                    self.project_dir, 'results', 'cluster_donor_composition.csv',
                )
                os.makedirs(os.path.dirname(composition_path), exist_ok=True)
                composition_frame.to_csv(composition_path, index=False)
                result_files.append({
                    'file_path': composition_path, 'file_type': 'csv', 'category': 'table',
                    'label': 'Cluster donor composition review',
                })
        else:
            if batch_is_invalid:
                self.progress(
                    80,
                    f"跳过批次组成图：'{batch_key}' 不是安全的分类列，"
                    f"{batch_info.get('reason', '')}",
                )
            # Do not leave a composition image from an earlier run visible
            # when the current input has no batch annotation (or the option is
            # disabled).
            try:
                self.remove_plot_artifacts(
                    'cluster_batch_composition.png',
                    'cluster_batch_composition.svg',
                )
            except OSError as exc:
                logger.warning("无法清理旧批次组成图：%s", exc)

        from modules.io_utils import batch_cluster_overlap
        batch_overlap = batch_cluster_overlap(
            adata, 'leiden', batch_key,
        ) if batch_values is not None and batch_info.get('valid') else None
        if batch_overlap and batch_overlap.get('warnings'):
            for warning in batch_overlap['warnings']:
                self.progress(-1, warning)

        # 主分辨率 UMAP：添加 cluster label，便于汇报和截图
        if self.params.get('show_labeled_umap', True) and 'X_umap' in adata.obsm and 'leiden' in adata.obs.columns:
            fig_label = umap_figure(
                adata, 'leiden',
                title=f'{clustering_method.title()} Clusters with Labels (primary res={selected_res})',
                point_size=5, opacity=0.68, label_categories=True,
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_label, plots_dir, 'cluster_umap_labeled.png',
                'umap', '带标签 Cluster UMAP', formats=('png', 'svg'), dpi=300,
            ))

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'clustering')

        self.progress(100, "Done")
        summary = {
            f'n_clusters_{res}': int(adata.obs[cluster_key_for_method(clustering_method, res)].nunique())
            for res in resolutions
            if cluster_key_for_method(clustering_method, res) in adata.obs.columns
        }
        summary['cluster_counts_by_resolution'] = {
            str(res): int(adata.obs[cluster_key_for_method(clustering_method, res)].nunique())
            for res in resolutions
            if cluster_key_for_method(clustering_method, res) in adata.obs.columns
        }
        summary['resolutions'] = resolutions
        summary['clustering_method'] = clustering_method
        summary['requested_use_corrected'] = bool(use_corrected)
        summary['resolved_use_rep'] = use_rep
        summary['neighbors_use_rep'] = neighbors_use_rep
        summary['use_rep'] = use_rep
        summary['requested_n_pcs'] = requested_n_pcs
        summary['n_pcs'] = used_n_pcs
        summary['pc_diagnostics'] = pc_diagnostics
        summary['umap_source_representation'] = use_rep
        summary['umap_recomputed_after_neighbors'] = True
        summary['metric_subsample'] = metric_subsample
        summary['resolution_scores'] = resolution_scores
        summary['resolution_selection_mode'] = resolution_selection_mode
        summary['primary_resolution'] = selected_res
        summary['primary_cluster_key'] = primary_key
        summary['primary_cluster_count'] = int(adata.obs[primary_key].nunique())
        summary['best_resolution'] = auto_best_res if auto_selection_valid else None
        summary['batch_overlap'] = batch_overlap
        summary['cluster_donor_composition'] = cluster_donor_composition
        summary['small_clusters_lt_50'] = int(
            (adata.obs['leiden'].astype(str).value_counts() < 50).sum()
        ) if 'leiden' in adata.obs.columns else 0
        summary['marker_selection'] = copy.deepcopy(marker_selection_payload)
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }
