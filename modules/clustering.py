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
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        from modules.native_figures import (
            bar_figure,
            grouped_bar_figure,
            umap_figure,
            umap_panel_figure,
        )
        from modules.figure_style import NATURE_PALETTE, NATURE_TEXT

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
        batch_key = self.params.get('batch_key', 'batch')
        primary_resolution = self.params.get('primary_resolution', None)
        if primary_resolution == '':
            primary_resolution = None

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

        self.progress(75, "Generating cluster UMAP plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        for res in resolutions:
            key = f'leiden_{res}'
            if key in adata.obs.columns:
                n_clusters = adata.obs[key].nunique()
                try:
                    fig_static = self.build_publication_umap(
                        adata, key, title=f'Leiden (res={res}, {n_clusters} clusters)'
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_static, plots_dir, f'cluster_umap_{res}.png',
                        'umap', f'Clusters (res={res})'
                    ))
                    import matplotlib.pyplot as plt
                    plt.close(fig_static)
                except Exception as exc:
                    self.progress(-1, f'聚类 UMAP 导出失败：{exc}')

        # 多分辨率 UMAP 比较图
        if 'X_umap' in adata.obsm:
            res_keys = [f'leiden_{res}' for res in resolutions
                        if f'leiden_{res}' in adata.obs.columns]
            fig_multi = umap_panel_figure(
                adata, res_keys,
                titles=[f'res={res}' for res in resolutions
                        if f'leiden_{res}' in adata.obs.columns],
                point_size=4.5, opacity=0.68,
                ncols=min(3, max(1, len(res_keys))),
            )
            if fig_multi is not None:
                result_files.extend(self.save_matplotlib_figure(
                    fig_multi, plots_dir, 'cluster_multi_res_umap.png',
                    'umap', '多分辨率聚类比较', formats=('png', 'svg'), dpi=300,
                    preserve_aspect=True,
                ))

        # Resolution Sankey：展示不同分辨率之间的簇分裂关系
        if self.params.get('show_resolution_sankey', True) and len(resolutions) > 1:
            sankey_res = sorted([res for res in resolutions if f'leiden_{res}' in adata.obs.columns])
            if len(sankey_res) > 1:
                node_labels = []
                node_index = {}
                for res in sankey_res:
                    key = f'leiden_{res}'
                    for cat in sorted(adata.obs[key].astype(str).unique(), key=lambda x: (len(x), x)):
                        label = f'res {res}: {cat}'
                        node_index[(res, cat)] = len(node_labels)
                        node_labels.append(label)
                sources, targets, values = [], [], []
                for left, right in zip(sankey_res[:-1], sankey_res[1:]):
                    left_key = f'leiden_{left}'
                    right_key = f'leiden_{right}'
                    flow = pd.crosstab(adata.obs[left_key].astype(str), adata.obs[right_key].astype(str))
                    for left_cat in flow.index:
                        for right_cat in flow.columns:
                            count = int(flow.loc[left_cat, right_cat])
                            if count <= 0:
                                continue
                            sources.append(node_index[(left, str(left_cat))])
                            targets.append(node_index[(right, str(right_cat))])
                            values.append(count)
                if values:
                    fig_sankey, ax_sankey = plt.subplots(
                        figsize=(max(8.0, 2.2 * len(sankey_res)),
                                 max(5.0, 0.32 * len(node_labels) + 2.5)), dpi=150,
                    )
                    max_value = max(values)
                    x_positions = np.linspace(0.08, 0.92, len(sankey_res))
                    node_y = {}
                    for res_index, res in enumerate(sankey_res):
                        cats = sorted(adata.obs[f'leiden_{res}'].astype(str).unique(),
                                      key=lambda x: (len(x), x))
                        y_positions = np.linspace(0.9, 0.1, max(1, len(cats)))
                        for y_value, cat in zip(y_positions, cats):
                            node_y[(res, cat)] = y_value
                            ax_sankey.scatter(x_positions[res_index], y_value, s=160,
                                              color=NATURE_PALETTE[res_index % len(NATURE_PALETTE)],
                                              edgecolor='white', linewidth=0.8, zorder=3)
                            ax_sankey.text(x_positions[res_index], y_value, cat,
                                           ha='center', va='center', fontsize=7,
                                           color='white', zorder=4)
                        ax_sankey.text(x_positions[res_index], 1.02, f'res {res}',
                                       ha='center', va='bottom', fontsize=9,
                                       color=NATURE_TEXT, fontweight='semibold')
                    for source, target, value in zip(sources, targets, values):
                        for res_index, res in enumerate(sankey_res[:-1]):
                            left_key = f'leiden_{res}'
                            right_key = f'leiden_{sankey_res[res_index + 1]}'
                            if node_labels[source].startswith(f'res {res}:') and node_labels[target].startswith(f'res {sankey_res[res_index + 1]}:'):
                                left_cat = node_labels[source].split(': ', 1)[1]
                                right_cat = node_labels[target].split(': ', 1)[1]
                                ax_sankey.plot(
                                    [x_positions[res_index], x_positions[res_index + 1]],
                                    [node_y[(res, left_cat)], node_y[(sankey_res[res_index + 1], right_cat)]],
                                    color=NATURE_PALETTE[3], alpha=0.28,
                                    linewidth=0.7 + 4.0 * value / max_value,
                                    zorder=1,
                                )
                                break
                    ax_sankey.set_xlim(0, 1)
                    ax_sankey.set_ylim(0, 1.08)
                    ax_sankey.axis('off')
                    ax_sankey.set_title('Cluster Resolution Flow', loc='left', pad=12,
                                        fontsize=10, fontweight='semibold', color=NATURE_TEXT)
                    result_files.extend(self.save_matplotlib_figure(
                        fig_sankey, plots_dir, 'cluster_resolution_sankey.png',
                        'sankey', '分辨率分群流向图', formats=('png', 'svg'), dpi=300,
                    ))

        # 聚类 marker dotplot：按每个 cluster 独立排名差异 marker，而不是
        # 从固定面板截取前几个基因。保留 Universal marker 作为少量人工
        # 验证锚点，且把最终基因和阈值写入 summary。
        marker_selection = {
            'cluster_markers': {}, 'genes': [], 'warnings': [],
            'n_data_driven': 0, 'n_classic_anchor': 0, 'n_relaxed': 0,
            'parameters': {},
        }
        compute_marker_preview = bool(self.params.get('compute_marker_preview', True))
        try:
            if not compute_marker_preview:
                marker_selection['warnings'].append(
                    'Marker preview disabled for proportions/expression-only batch run.'
                )
                raise StopIteration
            from modules.annotation import (
                DEFAULT_UNIVERSAL_MARKERS,
                select_cluster_marker_genes,
            )
            first_key = f'leiden_{resolutions[0]}'
            if first_key in adata.obs.columns:
                marker_selection = select_cluster_marker_genes(
                    adata,
                    first_key,
                    classic_markers=DEFAULT_UNIVERSAL_MARKERS,
                    method=self.params.get('marker_selection_method', 'wilcoxon'),
                    n_rank_genes=int(self.params.get('marker_rank_genes', 200)),
                    min_markers=int(self.params.get('marker_min_per_cluster', 2)),
                    max_markers=int(self.params.get('marker_max_per_cluster', 5)),
                    padj_cutoff=float(self.params.get('marker_padj_cutoff', 0.05)),
                    min_pct=float(self.params.get('marker_min_pct', 0.10)),
                    min_delta_pct=float(self.params.get('marker_min_delta_pct', 0.05)),
                )
                dotplot_genes = marker_selection.get('genes', [])[:40]
                if marker_selection.get('warnings'):
                    self.progress(-1, '；'.join(marker_selection['warnings'][:2]))
                if dotplot_genes:
                    import scanpy as sc2
                    sc2.tl.dendrogram(adata, groupby=first_key)
                    fig_dot = sc2.pl.dotplot(
                        adata, var_names=dotplot_genes, groupby=first_key,
                        return_fig=True,
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_dot, plots_dir, 'cluster_marker_dotplot.png', 'dotplot',
                        'Cluster-specific Marker Dotplot'
                    ))
                    import matplotlib.pyplot as plt
                    plt.close('all')
        except StopIteration:
            pass
        except Exception as e:
            logger.warning("生成 marker dotplot 失败（注释模块可能不可用）: %s", e)

        adata.uns['marker_selection'] = {
            'cluster_key': f'leiden_{resolutions[0]}',
            'cluster_markers': marker_selection.get('cluster_markers', {}),
            'n_data_driven': int(marker_selection.get('n_data_driven', 0)),
            'n_classic_anchor': int(marker_selection.get('n_classic_anchor', 0)),
            'n_relaxed': int(marker_selection.get('n_relaxed', 0)),
            'warnings': marker_selection.get('warnings', []),
            'parameters': marker_selection.get('parameters', {}),
        }

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

        if primary_resolution is not None:
            try:
                requested_res = float(primary_resolution)
                if requested_res in resolutions and f'leiden_{requested_res}' in adata.obs.columns:
                    best_res = requested_res
                else:
                    logger.warning("指定主分辨率 %s 不在已计算分辨率中，使用 %s", primary_resolution, best_res)
            except (TypeError, ValueError):
                logger.warning("无法解析主分辨率 %s，使用 %s", primary_resolution, best_res)

        primary_key = f'leiden_{best_res}'
        if primary_key in adata.obs.columns:
            adata.obs['leiden'] = adata.obs[primary_key].copy()
        elif 'leiden' not in adata.obs.columns:
            adata.obs['leiden'] = adata.obs[f'leiden_{resolutions[0]}'].copy()

        if self.params.get('show_cluster_size_bar', True) and 'leiden' in adata.obs.columns:
            cluster_labels = adata.obs['leiden'].astype(str)
            cluster_counts = cluster_labels.value_counts().sort_index(key=lambda idx: idx.map(lambda x: (len(x), x)))
            fig_size = bar_figure(
                cluster_counts.index.tolist(), cluster_counts.values,
                title=f'Cluster Size Distribution (primary res={best_res})',
                x_label='Cluster', y_label='Cell count',
                annotations=[f'{int(value):,}' for value in cluster_counts.values],
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_size, plots_dir, 'cluster_size_bar.png',
                'bar', 'Cluster Size Distribution', formats=('png', 'svg'), dpi=300,
            ))

        batch_values = adata.obs[batch_key] if batch_key in adata.obs.columns else None
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
            for stale_name in ('cluster_batch_composition.png',
                               'cluster_batch_composition.svg'):
                stale_path = os.path.join(plots_dir, stale_name)
                if os.path.exists(stale_path):
                    try:
                        os.remove(stale_path)
                    except OSError:
                        logger.warning("无法清理旧批次组成图: %s", stale_path)

        # 主分辨率 UMAP：添加 cluster label，便于汇报和截图
        if self.params.get('show_labeled_umap', True) and 'X_umap' in adata.obsm and 'leiden' in adata.obs.columns:
            fig_label = umap_figure(
                adata, 'leiden',
                title=f'Leiden Clusters with Labels (primary res={best_res})',
                point_size=5, opacity=0.68, label_categories=True,
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_label, plots_dir, 'cluster_umap_labeled.png',
                'umap', '带标签 Cluster UMAP', formats=('png', 'svg'), dpi=300,
            ))

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'clustering')

        self.progress(100, "Done")
        summary = {f'n_clusters_{res}': int(adata.obs[f'leiden_{res}'].nunique()) for res in resolutions if f'leiden_{res}' in adata.obs.columns}
        summary['resolutions'] = resolutions
        summary['best_resolution'] = best_res
        summary['primary_resolution'] = best_res
        summary['marker_selection'] = {
            'cluster_key': f'leiden_{resolutions[0]}',
            'cluster_markers': marker_selection.get('cluster_markers', {}),
            'n_data_driven': int(marker_selection.get('n_data_driven', 0)),
            'n_classic_anchor': int(marker_selection.get('n_classic_anchor', 0)),
            'n_relaxed': int(marker_selection.get('n_relaxed', 0)),
            'warnings': marker_selection.get('warnings', []),
            'parameters': marker_selection.get('parameters', {}),
        }
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }
