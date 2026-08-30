"""Targeted re-clustering and marker/enrichment analysis for one parent cluster."""

from modules.base import BaseAnalysis
from modules.io_utils import obs_grouping_info


class SubclusterAnalysis(BaseAnalysis):
    MODULE_NAME = "subcluster"
    DISPLAY_NAME = "子簇精细分析"
    DESCRIPTION = "选定一个已有簇进行重聚类、marker、热图和通路富集"
    INPUT_REQUIRES = ['leiden']

    @staticmethod
    def _matrix_to_array(matrix):
        import numpy as np
        return matrix.toarray() if hasattr(matrix, 'toarray') else np.asarray(matrix)

    def validate_input(self, adata):
        cluster_key = self.params.get('source_cluster_key', 'leiden')
        if cluster_key not in adata.obs.columns:
            return f"Cluster column '{cluster_key}' not found. Run clustering first or select an existing obs column."
        grouping = obs_grouping_info(
            adata, cluster_key, max_categories=50,
            max_numeric_categories=20, require_multiple=False,
        )
        if not grouping['valid']:
            return f"Cluster column '{cluster_key}' is not a valid categorical cluster column: {grouping['reason']}"
        target = str(self.params.get('target_cluster', '')).strip()
        if not target:
            return "Please select a target cluster to refine."
        if target not in set(adata.obs[cluster_key].astype(str)):
            return f"Target cluster '{target}' does not exist in '{cluster_key}'."
        return None

    def _enrich_markers(self, marker_df, results_dir, plots_dir, database, organism, top_n, pval_cutoff):
        """Run best-effort Enrichr independently for each subcluster.

        Enrichment is intentionally non-fatal: network-backed Enrichr must not
        discard successful re-clustering and differential-expression results.
        """
        import os
        import pandas as pd

        # logreg 等方法的 scanpy 输出只有 names/scores，没有 pvals_adj 与
        # logfoldchanges；直接取列会 KeyError。富集属 best-effort，缺列时
        # 跳过并保留聚类/DEG 结果（与模块的降级设计一致）。
        required_columns = {'pvals_adj', 'logfoldchanges'}
        missing_columns = required_columns - set(marker_df.columns)
        if missing_columns or marker_df.empty:
            self.progress(
                -1,
                f'marker 表缺少 {sorted(missing_columns)} 列'
                '（如 logreg 方法不提供 p 值），跳过通路富集。'
            )
            return [], 0

        try:
            import gseapy as gp
        except ImportError:
            self.progress(-1, 'gseapy 未安装，跳过通路富集；子簇聚类和差异结果已保留。')
            return [], 0

        result_files, all_results = [], []
        libraries = [database] if database else ['GO_Biological_Process_2023']
        for group, group_df in marker_df.groupby('subcluster'):
            genes = group_df.loc[
                (group_df['pvals_adj'] <= pval_cutoff) & (group_df['logfoldchanges'] > 0), 'names'
            ].astype(str).drop_duplicates().head(200).tolist()
            if len(genes) < 3:
                self.progress(-1, f'子簇 {group} 显著上调 marker 少于 3 个，跳过富集。')
                continue
            try:
                enr = gp.enrichr(
                    gene_list=genes, gene_sets=libraries, organism=organism,
                    outdir=None, no_plot=True,
                ).results
            except Exception as exc:
                self.progress(-1, f'子簇 {group} 通路富集不可用：{exc}')
                continue
            if enr is None or enr.empty:
                continue
            enr = enr.copy()
            enr['subcluster'] = str(group)
            all_results.append(enr)

        if not all_results:
            return result_files, 0

        combined = pd.concat(all_results, ignore_index=True)
        csv_path = os.path.join(results_dir, 'subcluster_enrichment_results.csv')
        combined.to_csv(csv_path, index=False)
        result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table',
                             'label': '子簇通路富集结果'})

        from modules.bulk_enrichment import _enrichment_figure

        # Render one publication-style chart per subcluster.  Keeping the
        # panels separate avoids mixing pathway labels from unrelated clusters
        # and matches the BP/MF/KEGG + overlap-count layout used by Bulk.
        term_col = 'Term' if 'Term' in combined.columns else combined.columns[0]
        p_col = 'Adjusted P-value' if 'Adjusted P-value' in combined.columns else 'P-value'
        plot_df = combined.sort_values(p_col).groupby('subcluster', group_keys=False).head(top_n)
        for group, group_df in plot_df.groupby('subcluster'):
            fig = _enrichment_figure(
                group_df,
                title=f'子簇 {group} enrich result',
                database=database,
            )
            if fig is None:
                continue
            safe_group = ''.join(ch if ch.isalnum() or ch in '._-' else '_' for ch in str(group)).strip('_') or 'cluster'
            result_files.extend(self.save_matplotlib_figure(
                fig,
                plots_dir,
                f'subcluster_enrichment_{safe_group}.png',
                'enrichment',
                f'子簇 {group} 通路富集图',
            ))
            import matplotlib.pyplot as plt
            plt.close(fig)
        return result_files, len(combined)

    def run(self, input_path):
        import os
        import numpy as np
        import pandas as pd
        import scanpy as sc
        from modules.native_figures import bar_figure, heatmap_figure, umap_figure

        self.progress(5, 'Loading clustered data...')
        adata = self.load_adata(input_path)
        adata = self.apply_scope(adata)
        source_key = self.params.get('source_cluster_key', 'leiden')
        target = str(self.params.get('target_cluster', '')).strip()
        if source_key not in adata.obs or not target:
            raise ValueError('必须提供存在的聚类列和目标簇。')
        grouping = obs_grouping_info(
            adata, source_key, max_categories=50,
            max_numeric_categories=20, require_multiple=False,
        )
        if not grouping['valid']:
            raise ValueError(f"source_cluster_key '{source_key}' 不是有效的分类聚类列：{grouping['reason']}")
        mask = adata.obs[source_key].astype(str) == target
        sub = adata[mask].copy()
        min_cells = int(self.params.get('min_cells', 30))
        if sub.n_obs < min_cells:
            raise ValueError(f"目标簇 {target} 仅有 {sub.n_obs} 个细胞，少于最小阈值 {min_cells}，不建议重聚类。")
        if sub.n_vars < 2:
            raise ValueError('目标簇的可用基因少于 2 个，无法进行子簇重聚类。')
        sub.obs['parent_cluster'] = target

        n_neighbors = min(int(self.params.get('n_neighbors', 15)), max(2, sub.n_obs - 1))
        # 与父级聚类保持一致：优先使用批次校正后的表示，而不是回退到
        # 未校正的 X_pca（否则子簇可能按 batch 分裂）。
        from modules.clustering import ClusteringAnalysis
        use_rep = next((
            key for key in ClusteringAnalysis.CORRECTED_REPRESENTATIONS
            if key in sub.obsm
        ), 'X_pca' if 'X_pca' in sub.obsm else None)
        if use_rep is None:
            self.progress(15, 'PCA missing; computing PCA for the selected cluster...')
            sc.pp.pca(sub, n_comps=min(50, sub.n_obs - 1, sub.n_vars - 1))
            use_rep = 'X_pca'
        from modules.pc_strategy import resolve_analysis_n_pcs
        used_n_pcs, pc_diagnostics = resolve_analysis_n_pcs(
            sub, int(self.params.get('n_pcs', 25)), representation_key=use_rep,
        )
        self.progress(25, f'Rebuilding neighbor graph for {sub.n_obs} cells ({use_rep}, {used_n_pcs} PCs)...')
        sc.pp.neighbors(sub, n_neighbors=n_neighbors, use_rep=use_rep,
                        n_pcs=used_n_pcs,
                        metric=self.params.get('distance_metric', 'euclidean'))
        resolution = float(self.params.get('resolution', 0.8))
        method = self.params.get('clustering_method', 'leiden')
        self.progress(40, f'Running {method} subclustering...')
        if method == 'louvain':
            sc.tl.louvain(sub, resolution=resolution, key_added='subcluster')
        else:
            sc.tl.leiden(sub, resolution=resolution, key_added='subcluster', flavor='igraph',
                         n_iterations=int(self.params.get('n_iterations', 2)))
        sub.obs['subcluster'] = sub.obs['subcluster'].astype('category')
        sc.tl.umap(sub, min_dist=float(self.params.get('umap_min_dist', 0.4)))

        plots_dir, results_dir = self.ensure_plots_dir(), os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        result_files = []
        fig_umap = umap_figure(
            sub, 'subcluster', title=f'{source_key}={target} 子簇重聚类',
            point_size=6, opacity=0.75,
        )
        result_files.extend(self.save_matplotlib_figure(
            fig_umap, plots_dir, 'subcluster_umap.png', 'umap',
            '目标簇子簇 UMAP', formats=('png', 'svg'), dpi=300,
        ))

        counts = sub.obs['subcluster'].astype(str).value_counts().sort_index()
        fig_counts = bar_figure(
            counts.index, counts.values, title='子簇细胞数',
            x_label='子簇', y_label='细胞数',
            annotations=[f'{int(value):,}' for value in counts.values],
        )
        result_files.extend(self.save_matplotlib_figure(
            fig_counts, plots_dir, 'subcluster_size_bar.png', 'bar',
            '子簇细胞数', formats=('png', 'svg'), dpi=300,
        ))

        self.progress(55, 'Computing subcluster marker genes...')
        groups = list(sub.obs['subcluster'].cat.categories)
        marker_df = pd.DataFrame()
        if len(groups) >= 2:
            sc.tl.rank_genes_groups(
                sub, groupby='subcluster', method=self.params.get('deg_method', 'wilcoxon'),
                n_genes=min(int(self.params.get('n_genes', 50)), sub.n_vars),
                use_raw=False,
            )
            marker_tables = []
            for group in groups:
                table = sc.get.rank_genes_groups_df(sub, group=group)
                table['subcluster'] = str(group)
                marker_tables.append(table)
            marker_df = pd.concat(marker_tables, ignore_index=True)
            marker_path = os.path.join(results_dir, 'subcluster_deg_results.csv')
            marker_df.to_csv(marker_path, index=False)
            result_files.append({'file_path': marker_path, 'file_type': 'csv', 'category': 'table',
                                 'label': '子簇差异表达结果（每簇 Top N）'})
        else:
            self.progress(-1, '当前参数仅得到 1 个子簇，已保留重聚类结果；请提高分辨率后再运行差异表达和富集。')

        heatmap_top_n = int(self.params.get('marker_heatmap_top_n', 5))
        marker_selection = {
            'cluster_markers': {}, 'genes': [], 'warnings': [],
            'n_data_driven': 0, 'n_classic_anchor': 0, 'n_relaxed': 0,
            'parameters': {},
        }
        try:
            from modules.annotation import select_cluster_marker_genes
            marker_selection = select_cluster_marker_genes(
                sub,
                'subcluster',
                method=self.params.get('deg_method', 'wilcoxon'),
                n_rank_genes=min(int(self.params.get('n_genes', 50)), sub.n_vars),
                min_markers=1,
                max_markers=max(1, heatmap_top_n),
                padj_cutoff=float(self.params.get('enrichment_pval_cutoff', 0.05)),
                min_pct=0.10,
                min_delta_pct=0.05,
            )
            if marker_selection.get('warnings'):
                self.progress(-1, '；'.join(marker_selection['warnings'][:2]))
        except Exception as exc:
            self.progress(-1, f'子簇 marker 选择失败，回退到 DEG 排名：{exc}')

        genes = marker_selection.get('genes', [])[:80]
        if not genes:
            # Keep the heatmap useful if a tiny/degenerate subcluster has no
            # gene passing the strict evidence filters.
            genes = []
            for group in groups:
                if not marker_df.empty:
                    genes.extend(marker_df[marker_df['subcluster'] == str(group)]['names'].head(heatmap_top_n).astype(str))
            genes = [gene for gene in dict.fromkeys(genes) if gene in sub.var_names][:80]
        if genes:
            expr = self._matrix_to_array(sub[:, genes].X).astype(float)
            group_labels = sub.obs['subcluster'].astype(str)
            means = np.vstack([expr[(group_labels == str(group)).values].mean(axis=0) for group in groups]).T
            z = (means - means.mean(axis=1, keepdims=True)) / (means.std(axis=1, keepdims=True) + 1e-10)
            fig_heatmap = heatmap_figure(
                np.clip(z, -3, 3), x_labels=[str(g) for g in groups], y_labels=genes,
                title='子簇 Marker 热图', x_label='子簇', y_label='Marker genes',
                colorbar_label='z-score',
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_heatmap, plots_dir, 'subcluster_marker_heatmap.png',
                'heatmap', '子簇 Marker 热图', formats=('png', 'svg'), dpi=300,
            ))

        enrichment_n = 0
        if self.params.get('run_enrichment', True) and not marker_df.empty:
            self.progress(75, 'Running pathway enrichment for subclusters...')
            enrichment_files, enrichment_n = self._enrich_markers(
                marker_df, results_dir, plots_dir, self.params.get('enrichment_database', 'GO_Biological_Process_2023'),
                self.params.get('organism', 'Human'), int(self.params.get('enrichment_top_n', 10)),
                float(self.params.get('enrichment_pval_cutoff', 0.05)),
            )
            result_files.extend(enrichment_files)

        self.progress(90, 'Saving subcluster dataset...')
        sub.uns['marker_selection'] = {
            'cluster_key': 'subcluster',
            'cluster_markers': marker_selection.get('cluster_markers', {}),
            'n_data_driven': int(marker_selection.get('n_data_driven', 0)),
            'n_classic_anchor': int(marker_selection.get('n_classic_anchor', 0)),
            'n_relaxed': int(marker_selection.get('n_relaxed', 0)),
            'warnings': marker_selection.get('warnings', []),
            'parameters': marker_selection.get('parameters', {}),
        }
        output_path = self.save_output(sub, 'subcluster')
        self.progress(100, 'Done')
        return {'output_adata': output_path, 'result_files': result_files,
                'summary': {'source_cluster_key': source_key, 'target_cluster': target,
                            'scope_key': str(self.params.get('scope_key', '') or '').strip() or None,
                            'scope_values': ([v.strip() for v in str(self.params.get('scope_values', '') or '').split(',') if v.strip()] or None),
                            'n_cells': int(sub.n_obs), 'n_subclusters': len(groups),
                            'subcluster_counts': {str(k): int(v) for k, v in counts.items()},
                            'marker_selection': {
                                'cluster_markers': marker_selection.get('cluster_markers', {}),
                                'n_data_driven': int(marker_selection.get('n_data_driven', 0)),
                                'n_classic_anchor': int(marker_selection.get('n_classic_anchor', 0)),
                                'n_relaxed': int(marker_selection.get('n_relaxed', 0)),
                                'warnings': marker_selection.get('warnings', []),
                                'parameters': marker_selection.get('parameters', {}),
                            },
                            'n_enrichment_terms': int(enrichment_n)}}
