"""Targeted re-clustering and marker/enrichment analysis for one parent cluster."""

from modules.base import BaseAnalysis


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
        import plotly.graph_objects as go

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

        term_col = 'Term' if 'Term' in combined.columns else combined.columns[0]
        p_col = 'Adjusted P-value' if 'Adjusted P-value' in combined.columns else 'P-value'
        plot_df = combined.sort_values(p_col).groupby('subcluster', group_keys=False).head(top_n)
        if not plot_df.empty:
            fig = go.Figure()
            for group, group_df in plot_df.groupby('subcluster'):
                fig.add_trace(go.Scatter(
                    x=-__import__('numpy').log10(group_df[p_col].clip(lower=1e-300)),
                    y=group_df[term_col], mode='markers', name=f'子簇 {group}',
                    marker=dict(size=10),
                    hovertemplate='%{y}<br>-log10(adjusted p): %{x:.2f}<extra></extra>',
                ))
            fig.update_layout(title='子簇 Top 通路富集', xaxis_title='-log10(Adjusted P-value)',
                              yaxis=dict(autorange='reversed'), plot_bgcolor='white', width=900,
                              height=max(450, 24 * len(plot_df)))
            result_files.append(self.save_plotly_json(
                fig, plots_dir, 'subcluster_enrichment_bubble.json', 'enrichment', '子簇通路富集气泡图'
            ))
        return result_files, len(combined)

    def run(self, input_path):
        import os
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import plotly.graph_objects as go
        from modules.visualization import umap_scatter

        self.progress(5, 'Loading clustered data...')
        adata = self.load_adata(input_path)
        source_key = self.params.get('source_cluster_key', 'leiden')
        target = str(self.params.get('target_cluster', '')).strip()
        if source_key not in adata.obs or not target:
            raise ValueError('必须提供存在的聚类列和目标簇。')
        mask = adata.obs[source_key].astype(str) == target
        sub = adata[mask].copy()
        min_cells = int(self.params.get('min_cells', 30))
        if sub.n_obs < min_cells:
            raise ValueError(f"目标簇 {target} 仅有 {sub.n_obs} 个细胞，少于最小阈值 {min_cells}，不建议重聚类。")
        if sub.n_vars < 2:
            raise ValueError('目标簇的可用基因少于 2 个，无法进行子簇重聚类。')
        sub.obs['parent_cluster'] = target

        n_neighbors = min(int(self.params.get('n_neighbors', 15)), max(2, sub.n_obs - 1))
        use_rep = 'X_pca' if 'X_pca' in sub.obsm else None
        if use_rep is None:
            self.progress(15, 'PCA missing; computing PCA for the selected cluster...')
            sc.pp.pca(sub, n_comps=min(50, sub.n_obs - 1, sub.n_vars - 1))
            use_rep = 'X_pca'
        self.progress(25, f'Rebuilding neighbor graph for {sub.n_obs} cells...')
        sc.pp.neighbors(sub, n_neighbors=n_neighbors, use_rep=use_rep,
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
        result_files.append({'file_path': os.path.join(plots_dir, 'subcluster_umap.json'),
                             'file_type': 'plotly_json', 'category': 'umap', 'label': '目标簇子簇 UMAP'})
        with open(result_files[-1]['file_path'], 'w') as handle:
            import json
            handle.write(json.dumps(umap_scatter(sub, 'subcluster', title=f'{source_key}={target} 子簇重聚类')))

        counts = sub.obs['subcluster'].astype(str).value_counts().sort_index()
        fig_counts = go.Figure(go.Bar(x=counts.index, y=counts.values, marker_color='#3949ab'))
        fig_counts.update_layout(title='子簇细胞数', xaxis_title='子簇', yaxis_title='细胞数', plot_bgcolor='white')
        result_files.append(self.save_plotly_json(fig_counts, plots_dir, 'subcluster_size_bar.json', 'bar', '子簇细胞数'))

        self.progress(55, 'Computing subcluster marker genes...')
        groups = list(sub.obs['subcluster'].cat.categories)
        marker_df = pd.DataFrame()
        if len(groups) >= 2:
            sc.tl.rank_genes_groups(sub, groupby='subcluster', method=self.params.get('deg_method', 'wilcoxon'),
                                    n_genes=min(int(self.params.get('n_genes', 50)), sub.n_vars))
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
            fig_heatmap = go.Figure(go.Heatmap(z=np.clip(z, -3, 3), x=[str(g) for g in groups], y=genes,
                                                colorscale='RdBu', zmid=0))
            fig_heatmap.update_layout(title='子簇 Marker 热图', xaxis_title='子簇', yaxis_title='Marker genes',
                                      width=max(650, 90 * len(groups)), height=max(450, 18 * len(genes)))
            result_files.append(self.save_plotly_json(fig_heatmap, plots_dir, 'subcluster_marker_heatmap.json',
                                                       'heatmap', '子簇 Marker 热图'))

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
        output_path = self.save_output(sub, 'subcluster')
        self.progress(100, 'Done')
        return {'output_adata': output_path, 'result_files': result_files,
                'summary': {'source_cluster_key': source_key, 'target_cluster': target,
                            'n_cells': int(sub.n_obs), 'n_subclusters': len(groups),
                            'subcluster_counts': {str(k): int(v) for k, v in counts.items()},
                            'n_enrichment_terms': int(enrichment_n)}}
