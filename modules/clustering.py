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
        from modules.visualization import umap_scatter, categorical_color_map
        import json
        import numpy as np
        import pandas as pd
        import plotly.graph_objects as go

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
                fig_json = json.dumps(umap_scatter(adata, key, title=f'Leiden (res={res}, {n_clusters} clusters)'))
                fpath = os.path.join(plots_dir, f'cluster_umap_{res}.json')
                with open(fpath, 'w') as f: f.write(fig_json)
                result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': f'Clusters (res={res})'})
                try:
                    fig_static = sc.pl.umap(
                        adata, color=key, title=f'Leiden (res={res}, {n_clusters} clusters)',
                        show=False, return_fig=True, frameon=False,
                        size=self.get_viz_params()['umap_point_size'], legend_loc='right margin',
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_static, plots_dir, f'cluster_umap_{res}.png',
                        'umap', f'Clusters (res={res})'
                    ))
                    import matplotlib.pyplot as plt
                    plt.close(fig_static)
                except Exception as exc:
                    self.progress(-1, f'静态聚类 UMAP 导出失败（不影响交互图）：{exc}')

        # 多分辨率 UMAP 比较图
        if 'X_umap' in adata.obsm:
            from plotly.subplots import make_subplots
            n_res = len(resolutions)
            fig_multi = make_subplots(rows=1, cols=n_res, subplot_titles=[f'res={r}' for r in resolutions])
            for ci, res in enumerate(resolutions, 1):
                key = f'leiden_{res}'
                if key in adata.obs.columns:
                    coords = adata.obsm['X_umap'][:, :2]
                    cats = adata.obs[key].values
                    color_map = categorical_color_map(adata, key)
                    for cat in sorted(set(cats)):
                        mask = cats == cat
                        fig_multi.add_trace(go.Scattergl(
                            x=coords[mask, 0], y=coords[mask, 1], mode='markers',
                            marker=dict(size=2, opacity=0.6, color=color_map.get(str(cat), '#bdbdbd')),
                            name=str(cat), showlegend=(ci==1)
                        ), row=1, col=ci)
            fig_multi.update_layout(height=400, width=350*n_res, title='多分辨率聚类比较')
            for i in range(1, n_res+1):
                fig_multi.update_xaxes(title_text='UMAP-1', row=1, col=i)
                fig_multi.update_yaxes(title_text='UMAP-2', row=1, col=i)
            fpath = os.path.join(plots_dir, 'cluster_multi_res_umap.json')
            with open(fpath, 'w') as f: f.write(fig_multi.to_json())
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': '多分辨率聚类比较'})

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
                    fig_sankey = go.Figure(data=[go.Sankey(
                        arrangement='snap',
                        node=dict(label=node_labels, pad=15, thickness=14),
                        link=dict(source=sources, target=targets, value=values),
                    )])
                    fig_sankey.update_layout(
                        title='Cluster Resolution Sankey',
                        width=900,
                        height=max(450, 120 * len(sankey_res)),
                        font=dict(size=11),
                    )
                    result_files.append(self.save_plotly_json(
                        fig_sankey, plots_dir, 'cluster_resolution_sankey.json',
                        'sankey', '分辨率分群流向图'
                    ))

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
                    result_files.extend(self.save_matplotlib_figure(
                        fig_dot, plots_dir, 'cluster_marker_dotplot.png', 'dotplot',
                        'Marker Dotplot'
                    ))
                    import matplotlib.pyplot as plt
                    plt.close('all')
        except Exception as e:
            logger.warning("生成 marker dotplot 失败（注释模块可能不可用）: %s", e)

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
            fig_size = go.Figure()
            fig_size.add_trace(go.Bar(
                x=cluster_counts.index.tolist(),
                y=cluster_counts.values.astype(int).tolist(),
                marker_color='#3949ab',
                text=cluster_counts.values.astype(int).tolist(),
                textposition='outside',
                hovertemplate='Cluster: %{x}<br>Cells: %{y}<extra></extra>',
            ))
            fig_size.update_layout(
                title=f'Cluster Size Distribution (primary res={best_res})',
                xaxis_title='Cluster',
                yaxis_title='Cell count',
                plot_bgcolor='white',
                width=max(650, 45 * max(1, len(cluster_counts))),
                height=430,
            )
            result_files.append(self.save_plotly_json(
                fig_size, plots_dir, 'cluster_size_bar.json',
                'bar', 'Cluster Size Distribution'
            ))

        if (
            self.params.get('show_cluster_batch_composition', True)
            and 'leiden' in adata.obs.columns
            and batch_key in adata.obs.columns
        ):
            batch_table = pd.crosstab(
                adata.obs['leiden'].astype(str),
                adata.obs[batch_key].astype(str),
                normalize='index',
            )
            fig_batch = go.Figure()
            for batch in batch_table.columns:
                fig_batch.add_trace(go.Bar(
                    x=batch_table.index.tolist(),
                    y=batch_table[batch].values,
                    name=str(batch),
                    hovertemplate='Cluster: %{x}<br>' + batch_key + ': ' + str(batch) + '<br>Fraction: %{y:.2%}<extra></extra>',
                ))
            fig_batch.update_layout(
                title=f'Batch Composition by Cluster ({batch_key})',
                xaxis_title='Cluster',
                yaxis_title='Fraction',
                barmode='stack',
                plot_bgcolor='white',
                width=max(700, 55 * max(1, len(batch_table.index))),
                height=460,
            )
            result_files.append(self.save_plotly_json(
                fig_batch, plots_dir, 'cluster_batch_composition.json',
                'bar', 'Cluster Batch Composition'
            ))

        # 主分辨率 UMAP：添加 cluster label，便于汇报和截图
        if self.params.get('show_labeled_umap', True) and 'X_umap' in adata.obsm and 'leiden' in adata.obs.columns:
            coords = adata.obsm['X_umap'][:, :2]
            labels = adata.obs['leiden'].astype(str)
            color_map = categorical_color_map(adata, 'leiden')
            fig_label = go.Figure()
            for cl in sorted(labels.unique(), key=lambda x: (len(x), x)):
                mask = labels == cl
                fig_label.add_trace(go.Scattergl(
                    x=coords[mask, 0],
                    y=coords[mask, 1],
                    mode='markers',
                    marker=dict(size=3, opacity=0.65, color=color_map.get(str(cl), '#bdbdbd')),
                    name=str(cl),
                    text=adata.obs_names[mask.values].tolist(),
                    hovertemplate='%{text}<br>Cluster: ' + str(cl) + '<extra></extra>',
                ))
                fig_label.add_annotation(
                    x=float(np.median(coords[mask, 0])),
                    y=float(np.median(coords[mask, 1])),
                    text=str(cl),
                    showarrow=False,
                    font=dict(size=14, color='black'),
                    bgcolor='rgba(255,255,255,0.75)',
                    bordercolor='rgba(0,0,0,0.25)',
                    borderwidth=1,
                )
            fig_label.update_layout(
                title=f'Leiden Clusters with Labels (primary res={best_res})',
                xaxis_title='UMAP-1',
                yaxis_title='UMAP-2',
                plot_bgcolor='white',
                width=720,
                height=540,
            )
            result_files.append(self.save_plotly_json(
                fig_label, plots_dir, 'cluster_umap_labeled.json',
                'umap', '带标签 Cluster UMAP'
            ))

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'clustering')

        self.progress(100, "Done")
        summary = {f'n_clusters_{res}': int(adata.obs[f'leiden_{res}'].nunique()) for res in resolutions if f'leiden_{res}' in adata.obs.columns}
        summary['resolutions'] = resolutions
        summary['best_resolution'] = best_res
        summary['primary_resolution'] = best_res
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }
