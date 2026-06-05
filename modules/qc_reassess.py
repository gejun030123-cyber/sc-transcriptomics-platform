import os
import json
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


class QCReassessAnalysis(BaseAnalysis):
    MODULE_NAME = "qc_reassess"
    DISPLAY_NAME = "QC 重新评估"
    DESCRIPTION = "聚类后检查 doublet 和 QC 指标，标记低质量簇"
    INPUT_REQUIRES = ['leiden']

    def validate_input(self, adata):
        cluster_key = self.params.get('cluster_key', 'leiden')
        if cluster_key not in adata.obs.columns:
            return f"Column '{cluster_key}' not found."
        return None

    def run(self, input_path):
        import scanpy as sc
        import plotly.graph_objects as go
        from modules.visualization import umap_scatter

        self.progress(5, "加载数据...")
        adata = sc.read_h5ad(input_path)

        cluster_key = self.params.get('cluster_key', 'leiden')
        doublet_threshold = float(self.params.get('doublet_threshold', 0.3))
        mt_threshold = float(self.params.get('mt_threshold', 15.0))

        self.progress(20, "计算各簇 QC 指标...")
        result_files = []
        plots_dir = os.path.join(self.project_dir, 'plots')
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)

        # 计算每个簇的统计
        clusters = adata.obs[cluster_key].unique()
        cluster_stats = []
        for c in sorted(clusters):
            mask = adata.obs[cluster_key] == c
            n_cells = mask.sum()
            doublet_frac = 0.0
            if 'predicted_doublet' in adata.obs.columns:
                doublet_frac = adata.obs.loc[mask, 'predicted_doublet'].mean()
            elif 'doublet_score' in adata.obs.columns:
                doublet_frac = (adata.obs.loc[mask, 'doublet_score'] > 0.5).mean()
            mt_mean = adata.obs.loc[mask, 'pct_counts_mt'].mean() if 'pct_counts_mt' in adata.obs.columns else 0
            counts_mean = adata.obs.loc[mask, 'total_counts'].mean() if 'total_counts' in adata.obs.columns else 0
            genes_mean = adata.obs.loc[mask, 'n_genes_by_counts'].mean() if 'n_genes_by_counts' in adata.obs.columns else 0
            is_low = doublet_frac > doublet_threshold or mt_mean > mt_threshold
            cluster_stats.append({
                'cluster': str(c), 'n_cells': n_cells,
                'doublet_fraction': round(doublet_frac, 4),
                'mean_pct_mt': round(mt_mean, 2),
                'mean_total_counts': round(counts_mean, 0),
                'mean_n_genes': round(genes_mean, 0),
                'low_quality': is_low
            })

        stats_df = pd.DataFrame(cluster_stats)
        csv_path = os.path.join(results_dir, 'qc_reassess_clusters.csv')
        stats_df.to_csv(csv_path, index=False)
        result_files.append({'file_path': csv_path, 'file_type': 'csv', 'category': 'table', 'label': '簇 QC 统计'})

        n_low = stats_df['low_quality'].sum()

        self.progress(50, "生成 UMAP 图...")
        # UMAP with doublet score
        if 'X_umap' in adata.obsm and 'doublet_score' in adata.obs.columns:
            fig = umap_scatter(adata, 'doublet_score', title='Doublet Score')
            fpath = os.path.join(plots_dir, 'qc_reassess_doublet_umap.json')
            with open(fpath, 'w') as f: f.write(json.dumps(fig))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': 'Doublet Score UMAP'})

        # UMAP with MT percentage
        if 'X_umap' in adata.obsm and 'pct_counts_mt' in adata.obs.columns:
            fig = umap_scatter(adata, 'pct_counts_mt', title='MT Percentage')
            fpath = os.path.join(plots_dir, 'qc_reassess_mt_umap.json')
            with open(fpath, 'w') as f: f.write(json.dumps(fig))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': 'MT% UMAP'})

        # UMAP with cluster highlighting
        if 'X_umap' in adata.obsm:
            fig = umap_scatter(adata, cluster_key, title=f'Clusters ({cluster_key})')
            fpath = os.path.join(plots_dir, 'qc_reassess_clusters_umap.json')
            with open(fpath, 'w') as f: f.write(json.dumps(fig))
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': '聚类 UMAP'})

        self.progress(90, "保存输出...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'qc_reassess_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_clusters': len(clusters),
                'n_low_quality': int(n_low),
                'low_quality_clusters': stats_df[stats_df['low_quality']]['cluster'].tolist(),
                'doublet_threshold': doublet_threshold,
                'mt_threshold': mt_threshold,
            }
        }
