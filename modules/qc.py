from modules.base import BaseAnalysis
from modules.constants import S_GENES, G2M_GENES


class QCAnalysis(BaseAnalysis):
    MODULE_NAME = "qc"
    DISPLAY_NAME = "质控"
    DESCRIPTION = "MT/ribo/hb 过滤 + Scrublet 双细胞 + 细胞周期评分 + 复杂度过滤"
    INPUT_REQUIRES = []

    def run(self, input_path):
        import scanpy as sc
        import omicverse as ov
        from modules.visualization import umap_scatter, violin_plot
        import plotly.graph_objects as go
        import numpy as np
        import os, json

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        # 应用自定义过滤规则
        adata = self.apply_filters(adata, 'qc')

        # 保存 counts 层
        save_counts = self.params.get('save_counts_layer', True)
        if save_counts:
            adata.layers["counts"] = adata.X.copy()

        # ── 1. 标记基因集 ─────────────────────────────────────────────────
        self.progress(10, "Flagging MT/ribo/hb genes...")
        if 'gene_name' in adata.var.columns:
            gene_names_str = adata.var['gene_name'].fillna('').astype(str)
        else:
            gene_names_str = adata.var_names.astype(str)
        adata.var["mt"] = gene_names_str.str.startswith("MT-")
        adata.var["ribo"] = gene_names_str.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = gene_names_str.str.contains("^HB[^(P)]")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo", "hb"],
            inplace=True, percent_top=[20], log1p=True
        )

        # ── 2. Novelty score（复杂度）──────────────────────────────────────
        self.progress(18, "Calculating novelty score...")
        adata.obs['novelty_score'] = (
            adata.obs['n_genes_by_counts'] / adata.obs['total_counts'].replace(0, np.nan)
        ).fillna(0)

        # ── 3. 细胞周期评分 ─────────────────────────────────────────────
        self.progress(25, "Scoring cell cycle phases...")
        # 筛选实际存在于数据中的基因（区分 Ensembl ID 和基因名）
        var_names_set = set(adata.var_names.astype(str))
        s_in = [g for g in S_GENES if g in var_names_set]
        g2m_in = [g for g in G2M_GENES if g in var_names_set]

        # 如果 var_names 是 Ensembl ID，尝试用 gene_name 列映射
        if len(s_in) < 5 and 'gene_name' in adata.var.columns:
            gene_name_to_idx = {}
            for idx, gn in zip(adata.var_names, adata.var['gene_name'].fillna('').astype(str)):
                gene_name_to_idx[gn] = idx
            s_in_ensembl = [gene_name_to_idx[g] for g in S_GENES if g in gene_name_to_idx]
            g2m_in_ensembl = [gene_name_to_idx[g] for g in G2M_GENES if g in gene_name_to_idx]
            if len(s_in_ensembl) > len(s_in):
                s_in = s_in_ensembl
            if len(g2m_in_ensembl) > len(g2m_in):
                g2m_in = g2m_in_ensembl

        if len(s_in) >= 5 and len(g2m_in) >= 5:
            # 需要先有 log1p normalized 数据用于打分
            adata_cc = adata.copy()
            if 'log1p_total_counts' not in adata_cc.obs.columns:
                sc.pp.normalize_total(adata_cc, target_sum=1e4)
                sc.pp.log1p(adata_cc)
            sc.tl.score_genes_cell_cycle(
                adata_cc, s_genes=s_in, g2m_genes=g2m_in
            )
            adata.obs['S_score'] = adata_cc.obs['S_score']
            adata.obs['G2M_score'] = adata_cc.obs['G2M_score']
            adata.obs['phase'] = adata_cc.obs['phase']
            cc_available = True
            del adata_cc
        else:
            cc_available = False

        # ── 4. Scrublet 双细胞检测 + 阈值过滤 ──────────────────────────
        self.progress(35, "Running Scrublet doublet detection + QC filtering...")
        n_before = adata.shape[0]
        requested_batch = self.params.get('batch_key', 'batch')
        batch_key = requested_batch if requested_batch in adata.obs.columns else None

        mito_perc = float(self.params.get('mito_perc', 0.2))
        nUMIs_min = int(self.params.get('nUMIs', 500))
        ngenes_min = int(self.params.get('detected_genes', 250))
        ngenes_max = int(self.params.get('max_detected_genes', 0))  # 0 = 不限制
        ribo_perc_max = float(self.params.get('ribo_perc', 0))  # 0 = 不过滤
        hb_perc_max = float(self.params.get('hb_perc', 0))  # 0 = 不过滤
        batch_adaptive = self.params.get('batch_adaptive_qc', False)
        mad_multiplier = float(self.params.get('mad_multiplier', 3.0))

        adata = ov.pp.qc(
            adata,
            thresh={
                'mito_perc': mito_perc,
                'nUMIs': nUMIs_min,
                'detected_genes': ngenes_min,
            },
            doublets_method='scrublet',
            batch_key=batch_key,
            filter_doublets=True,
        )
        n_after_scrublet = adata.shape[0]

        # ── 5. 额外过滤：基因数上限 ──────────────────────────────────────
        if ngenes_max > 0:
            mask = adata.obs['n_genes_by_counts'] <= ngenes_max
            n_before_upper = adata.shape[0]
            adata = adata[mask].copy()

        # ── 6. 额外过滤：核糖体比例上限 ─────────────────────────────────
        if ribo_perc_max > 0 and 'pct_counts_ribo' in adata.obs.columns:
            mask = adata.obs['pct_counts_ribo'] <= ribo_perc_max * 100
            adata = adata[mask].copy()

        # ── 6b. 额外过滤：血红蛋白比例上限 ──────────────────────────────
        if hb_perc_max > 0 and 'pct_counts_hb' in adata.obs.columns:
            mask = adata.obs['pct_counts_hb'] <= hb_perc_max * 100
            adata = adata[mask].copy()

        # ── 6c. 批次自适应 QC（MAD 方法）──────────────────────────────
        if batch_adaptive and batch_key and batch_key in adata.obs.columns:
            self.progress(55, "Batch-adaptive QC filtering (MAD)...")
            for qc_col in ['n_genes_by_counts', 'total_counts']:
                if qc_col not in adata.obs.columns:
                    continue
                keep_mask = np.ones(adata.n_obs, dtype=bool)
                for batch_val in adata.obs[batch_key].unique():
                    batch_mask = adata.obs[batch_key] == batch_val
                    vals = adata.obs.loc[batch_mask, qc_col]
                    median_val = vals.median()
                    mad_val = np.median(np.abs(vals - median_val))
                    if mad_val > 0:
                        lower = median_val - mad_multiplier * mad_val * 1.4826
                        upper = median_val + mad_multiplier * mad_val * 1.4826
                        keep_mask[batch_mask] = (vals >= lower) & (vals <= upper)
                adata = adata[keep_mask].copy()

        n_after = adata.shape[0]

        # ── 7. 生成图表 ─────────────────────────────────────────────────
        self.progress(60, "Generating QC plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        # QC Violin
        violin_keys = [k for k in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo']
                       if k in adata.obs.columns]
        if violin_keys:
            fig_json = json.dumps(violin_plot(adata, violin_keys, title='QC Metrics'))
            fpath = os.path.join(plots_dir, 'qc_violin.json')
            with open(fpath, 'w') as f: f.write(fig_json)
            result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'violin', 'label': 'QC Violin Plots'})

        # QC 散点图（Counts vs Genes，颜色 = MT%）
        if 'total_counts' in adata.obs.columns and 'n_genes_by_counts' in adata.obs.columns:
            fig_scatter = go.Figure()
            color_vals = adata.obs['pct_counts_mt'].values if 'pct_counts_mt' in adata.obs.columns else None
            fig_scatter.add_trace(go.Scattergl(
                x=adata.obs['total_counts'], y=adata.obs['n_genes_by_counts'],
                mode='markers', marker=dict(size=3, color=color_vals, colorscale='Reds',
                                            colorbar=dict(title='MT%'), opacity=0.6),
                text=adata.obs.index.tolist(),
                hovertemplate='%{text}<br>Counts: %{x:.0f}<br>Genes: %{y:.0f}<br>MT%: %{marker.color:.1f}'
            ))
            fig_scatter.update_layout(title='QC: Counts vs Genes', xaxis_title='Total Counts',
                                     yaxis_title='Detected Genes', plot_bgcolor='white', width=600, height=400)
            result_files.append(self.save_plotly_json(fig_scatter, plots_dir, 'qc_scatter.json', 'scatter', 'QC Scatter'))

        # Novelty score 散点图
        if 'novelty_score' in adata.obs.columns:
            fig_nov = go.Figure()
            color_vals = adata.obs['pct_counts_mt'].values if 'pct_counts_mt' in adata.obs.columns else None
            fig_nov.add_trace(go.Scattergl(
                x=adata.obs['total_counts'], y=adata.obs['novelty_score'],
                mode='markers', marker=dict(size=3, color=color_vals, colorscale='Reds',
                                            colorbar=dict(title='MT%'), opacity=0.6),
                text=adata.obs.index.tolist(),
                hovertemplate='%{text}<br>Counts: %{x:.0f}<br>Novelty: %{y:.3f}<br>MT%: %{marker.color:.1f}'
            ))
            fig_nov.update_layout(title='QC: Novelty Score vs Counts', xaxis_title='Total Counts',
                                  yaxis_title='Novelty Score (n_genes / total_counts)',
                                  plot_bgcolor='white', width=600, height=400)
            result_files.append(self.save_plotly_json(fig_nov, plots_dir, 'qc_novelty.json', 'scatter', 'Novelty Score'))

        # 细胞周期散点图（S_score vs G2M_score，颜色 = phase）
        if cc_available and 'S_score' in adata.obs.columns:
            phase_colors = {'G1': '#1f77b4', 'S': '#ff7f0e', 'G2M': '#2ca02c'}
            fig_cc = go.Figure()
            for phase, color in phase_colors.items():
                mask = adata.obs['phase'] == phase
                if mask.sum() == 0:
                    continue
                fig_cc.add_trace(go.Scattergl(
                    x=adata.obs.loc[mask, 'S_score'],
                    y=adata.obs.loc[mask, 'G2M_score'],
                    mode='markers', marker=dict(size=3, color=color, opacity=0.6),
                    name=phase,
                    text=adata.obs.index[mask].tolist(),
                    hovertemplate='%{text}<br>S_score: %{x:.3f}<br>G2M_score: %{y:.3f}<br>Phase: ' + phase
                ))
            fig_cc.update_layout(title='Cell Cycle Scoring', xaxis_title='S_score',
                                yaxis_title='G2M_score', plot_bgcolor='white',
                                width=600, height=400,
                                legend=dict(title='Phase'))
            result_files.append(self.save_plotly_json(fig_cc, plots_dir, 'qc_cell_cycle.json', 'scatter', 'Cell Cycle Scoring'))

        # UMAP（如有）
        if 'X_umap' in adata.obsm:
            for color_key in ['batch', 'leiden', 'phase']:
                if color_key in adata.obs.columns:
                    fig_json = json.dumps(umap_scatter(adata, color_key, title=f'UMAP by {color_key}'))
                    fpath = os.path.join(plots_dir, f'qc_umap_{color_key}.json')
                    with open(fpath, 'w') as f: f.write(fig_json)
                    result_files.append({'file_path': fpath, 'file_type': 'plotly_json', 'category': 'umap', 'label': f'UMAP by {color_key}'})

        # ── 8. 保存输出 ─────────────────────────────────────────────────
        self.progress(85, "Saving output...")
        output_path = self.save_output(adata, 'qc')

        # 汇总细胞周期比例
        phase_counts = {}
        if cc_available and 'phase' in adata.obs.columns:
            phase_counts = adata.obs['phase'].value_counts().to_dict()

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'cells_before': n_before,
                'cells_after': n_after,
                'cells_removed': n_before - n_after,
                'pct_removed': round((n_before - n_after) / max(n_before, 1) * 100, 1),
                'n_genes': adata.shape[1],
                'doublets_removed': n_before - n_after_scrublet,
                'novelty_median': round(float(adata.obs['novelty_score'].median()), 4) if 'novelty_score' in adata.obs.columns else None,
                'cell_cycle_available': cc_available,
                'phase_counts': phase_counts,
                's_genes_found': len(s_in),
                'g2m_genes_found': len(g2m_in),
            }
        }
