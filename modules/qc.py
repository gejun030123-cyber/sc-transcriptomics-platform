from modules.base import BaseAnalysis
from modules.constants import S_GENES, G2M_GENES
from modules.io_utils import resolve_obs_grouping


class QCAnalysis(BaseAnalysis):
    MODULE_NAME = "qc"
    DISPLAY_NAME = "质控"
    DESCRIPTION = "MT/ribo/hb 过滤 + Scrublet 双细胞 + 细胞周期评分 + 复杂度过滤"
    INPUT_REQUIRES = []

    def run(self, input_path):
        import scanpy as sc
        import omicverse as ov
        import numpy as np
        import os
        import pandas as pd
        import matplotlib.pyplot as plt
        from modules.figure_style import (
            NATURE_AXIS, NATURE_GRID, NATURE_MUTED, NATURE_PALETTE, NATURE_TEXT,
        )

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
        score_cell_cycle = bool(self.params.get('score_cell_cycle', True))
        self.progress(
            25,
            "Scoring cell cycle phases..." if score_cell_cycle
            else "Skipping optional cell-cycle scoring...",
        )
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

        if score_cell_cycle and len(s_in) >= 5 and len(g2m_in) >= 5:
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
        qc_before = adata.obs.copy()
        n_genes_before = adata.shape[1]
        requested_batch = self.params.get('batch_key', 'batch')
        batch_key, batch_info = resolve_obs_grouping(
            adata, requested_batch, max_categories=50,
            max_numeric_categories=20, require_multiple=True,
        )
        if batch_key and not hasattr(adata.obs[batch_key].dtype, 'categories'):
            adata.obs[batch_key] = adata.obs[batch_key].astype(str).astype('category')
        if requested_batch and requested_batch in adata.obs.columns and not batch_info.get('requested_valid', False):
            self.progress(-1, f"批次列已跳过：{batch_info.get('requested_reason', '不是有效分类列')}")

        mito_perc = float(self.params.get('mito_perc', 0.2))  # 0-1 scale, 0.2 = 20%
        nUMIs_min = int(self.params.get('nUMIs', 500))
        ngenes_min = int(self.params.get('detected_genes', 250))
        ngenes_max = int(self.params.get('max_detected_genes', 0))  # 0 = 不限制
        ribo_perc_max = float(self.params.get('ribo_perc', 0))  # 0 = 不过滤
        hb_perc_max = float(self.params.get('hb_perc', 0))  # 0 = 不过滤
        batch_adaptive = self.params.get('batch_adaptive_qc', False)
        mad_multiplier = float(self.params.get('mad_multiplier', 3.0))

        # OmicVerse keeps the historical argument spelling ``tresh``.
        # Passing ``thresh`` is silently accepted via **kwargs but ignored,
        # which would fall back to the library defaults (notably 15% MT).
        adata = ov.pp.qc(
            adata,
            tresh={
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

        # ── 7. 生成 Nature-style 图表 ────────────────────────────────────
        self.progress(60, "Generating QC plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        def _style_axis(axis, title, xlabel=None, ylabel=None):
            axis.set_title(title, loc='left', pad=7, fontsize=9,
                           fontweight='semibold', color=NATURE_TEXT)
            if xlabel:
                axis.set_xlabel(xlabel, color=NATURE_TEXT)
            if ylabel:
                axis.set_ylabel(ylabel, color=NATURE_TEXT)
            axis.tick_params(labelsize=8, length=3, width=0.7, colors=NATURE_AXIS)
            axis.grid(axis='y', color=NATURE_GRID, linewidth=0.55, alpha=0.72)
            axis.set_axisbelow(True)
            for spine_name, spine in axis.spines.items():
                spine.set_visible(spine_name in ('left', 'bottom'))
                spine.set_color('#98A2B3')
                spine.set_linewidth(0.7)

        def _finite_values(series):
            values = pd.to_numeric(series, errors='coerce').to_numpy(dtype=float)
            return values[np.isfinite(values)]

        def _save(fig, filename, category, label):
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, filename, category, label,
                formats=('png', 'svg'), dpi=300,
            ))
            plt.close(fig)

        # 过滤前后 QC 指标对比
        if self.params.get('show_qc_filter_summary', True):
            metrics = [
                ('n_cells', 'Cells', n_before, n_after),
                ('n_genes', 'Genes', n_genes_before, adata.shape[1]),
            ]
            for col, label in [
                ('n_genes_by_counts', 'Median detected genes'),
                ('total_counts', 'Median total counts'),
                ('pct_counts_mt', 'Median MT%'),
                ('pct_counts_ribo', 'Median ribo%'),
            ]:
                if col in qc_before.columns and col in adata.obs.columns:
                    metrics.append((col, label, float(qc_before[col].median()), float(adata.obs[col].median())))
            panel_specs = [
                (m[1], m[2], m[3], '{:,.0f}') if m[0] in ('n_cells', 'n_genes')
                else (m[1], m[2], m[3], '{:,.1f}')
                for m in metrics
            ]
            fig_filter, axes = plt.subplots(2, 2, figsize=(8.6, 6.0), squeeze=False)
            for axis, (label, before, after, number_format) in zip(axes.ravel(), panel_specs):
                bars = axis.bar(
                    [0, 1], [before, after], width=0.56,
                    color=[NATURE_PALETTE[6], NATURE_PALETTE[0]],
                    edgecolor='white', linewidth=0.5,
                )
                axis.set_xticks([0, 1], ['Before QC', 'After QC'])
                _style_axis(axis, label, ylabel='Value')
                ymax = max(abs(float(before)), abs(float(after)), 1.0)
                axis.set_ylim(0, ymax * 1.2)
                for bar, value in zip(bars, [before, after]):
                    axis.text(bar.get_x() + bar.get_width() / 2,
                              bar.get_height() + ymax * 0.035,
                              number_format.format(value), ha='center', va='bottom',
                              fontsize=8, color=NATURE_TEXT)
            fig_filter.suptitle('QC filtering overview', x=0.06, ha='left',
                                fontsize=13, fontweight='semibold', color=NATURE_TEXT)
            fig_filter.text(0.06, 0.01,
                            'Each metric uses its own scale so filtering effects remain readable.',
                            fontsize=7.5, color=NATURE_MUTED)
            _save(fig_filter, 'qc_filter_summary.png', 'qc', 'QC 过滤前后对比')

        # Scrublet doublet score 分布
        if self.params.get('show_doublet_histogram', True) and 'doublet_score' in adata.obs.columns:
            scores = _finite_values(adata.obs['doublet_score'])
            if len(scores):
                fig_doublet, axis = plt.subplots(figsize=(7.6, 4.8))
                axis.hist(scores, bins=42, color=NATURE_PALETTE[0], alpha=0.78,
                          edgecolor='white', linewidth=0.35)
                median = float(np.median(scores))
                axis.axvline(median, color=NATURE_PALETTE[3], lw=1.1,
                             ls=(0, (3, 2)), label=f'Median {median:.3f}')
                _style_axis(axis, 'Scrublet doublet score distribution',
                            xlabel='Doublet score', ylabel='Cells')
                axis.legend(loc='upper right', fontsize=8)
                _save(fig_doublet, 'qc_doublet_score_histogram.png', 'histogram',
                      'Scrublet Doublet Score 分布')

        # QC Violin
        violin_keys = [k for k in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo']
                       if k in adata.obs.columns]
        if violin_keys:
            fig_violin, axes = plt.subplots(1, len(violin_keys),
                                            figsize=(max(8.5, 2.7 * len(violin_keys)), 4.8),
                                            squeeze=False)
            metric_labels = {
                'n_genes_by_counts': 'Detected genes',
                'total_counts': 'Total counts',
                'pct_counts_mt': 'Mitochondrial reads (%)',
                'pct_counts_ribo': 'Ribosomal reads (%)',
            }
            for index, key in enumerate(violin_keys):
                axis = axes.ravel()[index]
                values = _finite_values(adata.obs[key])
                if len(values):
                    parts = axis.violinplot(values, positions=[1], widths=0.72,
                                            showmeans=False, showmedians=True,
                                            showextrema=False)
                    color = NATURE_PALETTE[index % len(NATURE_PALETTE)]
                    for body in parts['bodies']:
                        body.set_facecolor(color)
                        body.set_edgecolor(color)
                        body.set_alpha(0.62)
                    parts['cmedians'].set_color(NATURE_TEXT)
                    parts['cmedians'].set_linewidth(1.1)
                    median = float(np.median(values))
                    axis.text(1.0, median, f'  {median:,.1f}', va='center',
                              fontsize=7.5, color=NATURE_TEXT)
                    if key in ('n_genes_by_counts', 'total_counts') and np.nanmin(values) > 0:
                        axis.set_yscale('log')
                axis.set_xticks([1], [metric_labels.get(key, key)])
                axis.tick_params(axis='x', labelrotation=25)
                _style_axis(axis, metric_labels.get(key, key), ylabel='Value')
            fig_violin.suptitle('QC metric distributions', x=0.04, ha='left',
                                fontsize=13, fontweight='semibold', color=NATURE_TEXT)
            _save(fig_violin, 'qc_violin.png', 'violin', 'QC Violin Plots')

        # QC 散点图（Counts vs Genes，颜色 = MT%）
        if 'total_counts' in adata.obs.columns and 'n_genes_by_counts' in adata.obs.columns:
            x = _finite_values(adata.obs['total_counts'])
            y = _finite_values(adata.obs['n_genes_by_counts'])
            count = min(len(x), len(y))
            color_vals = (_finite_values(adata.obs['pct_counts_mt'])[:count]
                          if 'pct_counts_mt' in adata.obs.columns else None)
            fig_scatter, axis = plt.subplots(figsize=(8.8, 5.7))
            points = axis.scatter(x[:count], y[:count], c=color_vals, cmap='RdYlBu_r'
                                  if color_vals is not None else None,
                                  color=NATURE_PALETTE[0] if color_vals is None else None,
                                  s=11, alpha=0.58, linewidths=0, rasterized=True)
            axis.set_xscale('log')
            axis.set_yscale('log')
            _style_axis(axis, 'Library complexity and mitochondrial burden',
                        xlabel='Total counts per cell', ylabel='Detected genes per cell')
            if color_vals is not None:
                colorbar = fig_scatter.colorbar(points, ax=axis, fraction=0.032, pad=0.02)
                colorbar.set_label('Mitochondrial reads (%)', fontsize=8)
                colorbar.outline.set_visible(False)
            axis.text(0.015, 0.96, f'n = {count:,} cells', transform=axis.transAxes,
                      ha='left', va='top', color=NATURE_MUTED, fontsize=8)
            _save(fig_scatter, 'qc_scatter.png', 'scatter', 'QC Scatter')

        # Novelty score 散点图
        if 'novelty_score' in adata.obs.columns:
            x = _finite_values(adata.obs['total_counts'])
            y = _finite_values(adata.obs['novelty_score'])
            count = min(len(x), len(y))
            fig_nov, axis = plt.subplots(figsize=(8.0, 5.1))
            axis.scatter(x[:count], y[:count], color=NATURE_PALETTE[2], s=11,
                         alpha=0.56, linewidths=0, rasterized=True)
            axis.set_xscale('log')
            _style_axis(axis, 'Novelty score versus library size',
                        xlabel='Total counts per cell',
                        ylabel='Novelty score (genes / counts)')
            _save(fig_nov, 'qc_novelty.png', 'scatter', 'Novelty Score')

        # 细胞周期散点图（S_score vs G2M_score，颜色 = phase）
        if cc_available and 'S_score' in adata.obs.columns:
            phase_colors = {'G1': '#1f77b4', 'S': '#ff7f0e', 'G2M': '#2ca02c'}
            fig_cc, axis = plt.subplots(figsize=(7.6, 5.0))
            for phase, color in phase_colors.items():
                mask = adata.obs['phase'] == phase
                if mask.sum() == 0:
                    continue
                axis.scatter(adata.obs.loc[mask, 'S_score'],
                             adata.obs.loc[mask, 'G2M_score'], s=12, color=color,
                             alpha=0.64, linewidths=0, label=phase, rasterized=True)
            _style_axis(axis, 'Cell-cycle scoring', xlabel='S score', ylabel='G2M score')
            axis.legend(title='Phase', loc='best', fontsize=8, title_fontsize=8)
            _save(fig_cc, 'qc_cell_cycle.png', 'scatter', 'Cell Cycle Scoring')

        # UMAP（如有）
        if 'X_umap' in adata.obsm:
            for color_key in ['batch', 'leiden', 'phase']:
                if color_key in adata.obs.columns:
                    fig_umap = self.build_publication_umap(
                        adata, color_key, title=f'UMAP by {color_key}'
                    )
                    _save(fig_umap, f'qc_umap_{color_key}.png', 'umap',
                          f'UMAP by {color_key}')

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
                'cells_removed_by_qc_and_doublet': n_before - n_after_scrublet,
                'novelty_median': round(float(adata.obs['novelty_score'].median()), 4) if 'novelty_score' in adata.obs.columns else None,
                'cell_cycle_available': cc_available,
                'phase_counts': phase_counts,
                's_genes_found': len(s_in),
                'g2m_genes_found': len(g2m_in),
                'requested_batch_key': str(requested_batch or ''),
                'batch_key': batch_key,
            }
        }
