import os
import re
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis


# Human mitochondrial genes in Ensembl use stable IDs rather than an ``MT-``
# prefix.  This local set covers the 37 genes encoded by the human mitochondrial
# genome (protein-coding genes, rRNAs and tRNAs).  It is deliberately only a
# fallback: a supplied gene symbol or chromosome annotation remains the
# preferred and more general source of truth.
#
# The tRNAs are the fragile part.  Different GTF releases assign the same
# genomic tRNA to different stable IDs (e.g. MT-TP is ENSG00000210195 in some
# GENCODE builds and ENSG00000210196 in others), so the set is a *union* over
# releases rather than a single-release snapshot.  A chromosome column or a
# symbol column therefore always wins over this list.
_HUMAN_MITOCHONDRIAL_ENSEMBL_IDS = frozenset({
    'ENSG00000210049', 'ENSG00000211459', 'ENSG00000210077', 'ENSG00000210082',
    'ENSG00000209082', 'ENSG00000198888', 'ENSG00000210100', 'ENSG00000210107',
    'ENSG00000210112', 'ENSG00000198763', 'ENSG00000210117', 'ENSG00000210127',
    'ENSG00000210135', 'ENSG00000210140', 'ENSG00000210144', 'ENSG00000198804',
    'ENSG00000210151', 'ENSG00000210154', 'ENSG00000198712', 'ENSG00000210156',
    'ENSG00000228253', 'ENSG00000198899', 'ENSG00000198938', 'ENSG00000210164',
    'ENSG00000198840', 'ENSG00000210165', 'ENSG00000210196', 'ENSG00000212907',
    'ENSG00000198886', 'ENSG00000210174', 'ENSG00000210176', 'ENSG00000210184',
    'ENSG00000198786', 'ENSG00000198695', 'ENSG00000210191', 'ENSG00000198727',
    'ENSG00000210194', 'ENSG00000210195',
})

_GENE_SYMBOL_COLUMNS = frozenset({
    'gene_name', 'genename', 'gene_symbol', 'genesymbol', 'gene_symbols',
    'symbol', 'feature_name', 'gene',
})
_GENE_ID_COLUMNS = frozenset({
    'gene_id', 'geneid', 'gene_ids', 'ensembl_gene_id', 'ensembl_id', 'feature_id',
})
_CHROMOSOME_COLUMNS = frozenset({
    'chromosome', 'chrom', 'chr', 'seqname', 'seq_name', 'sequence_name',
})
_MITOCHONDRIAL_CHROMOSOMES = frozenset({
    'MT', 'M', 'CHRMT', 'CHRM', 'MITOCHONDRION', 'MITOCHONDRIAL',
})


class BulkQCAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_qc"
    DISPLAY_NAME = "Bulk RNA-seq 质控"
    DESCRIPTION = "计数矩阵质控：文库大小、基因检测、离群值过滤"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import matplotlib.pyplot as plt
        from modules.native_figures import (
            heatmap_figure, scatter_figure, bar_figure,
            correlation_heatmap_figure, correlation_pairwise_table,
            summarize_correlation_pairs,
        )
        from modules.figure_style import NATURE_PALETTE, NATURE_GRID, NATURE_TEXT

        self.progress(5, "加载表达矩阵...")
        from modules.io_utils import (
            read_expression_matrix, materialize_auto_sample_metadata,
            resolve_expression_measurement, run_bulk_sample_pca,
            validate_bulk_raw_counts,
        )
        adata = read_expression_matrix(input_path)
        input_measurement, measurement_info = resolve_expression_measurement(
            adata, input_path, self.params.get('input_measurement', 'auto'),
        )
        if measurement_info['requires_confirmation_for_log_transform']:
            raise ValueError(
                '自动检测到非整数连续值，但无法仅靠数值区分线性 FPKM/TPM 与已 log 的表达矩阵。'
                '请在“输入表达量尺度”中明确选择 continuous_expression 或 log_transformed 后再运行 QC。'
            )
        adata.uns['input_measurement'] = input_measurement
        adata.uns['input_measurement_provenance'] = measurement_info
        if input_measurement == 'raw_counts':
            validate_bulk_raw_counts(adata, context='Bulk QC')

        # 应用自定义过滤规则
        adata = self.apply_filters(adata, 'bulk_qc')

        min_counts = int(self.params.get('min_counts', 100000))
        min_genes = int(self.params.get('min_genes', 5000))
        max_mt_pct = float(self.params.get('max_mt_pct', 20.0))
        max_ribo_pct = float(self.params.get('max_ribo_pct', 40.0))
        max_gini = float(self.params.get('max_gini', 0))
        min_sample_expr = int(self.params.get('min_sample_expr', 0))
        min_count_threshold = float(self.params.get('min_count_threshold', 1))
        detect_outliers = self.params.get('detect_outliers', True)
        filter_strategy = self.params.get('filter_strategy', 'standard')

        # 过滤策略预设覆盖阈值
        if filter_strategy == 'strict':
            min_counts = max(min_counts, 200000)
            min_genes = max(min_genes, 8000)
            max_mt_pct = min(max_mt_pct, 15.0)
            max_ribo_pct = min(max_ribo_pct, 30.0)

        if input_measurement != 'raw_counts':
            self.progress(-1, '检测到连续或已标准化表达值：跳过 count 文库大小、MT/Ribo 百分比硬过滤。')
            min_counts, min_genes, max_mt_pct, max_ribo_pct = 0, 0, 100.0, 100.0

        if input_measurement == 'raw_counts':
            quantity_key = 'library_size'
            quantity_label = 'Library size (counts)'
            gene_metric_key = 'n_genes_detected'
            gene_metric_label = 'Detected genes'
            mt_metric_label = 'Mitochondrial fraction'
            ribo_metric_label = 'Ribosomal protein gene fraction'
        elif input_measurement == 'log_transformed':
            quantity_key = 'total_transformed_expression'
            quantity_label = 'Total transformed expression'
            gene_metric_key = 'n_expressed_genes'
            gene_metric_label = 'Expressed genes'
            mt_metric_label = 'Mitochondrial expression fraction'
            ribo_metric_label = 'Ribosomal protein gene expression fraction'
        else:
            quantity_key = 'total_expression'
            quantity_label = 'Total expression'
            gene_metric_key = 'n_expressed_genes'
            gene_metric_label = 'Expressed genes'
            mt_metric_label = 'Mitochondrial expression fraction'
            ribo_metric_label = 'Ribosomal protein gene expression fraction'
        materialize_auto_sample_metadata(adata, self.params.get('_auto_group_mapping', {}))

        self.progress(20, "计算质控指标...")
        # ``var_names`` may be Ensembl IDs, so MT- prefix matching alone is
        # not enough.  Prefer a supplied symbol, then chromosome annotation,
        # then the local canonical human Ensembl-MT IDs.  Never silently turn
        # an unidentifiable MT fraction into 0%.
        mt_detection = _detect_mitochondrial_genes(adata)
        adata.var['mt'] = mt_detection['mask']
        ribo_names = _preferred_gene_symbols(adata)
        adata.var['ribo'] = ribo_names.str.upper().str.startswith(('RPL', 'RPS'))
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
        # A fraction is only meaningful when matching mitochondrial features
        # actually carry signal.  ``calculate_qc_metrics`` returns a
        # deceptively reassuring 0.0 both when no MT feature was identified and
        # when the identified MT features are entirely absent from the matrix,
        # so keep those two failure modes explicit instead of reporting 0% as
        # evidence of clean samples.
        mt_evaluable = bool(mt_detection['n_mitochondrial_genes'] > 0)
        mt_feature_totals = (
            pd.to_numeric(adata.var.loc[adata.var['mt'], 'total_counts'], errors='coerce')
            if 'total_counts' in adata.var.columns else pd.Series(dtype=float)
        )
        mt_genes_with_counts = (
            int((mt_feature_totals.fillna(0.0) > 0).sum()) if mt_evaluable else 0
        )
        mt_total_counts = (
            float(np.nansum(adata.obs['total_counts_mt'].to_numpy(dtype=float)))
            if 'total_counts_mt' in adata.obs.columns else 0.0
        )
        mt_has_signal = bool(mt_evaluable and mt_total_counts > 0)
        # Hard filtering only makes sense for raw counts with usable MT signal;
        # comparing an all-zero column against max_mt_pct would "apply" a
        # threshold that can never reject a sample.
        mt_filter_enabled = bool(input_measurement == 'raw_counts' and mt_has_signal)
        if not mt_evaluable:
            adata.obs['pct_counts_mt'] = np.nan
            mt_detection.update({
                'status': 'not_identified',
                'figure_annotation': '',
                'message': (
                    '未识别到线粒体基因；MT% 无法计算，max_mt_pct 阈值未应用。'
                    '请提供 gene symbol、chromosome=MT/chrM 注释，或使用标准 human Ensembl gene ID。'
                ),
            })
        elif not mt_has_signal:
            mt_detection.update({
                'status': 'identified_no_reads',
                'figure_annotation': 'MT genes found, 0 MT reads',
                'message': (
                    f'识别到 {mt_detection["n_mitochondrial_genes"]} 个线粒体基因，但全部样本的 MT counts 均为 0；'
                    'MT% 恒为 0，不能作为样本质量证据，max_mt_pct 阈值未应用。'
                    '请检查上游定量是否遗漏 chrM/MT（例如 featureCounts 参考 GTF 未包含线粒体序列，'
                    '或表达矩阵在导入前已剔除 MT 基因）。'
                ),
            })
        else:
            mt_detection.update({
                'status': 'identified',
                'figure_annotation': (
                    f'{mt_genes_with_counts}/{mt_detection["n_mitochondrial_genes"]} MT genes with counts'
                ),
                'message': (
                    '' if mt_filter_enabled else
                    '已计算 MT expression fraction；连续/已变换表达输入不应用 max_mt_pct 硬过滤。'
                ),
            })
        mt_detection.update({
            'informative': mt_has_signal,
            'total_mt_counts': mt_total_counts,
            'mt_genes_with_counts': mt_genes_with_counts,
            'max_mt_pct': float(max_mt_pct),
            'threshold_applied': mt_filter_enabled,
        })
        if mt_detection['message']:
            self.progress(-1, mt_detection['message'])
        # Keep the source and threshold decision with the h5ad handoff so a
        # downstream DEG result can be audited without relying on task logs.
        adata.uns['bulk_qc_mitochondrial_detection'] = {
            key: value for key, value in mt_detection.items() if key != 'mask'
        }

        n_before = adata.n_obs
        lib_sizes = adata.obs['total_counts'].values
        n_genes_detected = adata.obs['n_genes_by_counts'].values
        mt_pct = adata.obs['pct_counts_mt'].values if 'pct_counts_mt' in adata.obs.columns else np.zeros(n_before)
        ribo_pct = adata.obs['pct_counts_ribo'].values if 'pct_counts_ribo' in adata.obs.columns else np.zeros(n_before)

        # 表达分布集中度 (Gini) 和新颖度
        raw_counts = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
        gini_values = np.array([_gini(raw_counts[i]) for i in range(n_before)])
        novelty_values = np.log10(n_genes_detected + 1) / np.log10(lib_sizes + 1)

        self.progress(40, "过滤样本...")
        # 分组推断。自动分组及自动因素已在上方正式写入 obs，后续模块可复用。
        sample_names = adata.obs.index.tolist()
        group_col = self.params.get('group_column', '').strip()
        if group_col == '_auto_group_':
            group_col = '_auto_group'
            groups = adata.obs[group_col].astype(str).tolist() if group_col in adata.obs.columns else _infer_groups(sample_names)
        elif group_col and group_col in adata.obs.columns:
            groups = adata.obs[group_col].astype(str).tolist()
        else:
            groups = adata.obs['_auto_group'].astype(str).tolist() if '_auto_group' in adata.obs.columns else _infer_groups(sample_names)
            group_col = '_auto_group'

        replicate_groups, replicate_group_info = _resolve_replicate_groups(
            adata, self.params.get('outlier_group_column', ''), groups,
        )

        # 样本过滤
        mt_mask = (mt_pct <= max_mt_pct) if mt_filter_enabled else np.ones(n_before, dtype=bool)
        mask = (lib_sizes >= min_counts) & (n_genes_detected >= min_genes) & mt_mask & (ribo_pct <= max_ribo_pct)
        if max_gini > 0:
            mask = mask & (gini_values <= max_gini)

        # 过滤日志
        filter_log_rows = []
        for i in range(n_before):
            fail_reasons = []
            if input_measurement == 'raw_counts' and lib_sizes[i] < min_counts:
                fail_reasons.append(f'library_size<{min_counts}')
            if input_measurement == 'raw_counts' and n_genes_detected[i] < min_genes:
                fail_reasons.append(f'n_genes<{min_genes}')
            if mt_filter_enabled and mt_pct[i] > max_mt_pct:
                fail_reasons.append(f'mt_pct>{max_mt_pct}')
            if ribo_pct[i] > max_ribo_pct:
                fail_reasons.append(f'ribo_pct>{max_ribo_pct}')
            if max_gini > 0 and gini_values[i] > max_gini:
                fail_reasons.append(f'gini>{max_gini}')
            filter_log_rows.append({
                'sample': sample_names[i],
                'group': groups[i],
                'replicate_group': replicate_groups[i],
                'passed': len(fail_reasons) == 0,
                quantity_key: float(lib_sizes[i]),
                gene_metric_key: int(n_genes_detected[i]),
                'mt_pct': round(float(mt_pct[i]), 2),
                'ribo_pct': round(float(ribo_pct[i]), 2),
                'gini': round(float(gini_values[i]), 4),
                'novelty': round(float(novelty_values[i]), 4),
                'fail_reasons': '; '.join(fail_reasons) if fail_reasons else '',
            })

        adata_filtered = adata[mask].copy()
        gini_filtered = gini_values[mask]
        n_after = adata_filtered.n_obs

        # 输出样本指标表和过滤日志
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        result_files = []

        metrics_df = pd.DataFrame(filter_log_rows)
        metrics_csv = os.path.join(results_dir, 'bulk_qc_sample_metrics.csv')
        metrics_df.to_csv(metrics_csv, index=False)

        filter_log_df = metrics_df[~metrics_df['passed']]
        filter_log_csv = os.path.join(results_dir, 'bulk_qc_filter_log.csv')
        filter_log_df.to_csv(filter_log_csv, index=False)

        result_files.append({'file_path': metrics_csv, 'file_type': 'csv', 'category': 'table', 'label': '样本 QC 指标'})
        result_files.append({'file_path': filter_log_csv, 'file_type': 'csv', 'category': 'table', 'label': '过滤日志'})

        # 基因层面过滤。连续表达量使用相同数值比较，但不再把阈值称为 count。
        genes_before_filter = adata_filtered.n_vars
        gene_filter_rows = []
        if min_sample_expr > 0:
            raw_filt = adata_filtered.X.toarray() if hasattr(adata_filtered.X, 'toarray') else np.asarray(adata_filtered.X)
            expr_count_per_gene = (raw_filt >= min_count_threshold).sum(axis=0)
            gene_mask = expr_count_per_gene >= min_sample_expr

            removed_indices = np.where(~gene_mask)[0]
            gene_filter_rows = [
                {
                    'gene': adata_filtered.var_names[idx],
                    'expressed_in_n_samples': int(expr_count_per_gene[idx]),
                    'action': 'removed',
                }
                for idx in removed_indices
            ]
            adata_filtered = adata_filtered[:, gene_mask].copy()

        # 管家基因稳定性检查
        housekeeping_genes = ['GAPDH', 'ACTB', 'B2M', 'HPRT1', 'TBP', 'UBC', 'YWHAZ', 'SDHA', 'HMBS', 'RPLP0']
        found_hk = [g for g in housekeeping_genes if g in adata_filtered.var_names]

        # 基因过滤日志
        if min_sample_expr > 0 and gene_filter_rows:
            gene_filter_csv = os.path.join(results_dir, 'bulk_qc_gene_filter.csv')
            pd.DataFrame(gene_filter_rows).to_csv(gene_filter_csv, index=False)
            result_files.append({'file_path': gene_filter_csv, 'file_type': 'csv', 'category': 'table', 'label': '基因过滤日志'})

        self.progress(60, "生成质控图表...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        from figure_engine import NatureFigureDirector, export_registered_figure
        director = NatureFigureDirector()
        nature_formats = ('svg', 'pdf', 'png')

        def _export_diagnostic(fig, stem, label, *, width='single', height_mm=None,
                               category='qc', plot_type='diagnostic'):
            spec = director.spec_from_params(
                plot_type, self.params, width=width, title=label,
            ).with_updates(formats=nature_formats, height_mm=height_mm)
            exported, report = export_registered_figure(
                fig, os.path.join(plots_dir, stem), spec,
                category=category, label=label,
                qa_path=os.path.join(results_dir, f'{stem}_nature_readiness.json'),
            )
            if not report.ready:
                self.progress(-1, f'{label} Nature readiness {report.score}/100；请查看 QA 报告。')
            plt.close(fig)
            return exported

        overview = [
            {'values': lib_sizes, 'title': quantity_label, 'ylabel': quantity_label},
            {'values': n_genes_detected, 'title': gene_metric_label, 'ylabel': gene_metric_label},
            {
                'values': mt_pct, 'title': mt_metric_label,
                'ylabel': 'MT expression %' if input_measurement != 'raw_counts' else 'MT%',
                # A true all-zero MT fraction is still scientifically useful;
                # render it rather than making the reader infer that the
                # metric was never calculated.  The annotation carries the
                # reason so an all-zero panel is never read as clean samples.
                'show_if_invariant': mt_evaluable,
                'annotation': mt_detection['figure_annotation'],
            },
            {'values': ribo_pct, 'title': ribo_metric_label, 'ylabel': 'RPL/RPS expression %' if input_measurement != 'raw_counts' else 'RPL/RPS %'},
            {'values': gini_values, 'title': 'Gini coefficient', 'ylabel': 'Gini'},
        ]
        invariant = []
        for item in overview:
            finite_values = np.asarray(item['values'], dtype=float)
            finite_values = finite_values[np.isfinite(finite_values)]
            if (finite_values.size == 0 or np.nanmax(finite_values) - np.nanmin(finite_values) <= 1e-12) and not item.get('show_if_invariant', False):
                invariant.append(item['title'])
        result_files.extend(_export_diagnostic(
            director.render(
                director.spec_from_params('diagnostic', self.params, width='double',
                                          title='Bulk RNA-seq QC overview').with_updates(
                                              extra={'kind': 'qc_overview'},
                                              formats=nature_formats, height_mm=130.0),
                {
                    'kind': 'qc_overview',
                    'sample_labels': sample_names,
                    'metrics': overview,
                    'library_size': lib_sizes,
                    'quantity_label': quantity_label,
                    'detected_genes': n_genes_detected,
                    'pass_colors': [('#4C78A8' if m else '#77808C') for m in mask],
                    'omitted_metrics': ', '.join(invariant),
                },
            ),
            'bulk_qc_overview', '质控总览', width='double', height_mm=130.0,
        ))

        # 保留经过筛选、但数值未变换的输入，用于后续模块输出。
        adata_raw_filtered = adata_filtered.copy()

        # 在独立副本上进行适合输入类型的 QC 相关性 / PCA 变换；不修改输出矩阵。
        self.progress(65, "标准化数据...")
        adata_normed = adata_filtered.copy()
        if input_measurement == 'raw_counts':
            sc.pp.normalize_total(adata_normed, target_sum=1e6)
            sc.pp.log1p(adata_normed)
            qc_transform = 'CPM (target 1e6) + log1p'
        elif input_measurement == 'continuous_expression':
            adata_normed.X = np.log2(np.maximum(adata_normed.X, 0) + 1)
            qc_transform = 'log2(max(expression, 0) + 1)'
        else:
            # Input provenance already declares a log-like transformation.
            # Applying log2 a second time compresses real sample distances.
            qc_transform = 'input log-transformed expression (no second log)'

        self.progress(70, "生成相关性热图...")
        corr_data = adata_normed.X if not hasattr(adata_normed.X, 'toarray') else adata_normed.X.toarray()
        corr_matrix = np.corrcoef(corr_data)
        sample_labels_corr = adata_normed.obs.index.tolist()

        # 分组条与样本排序均基于过滤后的样本；相关性原始矩阵不改动。
        filtered_groups = [groups[sample_names.index(s)] for s in sample_labels_corr]
        filtered_replicate_groups = [replicate_groups[sample_names.index(s)] for s in sample_labels_corr]
        unique_groups = sorted(set(filtered_groups))
        from figure_engine import NatureFigureDirector, export_registered_figure
        director = NatureFigureDirector()
        corr_spec = director.spec_from_params(
            'correlation', self.params, width='double',
            title='样本相关性热图 (Pearson)', correlation_method='pearson',
        )
        fig_corr = director.render(corr_spec, {
            'correlation_matrix': corr_matrix,
            'sample_labels': sample_labels_corr,
            'groups': filtered_groups,
        })
        corr_metadata = fig_corr._nature_correlation_metadata
        exported, corr_readiness = export_registered_figure(
            fig_corr, os.path.join(plots_dir, 'bulk_qc_corr'), corr_spec,
            category='heatmap', label='样本相关性热图',
            qa_path=os.path.join(results_dir, 'bulk_qc_corr_nature_readiness.json'),
        )
        result_files.extend(exported)
        if not corr_readiness.ready:
            self.progress(-1, f'Correlation Nature readiness {corr_readiness.score}/100；请查看 QA 报告。')
        plt.close(fig_corr)
        corr_pairs = correlation_pairwise_table(
            corr_matrix, sample_labels_corr, filtered_replicate_groups, method='pearson')
        corr_pairs = _add_factor_columns_to_correlation_pairs(corr_pairs, adata_normed)
        corr_pairs_csv = os.path.join(results_dir, 'bulk_qc_correlation_pairs.csv')
        corr_pairs.to_csv(corr_pairs_csv, index=False)
        corr_summary = _summarize_factor_aware_correlations(corr_pairs, method='pearson')
        corr_summary_csv = os.path.join(results_dir, 'bulk_qc_correlation_summary.csv')
        corr_summary.to_csv(corr_summary_csv, index=False)
        result_files.extend([
            {'file_path': corr_pairs_csv, 'file_type': 'csv', 'category': 'table',
             'label': '样本两两 Pearson 相关性'},
            {'file_path': corr_summary_csv, 'file_type': 'csv', 'category': 'table',
             'label': '样本相关性整体、重复组与因素摘要'},
        ])

        if n_before > n_after:
            fig_r = director.render(
                director.spec_from_params('diagnostic', self.params, width='single',
                                          title='Sample filtering').with_updates(
                                              extra={'kind': 'normalization_library'},
                                              formats=nature_formats, height_mm=60.0),
                {'kind': 'normalization_library', 'sample_labels': ['Before', 'After'],
                 'raw_values': [n_before, n_before], 'normalized_values': [n_after, n_after],
                 'raw_ylabel': 'Samples', 'normalized_ylabel': 'Samples'},
            )
            result_files.extend(_export_diagnostic(
                fig_r, 'bulk_qc_filter', '样本过滤结果', width='single', height_mm=60.0,
            ))

        self.progress(75, "PCA 离群检测...")
        outlier_samples = []
        outlier_diagnostics = pd.DataFrame()
        outlier_detection_details = {
            'method': 'within_replicate_group_agreement',
            'replicate_group_column': replicate_group_info['used'],
            'replicate_group_source': replicate_group_info['source'],
            'warning': '',
        }
        n_comps = min(10, n_after - 1, adata_normed.n_vars - 1)
        pca_preprocessing = None
        if n_comps >= 2:
            pca_preprocessing = run_bulk_sample_pca(adata_normed, n_comps=n_comps)
            pc = adata_normed.obsm['X_pca']
            filtered_sample_names = adata_normed.obs.index.tolist()
            filtered_groups_pca = [groups[sample_names.index(s)] for s in filtered_sample_names]
            filtered_replicate_groups_pca = [replicate_groups[sample_names.index(s)] for s in filtered_sample_names]
            unique_groups_pca = sorted(set(filtered_groups_pca))

            # Do not test samples against the global PCA centre: an entire,
            # coherent B/En stratum is a biological structure, not three
            # independent outliers.  A flag instead needs agreement between
            # within-replicate correlation and within-replicate PCA distance;
            # QC metrics are retained as an additional auditable signal.
            if detect_outliers:
                qc_metrics_for_outliers = pd.DataFrame({
                    quantity_key: adata_normed.obs['total_counts'].to_numpy(dtype=float),
                    gene_metric_key: adata_normed.obs['n_genes_by_counts'].to_numpy(dtype=float),
                    'mt_pct': adata_normed.obs['pct_counts_mt'].to_numpy(dtype=float),
                    'ribo_pct': adata_normed.obs['pct_counts_ribo'].to_numpy(dtype=float),
                    'gini': gini_filtered,
                }, index=filtered_sample_names)
                outlier_diagnostics, outlier_detection_details = _detect_within_replicate_outliers(
                    pc, corr_matrix, filtered_sample_names, filtered_replicate_groups_pca,
                    qc_metrics=qc_metrics_for_outliers,
                    replicate_group_column=replicate_group_info['used'],
                    replicate_group_source=replicate_group_info['source'],
                )
                outlier_samples = outlier_diagnostics.loc[
                    outlier_diagnostics['outlier_flag'], 'sample_id'
                ].astype(str).tolist()

                outlier_csv = os.path.join(results_dir, 'bulk_qc_within_group_outlier_diagnostics.csv')
                outlier_diagnostics.to_csv(outlier_csv, index=False)
                result_files.append({
                    'file_path': outlier_csv, 'file_type': 'csv', 'category': 'table',
                    'label': '重复组内离群诊断（相关性、PCA 与 QC 指标）',
                })

            pca_variance_qc = adata_normed.uns.get('pca', {}).get('variance_ratio', [])
            pca_color_groups, pca_marker_groups, pca_group_label, pca_batch_label = _resolve_qc_pca_encodings(
                adata_normed, filtered_groups_pca,
            )
            pca_spec = director.spec_from_params(
                'pca', self.params, width='single', title='QC-filtered sample PCA (replicate-aware flags)',
                show_legend=len(set(pca_color_groups)) <= 8,
            ).with_updates(formats=nature_formats, height_mm=82.0,
                           outlier_labels=tuple(outlier_samples))
            fig_pca = director.render(pca_spec, {
                'coordinates': pc[:, :2],
                'groups': pca_color_groups,
                'batches': pca_marker_groups,
                'group_label': pca_group_label,
                'batch_label': pca_batch_label,
                'samples': filtered_sample_names,
                'explained_variance': pca_variance_qc[:2],
            })
            result_files.extend(_export_diagnostic(
                fig_pca, 'bulk_qc_pca', '样本 PCA', width='single', height_mm=82.0, category='pca',
                plot_type='pca',
            ))

        # PCA 方差解释 elbow 图（使用 scanpy 存储的 variance_ratio）
        pca_variance = adata_normed.uns.get('pca', {}).get('variance_ratio', None)
        if pca_variance is None and n_comps >= 2:
            pca_var = np.var(pc, axis=0)
            pca_variance = pca_var / pca_var.sum()
        if pca_variance is not None and len(pca_variance) > 0:
            n_pcs = len(pca_variance)
            pc_labels = [f'PC{i+1}' for i in range(n_pcs)]
            fig_elbow = director.render(
                director.spec_from_params('diagnostic', self.params, width='single',
                                          title='PCA variance explained').with_updates(
                                              extra={'kind': 'variance'},
                                              formats=nature_formats, height_mm=70.0),
                {'kind': 'variance', 'variance': pca_variance, 'labels': pc_labels},
            )
            result_files.extend(_export_diagnostic(
                fig_elbow, 'bulk_qc_pca_elbow', 'PCA 方差解释', width='single', height_mm=70.0,
                category='pca',
            ))

        # These pairwise distances share samples, so they are descriptive QC
        # quantities only.  Split factor 2 where available instead of adding
        # an invalid independent-pairs p value to the title.
        if len(set(filtered_replicate_groups)) > 1 and n_after > 3:
            distance_sets = _descriptive_distance_sets(corr_pairs)
            box_data = [values for _, values in distance_sets if values]
            box_labels = [label for label, values in distance_sets if values]
            if box_data:
                fig_dist = director.render(
                    director.spec_from_params('diagnostic', self.params, width='single',
                                              title='Correlation distance (descriptive)').with_updates(
                                                  extra={'kind': 'boxplot'},
                                                  formats=nature_formats, height_mm=68.0),
                    {'kind': 'boxplot', 'groups': box_labels, 'values': box_data,
                     'ylabel': '1 - Pearson r'},
                )
                result_files.extend(_export_diagnostic(
                    fig_dist, 'bulk_qc_group_distance', '重复组与因素距离（描述性）', width='single', height_mm=68.0,
                ))

        # QC 指标散点矩阵 (Pairs Plot)
        obs_filtered = adata_filtered.obs
        pairs_groups = [groups[sample_names.index(s)] for s in obs_filtered.index.tolist()]
        unique_pg = sorted(set(pairs_groups))
        pg_color_map = ({g: NATURE_PALETTE[0] for g in unique_pg}
                        if len(unique_pg) > 3 else
                        {g: NATURE_PALETTE[i % len(NATURE_PALETTE)]
                         for i, g in enumerate(unique_pg)})
        from figure_engine.style import get_style
        pg_marker_map = get_style('nature').batch_markers(unique_pg)
        pg_colors = [pg_color_map[g] for g in pairs_groups]
        pair_values = [
            obs_filtered['total_counts'].values,
            obs_filtered['n_genes_by_counts'].values,
            obs_filtered['pct_counts_mt'].values if 'pct_counts_mt' in obs_filtered.columns else np.zeros(n_after),
            obs_filtered['pct_counts_ribo'].values if 'pct_counts_ribo' in obs_filtered.columns else np.zeros(n_after),
        ]
        pair_labels = [quantity_label, gene_metric_label,
                       'MT expression %' if input_measurement != 'raw_counts' else 'MT%',
                       'RPL/RPS expression %' if input_measurement != 'raw_counts' else 'Ribo%']
        fig_pairs = director.render(
            director.spec_from_params('diagnostic', self.params, width='double',
                                      title='QC metric pairs').with_updates(
                                          extra={'kind': 'pairs'}, formats=nature_formats,
                                          height_mm=112.0),
            {'kind': 'pairs', 'labels': pair_labels, 'values': pair_values,
             'groups': pairs_groups, 'group_colors': pg_color_map,
             'group_markers': pg_marker_map},
        )
        result_files.extend(_export_diagnostic(
            fig_pairs, 'bulk_qc_pairs_plot', 'QC 指标散点矩阵', width='double', height_mm=112.0,
        ))

        # n=3 replicate groups are shown as individual points with median/IQR,
        # never as a smooth density that suggests unobserved observations.
        if len(unique_groups) > 1:
            violin_data = []
            for s in obs_filtered.index.tolist():
                g = groups[sample_names.index(s)]
                violin_data.append({
                    'group': g,
                    'MT expression %' if input_measurement != 'raw_counts' else 'MT%': float(obs_filtered.loc[s, 'pct_counts_mt']) if 'pct_counts_mt' in obs_filtered.columns else 0,
                    'RPL/RPS expression %' if input_measurement != 'raw_counts' else 'Ribo%': float(obs_filtered.loc[s, 'pct_counts_ribo']) if 'pct_counts_ribo' in obs_filtered.columns else 0,
                    quantity_label: float(obs_filtered.loc[s, 'total_counts']),
                    gene_metric_label: float(obs_filtered.loc[s, 'n_genes_by_counts']),
                })
            violin_df = pd.DataFrame(violin_data)
            violin_metrics = []
            for metric in (
                'MT expression %' if input_measurement != 'raw_counts' else 'MT%',
                'RPL/RPS expression %' if input_measurement != 'raw_counts' else 'Ribo%',
                quantity_label, gene_metric_label,
            ):
                vals_by_group = [violin_df[violin_df['group'] == g][metric].values
                                 for g in unique_groups]
                if any((lambda finite: finite.size > 0 and
                        np.nanmax(finite) - np.nanmin(finite) > 1e-12)
                       (np.asarray(v, dtype=float)[np.isfinite(v)]) for v in vals_by_group):
                    violin_metrics.append({'label': metric, 'values': vals_by_group})
            if violin_metrics:
                fig_violin = director.render(
                    director.spec_from_params('diagnostic', self.params, width='double',
                                              title='QC metrics by group (individual samples)').with_updates(
                                                  extra={'kind': 'violin'}, formats=nature_formats,
                                                  height_mm=90.0),
                    {'kind': 'violin', 'groups': unique_groups, 'metrics': violin_metrics},
                )
                result_files.extend(_export_diagnostic(
                    fig_violin, 'bulk_qc_violin_by_group', '各组 QC 指标（样本点与中位数）', width='double', height_mm=90.0,
                ))

        # 管家基因稳定性热图
        if found_hk:
            # Reuse the exact QC transform; continuous/log inputs must never
            # be reinterpreted as raw counts for this display.
            norm_hk = adata_normed[:, found_hk].copy()
            hk_data = norm_hk.X if not hasattr(norm_hk.X, 'toarray') else norm_hk.X.toarray()
            hk_genes = norm_hk.var_names.tolist()
            hk_samples = norm_hk.obs.index.tolist()
            hk_cv = {}
            for j, g in enumerate(hk_genes):
                vals = hk_data[:, j]
                cv = float(np.std(vals) / (np.mean(vals) + 1e-10))
                hk_cv[g] = cv
            cv_labels = [f'{g} (CV={hk_cv[g]:.2f})' for g in hk_genes]
            hk_spec = director.spec_from_params(
                'heatmap', self.params, width='double', title='Housekeeping gene expression consistency',
                zscore='row', row_cluster=True, col_cluster=True,
            ).with_updates(formats=nature_formats, height_mm=100.0, max_row_labels=20,
                           max_col_labels=18)
            fig_hk = director.render(hk_spec, {
                'matrix': hk_data.T,
                'gene_labels': cv_labels,
                'sample_labels': hk_samples,
                # Auto-grouping can create one legend entry per sample.  A
                # long annotation legend is not informative and collides with
                # the heatmap/colorbar at publication size, so only retain it
                # when it remains compact.
                'annotations': ({'Group': [groups[sample_names.index(s)] for s in hk_samples]}
                                if len(set(groups[sample_names.index(s)] for s in hk_samples)) <= 6
                                else {}),
                'colorbar_label': 'Gene-wise z-score (relative expression)',
            })
            result_files.extend(_export_diagnostic(
                fig_hk, 'bulk_qc_housekeeping', '管家基因表达一致性', width='double', height_mm=100.0,
                category='heatmap', plot_type='heatmap',
            ))

        self.progress(90, "保存输出...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_qc_output.h5ad')
        adata_raw_filtered.write_h5ad(output_path)

        removed_samples = [r['sample'] for r in filter_log_rows if not r['passed']]
        removed_reasons = {r['sample']: r['fail_reasons'] for r in filter_log_rows if not r['passed']}

        self.progress(100, "完成")
        summary = {
            'samples_before': n_before,
            'samples_after': n_after,
            'samples_removed': n_before - n_after,
            'removed_samples': removed_samples,
            'removed_reasons': removed_reasons,
            'genes_before': genes_before_filter,
            'genes_after': adata_filtered.n_vars,
            'genes_removed': genes_before_filter - adata_filtered.n_vars,
            'gene_expression_filter': {
                'applied_in_qc': min_sample_expr > 0,
                'minimum_samples': int(min_sample_expr),
                'minimum_value': float(min_count_threshold),
                'value_scale': 'raw_counts' if input_measurement == 'raw_counts' else input_measurement,
                'n_genes_removed_in_qc': int(genes_before_filter - adata_filtered.n_vars),
                'required_before_deg': True,
                'message': (
                    'QC 阶段未启用低表达基因过滤（min_sample_expr=0）；'
                    'genes_removed=0 是预期结果，不代表无需过滤。'
                    '进入 Bulk DEG（DESeq2/edgeR/limma）前必须经过 Bulk normalization，'
                    '其默认按 CPM≥1 且至少 3 个样本表达过滤低表达基因，并记录实际过滤数量。'
                    if min_sample_expr == 0 else
                    '已在 QC 阶段应用低表达基因过滤。'
                ),
            },
            'outlier_samples': outlier_samples if detect_outliers else [],
            'outlier_detection': outlier_detection_details,
            'filter_strategy': filter_strategy,
            'input_measurement': input_measurement,
            'input_measurement_source': measurement_info['source'],
            'input_measurement_confidence': measurement_info['confidence'],
            'qc_transform_for_correlation_and_pca': qc_transform,
            'pca_preprocessing': pca_preprocessing,
            'count_qc_hard_filters_applied': input_measurement == 'raw_counts',
            'mitochondrial_qc': {
                key: value for key, value in mt_detection.items() if key != 'mask'
            },
            'median_mt_pct': (
                round(float(np.nanmedian(adata_filtered.obs['pct_counts_mt'])), 2)
                if mt_evaluable and adata_filtered.n_obs else None
            ),
            'median_genes': int(np.median(adata_filtered.obs['n_genes_by_counts'])),
            'median_ribo_pct': round(float(np.median(adata_filtered.obs['pct_counts_ribo'])), 2) if 'pct_counts_ribo' in adata_filtered.obs.columns else 0,
            'median_gini': round(float(np.median(gini_filtered)), 4),
            'correlation_medians': _correlation_medians_for_manifest(corr_summary),
            'correlation_off_diagonal': corr_metadata['off_diagonal'],
            'correlation_display_range': [
                corr_metadata['display_vmin'], corr_metadata['display_vmax'],
            ],
        }
        if pca_variance is not None and len(pca_variance) >= 2:
            summary['pc1_variance_pct'] = round(float(pca_variance[0] * 100), 2)
            summary['pc2_variance_pct'] = round(float(pca_variance[1] * 100), 2)
        if input_measurement == 'raw_counts':
            summary['median_library_size'] = int(np.median(adata_filtered.obs['total_counts']))
            summary['median_detected_genes'] = int(np.median(adata_filtered.obs['n_genes_by_counts']))
        else:
            summary[f'median_{quantity_key}'] = round(
                float(np.median(adata_filtered.obs['total_counts'])), 4,
            )
            summary['median_expressed_genes'] = int(np.median(adata_filtered.obs['n_genes_by_counts']))

        # Surface actionable QC caveats in one place so a manifest reader never
        # has to infer them from a null threshold or an unchanged gene count.
        qc_warnings = []
        if mt_detection['message']:
            qc_warnings.append(mt_detection['message'])
        if min_sample_expr == 0:
            qc_warnings.append(summary['gene_expression_filter']['message'])
        summary['warnings'] = qc_warnings

        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': summary,
        }


# --- 辅助函数 ---


def _normalise_annotation_column_name(column):
    """Normalise a ``var`` column name only for known annotation aliases."""
    return re.sub(r'[^a-z0-9]+', '_', str(column).strip().lower()).strip('_')


def _text_series(values, index):
    """Return a nullable-safe, stripped text Series for AnnData annotations."""
    return pd.Series(values, index=index, dtype='string').fillna('').str.strip()


def _annotation_columns(adata, accepted_names):
    """Find all matching annotation columns while preserving their order."""
    return [
        column for column in adata.var.columns
        if _normalise_annotation_column_name(column) in accepted_names
    ]


def _preferred_gene_symbols(adata):
    """Return the best available symbol/name per gene, falling back to var_names.

    This is intentionally a per-row fallback: imported matrices can have an
    incomplete ``gene_name`` column, and a valid secondary symbol column should
    still identify the remaining mitochondrial or ribosomal genes.
    """
    values = _text_series(adata.var_names, adata.var.index)
    for column in _annotation_columns(adata, _GENE_SYMBOL_COLUMNS):
        candidate = _text_series(adata.var[column], adata.var.index)
        values = candidate.where(candidate.ne(''), values)
    return values


def _ensembl_without_version(values):
    """Normalise Ensembl IDs without altering non-Ensembl identifiers."""
    return values.str.replace(r'^(?:gene:)?(ENSG\d+)(?:\.\d+)?$', r'\1', regex=True).str.upper()


def _detect_mitochondrial_genes(adata):
    """Identify MT features from symbols, chromosome metadata, or human IDs.

    The return value deliberately contains source-level counts.  A value such
    as ``pct_counts_mt=0`` is not evidence that an Ensembl-only input was
    handled correctly; the caller needs to know exactly how the MT mask was
    derived and whether the threshold was allowed to act.
    """
    index = adata.var.index
    mask = pd.Series(False, index=index, dtype=bool)
    source_counts = {}

    symbols = _preferred_gene_symbols(adata)
    symbol_mask = symbols.str.upper().str.startswith('MT-')
    if bool(symbol_mask.any()):
        mask |= symbol_mask
        source_counts['gene_symbol'] = int(symbol_mask.sum())

    chromosome_mask = pd.Series(False, index=index, dtype=bool)
    for column in _annotation_columns(adata, _CHROMOSOME_COLUMNS):
        chromosome = _text_series(adata.var[column], index).str.upper().str.replace(r'\s+', '', regex=True)
        hits = chromosome.isin(_MITOCHONDRIAL_CHROMOSOMES)
        if bool(hits.any()):
            chromosome_mask |= hits
    if bool(chromosome_mask.any()):
        mask |= chromosome_mask
        source_counts['chromosome'] = int(chromosome_mask.sum())

    ensembl_mask = pd.Series(False, index=index, dtype=bool)
    # ``var_names`` is included even when ``gene_id`` exists: many imported
    # h5ad files keep the raw stable ID only in the index.
    id_sources = [_text_series(adata.var_names, index)]
    id_sources.extend(_text_series(adata.var[column], index)
                      for column in _annotation_columns(adata, _GENE_ID_COLUMNS))
    for identifiers in id_sources:
        hits = _ensembl_without_version(identifiers).isin(_HUMAN_MITOCHONDRIAL_ENSEMBL_IDS)
        ensembl_mask |= hits
    if bool(ensembl_mask.any()):
        mask |= ensembl_mask
        source_counts['human_ensembl_id'] = int(ensembl_mask.sum())

    return {
        'mask': mask.to_numpy(dtype=bool),
        'n_mitochondrial_genes': int(mask.sum()),
        'sources': source_counts,
        'gene_symbol_columns': _annotation_columns(adata, _GENE_SYMBOL_COLUMNS),
        'chromosome_columns': _annotation_columns(adata, _CHROMOSOME_COLUMNS),
        'gene_id_columns': _annotation_columns(adata, _GENE_ID_COLUMNS),
    }


def _gini(values):
    """计算 Gini 系数。values 为原始 count 数组（非负）。"""
    vals = np.sort(values[values >= 0])
    n = len(vals)
    if n == 0:
        return 0.0
    if np.sum(vals) == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return (2.0 * np.sum(index * vals) / (n * np.sum(vals))) - (n + 1) / n


def _infer_groups(sample_names):
    """从样本名推断推荐分组，优先保留多因素联合组。"""
    from modules.io_utils import infer_sample_group_candidates

    candidates = infer_sample_group_candidates(sample_names)
    if candidates:
        mapping = candidates[0]['mapping']
        return [mapping[str(name)] for name in sample_names]

    groups = []
    for name in sample_names:
        name = str(name)
        assigned = False
        for sep in ['-', '_']:
            if sep in name:
                groups.append(name.split(sep)[0])
                assigned = True
                break
        if not assigned:
            groups.append(name)
    return groups


def _resolve_replicate_groups(adata, requested, fallback_groups):
    """Resolve the biological replicate group used for QC outlier checks.

    A combined auto group (for example ``Ctr_B``) is deliberately preferred
    over a treatment-only plotting group (``Ctr``).  This prevents a coherent
    second experimental factor from being measured against the wrong centre.
    """
    requested = str(requested or '').strip()
    if requested == '_auto_group_':
        requested = '_auto_group'
    if requested and requested in adata.obs.columns:
        return adata.obs[requested].astype(str).tolist(), {
            'requested': requested, 'used': requested, 'source': 'requested_obs_column',
        }
    if '_auto_group' in adata.obs.columns:
        return adata.obs['_auto_group'].astype(str).tolist(), {
            'requested': requested, 'used': '_auto_group', 'source': 'sample_name_combined_group',
        }
    return [str(value) for value in fallback_groups], {
        'requested': requested, 'used': 'group_column', 'source': 'display_group_fallback',
    }


def _add_factor_columns_to_correlation_pairs(pairwise_table, adata):
    """Attach reproducible auto-factor labels to each exported sample pair."""
    table = pairwise_table.copy()
    if table.empty:
        return table
    for number in range(1, 5):
        column = f'_auto_factor{number}'
        if column not in adata.obs.columns:
            continue
        values = adata.obs[column].astype(str).to_dict()
        name = f'factor{number}'
        table[f'{name}_1'] = table['sample_1'].map(values).fillna('')
        table[f'{name}_2'] = table['sample_2'].map(values).fillna('')
        table[f'{name}_relation'] = np.where(
            table[f'{name}_1'] == table[f'{name}_2'], 'same', 'different',
        )
    # Preserve an explicit alias for downstream use even when a project has
    # no parseable multi-factor sample naming scheme.
    table['replicate_group_1'] = table['group_1']
    table['replicate_group_2'] = table['group_2']
    return table


def _summarize_factor_aware_correlations(pairwise_table, method='pearson'):
    """Summarize overall, repeat-group and available factor-level correlations."""
    from modules.native_figures import summarize_correlation_pairs

    base = summarize_correlation_pairs(pairwise_table, method=method).copy()
    for column in ('summary_scope', 'factor', 'level_1', 'level_2'):
        base[column] = ''
    base['summary_scope'] = np.where(
        base['pair_set'] == 'all_pairs', 'overall', 'replicate_group',
    )
    value_column = f'{method}_r'
    extra_rows = []
    for number in range(1, 5):
        first, second = f'factor{number}_1', f'factor{number}_2'
        if first not in pairwise_table.columns or second not in pairwise_table.columns:
            continue
        levels = sorted({str(value) for value in pd.concat([
            pairwise_table[first], pairwise_table[second],
        ], ignore_index=True) if str(value)})
        for level in levels:
            subset = pairwise_table[(pairwise_table[first] == level) & (pairwise_table[second] == level)]
            values = pd.to_numeric(subset[value_column], errors='coerce').dropna()
            extra_rows.append({
                'pair_set': f'factor{number}_within:{level}',
                'n_pairs': int(values.shape[0]),
                'mean_r': float(values.mean()) if not values.empty else np.nan,
                'median_r': float(values.median()) if not values.empty else np.nan,
                'min_r': float(values.min()) if not values.empty else np.nan,
                'max_r': float(values.max()) if not values.empty else np.nan,
                'summary_scope': 'within_factor_level', 'factor': f'factor{number}',
                'level_1': level, 'level_2': level,
            })
        for index, first_level in enumerate(levels):
            for second_level in levels[index + 1:]:
                subset = pairwise_table[
                    ((pairwise_table[first] == first_level) & (pairwise_table[second] == second_level)) |
                    ((pairwise_table[first] == second_level) & (pairwise_table[second] == first_level))
                ]
                values = pd.to_numeric(subset[value_column], errors='coerce').dropna()
                extra_rows.append({
                    'pair_set': f'factor{number}_between:{first_level}_vs_{second_level}',
                    'n_pairs': int(values.shape[0]),
                    'mean_r': float(values.mean()) if not values.empty else np.nan,
                    'median_r': float(values.median()) if not values.empty else np.nan,
                    'min_r': float(values.min()) if not values.empty else np.nan,
                    'max_r': float(values.max()) if not values.empty else np.nan,
                    'summary_scope': 'between_factor_levels', 'factor': f'factor{number}',
                    'level_1': first_level, 'level_2': second_level,
                })
    if extra_rows:
        return pd.concat([base, pd.DataFrame(extra_rows)], ignore_index=True)
    return base


def _correlation_medians_for_manifest(summary_table):
    """Create a compact, JSON-safe correlation summary without hiding strata."""
    result = {}
    for _, row in summary_table.iterrows():
        value = pd.to_numeric(pd.Series([row.get('median_r')]), errors='coerce').iloc[0]
        result[str(row.get('pair_set'))] = float(value) if pd.notna(value) else None
    return result


def _robust_tail_flags(values, *, direction):
    """Return conservative one-sided robust flags and a transparent threshold."""
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    threshold = np.nan
    flags = np.zeros(values.shape[0], dtype=bool)
    if finite.sum() < 4:
        return flags, threshold
    valid = values[finite]
    median = float(np.median(valid))
    mad = float(np.median(np.abs(valid - median)))
    scale = 1.4826 * mad
    if scale <= 1e-12:
        q25, q75 = np.percentile(valid, [25, 75])
        scale = float((q75 - q25) / 1.349)
    if scale <= 1e-12:
        return flags, threshold
    if direction == 'low':
        threshold = median - 3.0 * scale
        flags[finite] = valid < threshold
    elif direction == 'high':
        threshold = median + 3.0 * scale
        flags[finite] = valid > threshold
    else:
        raise ValueError("direction must be 'low' or 'high'")
    return flags, float(threshold)


def _finite_number_or_none(value):
    """Return a JSON-safe float for manifest fields, or ``None`` if unavailable."""
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if np.isfinite(numeric) else None


def _invariant_metric_note(values):
    """Explain why a QC metric cannot yield a robust deviation threshold.

    An all-zero or single-valued metric (for example MT% when the input has no
    mitochondrial reads at all) produces ``threshold=None``.  Returning the
    reason keeps that ``null`` auditable instead of looking like a missing
    calculation.
    """
    finite = np.asarray(values, dtype=float).ravel()
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 'no_finite_values'
    if float(np.max(finite) - np.min(finite)) <= 1e-12:
        return 'invariant_metric_no_threshold'
    return ''


def _unique_within_group_extreme(values, groups, *, direction):
    """Identify a single directional extreme inside each replicate group.

    With exactly three replicates, the two normal samples both see the bad
    sample as a peer.  Requiring a *unique* local extreme prevents those two
    normal samples from inheriting the same warning merely because they share
    that bad peer.
    """
    values = np.asarray(values, dtype=float)
    result = np.zeros(values.shape[0], dtype=bool)
    for group in dict.fromkeys(groups):
        indices = np.asarray([index for index, value in enumerate(groups) if value == group], dtype=int)
        finite_indices = indices[np.isfinite(values[indices])]
        if finite_indices.size < 3:
            continue
        local = values[finite_indices]
        extreme_value = np.min(local) if direction == 'low' else np.max(local)
        candidates = finite_indices[np.isclose(local, extreme_value, rtol=1e-8, atol=1e-12)]
        if candidates.size != 1:
            continue
        other_values = local[~np.isclose(local, extreme_value, rtol=1e-8, atol=1e-12)]
        if other_values.size and ((direction == 'low' and extreme_value < np.min(other_values)) or
                                  (direction == 'high' and extreme_value > np.max(other_values))):
            result[candidates[0]] = True
    return result


def _detect_within_replicate_outliers(pca_coords, correlation_matrix, sample_names,
                                      replicate_groups, *, qc_metrics=None,
                                      replicate_group_column='_auto_group',
                                      replicate_group_source='sample_name_combined_group'):
    """Flag individual samples only relative to their own biological replicates.

    The diagnostic uses the median correlation to other members of the same
    replicate group, median standardized PCA distance to those peers, and
    group-centred QC metric deviations.  A final automatic flag requires at
    least two independent signals, so a whole, well-replicated experimental
    stratum cannot become a global-PCA false positive.
    """
    coords = np.asarray(pca_coords, dtype=float)
    corr = np.asarray(correlation_matrix, dtype=float)
    samples = [str(value) for value in sample_names]
    groups = [str(value) for value in replicate_groups]
    n_samples = len(samples)
    if coords.ndim != 2 or coords.shape[0] != n_samples:
        raise ValueError('PCA coordinates must have one row per sample.')
    if corr.shape != (n_samples, n_samples):
        raise ValueError('Correlation matrix must be square and match samples.')
    if len(groups) != n_samples:
        raise ValueError('Replicate groups must have one entry per sample.')

    n_components = min(10, coords.shape[1])
    usable = coords[:, :n_components]
    std = np.nanstd(usable, axis=0, ddof=1)
    std[~np.isfinite(std) | (std <= 1e-12)] = 1.0
    scaled = usable / std
    peer_count = np.zeros(n_samples, dtype=int)
    within_corr = np.full(n_samples, np.nan, dtype=float)
    within_pca_distance = np.full(n_samples, np.nan, dtype=float)
    group_metric_deviation = np.full(n_samples, np.nan, dtype=float)
    metric_names = []

    qc_frame = None
    if qc_metrics is not None:
        qc_frame = pd.DataFrame(qc_metrics).reindex(samples)
        qc_frame = qc_frame.select_dtypes(include=[np.number])
        metric_names = qc_frame.columns.tolist()

    for group in dict.fromkeys(groups):
        indices = np.asarray([index for index, value in enumerate(groups) if value == group], dtype=int)
        # With fewer than three samples, every sample has indistinguishable
        # one-peer distances.  Export the values but make no individual flag.
        if indices.size < 3:
            continue
        for index in indices:
            peers = indices[indices != index]
            peer_count[index] = peers.size
            correlations = corr[index, peers]
            correlations = correlations[np.isfinite(correlations)]
            if correlations.size:
                within_corr[index] = float(np.median(correlations))
            deltas = scaled[peers] - scaled[index]
            distances = np.sqrt(np.nansum(deltas ** 2, axis=1))
            distances = distances[np.isfinite(distances)]
            if distances.size:
                within_pca_distance[index] = float(np.median(distances))
        if qc_frame is not None and metric_names:
            values = qc_frame.iloc[indices]
            centres = values.median(axis=0, skipna=True)
            deviations = (values - centres).abs()
            # Different QC quantities are evaluated separately below.  The
            # row score counts how many group-centred metrics are unusual.
            for local_index, index in enumerate(indices):
                group_metric_deviation[index] = float(
                    np.nanmedian(deviations.iloc[local_index].to_numpy(dtype=float))
                )

    correlation_flag, correlation_threshold = _robust_tail_flags(within_corr, direction='low')
    correlation_flag &= _unique_within_group_extreme(within_corr, groups, direction='low')
    pca_flag, pca_threshold = _robust_tail_flags(within_pca_distance, direction='high')
    pca_flag &= _unique_within_group_extreme(within_pca_distance, groups, direction='high')
    qc_metric_flags = np.zeros(n_samples, dtype=bool)
    qc_thresholds = {}
    qc_threshold_notes = {}
    if qc_frame is not None and metric_names:
        for metric in metric_names:
            deviations = np.full(n_samples, np.nan, dtype=float)
            for group in dict.fromkeys(groups):
                indices = np.asarray([index for index, value in enumerate(groups) if value == group], dtype=int)
                if indices.size < 3:
                    continue
                values = qc_frame.iloc[indices][metric].to_numpy(dtype=float)
                deviations[indices] = np.abs(values - np.nanmedian(values))
            flags, threshold = _robust_tail_flags(deviations, direction='high')
            flags &= _unique_within_group_extreme(deviations, groups, direction='high')
            qc_metric_flags |= flags
            qc_thresholds[metric] = threshold
            note = _invariant_metric_note(qc_frame[metric].to_numpy(dtype=float))
            if note:
                qc_threshold_notes[metric] = note

    signal_count = correlation_flag.astype(int) + pca_flag.astype(int) + qc_metric_flags.astype(int)
    outlier_flag = (peer_count >= 2) & (signal_count >= 2)
    reasons = []
    for index in range(n_samples):
        labels = []
        if correlation_flag[index]:
            labels.append('low_within_group_correlation')
        if pca_flag[index]:
            labels.append('high_within_group_pca_distance')
        if qc_metric_flags[index]:
            labels.append('group_centred_qc_metric_deviation')
        reasons.append(';'.join(labels))
    table = pd.DataFrame({
        'sample_id': samples,
        'replicate_group': groups,
        'peer_count': peer_count,
        'within_group_median_correlation': within_corr,
        'within_group_median_pca_distance': within_pca_distance,
        'within_group_qc_deviation_median': group_metric_deviation,
        'low_correlation_flag': correlation_flag,
        'high_pca_distance_flag': pca_flag,
        'qc_metric_deviation_flag': qc_metric_flags,
        'outlier_signal_count': signal_count,
        'outlier_flag': outlier_flag,
        'outlier_reasons': reasons,
    })
    if qc_frame is not None:
        for metric in metric_names:
            table[f'qc_{metric}'] = qc_frame[metric].to_numpy(dtype=float)
    details = {
        'method': 'within_replicate_group_agreement',
        'replicate_group_column': replicate_group_column,
        'replicate_group_source': replicate_group_source,
        'pca_components_used': int(n_components),
        'minimum_peers_required': 2,
        'flag_rule': 'at_least_two_of_within_group_correlation_pca_distance_qc_deviation',
        'low_correlation_threshold': _finite_number_or_none(correlation_threshold),
        'high_pca_distance_threshold': _finite_number_or_none(pca_threshold),
        'qc_metric_deviation_thresholds': {
            metric: _finite_number_or_none(threshold)
            for metric, threshold in qc_thresholds.items()
        },
        # A null threshold means "no usable spread", not "metric was never
        # computed".  Record the reason so an invariant metric such as an
        # all-zero MT% is not read as a platform failure.
        'qc_metric_deviation_threshold_notes': qc_threshold_notes,
        'invariant_qc_metrics': sorted(qc_threshold_notes),
        'warning': ('Groups with fewer than 3 retained samples are exported but not automatically '
                    'flagged because they lack two independent peers.'),
    }
    return table, details


def _resolve_qc_pca_encodings(adata, fallback_groups):
    """Use formal factor metadata for QC PCA colour/marker encodings when present."""
    if '_auto_factor1' in adata.obs.columns:
        colors = adata.obs['_auto_factor1'].astype(str).tolist()
        color_label = 'Factor 1'
    else:
        colors = [str(value) for value in fallback_groups]
        color_label = 'Group'
    if '_auto_factor2' in adata.obs.columns:
        markers = adata.obs['_auto_factor2'].astype(str).tolist()
        marker_label = 'Factor 2'
    else:
        markers = None
        marker_label = ''
    return colors, markers, color_label, marker_label


def _descriptive_distance_sets(correlation_pairs):
    """Create non-inferential distance panels from an auditable pair table."""
    if correlation_pairs.empty:
        return []
    values = pd.to_numeric(correlation_pairs['pearson_r'], errors='coerce')
    distances = 1.0 - values
    sets = [('Within replicate', distances[correlation_pairs['pair_type'] == 'within_group'].dropna().tolist())]
    if 'factor2_relation' in correlation_pairs.columns:
        between = correlation_pairs['pair_type'] == 'between_group'
        factor_levels = sorted({str(value) for value in pd.concat([
            correlation_pairs['factor2_1'], correlation_pairs['factor2_2'],
        ], ignore_index=True) if str(value)})
        factor_contrast = (f'{factor_levels[0]} vs {factor_levels[1]}'
                           if len(factor_levels) == 2 else 'Between factor 2\nlevels')
        sets.append((
            'Between groups,\nsame factor 2',
            distances[between & (correlation_pairs['factor2_relation'] == 'same')].dropna().tolist(),
        ))
        sets.append((
            factor_contrast,
            distances[between & (correlation_pairs['factor2_relation'] == 'different')].dropna().tolist(),
        ))
    else:
        sets.append(('Between replicate groups', distances[correlation_pairs['pair_type'] == 'between_group'].dropna().tolist()))
    return sets


def _detect_outliers_mahal(pca_coords, sample_names, *, return_details=False):
    """Detect PCA outliers across the leading usable PCs.

    A PC3/PC4-only anomaly is invisible to the former PC1--PC3 calculation
    when its distance happens to be diluted by the first two axes.  Use up to
    ten components while keeping dimensionality below half the sample count;
    this avoids an unstable singular covariance matrix in small Bulk studies.
    ``MinCovDet`` supplies a robust centre/covariance when there are enough
    samples, with a deterministic empirical fallback for very small cohorts.
    """
    coords_all = np.asarray(pca_coords, dtype=float)
    n_samples = coords_all.shape[0] if coords_all.ndim == 2 else 0
    details = {
        'distances': np.full(n_samples, np.nan, dtype=float),
        'threshold': np.nan,
        'n_components': 0,
        'method': 'not_available',
    }
    if n_samples < 4 or coords_all.ndim != 2 or coords_all.shape[1] < 2:
        return ([], details) if return_details else []

    n_components = min(10, coords_all.shape[1], max(2, n_samples // 2))
    coords = coords_all[:, :n_components]
    finite = np.isfinite(coords).all(axis=1)
    if finite.sum() < 4:
        return ([], details) if return_details else []
    valid_coords = coords[finite]

    try:
        # Robust covariance needs appreciably more observations than features.
        if valid_coords.shape[0] >= valid_coords.shape[1] * 2 + 2:
            from sklearn.covariance import MinCovDet

            # A high support fraction avoids treating ordinary tail points as
            # contaminants in the small cohorts typical of this platform,
            # while still isolating a genuinely extreme sample.
            estimator = MinCovDet(support_fraction=0.9, random_state=42).fit(valid_coords)
            center = estimator.location_
            covariance = estimator.covariance_
            details['method'] = 'robust_mahalanobis_mcd'
        else:
            center = valid_coords.mean(axis=0)
            covariance = np.cov(valid_coords, rowvar=False)
            details['method'] = 'mahalanobis_empirical_small_n'
        covariance = np.atleast_2d(covariance)
        covariance_inv = np.linalg.pinv(covariance)
        differences = valid_coords - center
        valid_distances = np.sqrt(np.einsum(
            'ij,jk,ik->i', differences, covariance_inv, differences,
        ))
    except (ValueError, np.linalg.LinAlgError):
        return ([], details) if return_details else []

    distances = np.full(n_samples, np.nan, dtype=float)
    distances[finite] = valid_distances
    med = float(np.median(valid_distances))
    mad = float(np.median(np.abs(valid_distances - med)))
    threshold = np.nan
    if mad >= 1e-10:
        threshold = med + 3 * 1.4826 * mad
        outliers = [str(sample_names[i]) for i, distance in enumerate(distances)
                    if np.isfinite(distance) and distance > threshold]
    else:
        outliers = []
    details.update({
        'distances': distances,
        'threshold': threshold,
        'n_components': n_components,
    })
    return (outliers, details) if return_details else outliers
