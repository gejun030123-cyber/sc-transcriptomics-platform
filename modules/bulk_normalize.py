import os
import logging
import numpy as np
import pandas as pd
from modules.base import BaseAnalysis

logger = logging.getLogger(__name__)

class BulkNormalizeAnalysis(BaseAnalysis):
    MODULE_NAME = "bulk_normalize"
    DISPLAY_NAME = "Bulk 数据标准化"
    DESCRIPTION = "计数矩阵标准化：DESeq2 size factors、CPM、分位数标准化"
    INPUT_REQUIRES = []

    def validate_input(self, adata):
        return None

    def run(self, input_path):
        import scanpy as sc
        import matplotlib.pyplot as plt
        from modules.native_figures import bar_figure
        from modules.figure_style import NATURE_PALETTE, NATURE_TEXT, NATURE_GRID

        self.progress(5, "加载数据...")
        from modules.io_utils import read_expression_matrix, resolve_expression_measurement, run_bulk_sample_pca
        adata = read_expression_matrix(input_path)

        method = str(self.params.get('method', 'deseq2')).strip().lower()
        self.progress(20, f"标准化方法: {method}...")

        # 前置过滤参数
        min_expr_value = float(self.params.get('min_expr_value', 1))
        min_expr_samples = int(self.params.get('min_expr_samples', 3))
        max_zero_pct = float(self.params.get('max_zero_pct', 0))

        input_measurement, measurement_info = resolve_expression_measurement(
            adata, input_path, self.params.get('input_measurement', 'auto'),
        )
        # Persist source-scale provenance through the QC/normalization handoff.
        # ``normalization.is_log_transformed`` below remains the authoritative
        # description of the *output* scale.
        adata.uns['input_measurement'] = input_measurement
        adata.uns['input_measurement_provenance'] = measurement_info
        count_only_methods = {'deseq2', 'tmm', 'cpm', 'vst', 'rlog'}
        if input_measurement != 'raw_counts' and method in count_only_methods:
            raise ValueError(
                f"检测到 {input_measurement}，{method} 仅适用于原始整数 counts。"
                "线性 FPKM/TPM 请使用 log2；已经 log 变换的输入请选择 none。"
            )
        if method == 'log2' and input_measurement == 'raw_counts':
            raise ValueError(
                'raw_counts 不能只做 log2(x+1)。请使用 DESeq2、TMM、CPM、VST 或 rlog 的文库大小标准化。'
            )
        if method == 'log2' and input_measurement == 'log_transformed':
            raise ValueError(
                '输入已经是 log-transformed expression，不能再次执行 log2(x+1)。请选择 none 保留当前表达尺度。'
            )
        if (method in {'log2', 'log2_quantile'}
                and measurement_info['requires_confirmation_for_log_transform']):
            raise ValueError(
                '自动检测到非整数连续值，但无法仅靠数值区分线性 FPKM/TPM 与已 log 的表达矩阵。'
                '请在“输入表达量尺度”中明确选择 continuous_expression 或 log_transformed 后重试。'
            )
        if method == 'log2_quantile' and input_measurement != 'continuous_expression':
            raise ValueError('log2_quantile 仅适用于确认的线性连续表达值；已 log 数据请选择 none。')
        if method == 'none' and input_measurement != 'log_transformed':
            raise ValueError('none 仅用于已经 log 变换的表达矩阵；线性 FPKM/TPM 请使用 log2。')

        # 低表达基因前置过滤
        n_genes_before = adata.n_vars
        # CPM 是 count 数据的文库大小标准化单位；对 FPKM/TPM 直接用
        # CPM 过滤会扭曲其含义，因此连续表达值只按原始表达阈值过滤。
        if min_expr_samples > 0 or max_zero_pct > 0:
            from scipy import sparse as _sp
            raw_for_filter = adata.X.toarray() if _sp.issparse(adata.X) else np.asarray(adata.X)
            gene_mask = np.ones(adata.n_vars, dtype=bool)
            if min_expr_samples > 0:
                if input_measurement == 'raw_counts':
                    lib_for_filter = raw_for_filter.sum(axis=1, keepdims=True)
                    lib_for_filter[lib_for_filter == 0] = 1
                    expr_for_filter = raw_for_filter / lib_for_filter * 1e6
                else:
                    expr_for_filter = raw_for_filter
                n_expr = np.array((expr_for_filter >= min_expr_value).sum(axis=0)).flatten()
                gene_mask &= n_expr >= min_expr_samples
            if max_zero_pct > 0:
                zero_pct = np.array((raw_for_filter == 0).sum(axis=0)).flatten() / adata.n_obs * 100
                gene_mask &= zero_pct <= max_zero_pct
            adata = adata[:, gene_mask].copy()

        if adata.n_vars == 0:
            raise ValueError("所有基因均被过滤掉（低表达），请降低 min_expr_samples 参数")

        from scipy import sparse as _sp
        raw_counts = np.asarray(adata.X.toarray()) if _sp.issparse(adata.X) else np.asarray(adata.X)
        raw_counts = raw_counts.astype(float)
        adata.layers['raw'] = raw_counts.copy()
        raw_lib = raw_counts.sum(axis=1)

        if method == 'deseq2':
            size_factors = _estimate_size_factors(raw_counts)
            adata.obs['size_factor'] = size_factors
            norm_counts = raw_counts / size_factors[:, None]
            adata.layers['normalized'] = norm_counts
            adata.X = np.log2(norm_counts + 1)
            adata.X = _sanitize(adata.X, 'deseq2 X')
            adata.layers['normalized'] = _sanitize(adata.layers['normalized'], 'deseq2 normalized')

        elif method == 'tmm':
            tmm_factors = _tmm_normalize(raw_counts)
            tmm_factors = np.where(tmm_factors > 0, tmm_factors, 1.0)
            lib_sizes = raw_counts.sum(axis=1)
            effective_lib_sizes = lib_sizes * tmm_factors
            effective_lib_sizes = np.where(effective_lib_sizes > 0, effective_lib_sizes, 1.0)
            adata.obs['tmm_factor'] = tmm_factors
            adata.obs['effective_library_size'] = effective_lib_sizes
            cpm_tmm = raw_counts / effective_lib_sizes[:, None] * 1e6
            adata.layers['normalized'] = cpm_tmm
            adata.X = np.log2(cpm_tmm + 1)
            adata.X = _sanitize(adata.X, 'tmm X')
            adata.layers['normalized'] = _sanitize(adata.layers['normalized'], 'tmm CPM normalized')

        elif method == 'cpm':
            cpm_target = float(self.params.get('cpm_target', 1e6))
            lib_sizes = raw_counts.sum(axis=1, keepdims=True)
            cpm = np.where(lib_sizes > 0, raw_counts / lib_sizes * cpm_target, 0.0)
            adata.layers['normalized'] = cpm
            adata.X = np.log2(cpm + 1)
            adata.X = _sanitize(adata.X, 'cpm X')
            adata.layers['normalized'] = _sanitize(adata.layers['normalized'], 'cpm normalized')

        elif method == 'log2':
            # FPKM/TPM 等已按文库归一化的连续表达值：仅做 log2 变换，
            # 保留样本间真实的整体分布差异。
            adata.layers['normalized'] = raw_counts.copy()
            adata.X = np.log2(raw_counts + 1)
            adata.X = _sanitize(adata.X, 'log2 X')

        elif method == 'none':
            # This is deliberately restricted to explicitly recognised log
            # input above.  It gives downstream PCA/heatmap steps a normal
            # pipeline artifact without applying a silent second log.
            adata.layers['normalized'] = raw_counts.copy()
            adata.X = _sanitize(raw_counts.copy(), 'preserved log-expression X')

        elif method == 'log2_quantile':
            log_counts = np.log2(raw_counts + 1)
            # 按样本（行）排序，计算每个秩次的参考分布
            sorted_per_sample = np.sort(log_counts, axis=1)  # (n_samples, n_genes)
            ref_distribution = np.mean(sorted_per_sample, axis=0)  # (n_genes,)
            # 按秩次替换：每个基因的值替换为该秩次对应的参考值
            from scipy.stats import rankdata
            norm = np.zeros_like(log_counts)
            for i in range(log_counts.shape[0]):  # 遍历样本
                ranks = rankdata(log_counts[i], method='ordinal').astype(int) - 1
                ranks = np.clip(ranks, 0, len(ref_distribution) - 1)
                norm[i] = ref_distribution[ranks]
            adata.layers['normalized'] = 2**norm - 1
            adata.X = norm
            adata.X = _sanitize(adata.X, 'quantile X')
            adata.layers['normalized'] = _sanitize(adata.layers['normalized'], 'quantile normalized')

        elif method == 'vst':
            size_factors = _estimate_size_factors(raw_counts)
            adata.obs['size_factor'] = size_factors
            vst_vals = _vst_transform(raw_counts, size_factors)
            adata.X = vst_vals
            adata.layers['normalized'] = vst_vals

        elif method == 'rlog':
            size_factors = _estimate_size_factors(raw_counts)
            adata.obs['size_factor'] = size_factors
            rlog_vals = _rlog_transform(raw_counts, size_factors)
            adata.X = rlog_vals
            adata.layers['normalized'] = rlog_vals

        else:
            raise ValueError(f"未知标准化方法: {method}，支持: deseq2/tmm/cpm/log2/log2_quantile/vst/rlog/none")

        # 统一输出标记
        linear_layer_methods = ('deseq2', 'tmm', 'cpm')
        adata.uns['normalization'] = {
            'method': method,
            'is_log_transformed': method in ('deseq2', 'tmm', 'cpm', 'vst', 'log2', 'log2_quantile', 'rlog', 'none'),
            'X_scale': ('log2(TMM-CPM+1)' if method == 'tmm' else
                        'log2(CPM+1)' if method in linear_layer_methods else
                        'log2(input+1)' if method == 'log2' else
                        'input log-expression (preserved)' if method == 'none' else
                        ('log2(normed+0.5)' if method in ('vst', 'rlog') else 'quantile-normalized')),
            'normalized_layer_scale': ('linear' if method in linear_layer_methods or method == 'log2' else
                                       ('same as X' if method in ('vst', 'rlog', 'none') else 'linear(2^X-1)')),
            'note': '近似实现：log2(normed + 0.5)，非 DESeq2 原始 VST/rlog' if method in ('vst', 'rlog') else '',
            'input_measurement': input_measurement,
            'input_measurement_provenance': measurement_info,
        }

        self.progress(60, "生成标准化前后对比图...")
        plots_dir = os.path.join(self.project_dir, 'plots')
        results_dir = os.path.join(self.project_dir, 'results')
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(results_dir, exist_ok=True)
        result_files = []
        from figure_engine import NatureFigureDirector, export_registered_figure
        director = NatureFigureDirector()
        nature_formats = ('svg', 'pdf', 'png')

        def _export_diagnostic(fig, stem, label, *, width='single', height_mm=None):
            spec = director.spec_from_params(
                'diagnostic', self.params, width=width, title=label,
            ).with_updates(
                formats=nature_formats,
                height_mm=height_mm,
            )
            exported, report = export_registered_figure(
                fig, os.path.join(plots_dir, stem), spec,
                category='qc', label=label,
                qa_path=os.path.join(results_dir, f'{stem}_nature_readiness.json'),
            )
            if not report.ready:
                self.progress(-1, f'{label} Nature readiness {report.score}/100；请查看 QA 报告。')
            import matplotlib.pyplot as _plt
            _plt.close(fig)
            return exported

        # Log2 is a value transform, not library-size normalization.  A former
        # implementation summed the preserved linear ``normalized`` layer in
        # this branch, yielding an unchanged, misleading "after" library-size
        # panel.  Do not render that invalid comparison.
        if method in ('log2', 'none'):
            self.progress(-1, f'{method} 不进行文库大小标准化，已跳过不适用的 library-size 前后对比图。')
        else:
            if method in ('vst', 'rlog'):
                # VST/rlog 输出不是 counts，展示标准化前后每样本均值对比。
                raw_means = np.log2(raw_counts + 1).mean(axis=1)
                norm_means = np.asarray(adata.X).mean(axis=1)
                raw_lib_for_plot = np.asarray(raw_means, dtype=float).reshape(-1)
                norm_lib_for_plot = np.asarray(norm_means, dtype=float).reshape(-1)
                raw_ylabel = 'Mean log2(raw + 1)'
                normalized_ylabel = 'Mean transformed expression'
                library_title = 'Transformation mean-value check'
            else:
                norm_layer = adata.layers.get('normalized', adata.X)
                norm_lib = norm_layer.sum(axis=1) if hasattr(norm_layer, 'sum') else np.ones(adata.n_obs)
                raw_lib_for_plot = np.asarray(raw_lib, dtype=float).reshape(-1)
                norm_lib_for_plot = np.asarray(norm_lib, dtype=float).reshape(-1)
                if input_measurement == 'raw_counts':
                    raw_ylabel = 'Library size (counts)'
                    normalized_ylabel = 'Normalized expression total'
                    library_title = 'Normalization library-size check'
                else:
                    raw_ylabel = 'Total expression'
                    normalized_ylabel = 'Normalized expression total'
                    library_title = 'Total-expression normalization check'
            result_files.extend(_export_diagnostic(
                director.render(
                    director.spec_from_params('diagnostic', self.params, width='double',
                                              title=library_title).with_updates(
                                                  extra={'kind': 'normalization_library'},
                                                  formats=nature_formats, height_mm=60.0),
                    {
                        'kind': 'normalization_library',
                        'sample_labels': adata.obs.index.tolist(),
                        'raw_values': raw_lib_for_plot,
                        'normalized_values': norm_lib_for_plot,
                        'raw_ylabel': raw_ylabel,
                        'normalized_ylabel': normalized_ylabel,
                    },
                ),
                'bulk_norm_libsize', '文库大小对比', width='double', height_mm=60.0,
            ))

        if method in ('deseq2', 'tmm') and ('size_factor' in adata.obs.columns or 'tmm_factor' in adata.obs.columns):
            factor_col = 'tmm_factor' if method == 'tmm' else 'size_factor'
            # Reuse the fixed library-size diagnostic contract for factors;
            # this keeps the legacy bar chart from choosing its own fonts and
            # exports the same physical SVG/PDF/PNG set as other diagnostics.
            result_files.extend(_export_diagnostic(
                director.render(
                    director.spec_from_params('diagnostic', self.params, width='double',
                                              title=f'{method.upper()} normalization factors').with_updates(
                                                  extra={'kind': 'normalization_library'},
                                                  formats=nature_formats, height_mm=60.0),
                    {
                        'kind': 'normalization_library',
                        'sample_labels': adata.obs.index.tolist(),
                        'raw_values': np.ones(adata.n_obs, dtype=float),
                        'normalized_values': np.asarray(adata.obs[factor_col], dtype=float),
                        'raw_ylabel': 'Reference',
                        'normalized_ylabel': 'Factor',
                    },
                ),
                'bulk_norm_sizefactors', 'Size factors', width='double', height_mm=60.0,
            ))

        # 按样本展示表达分布，才能判断样本分布是否真正被对齐；将全部值
        # 压成一个箱线图会掩盖这一关键信息。
        sample_labels = adata.obs.index.tolist()
        norm_values = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
        if method == 'log2':
            box_title = 'Raw input vs log2-transformed expression'
            box_note = '左侧为未经变换的输入值；右侧为实际输出 log2(input + 1)。log2 只压缩数值范围，不强制各样本同分布。'
            box_panel_titles = ['Raw input', 'Actual log2-transformed output']
        elif method == 'none':
            box_title = 'Already log-transformed input (preserved)'
            box_note = '输入已被明确标记为 log 表达值；平台未执行第二次 log2 变换。两侧相同是预期行为。'
            box_panel_titles = ['Input log-expression', 'Preserved output']
        else:
            box_title = 'Raw input vs normalized expression distribution'
            box_note = '左侧始终为真实输入值，右侧为实际标准化/变换输出；不会再把预期 log2 值误标为变换前数据。'
            box_panel_titles = ['Raw input', 'Actual normalized output']
        result_files.extend(_export_diagnostic(
            director.render(
                director.spec_from_params('diagnostic', self.params, width='double',
                                          title=box_title).with_updates(
                                              extra={'kind': 'normalization_boxplot'},
                                              formats=nature_formats, height_mm=76.0),
                {
                    'kind': 'normalization_boxplot',
                    'sample_labels': sample_labels,
                    'raw_matrix': raw_counts,
                    'normalized_matrix': norm_values,
                    'note': box_note,
                    'ylabel': 'Expression',
                    'panel_titles': box_panel_titles,
                },
            ),
            'bulk_norm_boxplot_compare', '表达分布对比', width='double', height_mm=76.0,
        ))

        # PCA 前后对比
        n_pcs = min(10, adata.n_obs - 1, adata.n_vars - 1)
        pca_preprocessing = None
        pca_compare_variance = None
        if n_pcs >= 2:
            # Left panel uses the real input.  Both PCA calls use the shared
            # Bulk convention: PCA mean-centres genes but never Z-scores them.
            # This makes the right log2 panel directly comparable to QC PCA.
            adata_raw_pca = sc.AnnData(X=raw_counts.copy(), obs=adata.obs.copy())
            pca_preprocessing = run_bulk_sample_pca(adata_raw_pca, n_comps=n_pcs)
            pc_raw = adata_raw_pca.obsm['X_pca']
            # 标准化后 PCA
            adata_norm_pca = sc.AnnData(X=adata.X.copy(), obs=adata.obs.copy())
            run_bulk_sample_pca(adata_norm_pca, n_comps=n_pcs)
            pc_norm = adata_norm_pca.obsm['X_pca']
            if method == 'log2':
                pca_titles = ['Raw input PCA', 'Actual log2-transformed output PCA']
                pca_title = 'Raw-vs-log2 PCA comparison'
            elif method == 'none':
                pca_titles = ['Input log-expression PCA', 'Preserved output PCA']
                pca_title = 'Preserved log-expression PCA check'
            else:
                pca_titles = ['Raw input PCA', 'Actual normalized output PCA']
                pca_title = 'Raw-vs-normalized PCA comparison'
            # 解释方差必须来自实际 PCA，而不能由调用方临时填充。
            # 旧实现传入 NaN，最终会在投稿图轴标签中显示为 ``(nan%)``。
            raw_variance = np.asarray(
                adata_raw_pca.uns.get('pca', {}).get('variance_ratio', []),
                dtype=float,
            ).reshape(-1)[:2]
            norm_variance = np.asarray(
                adata_norm_pca.uns.get('pca', {}).get('variance_ratio', []),
                dtype=float,
            ).reshape(-1)[:2]
            explained_variance = np.concatenate((raw_variance, norm_variance))
            if explained_variance.size < 4:
                explained_variance = np.pad(
                    explained_variance, (0, 4 - explained_variance.size),
                    constant_values=np.nan,
                )
            pca_compare_variance = {
                'input': [float(value) for value in raw_variance],
                'output': [float(value) for value in norm_variance],
            }
            result_files.extend(_export_diagnostic(
                director.render(
                    director.spec_from_params('diagnostic', self.params, width='double',
                                              title=pca_title).with_updates(
                                                  extra={'kind': 'normalization_pca'},
                                                  formats=nature_formats, height_mm=70.0),
                    {
                        'kind': 'normalization_pca',
                        'coordinates': [pc_raw[:, :2], pc_norm[:, :2]],
                        'groups': ['All'] * adata.n_obs,
                        'sample_labels': sample_labels,
                        'panel_titles': pca_titles,
                        'explained_variance': explained_variance.tolist(),
                    },
                ),
                'bulk_norm_pca_compare', 'PCA 前后对比', width='double', height_mm=70.0,
            ))

        self.progress(85, "保存结果...")
        intermediate_dir = os.path.join(self.project_dir, 'intermediate')
        os.makedirs(intermediate_dir, exist_ok=True)
        output_path = os.path.join(intermediate_dir, 'bulk_normalize_output.h5ad')
        adata.write_h5ad(output_path)

        self.progress(100, "完成")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'method': method,
                'input_measurement': input_measurement,
                'n_samples': adata.n_obs,
                'n_genes_before_filter': n_genes_before,
                'n_genes_after_filter': adata.n_vars,
                'genes_filtered': n_genes_before - adata.n_vars,
                'median_size_factor': round(float(adata.obs['size_factor'].median()), 3) if 'size_factor' in adata.obs.columns else None,
                'median_tmm_factor': round(float(adata.obs['tmm_factor'].median()), 3) if 'tmm_factor' in adata.obs.columns else None,
                'is_log_transformed': adata.uns.get('normalization', {}).get('is_log_transformed', True),
                'input_measurement_source': measurement_info['source'],
                'input_measurement_confidence': measurement_info['confidence'],
                'pca_preprocessing': pca_preprocessing if n_pcs >= 2 else None,
                'pca_compare_variance_ratio': pca_compare_variance,
                'max_zero_pct': max_zero_pct,
                'max_zero_pct_filter_enabled': max_zero_pct > 0,
                'max_zero_pct_semantics': '0 disables this optional zero-value filter',
            }
        }


# --- 辅助函数 ---


def _infer_measurement_type(matrix, input_path=''):
    """Backward-compatible matrix-only classifier used by older callers/tests."""
    from anndata import AnnData
    from modules.io_utils import infer_expression_measurement
    return infer_expression_measurement(AnnData(X=matrix), input_path)


def _tmm_normalize(counts):
    """TMM 标准化（edgeR 风格）。counts: (n_samples, n_genes) ndarray of raw counts.
    返回 per-sample TMM 因子。"""
    n_samples, n_genes = counts.shape
    lib_sizes = counts.sum(axis=1)
    lib_sizes[lib_sizes == 0] = 1

    # 参考样本：文库大小最接近中位数
    median_lib = np.median(lib_sizes)
    ref_idx = int(np.argmin(np.abs(lib_sizes - median_lib)))

    factors = np.ones(n_samples)
    for i in range(n_samples):
        if i == ref_idx:
            continue
        mask = (counts[i] > 0) & (counts[ref_idx] > 0)
        if mask.sum() < 10:
            continue
        fi = lib_sizes[i]
        fr = lib_sizes[ref_idx]
        Mi = np.log2((counts[i, mask] / fi) / (counts[ref_idx, mask] / fr))
        Ai = 0.5 * (np.log2(counts[i, mask] / fi) + np.log2(counts[ref_idx, mask] / fr))

        m_lo, m_hi = np.percentile(Mi, [30, 70])
        a_hi = np.percentile(Ai, 95)
        keep = (Mi >= m_lo) & (Mi <= m_hi) & (Ai <= a_hi)

        if keep.sum() > 0:
            # edgeR delta method 近似方差：Var(M) ≈ 1/ci + 1/cr
            wi = 1.0 / counts[i, mask][keep]
            wr = 1.0 / counts[ref_idx, mask][keep]
            weights = 1.0 / (wi + wr + 1e-30)
            factors[i] = 2 ** (-np.average(Mi[keep], weights=weights))

    # 归一化使几何均值为 1
    log_factors = np.log(factors[factors > 0])
    if len(log_factors) > 0:
        geo_mean = np.exp(np.mean(log_factors))
        factors = factors / geo_mean

    return factors


def _vst_transform(counts, size_factors):
    """近似 VST。counts: raw counts (n_samples, n_genes), size_factors: (n_samples,)."""
    normed = counts / size_factors[:, None]
    normed = np.maximum(normed, 0)
    vst = np.log2(normed + 0.5)
    return np.nan_to_num(vst, nan=0.0, posinf=0.0, neginf=0.0)


def _rlog_transform(counts, size_factors, prior_mean=None):
    """近似 rlog。小样本时通过正则化收缩基因效应。"""
    normed = counts / size_factors[:, None]
    normed = np.maximum(normed, 0)
    log_normed = np.log2(normed + 0.5)

    n_samples = counts.shape[0]
    if prior_mean is None:
        prior_mean = float(np.mean(log_normed))

    gene_effects = log_normed.mean(axis=0) - prior_mean
    shrinkage = min(1.0, n_samples / 30.0)

    sample_residuals = log_normed - log_normed.mean(axis=0, keepdims=True)
    rlog = prior_mean + shrinkage * gene_effects[None, :] + sample_residuals

    return np.nan_to_num(rlog, nan=0.0, posinf=0.0, neginf=0.0)


def _sanitize(arr, label='data'):
    """替换 inf/NaN 为 0 并记录警告。"""
    n_bad = int(np.sum(~np.isfinite(arr)))
    if n_bad > 0:
        logger.warning(f"[bulk_normalize] {label}: {n_bad} 个 inf/NaN 值被替换为 0")
    return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)


def _estimate_size_factors(raw_counts):
    """DESeq2 中位比率法计算 size factors。raw_counts: (n_samples, n_genes)."""
    from scipy.stats import gmean
    # 仅使用所有样本中都非零的基因计算几何均值
    nonzero_all = (raw_counts > 0).all(axis=0)
    n_usable = int(nonzero_all.sum())
    if n_usable == 0:
        # 无全非零基因时，使用 log-几何均值近似（加 1 避免 log(0)）
        geo_means = gmean(raw_counts + 1, axis=0)
        logger.warning(f"[bulk_normalize] 无全非零基因，使用 log-几何均值近似计算 size factors")
    else:
        geo_means = np.ones(raw_counts.shape[1])
        geo_means[nonzero_all] = gmean(raw_counts[:, nonzero_all], axis=0)
    ratios = raw_counts / (geo_means + 1e-10)
    size_factors = np.median(ratios, axis=1)
    size_factors = np.where(size_factors > 0, size_factors, 1.0)
    return size_factors
