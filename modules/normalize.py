from modules.base import BaseAnalysis


def _matrix_values(matrix):
    """Return stored/dense matrix values without densifying sparse counts."""
    import numpy as np

    if hasattr(matrix, 'data') and hasattr(matrix, 'tocsr'):
        return np.asarray(matrix.data).ravel()
    return np.asarray(matrix).ravel()


def _is_raw_count_matrix(matrix):
    """Return whether a matrix is finite, non-negative and integer-valued."""
    import numpy as np

    values = _matrix_values(matrix)
    if values.size == 0:
        return True
    return bool(
        np.isfinite(values).all()
        and (values >= 0).all()
        and np.allclose(values, np.rint(values), rtol=0.0, atol=1e-6)
    )


def _prepare_counts_input(adata):
    """Select trustworthy raw counts and make them the normalization input.

    An existing ``counts`` layer is authoritative and is never overwritten.
    Without one, only an integer, non-negative ``X`` can be promoted to counts;
    this prevents an already-normalized matrix from being silently relabeled as
    raw counts.
    """
    if 'counts' in adata.layers:
        counts = adata.layers['counts']
        counts_source = 'layers[counts]'
        if counts.shape != adata.shape:
            raise ValueError(
                "layers['counts'] 的形状与 adata 不一致，无法安全标准化。"
            )
        if not _is_raw_count_matrix(counts):
            raise ValueError(
                "layers['counts'] 不是非负整数原始计数；为避免静默数据污染，"
                "标准化已停止。"
            )
    else:
        if not _is_raw_count_matrix(adata.X):
            raise ValueError(
                "未找到 layers['counts']，且 adata.X 不是非负整数原始计数。"
                "请提供真实 counts 后再运行标准化。"
            )
        adata.layers['counts'] = adata.X.copy()
        counts = adata.layers['counts']
        counts_source = 'adata.X'

    # Always normalize from raw counts, even when input X is already transformed.
    adata.X = counts.copy()
    return counts_source


def _require_pearson_normalizer(scanpy_module):
    """Return Scanpy's Pearson implementation or fail without a data fallback."""
    experimental = getattr(scanpy_module, 'experimental', None)
    experimental_pp = getattr(experimental, 'pp', None)
    normalizer = getattr(experimental_pp, 'normalize_pearson_residuals', None)
    if not callable(normalizer):
        raise RuntimeError(
            "当前 Scanpy 不提供 experimental.pp.normalize_pearson_residuals；"
            "不能把原始 counts 静默写成 normalized。请升级 Scanpy 或改用 log1p。"
        )
    return normalizer


def _pearson_clip_setting(value):
    """Translate the checkbox contract to Scanpy's clip argument."""
    import numpy as np

    if isinstance(value, str):
        enabled = value.strip().lower() not in {'', '0', 'false', 'no', 'off'}
    else:
        enabled = bool(value)
    return enabled, None if enabled else np.inf


def _row_sums(matrix):
    import numpy as np

    return np.asarray(matrix.sum(axis=1), dtype=float).reshape(-1)


def _normalize_total_matrix(matrix, target_sum):
    """Return per-cell total-count normalized data without densifying it."""
    import numpy as np

    totals = _row_sums(matrix)
    scale = np.divide(
        float(target_sum), totals,
        out=np.zeros_like(totals, dtype=float), where=totals > 0,
    )
    if hasattr(matrix, 'tocsr'):
        dtype = np.float64 if matrix.dtype == np.float64 else np.float32
        normalized = matrix.astype(dtype, copy=True).multiply(
            scale.astype(dtype)[:, None]
        ).tocsr()
        normalized.eliminate_zeros()
        return normalized
    dtype = np.float64 if np.asarray(matrix).dtype == np.float64 else np.float32
    return np.asarray(matrix, dtype=dtype) * scale.astype(dtype)[:, None]


def _log1p_matrix(matrix):
    """Apply natural log1p while retaining sparse structure."""
    import numpy as np

    if hasattr(matrix, 'tocsr'):
        transformed = matrix.tocsr(copy=True)
        transformed.data = np.log1p(transformed.data)
        return transformed
    return np.log1p(np.asarray(matrix))


def _sample_expression_values(matrix, *, nonzero_only=False, max_values=100000):
    """Deterministically sample finite expression values for diagnostics."""
    import numpy as np

    values = _matrix_values(matrix)
    values = values[np.isfinite(values)]
    if nonzero_only:
        values = values[values != 0]
    if values.size > max_values:
        rng = np.random.default_rng(0)
        values = values[rng.choice(values.size, max_values, replace=False)]
    return values


def _expression_distribution_payload(method, matrix):
    """Return one internally consistent output-scale expression distribution."""
    if method == 'pearson_residuals':
        values = _sample_expression_values(matrix, nonzero_only=False)
        return values, 'Pearson residuals', 'Pearson Residual Distribution'
    values = _sample_expression_values(matrix, nonzero_only=True)
    return (
        values,
        'log1p normalized expression (non-zero)',
        'Normalized Non-zero Expression Distribution',
    )


class NormalizeAnalysis(BaseAnalysis):
    MODULE_NAME = "normalize"
    DISPLAY_NAME = "标准化"
    DESCRIPTION = "数据标准化（log1p / Pearson 残差）"
    INPUT_REQUIRES = []

    def run(self, input_path):
        import numpy as np
        from modules.native_figures import histogram_figure

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        input_shape = tuple(adata.shape)
        input_obs_names = adata.obs_names.copy()
        input_var_names = adata.var_names.copy()
        counts_source = _prepare_counts_input(adata)

        method = str(self.params.get('method', 'log1p') or '').strip().lower()
        if method not in {'log1p', 'pearson_residuals'}:
            raise ValueError(f"不支持的标准化方法: {method!r}")
        target_sum = None
        if method == 'log1p':
            target_sum = float(self.params.get('target_sum', 10000))
            if not np.isfinite(target_sum) or target_sum <= 0:
                raise ValueError("target_sum 必须是大于 0 的有限数值。")

        raw_library_sizes = _row_sums(adata.layers['counts'])
        linear_normalized_library_sizes = None
        normalized_layer = None
        zero_count_genes = 0

        self.progress(20, f"Normalizing ({method})...")
        if method == 'pearson_residuals':
            import scanpy as sc

            self.progress(30, "Computing Pearson residuals...")
            if np.any(raw_library_sizes <= 0):
                raise ValueError(
                    "Pearson 残差要求每个细胞至少有一个 count；"
                    "请先移除总计数为 0 的细胞。"
                )
            gene_totals = np.asarray(
                adata.layers['counts'].sum(axis=0), dtype=float
            ).reshape(-1)
            nonzero_gene_mask = gene_totals > 0
            zero_count_genes = int((~nonzero_gene_mask).sum())
            clip_values, clip_argument = _pearson_clip_setting(
                self.params.get('clip_values', True)
            )
            pearson_normalizer = _require_pearson_normalizer(sc)
            # Scanpy uses clip=None for its documented ±sqrt(n_obs) default.
            # np.inf explicitly disables clipping.
            if zero_count_genes:
                # Residuals are undefined for genes with zero observations.
                # Compute on expressed genes, then restore zero columns so the
                # output keeps the exact input gene set and order.
                pearson_input = adata[:, nonzero_gene_mask].copy()
                pearson_normalizer(
                    pearson_input, clip=clip_argument, inplace=True,
                )
                # 保持稀疏性：np.asarray(csr) 会得到 0 维 object 数组，
                # 在保持稀疏输出的 scanpy 版本上 slice 赋值会失败。
                from scipy import sparse
                if sparse.issparse(pearson_input.X):
                    residuals = sparse.csr_matrix(
                        adata.shape, dtype=pearson_input.X.dtype,
                    )
                    residuals[:, nonzero_gene_mask] = pearson_input.X
                else:
                    residuals = np.zeros(
                        adata.shape, dtype=pearson_input.X.dtype,
                    )
                    residuals[:, nonzero_gene_mask] = np.asarray(
                        pearson_input.X,
                    )
                adata.X = residuals
                adata.uns['pearson_residuals_normalization'] = dict(
                    pearson_input.uns['pearson_residuals_normalization']
                )
                del pearson_input
            else:
                pearson_normalizer(
                    adata, clip=clip_argument, inplace=True,
                )
        else:
            adata.X = _normalize_total_matrix(adata.X, target_sum)
            # Capture the meaningful linear-scale totals before log1p. Summing
            # log-transformed X is not a library size.
            linear_normalized_library_sizes = _row_sums(adata.X)
            normalized_layer = adata.X.copy()
            adata.X = _log1p_matrix(adata.X)
            adata.uns['log1p'] = {'base': None}
            clip_values = None

        if tuple(adata.shape) != input_shape:
            raise RuntimeError("标准化不应改变细胞或基因数量。")
        if not adata.obs_names.equals(input_obs_names):
            raise RuntimeError("标准化不应改变细胞名称或顺序。")
        if not adata.var_names.equals(input_var_names):
            raise RuntimeError("标准化不应改变基因名称或顺序。")

        normalized_values = _matrix_values(adata.X)
        if normalized_values.size and not np.isfinite(normalized_values).all():
            raise RuntimeError("标准化结果包含 NaN 或 Inf。")
        # Keep the linear normalized matrix separate from log1p X.  For
        # Pearson residuals both are intentionally the same residual scale.
        adata.layers['normalized'] = (
            normalized_layer if normalized_layer is not None else adata.X.copy()
        )
        normalization_metadata = {
            'method': method,
            'counts_source': counts_source,
            'x_contains': method,
            'normalized_layer_matches_x': bool(method == 'pearson_residuals'),
            'layer_semantics': {
                'counts': 'raw non-negative integer counts before normalization',
                'normalized': (
                    'library-size normalized linear expression before log1p'
                    if method == 'log1p' else
                    'Pearson residuals on the same scale as X'
                ),
                'X': (
                    'log1p(library-size normalized expression)'
                    if method == 'log1p' else
                    'Pearson residuals'
                ),
            },
        }
        if method == 'log1p':
            normalization_metadata['target_sum'] = float(target_sum)
        else:
            normalization_metadata['pearson_clipped'] = bool(clip_values)
            normalization_metadata['zero_count_genes_set_to_zero'] = zero_count_genes
        adata.uns['normalization'] = normalization_metadata

        self.progress(70, "Generating normalization plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []
        donor_key = None
        normalization_by_batch = {}
        if method == 'log1p':
            from modules.sc_figure_diagnostics import (
                library_size_by_batch_figure,
                normalization_batch_medians_figure,
                resolve_donor_key,
            )
            donor_key, _ = resolve_donor_key(
                adata, self.params.get('batch_key', ''), require_multiple=True,
            )

        # Total-count normalization has a meaningful before/after library size.
        # Pearson residuals do not preserve or target per-cell library sums, so
        # deliberately omit this plot for that method.
        if method == 'log1p':
            raw_library_sizes_log10 = np.log10(np.clip(raw_library_sizes, 0, None) + 1.0)
            fig = histogram_figure(
                [raw_library_sizes_log10], labels=['Raw counts'],
                bins=50,
                title='Input Library Size and Normalization Target',
                x_label='log10(raw total counts + 1) per cell',
            )
            axis = fig.axes[0]
            axis.axvline(
                np.log10(float(target_sum) + 1.0), color='#B64342',
                linestyle='--', linewidth=1.4,
                label=f'Normalization target ({target_sum:g})',
            )
            nonempty_mask = raw_library_sizes > 0
            if np.any(nonempty_mask):
                normalized_nonempty = linear_normalized_library_sizes[nonempty_mask]
                max_deviation = float(
                    np.max(np.abs(normalized_nonempty - target_sum))
                )
                axis.text(
                    0.98, 0.94,
                    'After total-count normalization\n'
                    f'non-empty cells: {int(nonempty_mask.sum()):,}\n'
                    f'max |total - target|: {max_deviation:.3g}',
                    transform=axis.transAxes, ha='right', va='top', fontsize=7.5,
                    color='#475467',
                    bbox={
                        'boxstyle': 'round,pad=0.3', 'facecolor': 'white',
                        'edgecolor': '#D0D5DD', 'alpha': 0.9, 'linewidth': 0.6,
                    },
                )
            axis.legend(frameon=False, fontsize=8)
            fig.tight_layout(pad=1.1)
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, 'normalize_libsize.png', 'histogram',
                'Input Library Size and Normalization Target',
                formats=('png', 'svg'), dpi=300,
            ))

        # Plot one well-defined output scale. Comparing raw counts with Pearson
        # residuals or normalized log expression on one axis is not meaningful.
        if self.params.get('show_expression_distribution', True):
            expression_values, expression_label, expression_title = (
                _expression_distribution_payload(method, adata.X)
            )
            if expression_values.size:
                fig_expr = histogram_figure(
                    [expression_values], labels=[expression_label], bins=80,
                    title=expression_title,
                    x_label=expression_label,
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_expr, plots_dir, 'normalize_expression_distribution.png',
                    'histogram', expression_title,
                    formats=('png', 'svg'), dpi=300,
                ))

        if method == 'log1p' and donor_key and donor_key in adata.obs.columns:
            fig_by_batch = library_size_by_batch_figure(
                raw_library_sizes,
                linear_normalized_library_sizes,
                adata.obs,
                donor_key,
                target_sum,
            )
            if fig_by_batch is not None:
                result_files.extend(self.save_matplotlib_figure(
                    fig_by_batch, plots_dir, 'normalize_libsize_by_donor.png',
                    'violin', 'Library size before and after normalization by donor',
                    formats=('png', 'svg'), dpi=300,
                ))
            fig_medians = normalization_batch_medians_figure(
                raw_library_sizes,
                linear_normalized_library_sizes,
                adata.obs,
                donor_key,
            )
            if fig_medians is not None:
                result_files.extend(self.save_matplotlib_figure(
                    fig_medians, plots_dir, 'normalize_donor_medians.png',
                    'bar', 'Donor-level normalization diagnostics',
                    formats=('png', 'svg'), dpi=300,
                ))
            donor_labels = adata.obs[donor_key].astype(str)
            for donor in sorted(donor_labels.unique(), key=lambda value: (len(value), value)):
                mask = donor_labels.to_numpy() == donor
                raw_values = raw_library_sizes[mask]
                normalized_values = linear_normalized_library_sizes[mask]
                normalized_values = normalized_values[np.isfinite(normalized_values)]
                deviations = np.abs(normalized_values - float(target_sum))
                normalization_by_batch[donor] = {
                    'n_cells': int(mask.sum()),
                    'raw_total_median': float(np.nanmedian(raw_values)),
                    'normalized_total_median': float(np.nanmedian(normalized_values)) if len(normalized_values) else None,
                    'max_abs_total_target_deviation': float(np.nanmax(deviations)) if len(deviations) else None,
                    'median_abs_total_target_deviation': float(np.nanmedian(deviations)) if len(deviations) else None,
                }

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'normalize')

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_cells': int(adata.n_obs),
                'n_genes': int(adata.n_vars),
                'method': method,
                'target_sum': float(target_sum) if method == 'log1p' else None,
                'counts_source': counts_source,
                'normalized_layer_matches_x': normalization_metadata['normalized_layer_matches_x'],
                'layer_semantics': normalization_metadata['layer_semantics'],
                'library_size_plot': method == 'log1p',
                'library_size_plot_scale': 'log10_total_counts_plus_1' if method == 'log1p' else None,
                'batch_key': donor_key,
                'normalization_by_batch': normalization_by_batch,
                'zero_count_genes_set_to_zero': zero_count_genes,
            }
        }
