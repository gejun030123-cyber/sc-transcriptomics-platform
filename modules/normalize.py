from modules.base import BaseAnalysis

class NormalizeAnalysis(BaseAnalysis):
    MODULE_NAME = "normalize"
    DISPLAY_NAME = "标准化"
    DESCRIPTION = "数据标准化（log1p / Pearson 残差）"
    INPUT_REQUIRES = []

    def run(self, input_path):
        import scanpy as sc
        import omicverse as ov
        import numpy as np
        from modules.native_figures import histogram_figure

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        adata.layers["counts"] = adata.X.copy()

        method = self.params.get('method', 'log1p')
        target_sum = float(self.params.get('target_sum', 10000))

        self.progress(20, f"Normalizing ({method})...")
        if method == 'pearson_residuals':
            self.progress(30, "Computing Pearson residuals...")
            clip_values = self.params.get('clip_values', True)
            adata.layers["normalized"] = sc.experimental.pp.normalize_pearson_residuals(
                adata, clip=clip_values, copy=True
            ).X if hasattr(sc.experimental.pp, 'normalize_pearson_residuals') else adata.X
            # Fallback: use omicverse pearson mode
            try:
                adata = ov.pp.preprocess(adata, mode='pearson', target_sum=target_sum)
            except Exception:
                sc.pp.normalize_total(adata, target_sum=target_sum)
                sc.pp.log1p(adata)
        else:
            try:
                adata = ov.pp.preprocess(adata, mode='shiftlog|seurat', target_sum=target_sum)
            except Exception:
                # omicverse 内部 HVG 绑定 bug 的回退方案
                sc.pp.normalize_total(adata, target_sum=target_sum)
                sc.pp.log1p(adata)

        self.progress(70, "Generating normalization plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        # Library size distribution before/after (static display figure)
        series = []
        labels = []
        if 'counts' in adata.layers:
            counts_lib = np.array(adata.layers["counts"].sum(axis=1)).flatten()
            series.append(counts_lib)
            labels.append('Before')
        norm_lib = np.array(adata.X.sum(axis=1)).flatten() if hasattr(adata.X, 'sum') else None
        if norm_lib is not None:
            series.append(norm_lib)
            labels.append('After')
        fig = histogram_figure(series, labels=labels, bins=50,
                               title='Library Size Distribution',
                               x_label='Total Counts')
        result_files.extend(self.save_matplotlib_figure(
            fig, plots_dir, 'normalize_libsize.png', 'histogram',
            'Library Size Distribution', formats=('png', 'svg'), dpi=300,
        ))

        # Expression value distribution after normalization
        if self.params.get('show_expression_distribution', True):
            def _sample_values(matrix, max_values=100000):
                if hasattr(matrix, 'toarray') and hasattr(matrix, 'data'):
                    vals = np.asarray(matrix.data).ravel()
                else:
                    vals = np.asarray(matrix).ravel()
                vals = vals[np.isfinite(vals)]
                if vals.size > max_values:
                    rng = np.random.default_rng(0)
                    vals = vals[rng.choice(vals.size, max_values, replace=False)]
                return vals

            expr_series = []
            expr_labels = []
            if 'counts' in adata.layers:
                raw_vals = _sample_values(adata.layers['counts'])
                if raw_vals.size:
                    expr_series.append(np.log1p(raw_vals))
                    expr_labels.append('log1p(raw counts)')
            norm_vals = _sample_values(adata.X)
            if norm_vals.size:
                expr_series.append(norm_vals)
                expr_labels.append('normalized X')
            fig_expr = histogram_figure(
                expr_series, labels=expr_labels, bins=80,
                title='Expression Value Distribution',
                x_label='Expression value',
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_expr, plots_dir, 'normalize_expression_distribution.png',
                'histogram', 'Expression Value Distribution',
                formats=('png', 'svg'), dpi=300,
            ))

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'normalize')

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'method': method,
                'target_sum': target_sum,
            }
        }
