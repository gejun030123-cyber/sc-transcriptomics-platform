import logging

from modules.base import BaseAnalysis
from modules.io_utils import resolve_obs_grouping
from modules.pc_strategy import DEFAULT_FINAL_N_PCS


logger = logging.getLogger(__name__)


def _embedding_color_keys(batch_key, *fallback_keys):
    """Return ordered, unique plot keys with the configured batch key first."""
    return list(dict.fromkeys(
        key for key in [batch_key, *fallback_keys] if key
    ))


def _invalidate_pca_results(adata):
    """Remove PCA artifacts before a fresh run so old coordinates cannot leak in."""
    removed = []
    for mapping, key in ((adata.obsm, 'X_pca'), (adata.varm, 'PCs'), (adata.uns, 'pca')):
        if key in mapping:
            del mapping[key]
            removed.append(key)
    return removed


def _pca_fingerprint(rep, max_components=10):
    """Return a reproducibility fingerprint for the first PCA coordinates."""
    import hashlib
    import numpy as np

    array = np.asarray(rep)
    if array.ndim != 2:
        raise ValueError(f'PCA representation must be 2D, got shape {array.shape}')
    n_components = min(int(max_components), array.shape[1])
    payload = np.ascontiguousarray(array[:, :n_components]).tobytes()
    return hashlib.md5(payload).hexdigest()


def _select_elbow_n_comps(variance_ratio, max_candidate):
    """Select an elbow PC count without using a cumulative-variance target.

    The old implementation took the argmax of the second difference of the
    cumulative curve.  For a normal monotonically decreasing PCA spectrum that
    quantity is usually largest in the long tail, which can select 45 PCs even
    when the visible elbow is near the first few dimensions.  We instead find
    the greatest distance between the normalized variance-ratio curve and its
    first-to-last chord.  This is a bounded geometric elbow rule and does not
    imply that single-cell PCA must reach a fixed cumulative percentage.
    """
    import numpy as np

    ratios = np.asarray(variance_ratio, dtype=float).reshape(-1)
    available = min(int(max_candidate), len(ratios))
    if available < 2:
        return available
    if available < 3:
        return available

    ratios = ratios[:available]
    if not np.isfinite(ratios).all():
        ratios = np.nan_to_num(ratios, nan=0.0, posinf=0.0, neginf=0.0)
    span = float(ratios[0] - ratios[-1])
    if span <= 0:
        return min(5, available)

    x = np.linspace(0.0, 1.0, available)
    curve = (ratios - ratios[-1]) / span
    chord = 1.0 - x
    distances = chord - curve
    candidate = int(np.argmax(distances) + 1)
    minimum = min(5, available)
    return max(minimum, min(candidate, available))


def _analytical_qc_summary(adata, *, pca_hvg_only, pca_n_genes, n_comps,
                           batch_overlap=None):
    """Build analysis-QC evidence separately from figure-quality scoring."""
    import numpy as np

    variance = np.asarray(
        adata.uns.get('pca', {}).get('variance_ratio', []), dtype=float,
    )
    pc1_pct = round(float(variance[0] * 100), 2) if len(variance) else None
    cumulative_pct = round(float(variance[:n_comps].sum() * 100), 2) if len(variance) else None
    warnings = []
    checks = []

    checks.append({
        'name': 'pca_input',
        'status': 'pass' if pca_hvg_only else 'review',
        'message': (
            f'PCA 使用 {pca_n_genes} 个 HVG。' if pca_hvg_only else
            f'PCA 使用全部 {pca_n_genes} 个基因；建议优先使用 HVG。'
        ),
    })
    if not pca_hvg_only:
        warnings.append(
            f'当前为全基因 PCA（{pca_n_genes} genes）；解释率通常更分散，'
            '这本身不等于数据质量差。'
        )

    if pc1_pct is not None:
        checks.append({
            'name': 'pc1_variance',
            'status': 'info',
            'value_pct': pc1_pct,
            'message': f'PC1 解释率 {pc1_pct:.2f}%。',
        })
    if cumulative_pct is not None:
        checks.append({
            'name': 'cumulative_variance',
            'status': 'info',
            'value_pct': cumulative_pct,
            'message': f'前 {min(n_comps, len(variance))} 个 PC 累计解释率 {cumulative_pct:.2f}%。',
        })
        if not pca_hvg_only and cumulative_pct < 20.0:
            warnings.append(
                f'前 {min(n_comps, len(variance))} 个 PC 累计解释率仅 {cumulative_pct:.2f}%；'
                '这是全基因 PCA 中常见的方差分散表现，不应单独作为低质量判据。'
            )

    overlap = batch_overlap or {}
    if overlap.get('warnings'):
        warnings.extend(overlap['warnings'])
        checks.append({
            'name': 'batch_cluster_overlap',
            'status': 'warning',
            'message': '；'.join(overlap['warnings']),
        })
    elif overlap.get('valid'):
        checks.append({
            'name': 'batch_cluster_overlap',
            'status': 'pass',
            'message': '未发现达到阈值的单一 batch 支配 cluster。',
        })

    return {
        'status': 'warning' if warnings else 'pass',
        'warnings': list(dict.fromkeys(warnings)),
        'checks': checks,
        'pca_hvg_only': bool(pca_hvg_only),
        'pca_n_genes': int(pca_n_genes),
        'selected_n_pcs': int(n_comps),
        'pc1_variance_pct': pc1_pct,
        'cumulative_variance_pct': cumulative_pct,
        'figure_quality_score_source': 'per-figure Nature readiness report',
    }


class DimredAnalysis(BaseAnalysis):
    MODULE_NAME = "dimred"
    DISPLAY_NAME = "降维分析"
    DESCRIPTION = "PCA、UMAP/t-SNE 降维，支持自动选 PC"
    INPUT_REQUIRES = []

    def run(self, input_path):
        import copy
        import os
        import scanpy as sc
        import omicverse as ov
        import numpy as np
        from figure_engine import NatureFigureDirector

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        requested_batch_key = str(self.params.get('batch_key', '') or '').strip()
        batch_key = None
        if requested_batch_key:
            batch_key, batch_info = resolve_obs_grouping(
                adata, requested_batch_key, max_categories=50,
                max_numeric_categories=20, require_multiple=True,
            )
            if not batch_info.get('requested_valid', False):
                self.progress(
                    -1,
                    f"批次列已跳过：{batch_info.get('requested_reason', '不是有效分类列')}",
                )
        requested_max_pcs = int(self.params.get('n_comps', DEFAULT_FINAL_N_PCS))
        if requested_max_pcs < 2:
            raise ValueError('最大候选 PC 数必须至少为 2')
        use_mde = self.params.get('use_mde', False)
        auto_n_comps = str(self.params.get('auto_n_comps', 'elbow') or 'elbow').lower()
        if auto_n_comps not in {'none', 'elbow', 'kneedle'}:
            raise ValueError(f'不支持的 PC 选择方法: {auto_n_comps}')
        pca_hvg_only = bool(self.params.get('pca_hvg_only', True))

        # UMAP parameters
        umap_n_neighbors = int(self.params.get('umap_n_neighbors', 15))
        umap_min_dist = float(self.params.get('umap_min_dist', 0.5))
        umap_metric = self.params.get('umap_metric', 'euclidean')
        umap_spread = float(self.params.get('umap_spread', 1.0))

        # t-SNE parameters
        enable_tsne = self.params.get('enable_tsne', False)
        tsne_perplexity = float(self.params.get('tsne_perplexity', 30))
        tsne_learning_rate = float(self.params.get('tsne_learning_rate', 1000))

        pca_n_genes = int(adata.n_vars)
        hvg_mask = None
        if pca_hvg_only and 'highly_variable' in adata.var.columns:
            hvg_mask = np.asarray(adata.var['highly_variable'].fillna(False), dtype=bool)
            if int(hvg_mask.sum()) < 2:
                hvg_mask = None
                self.progress(-1, "HVG 数量不足，PCA 回退到全基因矩阵")

        if hvg_mask is not None:
            pca_n_genes = int(hvg_mask.sum())
        max_available_comps = min(adata.n_obs - 1, pca_n_genes - 1)
        if max_available_comps < 2:
            raise ValueError(
                f'PCA 至少需要 3 个细胞和 3 个输入基因，当前为 '
                f'{adata.n_obs} cells × {pca_n_genes} genes'
            )
        candidate_n_comps = min(requested_max_pcs, max_available_comps)
        if requested_max_pcs > max_available_comps:
            self.progress(
                -1,
                f'请求最大 {requested_max_pcs} PCs 超过数据上限，已限制为 {max_available_comps} PCs',
            )

        stale_pca = _invalidate_pca_results(adata)
        if stale_pca:
            logger.info('Removed stale PCA artifacts before recomputation: %s', stale_pca)

        if hvg_mask is not None:
            self.progress(20, f"Scaling {pca_n_genes} HVGs for memory-safe PCA...")
            pca_adata = adata[:, hvg_mask].copy()
            _invalidate_pca_results(pca_adata)
            ov.pp.scale(pca_adata, max_value=10)
            self.progress(35, f"Running PCA (maximum {candidate_n_comps} components, HVG-only)...")
            try:
                # The object still contains ``highly_variable`` metadata.  Make
                # the intended input explicit so Scanpy cannot apply a second,
                # implicit mask or reuse the old representation.
                sc.pp.pca(
                    pca_adata, n_comps=candidate_n_comps,
                    layer='scaled', mask_var=None,
                )
            except TypeError:
                if 'scaled' in pca_adata.layers:
                    pca_adata.X = pca_adata.layers['scaled']
                sc.pp.pca(
                    pca_adata, n_comps=candidate_n_comps,
                    use_highly_variable=False,
                )
            adata.obsm['X_pca'] = np.asarray(pca_adata.obsm['X_pca'])
            adata.uns['pca'] = copy.deepcopy(pca_adata.uns['pca'])
            if 'PCs' in pca_adata.varm:
                full_loadings = np.zeros(
                    (adata.n_vars, pca_adata.varm['PCs'].shape[1]), dtype=np.float32,
                )
                full_loadings[hvg_mask, :] = np.asarray(
                    pca_adata.varm['PCs'], dtype=np.float32,
                )
                adata.varm['PCs'] = full_loadings
            del pca_adata
        else:
            self.progress(20, "Scaling data...")
            ov.pp.scale(adata, max_value=10)

            self.progress(35, f"Running PCA (maximum {candidate_n_comps} components)...")
            try:
                # Scanpy otherwise auto-detects adata.var['highly_variable']
                # even when the user explicitly requested full-gene PCA.
                sc.pp.pca(
                    adata, n_comps=candidate_n_comps,
                    layer='scaled', mask_var=None,
                )
            except TypeError:
                # 旧版 scanpy 没有 mask_var 参数：用 layer='scaled' +
                # use_highly_variable=False 回退。切勿把 adata.X 永久替换为
                # scaled 值——那会把零中心化（含负值）的矩阵写进输出
                # h5ad，污染下游 marker 检测率与 DEG 表达。
                sc.pp.pca(
                    adata, n_comps=candidate_n_comps,
                    layer='scaled', use_highly_variable=False,
                )

        pca_result = np.asarray(adata.obsm['X_pca'])
        pca_fingerprint = _pca_fingerprint(pca_result)
        pca_variance_ratio = np.asarray(
            adata.uns.get('pca', {}).get('variance_ratio', []), dtype=float,
        )
        pca_diagnostics = {
            'input_n_cells': int(adata.n_obs),
            'input_n_genes': int(pca_n_genes),
            'hvg_only': bool(hvg_mask is not None),
            'requested_max_pcs': int(requested_max_pcs),
            'candidate_n_pcs': int(candidate_n_comps),
            'pc_selection_method': auto_n_comps,
            'selected_n_pcs': int(candidate_n_comps),
            'result_shape': [int(value) for value in pca_result.shape],
            'pc1_variance_ratio': (
                float(pca_variance_ratio[0]) if len(pca_variance_ratio) else None
            ),
            'fingerprint_first_10_pcs': pca_fingerprint,
        }
        adata.uns['dimred_pca'] = copy.deepcopy(pca_diagnostics)
        logger.info(
            'PCA input: n_cells=%d, n_genes=%d, hvg_only=%s; '
            'result shape=%s, PC1 variance=%.6f, fingerprint=%s',
            adata.n_obs, pca_n_genes, bool(hvg_mask is not None),
            tuple(pca_result.shape),
            float(pca_variance_ratio[0]) if len(pca_variance_ratio) else float('nan'),
            pca_fingerprint,
        )
        self.progress(
            38,
            f"PCA 完成：输入 {pca_n_genes} genes，结果 {tuple(pca_result.shape)}，"
            f"fingerprint={pca_fingerprint[:12]}",
        )

        # Auto-select number of PCs
        selected_n_comps = candidate_n_comps
        if auto_n_comps != 'none' and 'pca' in adata.uns:
            variance_ratio = adata.uns['pca']['variance_ratio']
            if auto_n_comps == 'elbow':
                n_comps_auto = _select_elbow_n_comps(variance_ratio, candidate_n_comps)
                self.progress(40, f"Auto-selected {n_comps_auto} PCs (elbow)")
                selected_n_comps = n_comps_auto
            elif auto_n_comps == 'kneedle':
                try:
                    from kneed import KneeLocator
                    available = min(candidate_n_comps, len(variance_ratio))
                    x = range(1, available + 1)
                    kl = KneeLocator(
                        x, variance_ratio[:available], curve='convex', direction='decreasing'
                    )
                    if kl.knee:
                        n_comps_auto = max(
                            min(5, available), min(int(kl.knee), available)
                        )
                        self.progress(40, f"Auto-selected {n_comps_auto} PCs (kneedle)")
                        selected_n_comps = n_comps_auto
                except ImportError:
                    self.progress(40, "kneed not installed, using maximum candidate PCs")

        pca_diagnostics['selected_n_pcs'] = int(selected_n_comps)
        adata.uns['dimred_pca'] = copy.deepcopy(pca_diagnostics)

        self.progress(50, "Computing neighbors...")
        sc.pp.neighbors(
            adata, n_pcs=selected_n_comps,
            n_neighbors=umap_n_neighbors, metric=umap_metric,
        )

        self.progress(65, "Computing UMAP...")
        embedding_method = 'umap'
        if use_mde:
            try:
                import pymde
                mde_input = np.asarray(adata.obsm['X_pca'])[:, :selected_n_comps]
                mde = pymde.preserve_neighbors(mde_input, embedding_dim=2, device='cpu')
                embedding = mde.embed(verbose=False)
                adata.obsm['X_mde'] = embedding.cpu().numpy() if hasattr(embedding, 'cpu') else embedding.numpy()
                adata.obsm['X_umap'] = adata.obsm['X_mde']
                embedding_method = 'mde'
            except ImportError:
                self.progress(66, "pymde not installed, falling back to UMAP...")
                sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread)
        else:
            sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread)

        # Optional t-SNE
        if enable_tsne:
            self.progress(72, "Computing t-SNE...")
            sc.tl.tsne(
                adata, perplexity=tsne_perplexity,
                learning_rate=int(tsne_learning_rate), n_pcs=selected_n_comps,
            )

        self.progress(80, "Generating embedding plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []

        embedding_color_keys = _embedding_color_keys(
            batch_key, 'total_counts', 'n_genes_by_counts', 'pct_counts_mt',
            'phase', 'doublet_score', 'leiden',
        )
        for color_key in embedding_color_keys:
            if color_key in adata.obs.columns:
                try:
                    fig_static = self.build_publication_umap(
                        adata, color_key, title=f'UMAP colored by {color_key}'
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig_static, plots_dir, f'dimred_umap_{color_key}.png',
                        'umap', f'UMAP by {color_key}'
                    ))
                    import matplotlib.pyplot as plt
                    plt.close(fig_static)
                except Exception as exc:
                    self.progress(-1, f'UMAP 静态导出失败：{exc}')

        # PCA scatter for early detection of outliers and batch/sample structure
        if self.params.get('show_pca_scatter', True) and 'X_pca' in adata.obsm:
            color_key = next((
                key for key in [batch_key, 'phase', 'n_genes_by_counts', 'total_counts']
                if key and key in adata.obs.columns
            ), None)
            pca_variance = np.asarray(adata.uns.get('pca', {}).get('variance_ratio', []), dtype=float)
            x_suffix = f' ({pca_variance[0] * 100:.1f}%)' if len(pca_variance) > 0 else ''
            y_suffix = f' ({pca_variance[1] * 100:.1f}%)' if len(pca_variance) > 1 else ''
            fig_pca = self.build_publication_embedding(
                adata, color_key or '', title='Cell-level PCA', basis='X_pca',
                x_label=f'PC1{x_suffix}', y_label=f'PC2{y_suffix}',
            )
            result_files.extend(self.save_matplotlib_figure(
                fig_pca, plots_dir, 'dimred_pca_scatter.png', 'pca',
                'Cell-level PCA',
            ))

        # t-SNE plot
        if enable_tsne and 'X_tsne' in adata.obsm:
            tsne_color_keys = _embedding_color_keys(batch_key, 'leiden')
            for color_key in tsne_color_keys:
                if color_key in adata.obs.columns:
                    fig = self.build_publication_embedding(
                        adata, color_key, title=f't-SNE by {color_key}', basis='X_tsne',
                        x_label='t-SNE 1', y_label='t-SNE 2',
                    )
                    result_files.extend(self.save_matplotlib_figure(
                        fig, plots_dir, f'dimred_tsne_{color_key}.png', 'tsne',
                        f't-SNE by {color_key}',
                    ))

        # Variance ratio plot
        if 'pca' in adata.uns:
            vr = adata.uns['pca']['variance_ratio'][:min(50, candidate_n_comps)]
            director = NatureFigureDirector()
            variance_spec = director.spec_from_params(
                'variance', self.params, width='single',
                title='PCA variance explained', evidence_role='quality_control',
            ).with_updates(height_mm=78.0)
            fig = director.render(variance_spec, {
                'kind': 'variance',
                'variance': vr,
                'labels': [f'PC{i + 1}' for i in range(len(vr))],
            })
            result_files.extend(self.save_matplotlib_figure(
                fig, plots_dir, 'dimred_pca_variance.png', 'pca',
                'PCA Variance Ratio',
            ))

            # The standard variance panel is retained for compatibility.  A
            # second diagnostic panel makes the selected PC count explicit and
            # shows cumulative variance on its own axis.
            from modules.sc_figure_diagnostics import pca_variance_diagnostic_figure
            fig_variance_diagnostic = pca_variance_diagnostic_figure(
                adata.uns['pca']['variance_ratio'], selected_n_comps,
            )
            if fig_variance_diagnostic is not None:
                result_files.extend(self.save_matplotlib_figure(
                    fig_variance_diagnostic, plots_dir,
                    'dimred_pca_variance_diagnostic.png', 'pca',
                    'PCA variance explained and selected PCs',
                    formats=('png', 'svg'), dpi=300,
                ))

        # Keep PCA donor and QC views separate.  A single automatic colour
        # choice can hide whether a separation is technical, biological, or a
        # quality gradient.
        if self.params.get('show_pca_scatter', True) and 'X_pca' in adata.obsm:
            pca_views = []
            if batch_key and batch_key in adata.obs.columns:
                pca_views.append((batch_key, f'PCA by {batch_key}', 'dimred_pca_by_donor.png'))
            qc_key = next((key for key in ('total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'phase')
                           if key in adata.obs.columns), None)
            if qc_key:
                pca_views.append((qc_key, f'PCA by {qc_key}', 'dimred_pca_by_qc.png'))
            for color_key, title, filename in pca_views:
                fig_view = self.build_publication_embedding(
                    adata, color_key, title=title, basis='X_pca',
                    x_label='PC1', y_label='PC2',
                )
                result_files.extend(self.save_matplotlib_figure(
                    fig_view, plots_dir, filename, 'pca', title,
                ))

        if 'X_pca' in adata.obsm and 'PCs' in adata.varm:
            from modules.sc_figure_diagnostics import pca_loading_figure
            loadings = np.asarray(adata.varm['PCs'])
            loading_genes = np.asarray(adata.var_names.astype(str))
            if hvg_mask is not None and loadings.shape[0] == adata.n_vars:
                loadings = loadings[hvg_mask]
                loading_genes = loading_genes[hvg_mask]
            fig_loadings = pca_loading_figure(
                loadings, loading_genes,
                variance_ratio=adata.uns.get('pca', {}).get('variance_ratio', []),
                n_pcs=min(5, selected_n_comps), top_n=8,
            )
            if fig_loadings is not None:
                result_files.extend(self.save_matplotlib_figure(
                    fig_loadings, plots_dir, 'dimred_pca_loading_genes.png',
                    'pca', 'Top PCA loading genes', formats=('png', 'svg'), dpi=300,
                ))

        from modules.io_utils import batch_cluster_overlap
        batch_overlap = batch_cluster_overlap(
            adata, 'leiden', batch_key,
        ) if batch_key else None
        if batch_overlap and batch_overlap.get('warnings'):
            for warning in batch_overlap['warnings']:
                self.progress(-1, warning)
        analytical_qc = _analytical_qc_summary(
            adata,
            pca_hvg_only=bool(hvg_mask is not None),
            pca_n_genes=pca_n_genes,
            n_comps=selected_n_comps,
            batch_overlap=batch_overlap,
        )
        batch_mixing = None
        if batch_key:
            # A pre-correction batch-coloured embedding is diagnostic evidence,
            # not a successful mixing assessment.  Keep the two states explicit
            # so a generated figure cannot make a strongly batch-separated PCA
            # look like a PASS.
            batch_mixing = {
                'status': 'review',
                'message': (
                    '已生成 Batch UMAP，但降维阶段未完成批次混合质量评价；'
                    '请在批次校正后查看共享抽样的 ASW、邻居混合和 cluster 组成。'
                ),
                'batch_key': batch_key,
                'source': 'pre_correction_embedding_diagnostic',
            }
            analytical_qc['batch_mixing'] = copy.deepcopy(batch_mixing)
        for warning in analytical_qc['warnings']:
            if not batch_overlap or warning not in batch_overlap.get('warnings', []):
                self.progress(-1, warning)

        self.progress(90, "Saving output...")
        output_path = self.save_output(adata, 'dimred')

        self.progress(100, "Done")
        return {
            'output_adata': output_path,
            'result_files': result_files,
            'summary': {
                'n_cells': adata.n_obs,
                'requested_max_pcs': requested_max_pcs,
                'candidate_n_pcs': candidate_n_comps,
                'pc_selection_method': auto_n_comps,
                'selected_n_pcs': selected_n_comps,
                # Backward-compatible alias; new consumers should use selected_n_pcs.
                'n_pcs': selected_n_comps,
                'pca_variance_ratio_top5': round(float(adata.uns['pca']['variance_ratio'][:5].sum()), 6) if 'pca' in adata.uns else None,
                'pca_pc1_variance_ratio': (
                    round(float(pca_variance_ratio[0]), 6)
                    if len(pca_variance_ratio) else None
                ),
                'pca_fingerprint': pca_fingerprint,
                'pca_result_shape': pca_diagnostics['result_shape'],
                'embedding_method': embedding_method,
                'tsne_enabled': enable_tsne,
                'pca_hvg_only': bool(hvg_mask is not None),
                'pca_n_genes': pca_n_genes,
                'requested_batch_key': requested_batch_key,
                'batch_key': batch_key,
                'batch_mixing': batch_mixing,
                'analytical_qc': analytical_qc,
            }
        }
