from modules.base import BaseAnalysis
from modules.constants import S_GENES, G2M_GENES
from modules.io_utils import resolve_obs_grouping, restore_scanpy_qc_percentages
from modules.normalize import _is_raw_count_matrix


QC_OBS_FIELD_ALIASES = {
    'total_counts': ['nUMIs', 'n_counts'],
    'n_genes_by_counts': ['detected_genes', 'n_genes'],
    'pct_counts_mt': ['mito_perc'],
    'pct_counts_ribo': ['ribo_perc'],
    'pct_counts_hb': ['hb_perc'],
}


def canonicalize_qc_obs_fields(adata):
    """Keep one canonical QC field while preserving the alias contract in metadata."""
    for canonical, aliases in QC_OBS_FIELD_ALIASES.items():
        if canonical not in adata.obs.columns:
            source = next((alias for alias in aliases if alias in adata.obs.columns), None)
            if source is not None:
                adata.obs[canonical] = adata.obs[source]
        drop_columns = [
            alias for alias in aliases
            if alias != canonical and alias in adata.obs.columns
        ]
        if drop_columns:
            adata.obs.drop(columns=drop_columns, inplace=True)
    return adata


def _as_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {'1', 'true', 'yes', 'on'}
    return bool(value)


def _ribosomal_gene_mask(gene_names):
    """Identify cytosolic ribosomal-protein symbols without kinase prefixes.

    A plain ``startswith(('RPS', 'RPL'))`` also matches symbols such as
    ``RPS6KA1``/``RPS6KB1`` (ribosomal S6 kinases).  Match the established RPS
    and RPL symbol forms instead, including RPSA, RPLP0-2 and common paralogs.
    """
    import re

    pattern = re.compile(
        r'^(?:RPS(?:A|\d+(?:AL|A|L\d*|X|Y[12]?)?)|'
        r'RPL(?:P[0-2]|\d+(?:AL|A|L\d*)?))$',
        flags=re.IGNORECASE,
    )
    return [bool(pattern.fullmatch(str(name).strip())) for name in gene_names]


def _first_gene_name_to_var_name(var_names, gene_names):
    """Map each gene symbol deterministically without overwriting duplicates."""
    mapping = {}
    for var_name, gene_name in zip(var_names, gene_names):
        symbol = str(gene_name).strip()
        if not symbol or symbol.lower() == 'nan':
            continue
        mapping.setdefault(symbol, var_name)
    return mapping


def _qc_overview_grid(n_metrics):
    """Return a compact grid that has room for every QC overview metric."""
    n_metrics = max(int(n_metrics), 1)
    ncols = 1 if n_metrics == 1 else 2
    nrows = (n_metrics + ncols - 1) // ncols
    return nrows, ncols


def _qc_removal_counts(n_before, n_after_initial, n_after, stage_removed):
    """Build complete and backward-compatible QC removal statistics."""
    initial = max(int(n_before) - int(n_after_initial), 0)
    total = max(int(n_before) - int(n_after), 0)
    additional = max(int(n_after_initial) - int(n_after), 0)
    return {
        # Backward-compatible field: now covers the complete QC pipeline.
        'cells_removed_by_qc_and_doublet': total,
        'cells_removed_by_initial_qc_and_doublet': initial,
        'cells_removed_by_additional_qc': additional,
        'cells_removed_by_max_genes': int(stage_removed.get('max_genes', 0)),
        'cells_removed_by_ribo': int(stage_removed.get('ribo', 0)),
        'cells_removed_by_hb': int(stage_removed.get('hb', 0)),
        'cells_removed_by_batch_adaptive_qc': int(
            stage_removed.get('batch_adaptive_qc', 0)
        ),
    }


def _finite_scores(values):
    """Flatten to finite float scores (tolerates lists and numpy arrays)."""
    import numpy as np

    arr = np.asarray(values, dtype=float).ravel()
    return arr[np.isfinite(arr)]


def _plain_scalar(value):
    """Convert numpy scalars to plain Python scalars for JSON manifests."""
    if isinstance(value, (int, float, str)) and not isinstance(value, bool):
        return value
    if hasattr(value, 'item'):
        return value.item()
    return value


SCDBLFINDER_PACKAGE = 'pyscdblfinder'
SCDBLFINDER_MIN_VERSION = '0.2.0'


def _scdblfinder_dependency_status(version_getter=None):
    """Return whether the supported pyscdblfinder backend is available.

    OmicVerse deliberately falls back to Scrublet when this optional package
    cannot be imported.  That is useful in a library, but is unsafe
    provenance for an analysis workbench where the user explicitly selected
    scDblFinder.  Keep the check at the platform boundary and make its result
    serializable for the QC manifest.
    """
    from importlib.metadata import PackageNotFoundError, version
    from packaging.version import InvalidVersion, Version

    version_getter = version_getter or version
    try:
        installed_version = str(version_getter(SCDBLFINDER_PACKAGE))
    except PackageNotFoundError:
        return {
            'available': False,
            'package': SCDBLFINDER_PACKAGE,
            'minimum_version': SCDBLFINDER_MIN_VERSION,
            'installed_version': None,
            'reason': 'pyscdblfinder_missing_or_version_below_0.2.0',
        }
    try:
        available = Version(installed_version) >= Version(SCDBLFINDER_MIN_VERSION)
    except InvalidVersion:
        available = False
    return {
        'available': bool(available),
        'package': SCDBLFINDER_PACKAGE,
        'minimum_version': SCDBLFINDER_MIN_VERSION,
        'installed_version': installed_version,
        'reason': None if available else 'pyscdblfinder_missing_or_version_below_0.2.0',
    }


def _effective_omicverse_doublets_method(adata):
    """Read the method OmicVerse says it actually resolved and ran."""
    status_args = adata.uns.get('status_args', {}) if hasattr(adata, 'uns') else {}
    qc_args = status_args.get('qc', {}) if isinstance(status_args, dict) else {}
    method = qc_args.get('doublets_method') if isinstance(qc_args, dict) else None
    if method is None:
        return None
    method = str(method).strip().lower()
    return method or None


def _doublet_method_provenance(requested, effective, *, verification_source,
                               dependency=None):
    """Build explicit requested-versus-effective doublet caller provenance."""
    requested = str(requested or 'none').strip().lower() or 'none'
    effective = str(effective or 'unknown').strip().lower() or 'unknown'
    fallback = requested != effective
    return {
        'requested_doublets_method': requested,
        'effective_doublets_method': effective,
        'fallback': fallback,
        'fallback_reason': (
            'omicverse_resolved_a_different_doublet_method' if fallback else None
        ),
        'verification_source': verification_source,
        'dependency': dependency,
    }


def _doublet_method_display(method):
    """Use caller names that are clear in user-facing figure labels."""
    labels = {
        'scdblfinder': 'scDblFinder',
        'scrublet': 'Scrublet',
        'sccomposite': 'scComposite',
        'doubletfinder': 'DoubletFinder',
        'none': 'No doublet caller',
    }
    return labels.get(str(method or '').lower(), str(method or 'Doublet caller'))


def _sim_score_distribution_summary(sim_scores, threshold=None):
    """JSON-safe summary of the simulated-doublet score distribution.

    Scrublet chooses the automatic threshold from the simulated doublet score
    histogram, so the simulated scores are the decisive evidence for judging
    that threshold.  Full arrays stay in adata.uns['scrublet']; the manifest
    receives quantiles + a fixed-bin histogram instead of dumping every score.
    """
    import numpy as np

    arr = _finite_scores(sim_scores)
    empty = {
        'n_simulated_doublets': 0,
        'sim_score_quantiles': None,
        'sim_score_histogram': None,
        'detectable_doublet_fraction': None,
    }
    if arr.size == 0:
        return empty
    try:
        threshold_value = float(threshold) if threshold is not None else None
    except (TypeError, ValueError):
        threshold_value = None
    quantile_names = ['min', 'p01', 'p05', 'p10', 'p25', 'median',
                      'p75', 'p90', 'p95', 'p99', 'max']
    quantiles = np.percentile(arr, [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100])
    counts, edges = np.histogram(arr, bins=30)
    detectable = None
    if threshold_value is not None:
        detectable = float(np.mean(arr > threshold_value))
    return {
        'n_simulated_doublets': int(arr.size),
        'sim_score_quantiles': {
            name: round(float(value), 6)
            for name, value in zip(quantile_names, quantiles)
        },
        'sim_score_histogram': {
            'edges': [round(float(edge), 6) for edge in edges],
            'counts': [int(count) for count in counts],
        },
        'detectable_doublet_fraction': (
            round(detectable, 6) if detectable is not None else None
        ),
    }


def _estimated_doublet_rate(detected_rate, detectable_fraction):
    """Estimate overall multiplet rate: detected_rate / detectable_fraction.

    Scrublet cannot recognise every multiplet (homotypic doublets can look like
    singletons), so the detected rate alone understates the real multiplet
    rate.  None means the estimate is unavailable (no simulated scores, or no
    simulated score reached the threshold).
    """
    if detected_rate is None or detectable_fraction is None:
        return None
    if detectable_fraction <= 0:
        return None
    return round(float(detected_rate) / float(detectable_fraction), 6)


def _doublet_review_status(detected_rate):
    """Flag unusually low detected doublet rates for manual review.

    A very low detected rate is a review signal, not proof the caller failed:
    the threshold and the simulated score distribution must be inspected before
    deciding whether the call is reasonable.
    """
    if detected_rate is None:
        return 'review'
    return 'review' if detected_rate < 0.005 else 'pass'


def _scrublet_summary(adata, batch_key, doublet_threshold=None,
                      scrublet_payload=None, run_params=None):
    """Return JSON-safe global and per-batch Scrublet evidence.

    run_params carries the parameters the platform itself passed to
    ov.pp.scrublet (expected_doublet_rate, sim_doublet_ratio, n_prin_comps,
    random_state, threshold_method), so the manifest records exactly what was
    controlled instead of inheriting library defaults implicitly.
    """
    labels = adata.obs.get('predicted_doublet')
    scrublet_payload = scrublet_payload if isinstance(scrublet_payload, dict) else {}
    batch_payload = scrublet_payload.get('batches')
    if not isinstance(batch_payload, dict):
        batch_payload = {}
    run_params = run_params if isinstance(run_params, dict) else {}
    try:
        threshold = float(doublet_threshold) if doublet_threshold is not None else None
    except (TypeError, ValueError):
        threshold = None
    if labels is None:
        return {
            'available': False,
            'method': 'scrublet',
            'threshold': threshold,
            'threshold_method': run_params.get('threshold_method'),
            'predicted_doublets': 0,
            'n_cells_evaluated': int(adata.n_obs),
            'doublet_rate': 0.0,
            'detected_doublet_rate': 0.0,
            'expected_doublet_rate': run_params.get('expected_doublet_rate'),
            'sim_doublet_ratio': run_params.get('sim_doublet_ratio'),
            'n_neighbors': run_params.get('n_neighbors'),
            'n_prin_comps': run_params.get('n_prin_comps'),
            'random_state': run_params.get('random_state'),
            'status': 'review',
            'by_batch': {},
        }
    labels = labels.fillna(False).astype(bool)
    n_cells = int(len(labels))
    detected_rate = round(float(labels.mean()) if n_cells else 0.0, 6)
    summary = {
        'available': True,
        'method': 'scrublet',
        'threshold': threshold,
        'threshold_method': run_params.get('threshold_method'),
        'predicted_doublets': int(labels.sum()),
        'n_cells_evaluated': n_cells,
        # Backward-compatible alias kept for existing consumers; the canonical
        # name is detected_doublet_rate (it is a detection rate, not the true
        # multiplet rate of the experiment).
        'doublet_rate': detected_rate,
        'detected_doublet_rate': detected_rate,
        'expected_doublet_rate': run_params.get('expected_doublet_rate'),
        'sim_doublet_ratio': run_params.get('sim_doublet_ratio'),
        'n_neighbors': run_params.get('n_neighbors'),
        'n_prin_comps': run_params.get('n_prin_comps'),
        'random_state': run_params.get('random_state'),
        'by_batch': {},
        'status': _doublet_review_status(detected_rate),
    }

    def _batch_parameters(metadata):
        parameters = metadata.get('parameters')
        return parameters if isinstance(parameters, dict) else {}

    if batch_key and batch_key in adata.obs.columns:
        for batch, indices in adata.obs.groupby(batch_key, observed=True).groups.items():
            batch_labels = labels.loc[indices]
            n_batch = int(len(batch_labels))
            batch_metadata = batch_payload.get(batch)
            if batch_metadata is None:
                batch_metadata = batch_payload.get(str(batch), {})
            if not isinstance(batch_metadata, dict):
                batch_metadata = {}
            parameters = _batch_parameters(batch_metadata)
            batch_threshold = batch_metadata.get('threshold')
            if batch_threshold is None:
                batch_threshold = parameters.get('threshold')
            try:
                batch_threshold = float(batch_threshold) if batch_threshold is not None else None
            except (TypeError, ValueError):
                batch_threshold = None
            batch_detected_rate = round(
                float(batch_labels.mean()) if n_batch else 0.0, 6,
            )
            sim_summary = _sim_score_distribution_summary(
                batch_metadata.get('doublet_scores_sim'), batch_threshold,
            )
            detectable = sim_summary.get('detectable_doublet_fraction')
            estimated = _estimated_doublet_rate(batch_detected_rate, detectable)
            summary['by_batch'][str(batch)] = {
                'n_cells_evaluated': n_batch,
                'predicted_doublets': int(batch_labels.sum()),
                'doublet_rate': batch_detected_rate,
                'detected_doublet_rate': batch_detected_rate,
                'threshold': batch_threshold if batch_threshold is not None else threshold,
                'threshold_method': run_params.get('threshold_method'),
                'expected_doublet_rate': _plain_scalar(
                    parameters.get('expected_doublet_rate')
                    if parameters.get('expected_doublet_rate') is not None
                    else run_params.get('expected_doublet_rate')
                ),
                'sim_doublet_ratio': _plain_scalar(
                    parameters.get('sim_doublet_ratio')
                    if parameters.get('sim_doublet_ratio') is not None
                    else run_params.get('sim_doublet_ratio')
                ),
                'n_neighbors': _plain_scalar(
                    parameters.get('n_neighbors')
                    if parameters.get('n_neighbors') is not None
                    else run_params.get('n_neighbors')
                ),
                'random_state': _plain_scalar(
                    parameters.get('random_state')
                    if parameters.get('random_state') is not None
                    else run_params.get('random_state')
                ),
                'status': _doublet_review_status(batch_detected_rate),
                **sim_summary,
                'estimated_overall_doublet_rate': estimated,
            }
    if summary['by_batch']:
        summary['threshold_by_batch'] = {
            batch: values.get('threshold')
            for batch, values in summary['by_batch'].items()
        }
        # Global simulated evidence: pooled quantiles are meaningless across
        # donor-specific thresholds, so aggregate the per-batch summary by the
        # number of simulated doublets actually generated.
        detect_entries = [
            (values.get('n_simulated_doublets', 0), values.get('detectable_doublet_fraction'))
            for values in summary['by_batch'].values()
            if values.get('detectable_doublet_fraction') is not None
        ]
        if detect_entries:
            sim_weight = sum(weight for weight, _ in detect_entries)
            if sim_weight:
                summary['n_simulated_doublets'] = sim_weight
                summary['detectable_doublet_fraction'] = round(
                    sum(weight * fraction for weight, fraction in detect_entries) / sim_weight, 6,
                )
                summary['estimated_overall_doublet_rate'] = _estimated_doublet_rate(
                    detected_rate, summary['detectable_doublet_fraction'],
                )
    else:
        # Unbatched run: the payload stores the simulated scores at top level.
        sim_summary = _sim_score_distribution_summary(
            scrublet_payload.get('doublet_scores_sim'), threshold,
        )
        summary.update(sim_summary)
        summary['estimated_overall_doublet_rate'] = _estimated_doublet_rate(
            detected_rate, summary.get('detectable_doublet_fraction'),
        )
        top_parameters = scrublet_payload.get('parameters')
        if isinstance(top_parameters, dict):
            for key in ('expected_doublet_rate', 'sim_doublet_ratio',
                        'n_neighbors', 'random_state'):
                value = top_parameters.get(key)
                if value is not None and summary.get(key) is None:
                    summary[key] = _plain_scalar(value)
    return summary


def _caller_doublet_summary(adata, batch_key, method, *, score_column='doublet_score',
                            class_column='predicted_doublet', run_params=None):
    """Summarize a non-Scrublet caller without inventing Scrublet evidence.

    scDblFinder exposes a classifier score and class, not Scrublet's simulated
    score distribution or an automatic Scrublet threshold.  Keeping these
    summaries separate prevents downstream reports from attaching a simulated
    threshold interpretation to a different algorithm.
    """
    import numpy as np
    import pandas as pd

    run_params = run_params if isinstance(run_params, dict) else {}
    labels = adata.obs.get(class_column)
    if labels is None and class_column != 'predicted_doublet':
        labels = adata.obs.get('predicted_doublet')
    if labels is None:
        labels = pd.Series(False, index=adata.obs_names)
    labels = labels.fillna(False).astype(bool)
    scores = adata.obs.get(score_column)
    if scores is None and score_column != 'doublet_score':
        scores = adata.obs.get('doublet_score')
    n_cells = int(len(labels))
    detected_rate = round(float(labels.mean()) if n_cells else 0.0, 6)
    summary = {
        'available': score_column in adata.obs.columns or class_column in adata.obs.columns,
        'method': str(method),
        'score_column': score_column if score_column in adata.obs.columns else None,
        'class_column': class_column if class_column in adata.obs.columns else 'predicted_doublet',
        'threshold': None,
        'score_threshold': None,
        'threshold_method': run_params.get('threshold_method', 'model_classification'),
        'predicted_doublets': int(labels.sum()),
        'n_cells_evaluated': n_cells,
        'doublet_rate': detected_rate,
        'detected_doublet_rate': detected_rate,
        'expected_doublet_rate': run_params.get('expected_doublet_rate'),
        'estimated_doublet_rate': run_params.get('estimated_doublet_rate'),
        'by_batch': {},
        'status': _doublet_review_status(detected_rate),
    }
    if batch_key and batch_key in adata.obs.columns:
        for batch, indices in adata.obs.groupby(batch_key, observed=True).groups.items():
            batch_labels = labels.loc[indices]
            n_batch = int(len(batch_labels))
            batch_rate = round(float(batch_labels.mean()) if n_batch else 0.0, 6)
            item = {
                'n_cells_evaluated': n_batch,
                'predicted_doublets': int(batch_labels.sum()),
                'doublet_rate': batch_rate,
                'detected_doublet_rate': batch_rate,
                'threshold': None,
                'score_threshold': None,
                'threshold_method': summary['threshold_method'],
                'status': _doublet_review_status(batch_rate),
            }
            if scores is not None:
                batch_scores = _finite_scores(scores.loc[indices])
                if batch_scores.size:
                    item['score_quantiles'] = {
                        'min': round(float(np.min(batch_scores)), 6),
                        'median': round(float(np.median(batch_scores)), 6),
                        'p95': round(float(np.percentile(batch_scores, 95)), 6),
                        'p99': round(float(np.percentile(batch_scores, 99)), 6),
                        'max': round(float(np.max(batch_scores)), 6),
                    }
            summary['by_batch'][str(batch)] = item
    return summary


def _scdblfinder_summary(adata, batch_key):
    """Return scDblFinder-specific evidence using its classifier contract."""
    return _caller_doublet_summary(
        adata,
        batch_key,
        'scdblfinder',
        score_column='scdblfinder_score',
        class_column='scdblfinder_doublet',
        run_params={
            # OmicVerse's wrapper exposes classifier calls but not a single
            # global score cutoff or a rate estimator.  Record that absence
            # rather than presenting a fabricated Scrublet-style threshold.
            'threshold_method': 'scdblfinder_model_classification',
        },
    )


def _gene_detection_counts(matrix, cell_mask):
    """Count non-zero cells per gene without densifying sparse matrices."""
    import numpy as np
    from scipy import sparse

    selected = matrix[np.asarray(cell_mask, dtype=bool)]
    if sparse.issparse(selected):
        return np.asarray(selected.getnnz(axis=0)).ravel().astype(int)
    return np.count_nonzero(np.asarray(selected), axis=0).astype(int)


def _qc_removal_reason_summary(qc_before, after_initial, after_final, *,
                               mito_perc, n_umis, n_genes_min,
                               n_genes_max, ribo_perc_max, hb_perc_max,
                               doublet_index=None, stage_removed=None):
    """Summarize non-exclusive removal reasons and their overlaps.

    Threshold reasons are flags, so their sum is intentionally allowed to be
    larger than the unique removed-cell count.  The overlap field makes that
    distinction explicit instead of implying that reasons partition cells.
    """
    import pandas as pd

    before_index = pd.Index(qc_before.index)
    initial_index = pd.Index(
        after_initial.obs_names if hasattr(after_initial, 'obs_names') else after_initial
    )
    final_index = pd.Index(after_final.obs_names)
    initial_removed = before_index.difference(initial_index)
    additional_removed = initial_index.difference(final_index)

    def _numeric(column):
        if column not in qc_before.columns:
            return pd.Series(float('nan'), index=qc_before.index)
        return pd.to_numeric(qc_before[column], errors='coerce')

    total_counts = _numeric('total_counts')
    detected_genes = _numeric('n_genes_by_counts')
    mt_pct = _numeric('pct_counts_mt')
    ribo_pct = _numeric('pct_counts_ribo')
    hb_pct = _numeric('pct_counts_hb')
    initial_mask = pd.Series(qc_before.index.isin(initial_removed), index=qc_before.index)
    doublet_index = pd.Index(doublet_index) if doublet_index is not None else pd.Index([])
    reasons = {
        'removed_low_umi': int((initial_mask & (total_counts < float(n_umis))).sum()),
        'removed_low_genes': int((initial_mask & (detected_genes < float(n_genes_min))).sum()),
        'removed_high_mt': int((initial_mask & (mt_pct > float(mito_perc) * 100.0)).sum()),
        'removed_high_genes': int((
            int((detected_genes.loc[additional_removed] > float(n_genes_max)).sum())
            if n_genes_max > 0 else 0
        )),
        'removed_high_ribo': int((
            int((ribo_pct.loc[additional_removed] > float(ribo_perc_max)).sum())
            if ribo_perc_max > 0 else 0
        )),
        'removed_high_hb': int((
            int((hb_pct.loc[additional_removed] > float(hb_perc_max)).sum())
            if hb_perc_max > 0 else 0
        )),
        'predicted_doublets': int(len(doublet_index)),
        'removed_doublets': int(len(doublet_index.intersection(initial_removed))),
        'removed_batch_adaptive_qc': int((stage_removed or {}).get('batch_adaptive_qc', 0)),
    }
    reasons['unique_removed_total'] = int(len(before_index.difference(final_index)))
    reason_flag_sum = sum(reasons[key] for key in (
        'removed_low_umi', 'removed_low_genes', 'removed_high_mt',
        'removed_high_genes', 'removed_high_ribo', 'removed_high_hb',
        'removed_doublets', 'removed_batch_adaptive_qc',
    ))
    reasons['reason_flag_sum'] = int(reason_flag_sum)
    reasons['reason_overlap_cells'] = max(
        int(reason_flag_sum - reasons['unique_removed_total']), 0,
    )
    reasons['initial_removed_total'] = int(len(initial_removed))
    reasons['additional_removed_total'] = int(len(additional_removed))
    return reasons


def _normalized_cell_cycle_copy(adata, scanpy_module):
    """Return a normalized/log1p copy suitable for cell-cycle scoring.

    ``scanpy.pp.calculate_qc_metrics(log1p=True)`` creates observation columns
    such as ``log1p_total_counts`` but does not transform ``adata.X``. The
    presence of those columns therefore cannot be used to infer expression
    scale.
    """
    adata_cc = adata.copy()
    # X 不是整数计数时（重跑 QC / 已标准化输入）从 counts 层重建原始计数，
    # 避免对 log1p 值二次 normalize+log1p。
    if not _is_raw_count_matrix(adata_cc.X) and 'counts' in adata_cc.layers:
        if _is_raw_count_matrix(adata_cc.layers['counts']):
            adata_cc.X = adata_cc.layers['counts'].copy()
    scanpy_module.pp.normalize_total(adata_cc, target_sum=1e4)
    scanpy_module.pp.log1p(adata_cc)
    return adata_cc


class QCAnalysis(BaseAnalysis):
    MODULE_NAME = "qc"
    DISPLAY_NAME = "质控"
    DESCRIPTION = "MT/ribo/hb 过滤 + 双细胞检测 + 细胞周期评分 + 复杂度过滤"
    INPUT_REQUIRES = []

    def run(self, input_path):
        import numpy as np
        import os
        import pandas as pd
        import matplotlib.pyplot as plt
        from figure_engine import NatureFigureDirector
        from modules.figure_style import (
            NATURE_AXIS, NATURE_GRID, NATURE_MUTED, NATURE_PALETTE, NATURE_TEXT,
        )

        self.progress(5, "Loading data...")
        adata = self.load_adata(input_path)
        # 应用自定义过滤规则
        adata = self.apply_filters(adata, 'qc')

        # Check this before importing/running the analysis stack.  OmicVerse
        # otherwise resolves a missing scDblFinder backend to Scrublet only
        # after other QC work has already started.
        requested_doublets_method = str(
            self.params.get('doublets_method', 'scdblfinder') or 'scdblfinder'
        ).strip().lower()
        doublets_disabled = requested_doublets_method in {
            '', 'none', 'off', 'false', '0', 'disabled',
        }
        doublets_method = 'scrublet' if doublets_disabled else requested_doublets_method
        if doublets_method not in {'scrublet', 'sccomposite', 'doubletfinder', 'scdblfinder'}:
            raise ValueError(
                f"不支持的 doublets_method '{requested_doublets_method}'；"
                "请选择 scrublet、sccomposite、doubletfinder、scdblfinder 或 none。"
            )
        scdblfinder_dependency = None
        if not doublets_disabled and doublets_method == 'scdblfinder':
            scdblfinder_dependency = _scdblfinder_dependency_status()
            if not scdblfinder_dependency['available']:
                installed = scdblfinder_dependency.get('installed_version')
                installed_text = (
                    f"当前检测到 {SCDBLFINDER_PACKAGE}={installed}" if installed
                    else f"当前未安装 {SCDBLFINDER_PACKAGE}"
                )
                raise ValueError(
                    "当前环境无法运行 scDblFinder：需要 "
                    f"{SCDBLFINDER_PACKAGE}>={SCDBLFINDER_MIN_VERSION}（{installed_text}）。"
                    "请安装或升级该依赖，或主动将双细胞检测方法切换为 Scrublet。"
                )

        import scanpy as sc
        import omicverse as ov

        # QC metrics and doublet callers are defined on un-normalized UMI counts.
        # Keep the incoming expression representation so it can be restored
        # after filtering; otherwise rerunning QC on a normalized object would
        # silently replace its X matrix with counts.
        original_x = adata.X.copy() if hasattr(adata.X, 'copy') else adata.X
        original_obs_names = adata.obs_names.copy()
        original_var_names = adata.var_names.copy()
        if 'counts' in adata.layers:
            if not _is_raw_count_matrix(adata.layers['counts']):
                raise ValueError(
                    "layers['counts'] 不是非负整数原始 UMI counts；"
                    "QC/双细胞检测不能在该层上运行，请从保留原始 counts 的输入重新开始。"
                )
            qc_expression_source = 'layers[counts]'
            adata.X = adata.layers['counts'].copy()
        elif _is_raw_count_matrix(adata.X):
            qc_expression_source = 'X (validated raw counts)'
        else:
            raise ValueError(
                "QC/双细胞检测需要原始 UMI counts，但输入既没有 layers['counts']，"
                "且 X 不是非负整数 counts；请提供带 counts 层的 AnnData。"
            )

        # 保存 counts 层：已有 counts 层且 X 不是整数计数（例如在已标准化
        # 的中间文件上重跑 QC）时保留原 counts，避免被 log1p 值覆盖。
        save_counts = self.params.get('save_counts_layer', True)
        if save_counts and 'counts' not in adata.layers:
            adata.layers["counts"] = adata.X.copy()

        # ── 1. 标记基因集 ─────────────────────────────────────────────────
        self.progress(10, "Flagging MT/ribo/hb genes...")
        var_names_str = adata.var_names.astype(str).to_series(index=adata.var_names)
        annotation_names = None
        if 'gene_name' in adata.var.columns:
            annotation_names = adata.var['gene_name'].fillna('').astype(str)
        # Use both the feature index and an optional gene_name annotation.  A
        # number of 10x exports retain a stale/empty gene_name column even
        # though var_names already contain standard symbols; relying on only
        # that column silently produced pct_counts_ribo == 0 for valid data.
        mt_mask = var_names_str.str.upper().str.startswith('MT-')
        ribo_var_mask = _ribosomal_gene_mask(var_names_str)
        # 只把真正的地中海贫血珠蛋白基因记为 hemoglobin，避免 HBEGF 等
        # 以 HB 开头但非珠蛋白的基因（肝素结合 EGF 样生长因子）污染 hb%。
        _HB_GLOBINS = {
            'HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2',
            'HBM', 'HBQ1', 'HBZ', 'HBBP1', 'HBZP1',
        }
        hb_mask = var_names_str.str.upper().isin(_HB_GLOBINS)
        gene_name_source = 'var_names'
        if annotation_names is not None:
            annotation_upper = annotation_names.str.upper()
            mt_mask = mt_mask | annotation_upper.str.startswith('MT-')
            ribo_annotation_mask = _ribosomal_gene_mask(annotation_names)
            ribo_var_mask = [left or right for left, right in zip(
                ribo_var_mask, ribo_annotation_mask,
            )]
            hb_mask = hb_mask | annotation_upper.isin(_HB_GLOBINS)
            has_annotation = any(
                str(value).strip().lower() not in {'', 'nan', 'none'}
                for value in annotation_names
            )
            if has_annotation:
                gene_name_source = 'var_names_and_gene_name'
        adata.var["mt"] = mt_mask.to_numpy(dtype=bool)
        adata.var["ribo"] = ribo_var_mask
        adata.var["hb"] = hb_mask.to_numpy(dtype=bool)
        mt_genes = adata.var_names[adata.var["mt"]].astype(str).tolist()
        ribo_genes = adata.var_names[adata.var["ribo"]].astype(str).tolist()
        hb_genes = adata.var_names[adata.var["hb"]].astype(str).tolist()
        n_mt_genes_found = len(mt_genes)
        n_ribo_genes_found = len(ribo_genes)
        n_hb_genes_found = len(hb_genes)
        qc_metric_warnings = []
        if n_ribo_genes_found == 0:
            qc_metric_warnings.append(
                '未识别到 RPL/RPS 核糖体蛋白基因；pct_counts_ribo=0 不能解释为真实低核糖体比例，'
                '请检查 var_names/gene_name 注释。'
            )
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
            gene_name_to_idx = _first_gene_name_to_var_name(
                adata.var_names,
                adata.var['gene_name'].fillna('').astype(str),
            )
            s_in_ensembl = [gene_name_to_idx[g] for g in S_GENES if g in gene_name_to_idx]
            g2m_in_ensembl = [gene_name_to_idx[g] for g in G2M_GENES if g in gene_name_to_idx]
            if len(s_in_ensembl) > len(s_in):
                s_in = s_in_ensembl
            if len(g2m_in_ensembl) > len(g2m_in):
                g2m_in = g2m_in_ensembl

        if score_cell_cycle and len(s_in) >= 5 and len(g2m_in) >= 5:
            # QC 指标列不代表 X 已标准化；始终在副本上进行 normalize + log1p。
            adata_cc = _normalized_cell_cycle_copy(adata, sc)
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

        # ── 4. 双细胞检测 + 阈值过滤 ───────────────────────────────────
        self.progress(35, "Running doublet detection + QC filtering...")
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
        # Scrublet 不再由 ov.pp.qc 内部调用：基础 QC 与双细胞检测解耦，
        # expected_doublet_rate / sim_doublet_ratio / threshold / n_prin_comps
        # 等参数由平台显式控制、显式记录，而不是继承 OmicVerse 的隐式默认值。
        scrublet_decoupled = (not doublets_disabled) and doublets_method == 'scrublet'
        expected_doublet_rate = float(self.params.get('expected_doublet_rate', 0.05))
        sim_doublet_ratio = float(self.params.get('sim_doublet_ratio', 2.0))
        n_prin_comps = int(self.params.get('n_prin_comps', 30))
        scrublet_threshold = self.params.get('scrublet_threshold')
        if isinstance(scrublet_threshold, str):
            scrublet_threshold = scrublet_threshold.strip()
        if scrublet_threshold in (None, '', 'auto', 'none'):
            # threshold=None 表示让 Scrublet 根据 simulated doublet score
            # 分布自动找阈值；这是唯一能保留自动阈值语义的传法。
            scrublet_threshold_value = None
        else:
            scrublet_threshold_value = float(scrublet_threshold)
        scrublet_random_state = int(self.params.get('scrublet_random_state', 1234))
        scrublet_use_gpu = _as_bool(self.params.get('scrublet_use_gpu', False), False)

        # OmicVerse 的 min_genes 是“细胞至少检测到的基因数”过滤，会把
        # 低于该值的细胞再滤一轮。此前硬编码 200，会静默覆盖用户
        # detected_genes 阈值（用户填 50 时仍按 200 过滤且归因无法闭合）。
        # 改为与用户阈值一致，并保留支持计数用于归因。
        gene_filter_min_cells = 3
        gene_filter_max_cells_ratio = 1.0
        gene_filter_max_genes_ratio = 1.0
        gene_filter_min_genes = max(1, ngenes_min)
        qc_input_var_names = pd.Index(adata.var_names)
        gene_filter_cell_mask = np.ones(n_before, dtype=bool)
        if 'total_counts' in qc_before.columns:
            gene_filter_cell_mask &= (
                pd.to_numeric(qc_before['total_counts'], errors='coerce').to_numpy()
                > nUMIs_min
            )
        if 'n_genes_by_counts' in qc_before.columns:
            gene_filter_cell_mask &= (
                pd.to_numeric(qc_before['n_genes_by_counts'], errors='coerce').to_numpy()
                > gene_filter_min_genes
            )
        if 'pct_counts_mt' in qc_before.columns:
            gene_filter_cell_mask &= (
                pd.to_numeric(qc_before['pct_counts_mt'], errors='coerce').to_numpy()
                / 100.0 < mito_perc
            )
        elif 'mito_perc' in qc_before.columns:
            gene_filter_cell_mask &= (
                pd.to_numeric(qc_before['mito_perc'], errors='coerce').to_numpy()
                < mito_perc
            )
        try:
            gene_support_before_omicverse = _gene_detection_counts(
                adata.X, gene_filter_cell_mask,
            )
        except Exception as exc:
            gene_support_before_omicverse = None
            self.progress(-1, f'无法计算 OmicVerse gene filter provenance：{exc}')

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
            min_cells=gene_filter_min_cells,
            min_genes=gene_filter_min_genes,
            max_cells_ratio=gene_filter_max_cells_ratio,
            max_genes_ratio=gene_filter_max_genes_ratio,
            # Scrublet 在基础 QC 之后单独显式运行（scrublet_decoupled），
            # 因此这里不再让 ov.pp.qc 内部处理 doublet。
            doublets=(not doublets_disabled) and not scrublet_decoupled,
            doublets_method=doublets_method,
            batch_key=batch_key,
            # Keep the full predicted_doublet mask long enough to report
            # per-batch rates and to distinguish doublets from QC thresholds.
            filter_doublets=False,
            mt_genes=mt_genes,
            ribo_genes=ribo_genes,
            hb_genes=hb_genes,
        )
        if doublets_disabled:
            effective_doublets_method = 'none'
            doublet_verification_source = 'platform_disabled'
        elif scrublet_decoupled:
            # The preceding ov.pp.qc call deliberately had doublets=False;
            # the explicit call below is the authoritative Scrublet execution.
            effective_doublets_method = 'scrublet'
            doublet_verification_source = 'platform_explicit_scrublet'
        else:
            effective_doublets_method = _effective_omicverse_doublets_method(adata)
            doublet_verification_source = 'omicverse.status_args.qc.doublets_method'
            if effective_doublets_method is None:
                raise RuntimeError(
                    "无法验证 OmicVerse 实际使用的双细胞检测方法；"
                    "为避免不透明的 caller fallback，本次 QC 已停止。"
                )
        doublet_method_provenance = _doublet_method_provenance(
            'none' if doublets_disabled else doublets_method,
            effective_doublets_method,
            verification_source=doublet_verification_source,
            dependency=scdblfinder_dependency,
        )
        if doublet_method_provenance['fallback']:
            raise RuntimeError(
                "双细胞检测方法发生了未允许的回退："
                f"请求 {doublet_method_provenance['requested_doublets_method']}，"
                f"实际为 {doublet_method_provenance['effective_doublets_method']}。"
                "请修复对应依赖后重试，或在页面中主动选择实际要使用的方法。"
            )
        # OmicVerse 的 *_perc 别名为 0-1 fraction，但 pct_counts_* 应按
        # Scanpy 约定使用 0-100。在额外过滤、绘图和下游分析前统一。
        restore_scanpy_qc_percentages(adata)
        # Scrublet 双细胞检测：基础 QC 完成后再单独运行，参数全部显式传递。
        # batch_key 存在时 OmicVerse 内部对每个 batch 独立建模并保存
        # uns['scrublet']['batches'][batch] 下的 doublet_scores_sim /
        # threshold / parameters，供绘图与 manifest 使用。
        if scrublet_decoupled:
            self.progress(40, "Running Scrublet doublet detection (explicit parameters)...")
            scrublet_batch_key = (
                batch_key if (batch_key and batch_key in adata.obs.columns) else None
            )
            ov.pp.scrublet(
                adata,
                batch_key=scrublet_batch_key,
                expected_doublet_rate=expected_doublet_rate,
                sim_doublet_ratio=sim_doublet_ratio,
                threshold=scrublet_threshold_value,
                n_prin_comps=n_prin_comps,
                random_state=scrublet_random_state,
                use_gpu=scrublet_use_gpu,
                verbose=False,
            )
            self.progress(45, "Scrublet doublet detection finished.")
        n_genes_after_omicverse = int(adata.n_vars)
        removed_gene_mask = ~qc_input_var_names.isin(pd.Index(adata.var_names))
        n_removed_genes = int(removed_gene_mask.sum())
        gene_filter_provenance = {
            'operation': 'omicverse.pp.qc final gene filtering',
            'min_cells': gene_filter_min_cells,
            'min_genes': gene_filter_min_genes,
            'max_cells_ratio': gene_filter_max_cells_ratio,
            'max_genes_ratio': gene_filter_max_genes_ratio,
            'n_cells_used_for_support_estimate': int(gene_filter_cell_mask.sum()),
            'n_genes_removed': n_removed_genes,
            'status': 'review',
        }
        if gene_support_before_omicverse is not None:
            removed_support = gene_support_before_omicverse[removed_gene_mask]
            below_min_cells = removed_support < gene_filter_min_cells
            above_max_cells = (
                removed_support > gene_filter_max_cells_ratio * gene_filter_cell_mask.sum()
            )
            explained_mask = below_min_cells | above_max_cells
            gene_filter_provenance.update({
                'n_removed_below_min_cells': int(below_min_cells.sum()),
                'n_removed_above_max_cells_ratio': int(above_max_cells.sum()),
                'n_removed_unexplained_by_recorded_rules': int(
                    max(n_removed_genes - int(explained_mask.sum()), 0)
                ),
                'status': 'pass' if int(explained_mask.sum()) == n_removed_genes else 'review',
            })
        else:
            gene_filter_provenance['status'] = 'review'
        if doublets_disabled:
            # A previous QC result may carry this column.  Selecting "none"
            # must neither inherit nor filter stale doublet calls.
            doublet_labels = pd.Series(False, index=adata.obs_names)
        else:
            doublet_labels = adata.obs.get('predicted_doublet')
            if doublet_labels is None and 'predicted_doublets' in adata.obs.columns:
                doublet_labels = adata.obs['predicted_doublets']
            if doublet_labels is None:
                doublet_labels = pd.Series(False, index=adata.obs_names)
        doublet_labels = doublet_labels.fillna(False).astype(bool)
        if 'predicted_doublet' not in adata.obs.columns:
            adata.obs['predicted_doublet'] = doublet_labels
        doublet_index = adata.obs_names[doublet_labels.to_numpy()].tolist()
        scrublet_payload = (
            adata.uns.get('scrublet', {}) or {}
            if effective_doublets_method == 'scrublet' else {}
        )
        doublet_threshold = None
        if effective_doublets_method == 'scrublet':
            doublet_threshold = scrublet_payload.get('threshold')
            if doublet_threshold is None and isinstance(scrublet_payload.get('parameters'), dict):
                doublet_threshold = scrublet_payload['parameters'].get('threshold')
            run_params = {
                'expected_doublet_rate': expected_doublet_rate,
                'sim_doublet_ratio': sim_doublet_ratio,
                'n_prin_comps': n_prin_comps,
                'random_state': scrublet_random_state,
                'threshold_method': (
                    'manual' if scrublet_threshold_value is not None
                    else 'automatic_simulated_distribution'
                ),
                'n_neighbors': None,
            }
            doublet_summary = _scrublet_summary(
                adata, batch_key, doublet_threshold, scrublet_payload,
                run_params=run_params,
            )
        elif effective_doublets_method == 'scdblfinder':
            # OmicVerse already creates these backend-specific aliases.  Keep
            # a stable class field in the platform output even if a future
            # wrapper exposes only the generic predicted_doublet mask.
            if 'scdblfinder_doublet' not in adata.obs.columns:
                adata.obs['scdblfinder_doublet'] = doublet_labels.to_numpy()
            if 'scdblfinder_score' not in adata.obs.columns and 'doublet_score' in adata.obs.columns:
                adata.obs['scdblfinder_score'] = adata.obs['doublet_score'].to_numpy()
            adata.obs['scdblfinder_class'] = np.where(
                doublet_labels.to_numpy(), 'doublet', 'singlet',
            )
            doublet_summary = _scdblfinder_summary(adata, batch_key)
        else:
            doublet_summary = _caller_doublet_summary(
                adata, batch_key, effective_doublets_method,
            )
            if doublets_disabled:
                doublet_summary['available'] = False
                doublet_summary['disabled'] = True
        # Simulated score distributions are Scrublet-only evidence.  Do not
        # reuse stale arrays from an input AnnData after another caller ran.
        doublet_sim_scores_by_batch = {}
        if effective_doublets_method == 'scrublet':
            batches_payload = scrublet_payload.get('batches')
            if isinstance(batches_payload, dict):
                for batch_key_value, batch_metadata in batches_payload.items():
                    if not isinstance(batch_metadata, dict):
                        continue
                    sim_scores = batch_metadata.get('doublet_scores_sim')
                    if sim_scores is not None:
                        doublet_sim_scores_by_batch[str(batch_key_value)] = np.asarray(
                            sim_scores, dtype=float,
                        )
            if not doublet_sim_scores_by_batch:
                top_sim = scrublet_payload.get('doublet_scores_sim')
                if top_sim is not None:
                    doublet_sim_scores_by_batch['all'] = np.asarray(top_sim, dtype=float)
        remove_doublets = _as_bool(self.params.get('filter_doublets', True), True)
        # 双细胞分数直方图应展示“过滤前”的完整分布；剔除后再画会得到
        # 被截断的分布，误导阈值判断。
        doublet_scores_before_filter = (
            adata.obs['doublet_score'].copy()
            if 'doublet_score' in adata.obs.columns else None
        )
        # Preserve donor labels with the unfiltered score distribution.  Once
        # predicted doublets are removed, rebuilding this view from ``adata``
        # would hide the very observations used to determine each threshold.
        doublet_scores_by_batch = {}
        if doublet_scores_before_filter is not None and batch_key and batch_key in adata.obs.columns:
            score_frame = pd.DataFrame({
                'batch': adata.obs[batch_key].astype(str),
                'score': pd.to_numeric(doublet_scores_before_filter, errors='coerce'),
            }, index=adata.obs_names)
            for batch, values in score_frame.groupby('batch', observed=True)['score']:
                doublet_scores_by_batch[str(batch)] = values.to_numpy(dtype=float)
        if remove_doublets:
            adata = adata[~doublet_labels].copy()
        initial_obs_names = adata.obs_names.copy()
        n_after_doublet = adata.shape[0]
        doublet_summary['filter_doublets'] = bool(remove_doublets)
        doublet_summary['removed_doublets'] = int(len(doublet_index)) if remove_doublets else 0
        adata.uns['qc_doublets'] = doublet_summary
        adata.uns['qc_doublet_method'] = doublet_method_provenance
        if effective_doublets_method == 'scrublet':
            # Backward-compatible alias for historical Scrublet outputs.
            adata.uns['qc_scrublet'] = doublet_summary
        else:
            adata.uns.pop('qc_scrublet', None)
        stage_removed = {
            'max_genes': 0,
            'ribo': 0,
            'hb': 0,
            'batch_adaptive_qc': 0,
        }

        # ── 5. 额外过滤：基因数上限 ──────────────────────────────────────
        if ngenes_max > 0:
            mask = adata.obs['n_genes_by_counts'] <= ngenes_max
            n_before_upper = adata.shape[0]
            adata = adata[mask].copy()
            stage_removed['max_genes'] = n_before_upper - adata.shape[0]

        # ── 6. 额外过滤：核糖体比例上限 ─────────────────────────────────
        if ribo_perc_max > 0 and 'pct_counts_ribo' in adata.obs.columns:
            n_before_ribo = adata.shape[0]
            mask = adata.obs['pct_counts_ribo'] <= ribo_perc_max
            adata = adata[mask].copy()
            stage_removed['ribo'] = n_before_ribo - adata.shape[0]

        # ── 6b. 额外过滤：血红蛋白比例上限 ──────────────────────────────
        if hb_perc_max > 0 and 'pct_counts_hb' in adata.obs.columns:
            n_before_hb = adata.shape[0]
            mask = adata.obs['pct_counts_hb'] <= hb_perc_max
            adata = adata[mask].copy()
            stage_removed['hb'] = n_before_hb - adata.shape[0]

        # ── 6c. 批次自适应 QC（MAD 方法）──────────────────────────────
        if batch_adaptive and batch_key and batch_key in adata.obs.columns:
            self.progress(55, "Batch-adaptive QC filtering (MAD)...")
            n_before_adaptive = adata.shape[0]
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
                        # 只过滤低侧：计数/基因数的高侧离群是真实大细胞
                        # （高深度、大转录组），双侧 MAD 会把它们误删。
                        keep_mask[batch_mask] = vals >= lower
                adata = adata[keep_mask].copy()
            stage_removed['batch_adaptive_qc'] = n_before_adaptive - adata.shape[0]

        n_after = adata.shape[0]
        # 空数据守卫：全部细胞（或全部基因）被过滤时给出可操作错误，
        # 而不是让后续绘图以 NaN 坐标轴晦涩崩溃。
        if adata.n_obs == 0:
            raise ValueError(
                'QC 过滤后没有剩余细胞。请放宽 mito_perc / nUMIs / '
                'detected_genes / ribo_perc 等阈值后重试。'
            )
        if adata.n_vars == 0:
            raise ValueError(
                'QC 过滤后没有剩余基因。请检查基因过滤参数（min_cells 等）。'
            )

        # ── 7. 生成 Nature-style 图表 ────────────────────────────────────
        self.progress(60, "Generating QC plots...")
        plots_dir = self.ensure_plots_dir()
        result_files = []
        doublet_method_label = _doublet_method_display(effective_doublets_method)
        from modules.sc_figure_diagnostics import (
            doublet_by_batch_figure,
            qc_by_batch_figure,
            summarize_qc_by_batch,
        )

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
            nrows, ncols = _qc_overview_grid(len(panel_specs))
            fig_filter, axes = plt.subplots(
                nrows, ncols,
                figsize=(8.6, max(3.2, 2.65 * nrows)),
                squeeze=False,
            )
            flat_axes = axes.ravel()
            for axis, (label, before, after, number_format) in zip(flat_axes, panel_specs):
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
            for axis in flat_axes[len(panel_specs):]:
                axis.set_visible(False)
            fig_filter.suptitle('QC filtering overview', x=0.06, ha='left',
                                fontsize=13, fontweight='semibold', color=NATURE_TEXT)
            fig_filter.text(0.06, 0.01,
                            'Each metric uses its own scale so filtering effects remain readable.',
                            fontsize=7.5, color=NATURE_MUTED)
            _save(fig_filter, 'qc_filter_summary.png', 'qc', 'QC 过滤前后对比')

        # Caller-specific score distribution, using the unfiltered cells.
        if self.params.get('show_doublet_histogram', True) and doublet_scores_before_filter is not None:
            scores = _finite_values(doublet_scores_before_filter)
            if len(scores):
                fig_doublet, axis = plt.subplots(figsize=(7.6, 4.8))
                axis.hist(scores, bins=42, color=NATURE_PALETTE[0], alpha=0.78,
                          edgecolor='white', linewidth=0.35)
                median = float(np.median(scores))
                axis.axvline(median, color=NATURE_PALETTE[3], lw=1.1,
                             ls=(0, (3, 2)), label=f'Median {median:.3f}')
                if effective_doublets_method == 'scrublet' and doublet_threshold is not None:
                    try:
                        axis.axvline(float(doublet_threshold), color=NATURE_PALETTE[1],
                                     lw=1.0, ls=':',
                                     label=f'Threshold {float(doublet_threshold):.3f}')
                    except (TypeError, ValueError):
                        pass
                _style_axis(axis, f'{doublet_method_label} doublet score distribution',
                            xlabel='Doublet score', ylabel='Cells')
                axis.legend(loc='upper right', fontsize=8)
                _save(fig_doublet, 'qc_doublet_score_histogram.png', 'histogram',
                      f'{doublet_method_label} Doublet Score 分布')

        # QC Violin：复用 Bulk diagnostics 的最终物理尺寸与 QA 合同，
        # 但输入改为单细胞过滤前后分布，不把细胞误作生物学重复。
        violin_keys = [k for k in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo']
                       if k in qc_before.columns and k in adata.obs.columns]
        if violin_keys:
            metric_labels = {
                'n_genes_by_counts': 'log10(detected genes + 1)',
                'total_counts': 'log10(total counts + 1)',
                'pct_counts_mt': 'Mitochondrial reads (%)',
                'pct_counts_ribo': 'Ribosomal reads (%)',
            }
            violin_metrics = []
            for key in violin_keys:
                before_values = _finite_values(qc_before[key])
                after_values = _finite_values(adata.obs[key])
                if key in ('n_genes_by_counts', 'total_counts'):
                    before_values = np.log10(np.clip(before_values, 0, None) + 1.0)
                    after_values = np.log10(np.clip(after_values, 0, None) + 1.0)
                # Very large cell sets do not improve the violin density after
                # this deterministic display sample, but can exhaust memory.
                if len(before_values) > 50000:
                    before_values = before_values[np.linspace(
                        0, len(before_values) - 1, 50000, dtype=int
                    )]
                if len(after_values) > 50000:
                    after_values = after_values[np.linspace(
                        0, len(after_values) - 1, 50000, dtype=int
                    )]
                violin_metrics.append({
                    'label': metric_labels.get(key, key),
                    'values': [before_values, after_values],
                })
            director = NatureFigureDirector()
            violin_spec = director.spec_from_params(
                'qc_violin', self.params, width='double',
                title='Single-cell QC distributions before and after filtering',
                evidence_role='quality_control',
            ).with_updates(height_mm=105.0)
            fig_violin = director.render(violin_spec, {
                'kind': 'violin',
                'groups': ['Before QC', 'After QC'],
                'metrics': violin_metrics,
            })
            _save(fig_violin, 'qc_violin.png', 'violin', 'QC Violin Plots')

        # Batch-stratified QC diagnostics are intentionally separate from the
        # all-cell violin: they make library-complexity shifts visible without
        # treating batch as a biological grouping.
        if (
            batch_key and batch_key in adata.obs.columns
            and _as_bool(self.params.get('show_qc_batch_metrics', True), True)
        ):
            batch_series = adata.obs[batch_key].astype(str)
            batch_order = batch_series.value_counts().head(12).index.tolist()
            batch_metric_keys = [
                key for key in [
                    'n_genes_by_counts', 'total_counts',
                    'pct_counts_mt', 'pct_counts_ribo', 'novelty_score',
                    'doublet_score',
                ] if key in adata.obs.columns
            ]
            if len(batch_order) >= 2 and batch_metric_keys:
                ncols = 2
                nrows = (len(batch_metric_keys) + ncols - 1) // ncols
                fig_batch, axes_batch = plt.subplots(
                    nrows, ncols, figsize=(10.0, max(3.5, 3.0 * nrows)),
                    squeeze=False,
                )
                batch_labels = {
                    'n_genes_by_counts': 'Detected genes',
                    'total_counts': 'Total counts',
                    'pct_counts_mt': 'MT reads (%)',
                    'pct_counts_ribo': 'Ribosomal reads (%)',
                    'novelty_score': 'Novelty score',
                    'doublet_score': 'Doublet score',
                }
                for axis, metric in zip(axes_batch.ravel(), batch_metric_keys):
                    valid_batches = []
                    values = []
                    for batch in batch_order:
                        batch_values = _finite_values(
                            adata.obs.loc[batch_series == batch, metric]
                        )
                        if len(batch_values):
                            valid_batches.append(batch)
                            values.append(batch_values)
                    if len(valid_batches) < 2:
                        axis.set_visible(False)
                        continue
                    if metric in {'n_genes_by_counts', 'total_counts'}:
                        values = [np.log10(np.clip(item, 0, None) + 1.0) for item in values]
                    axis.violinplot(values, showmeans=True, showmedians=True, widths=0.82)
                    axis.set_xticks(
                        range(1, len(valid_batches) + 1), valid_batches,
                        rotation=35, ha='right',
                    )
                    _style_axis(axis, batch_labels.get(metric, metric), ylabel='Value')
                for axis in axes_batch.ravel()[len(batch_metric_keys):]:
                    axis.set_visible(False)
                fig_batch.suptitle(
                    f'QC metrics by batch ({batch_key})', x=0.06, ha='left',
                    fontsize=13, fontweight='semibold', color=NATURE_TEXT,
                )
                _save(fig_batch, 'qc_batch_metrics.png', 'violin', 'QC Metrics by Batch')

        # Donor-level retention and non-exclusive reason flags.  The reason
        # panel is explicitly scaled per 1,000 input cells, so overlapping
        # threshold violations are not misread as a partition of removals.
        qc_by_batch = {}
        if batch_key and batch_key in qc_before.columns:
            qc_by_batch = summarize_qc_by_batch(
                qc_before,
                adata.obs_names,
                batch_key,
                mito_perc=mito_perc,
                n_umis=nUMIs_min,
                n_genes_min=ngenes_min,
                n_genes_max=ngenes_max,
                ribo_perc_max=ribo_perc_max,
                hb_perc_max=hb_perc_max,
                doublet_index=doublet_index,
                doublet_summary=doublet_summary,
            )
            if qc_by_batch:
                fig_qc_batch = qc_by_batch_figure(qc_by_batch)
                if fig_qc_batch is not None:
                    result_files.extend(self.save_matplotlib_figure(
                        fig_qc_batch, plots_dir, 'qc_by_donor.png', 'qc',
                        'QC retention and non-exclusive flags by donor',
                        formats=('png', 'svg'), dpi=300,
                    ))
                if doublet_scores_by_batch:
                    fig_doublet_batch = doublet_by_batch_figure(
                        doublet_scores_by_batch,
                        doublet_summary,
                        method=effective_doublets_method,
                    )
                    if fig_doublet_batch is not None:
                        result_files.extend(self.save_matplotlib_figure(
                            fig_doublet_batch, plots_dir, 'qc_doublet_by_donor.png',
                            'histogram', f'{doublet_method_label} score and rate by donor',
                            formats=('png', 'svg'), dpi=300,
                        ))
                    # 关键 Scrublet QC 图：observed vs simulated + 自动阈值。
                    # 自动阈值是从 simulated doublet score 分布中找出来的，
                    # 所以只画 observed 分布无法判断阈值是否合理。
                    if effective_doublets_method == 'scrublet' and doublet_sim_scores_by_batch:
                        from modules.sc_figure_diagnostics import (
                            doublet_obs_vs_sim_figure,
                        )
                        fig_obs_sim = doublet_obs_vs_sim_figure(
                            doublet_scores_by_batch, doublet_sim_scores_by_batch,
                            doublet_summary,
                        )
                        if fig_obs_sim is not None:
                            result_files.extend(self.save_matplotlib_figure(
                                fig_obs_sim, plots_dir, 'qc_doublet_obs_vs_sim.png',
                                'histogram',
                                'Scrublet observed vs simulated doublet score distribution',
                                formats=('png', 'svg'), dpi=300,
                            ))

        # QC 散点图（Counts vs Genes，颜色 = MT%）
        if 'total_counts' in adata.obs.columns and 'n_genes_by_counts' in adata.obs.columns:
            # 三列必须来自同一批细胞：逐列过滤有限值后再按下标拼接会把
            # A 细胞的 counts 与 B 细胞的 genes/MT% 画在一起。
            scatter_columns = ['total_counts', 'n_genes_by_counts']
            if 'pct_counts_mt' in adata.obs.columns:
                scatter_columns.append('pct_counts_mt')
            scatter_frame = adata.obs[scatter_columns].apply(
                pd.to_numeric, errors='coerce'
            ).replace([np.inf, -np.inf], np.nan).dropna()
            x = scatter_frame['total_counts'].to_numpy(dtype=float)
            y = scatter_frame['n_genes_by_counts'].to_numpy(dtype=float)
            count = int(len(scatter_frame))
            color_vals = (
                scatter_frame['pct_counts_mt'].to_numpy(dtype=float)
                if 'pct_counts_mt' in scatter_frame.columns else None
            )
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
            # 同样成对过滤：novelty 更容易出现 NaN，独立过滤会错位。
            novelty_frame = adata.obs[['total_counts', 'novelty_score']].apply(
                pd.to_numeric, errors='coerce'
            ).replace([np.inf, -np.inf], np.nan).dropna()
            x = novelty_frame['total_counts'].to_numpy(dtype=float)
            y = novelty_frame['novelty_score'].to_numpy(dtype=float)
            count = int(len(novelty_frame))
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
            color_keys = list(dict.fromkeys(
                key for key in [batch_key, 'leiden', 'phase'] if key
            ))
            for color_key in color_keys:
                if color_key in adata.obs.columns:
                    fig_umap = self.build_publication_umap(
                        adata, color_key, title=f'UMAP by {color_key}'
                    )
                    _save(fig_umap, f'qc_umap_{color_key}.png', 'umap',
                          f'UMAP by {color_key}')

        # 汇总细胞周期比例
        phase_counts = {}
        if cc_available and 'phase' in adata.obs.columns:
            phase_counts = adata.obs['phase'].value_counts().to_dict()
        removal_counts = _qc_removal_counts(
            n_before, n_after_doublet, n_after, stage_removed,
        )
        removal_reasons = _qc_removal_reason_summary(
            qc_before,
            initial_obs_names,
            adata,
            mito_perc=mito_perc,
            n_umis=nUMIs_min,
            n_genes_min=ngenes_min,
            n_genes_max=ngenes_max,
            ribo_perc_max=ribo_perc_max,
            hb_perc_max=hb_perc_max,
            doublet_index=doublet_index,
            stage_removed=stage_removed,
        )
        gene_count_change = {
            'n_genes_before_qc': int(n_genes_before),
            'n_genes_after_omicverse_qc': int(n_genes_after_omicverse),
            'n_genes_after_all_qc': int(adata.n_vars),
            'genes_removed_by_omicverse_qc': int(n_genes_before - n_genes_after_omicverse),
            'genes_removed_by_additional_qc': int(n_genes_after_omicverse - adata.n_vars),
            'gene_filter_provenance': gene_filter_provenance,
            'explanation': (
                '基因数变化已归因到 OmicVerse QC 的 final gene filtering；'
                f"记录规则为 min_cells={gene_filter_min_cells}, "
                f"max_cells_ratio={gene_filter_max_cells_ratio}。"
                if n_genes_before != adata.n_vars else
                'QC 阶段未改变 adata.n_vars；页面中的 gene 数应解释为矩阵 var 数。'
            ),
        }
        adata.uns['qc_removal_reasons'] = removal_reasons
        adata.uns['qc_gene_count_change'] = gene_count_change
        adata.uns['qc_gene_filter_provenance'] = gene_filter_provenance
        adata.uns['qc_obs_field_aliases'] = QC_OBS_FIELD_ALIASES
        adata.uns['qc_expression_source'] = qc_expression_source
        canonicalize_qc_obs_fields(adata)

        # Preserve the representation supplied by the caller while retaining
        # all QC metrics and filters computed from counts.  AnnData indexing
        # handles both dense and sparse matrices without densifying the data.
        # ``Index.get_indexer`` 要求索引唯一，重复基因名会抛
        # InvalidIndexError；这里显式守卫，重复时保留 counts 表示而不是崩溃。
        if original_obs_names.is_unique and original_var_names.is_unique:
            obs_pos = original_obs_names.get_indexer(adata.obs_names)
            var_pos = original_var_names.get_indexer(adata.var_names)
            if (obs_pos >= 0).all() and (var_pos >= 0).all():
                adata.X = original_x[obs_pos][:, var_pos]
        else:
            self.progress(
                -1,
                "输入包含重复细胞或基因名，已跳过原始表达表示回填；"
                "输出保留基于 counts 的 QC 结果，建议先用 var_names_make_unique 清理输入。",
            )

        # Save only after all provenance metadata and canonical obs fields have
        # been attached; previously these uns entries were created after the
        # H5AD write and therefore were absent from the saved result.
        self.progress(85, "Saving output...")
        output_path = self.save_output(adata, 'qc')

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
                'n_genes_before_qc': n_genes_before,
                'n_genes_after_omicverse_qc': n_genes_after_omicverse,
                'genes_removed_by_omicverse_qc': n_genes_before - n_genes_after_omicverse,
                'genes_removed_by_additional_qc': n_genes_after_omicverse - adata.n_vars,
                'gene_count_change': gene_count_change,
                **removal_counts,
                **removal_reasons,
                'removal_reasons': removal_reasons,
                'doublets': doublet_summary,
                'doublet_method_provenance': doublet_method_provenance,
                'requested_doublets_method': doublet_method_provenance['requested_doublets_method'],
                'effective_doublets_method': doublet_method_provenance['effective_doublets_method'],
                'doublet_method_fallback': doublet_method_provenance['fallback'],
                'doublet_method_fallback_reason': doublet_method_provenance['fallback_reason'],
                **(
                    {'scrublet': doublet_summary}
                    if effective_doublets_method == 'scrublet' else {}
                ),
                'novelty_median': round(float(adata.obs['novelty_score'].median()), 4) if 'novelty_score' in adata.obs.columns else None,
                'cell_cycle_available': cc_available,
                'phase_counts': phase_counts,
                's_genes_found': len(s_in),
                'g2m_genes_found': len(g2m_in),
                'n_mt_genes_found': n_mt_genes_found,
                'n_ribo_genes_found': n_ribo_genes_found,
                'n_hb_genes_found': n_hb_genes_found,
                'gene_name_source': gene_name_source,
                'ribo_detection_status': 'pass' if n_ribo_genes_found else 'warning',
                'qc_metric_warnings': qc_metric_warnings,
                'qc_expression_source': qc_expression_source,
                # Retain the historic key, but make it the *effective* method
                # so it cannot contradict the detailed provenance above.
                'doublets_method': effective_doublets_method,
                'obs_field_aliases': QC_OBS_FIELD_ALIASES,
                'requested_batch_key': str(requested_batch or ''),
                'batch_key': batch_key,
                'qc_by_batch': qc_by_batch,
            }
        }
