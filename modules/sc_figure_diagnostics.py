"""Review-facing diagnostic figures for the single-cell pipeline.

The regular analysis modules own the scientific calculations.  This module
only turns already available AnnData/summary values into bounded, static
figures and small JSON-safe tables.  In particular, donor/sample is treated
as metadata for review; no cell-level value is presented as an inferential
replicate.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from modules.figure_style import (
    NATURE_AXIS,
    NATURE_GRID,
    NATURE_MUTED,
    NATURE_PALETTE,
    NATURE_TEXT,
    nature_continuous_cmap,
    stable_category_colors,
)


def resolve_donor_key(adata, requested=None, *, require_multiple=True):
    """Resolve a donor/sample/batch column without accepting cell-like IDs."""
    from modules.io_utils import resolve_obs_grouping

    requested = str(requested or '').strip()
    fallbacks = [
        'donor_id', 'donor', 'subject_id', 'subject', 'patient_id',
        'sample_id', 'sample', 'library_id', 'batch', 'technical_batch',
    ]
    if requested:
        key, info = resolve_obs_grouping(
            adata, requested, fallbacks=[], max_categories=50,
            max_numeric_categories=20, require_multiple=require_multiple,
        )
        return key, info
    return resolve_obs_grouping(
        adata, '', fallbacks=fallbacks, max_categories=50,
        max_numeric_categories=20, require_multiple=require_multiple,
    )


def _finite(values):
    values = pd.to_numeric(pd.Series(values), errors='coerce').to_numpy(dtype=float)
    return values[np.isfinite(values)]


def _natural_labels(values):
    return sorted(
        [str(value) for value in values],
        key=lambda value: (len(value), value),
    )


def _style_axis(ax, *, grid=True):
    ax.tick_params(labelsize=8, length=3, width=0.7, colors=NATURE_AXIS)
    ax.grid(axis='y' if grid else 'both', color=NATURE_GRID, linewidth=0.55, alpha=0.72)
    ax.set_axisbelow(True)
    for name, spine in ax.spines.items():
        spine.set_visible(name in ('left', 'bottom'))
        spine.set_color('#98A2B3')
        spine.set_linewidth(0.7)


def _donor_figure_width(labels, *, minimum):
    """Return enough horizontal space for two donor-labelled panels.

    Each of these figures has two axes, so a donor needs enough room in half
    the total figure width for both its rotated tick label and any value label
    above a bar.  Keeping that allocation proportional to donor count avoids
    reintroducing label collisions for larger projects.
    """
    return max(float(minimum), 7.0 + 1.0 * len(labels))


def _entropy(values):
    values = np.asarray(values, dtype=float)
    total = float(values.sum())
    if total <= 0 or len(values) <= 1:
        return 0.0
    probs = values[values > 0] / total
    return float(-(probs * np.log(probs)).sum() / np.log(len(values)))


def summarize_qc_by_batch(
    qc_before,
    final_obs_names,
    batch_key,
    *,
    mito_perc=0.2,
    n_umis=500,
    n_genes_min=250,
    n_genes_max=0,
    ribo_perc_max=0,
    hb_perc_max=0,
    doublet_index=None,
    doublet_summary=None,
    scrublet_summary=None,
):
    """Return donor-level QC retention and non-exclusive reason flags.

    Reason flags are calculated per 1,000 input cells and are intentionally
    not forced to sum to the removed count: one cell may violate more than one
    threshold.  This is the representation used by the companion figure.
    """
    if qc_before is None or batch_key not in qc_before.columns:
        return {}
    frame = qc_before.copy()
    labels = frame[batch_key].astype(str)
    final_index = pd.Index(final_obs_names)
    doublet_index = pd.Index([] if doublet_index is None else doublet_index)
    # ``scrublet_summary`` is accepted for backwards compatibility with
    # existing callers; new QC runs pass caller-neutral evidence here.
    doublet_by_batch = (doublet_summary or scrublet_summary or {}).get('by_batch', {})

    def numeric(column):
        if column not in frame.columns:
            return pd.Series(np.nan, index=frame.index)
        return pd.to_numeric(frame[column], errors='coerce')

    counts = numeric('total_counts')
    genes = numeric('n_genes_by_counts')
    mt = numeric('pct_counts_mt')
    ribo = numeric('pct_counts_ribo')
    hb = numeric('pct_counts_hb')
    flags = {
        'low_umi': counts < float(n_umis),
        'low_genes': genes < float(n_genes_min),
        'high_mt': mt > float(mito_perc) * 100.0,
        'high_genes': (genes > float(n_genes_max)) if float(n_genes_max) > 0 else pd.Series(False, index=frame.index),
        'high_ribo': (ribo > float(ribo_perc_max)) if float(ribo_perc_max) > 0 else pd.Series(False, index=frame.index),
        'high_hb': (hb > float(hb_perc_max)) if float(hb_perc_max) > 0 else pd.Series(False, index=frame.index),
        'doublet': frame.index.isin(doublet_index),
    }
    rows = {}
    for batch in _natural_labels(labels.unique()):
        mask = labels == batch
        indices = frame.index[mask]
        n_input = int(mask.sum())
        n_retained = int(len(final_index.intersection(indices)))
        payload = {
            'n_before': n_input,
            'n_after': n_retained,
            'n_removed': max(n_input - n_retained, 0),
            'pct_removed': round((n_input - n_retained) / max(n_input, 1) * 100.0, 3),
            'retained_fraction': round(n_retained / max(n_input, 1), 6),
            'flags_per_1000': {},
            'flagged_cells': {},
        }
        for name, flag in flags.items():
            number = int(np.asarray(flag.loc[indices] if hasattr(flag, 'loc') else flag[mask]).sum())
            payload['flagged_cells'][name] = number
            payload['flags_per_1000'][name] = round(number / max(n_input, 1) * 1000.0, 3)
        batch_doublet = doublet_by_batch.get(batch, doublet_by_batch.get(str(batch), {}))
        if isinstance(batch_doublet, dict):
            payload['doublet_rate'] = batch_doublet.get(
                'detected_doublet_rate', batch_doublet.get('doublet_rate'),
            )
            payload['doublet_threshold'] = batch_doublet.get('threshold')
            payload['n_cells_evaluated'] = batch_doublet.get('n_cells_evaluated')
            payload['detectable_doublet_fraction'] = batch_doublet.get(
                'detectable_doublet_fraction',
            )
            payload['estimated_overall_doublet_rate'] = batch_doublet.get(
                'estimated_overall_doublet_rate',
            )
        rows[batch] = payload
    return rows


def qc_by_batch_figure(summary, title='QC retention and non-exclusive flags by donor'):
    """Draw retention/removal and per-1,000 QC flags in two panels."""
    import matplotlib.pyplot as plt

    labels = list(summary)
    if not labels:
        return None
    x = np.arange(len(labels))
    fig, axes = plt.subplots(
        1, 2,
        figsize=(_donor_figure_width(labels, minimum=11.0), 5.5),
        dpi=150,
    )
    retained = np.asarray([summary[key]['n_after'] for key in labels], dtype=float)
    removed = np.asarray([summary[key]['n_removed'] for key in labels], dtype=float)
    axes[0].bar(x, retained, color=NATURE_PALETTE[0], label='Retained')
    axes[0].bar(x, removed, bottom=retained, color=NATURE_PALETTE[3], alpha=0.82, label='Removed')
    axes[0].set_xticks(x, labels, rotation=35, ha='right')
    axes[0].set_xlabel('Donor / batch')
    axes[0].set_ylabel('Input cells')
    axes[0].set_title('Retention by donor', loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    axes[0].legend(frameon=False, fontsize=8)
    for index, key in enumerate(labels):
        axes[0].text(index, retained[index] + removed[index] + max(1.0, (retained[index] + removed[index]) * 0.02),
                     f"{summary[key]['pct_removed']:.1f}% removed", ha='center', fontsize=7, color=NATURE_MUTED)
    total_cells = retained + removed
    if total_cells.size:
        axes[0].set_ylim(top=max(float(total_cells.max()) * 1.12, 1.0))

    reason_names = ['low_genes', 'low_umi', 'high_mt', 'doublet', 'high_genes', 'high_ribo', 'high_hb']
    reason_labels = ['Low genes', 'Low UMI', 'High MT', 'Doublet', 'High genes', 'High ribo', 'High HB']
    width = min(0.78 / len(reason_names), 0.14)
    for index, (reason, label) in enumerate(zip(reason_names, reason_labels)):
        values = [summary[key]['flags_per_1000'].get(reason, 0.0) for key in labels]
        axes[1].bar(x + (index - (len(reason_names) - 1) / 2) * width, values,
                    width=width, color=NATURE_PALETTE[index % len(NATURE_PALETTE)], label=label)
    axes[1].set_xticks(x, labels, rotation=35, ha='right')
    axes[1].set_xlabel('Donor / batch')
    axes[1].set_ylabel('Flagged cells / 1,000 input cells')
    axes[1].set_title('Removal reason flags', loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    axes[1].legend(frameon=False, fontsize=7, ncol=2)
    for ax in axes:
        _style_axis(ax)
    fig.suptitle(title, x=0.06, ha='left', fontsize=12, fontweight='semibold', color=NATURE_TEXT)
    # A figure-level footer has dedicated space below both panels.  Placing it
    # in axis coordinates made it overlap the rotated donor labels.
    fig.text(0.5, 0.025, 'Flags are non-exclusive; their sum may exceed removed cells.',
             ha='center', va='bottom', fontsize=7.5, color=NATURE_MUTED)
    layout_rect = (0, 0.10, 1, 0.93)
    fig.tight_layout(pad=1.2, rect=layout_rect)
    # The shared exporter applies a final tight layout.  Preserve this exact
    # reservation so the saved PNG/SVG keeps the footer clear of donor ticks.
    fig._native_layout_rect = layout_rect
    return fig


def doublet_by_batch_figure(scores_by_batch, scrublet_summary=None, title=None, *,
                            method='scrublet'):
    """Draw caller-specific donor score distributions and predicted rates.

    ``scrublet_summary`` remains the argument name for callers outside this
    module.  It accepts the generic QC doublet summary too; only Scrublet
    summaries are allowed to render simulated-score thresholds.
    """
    import matplotlib.pyplot as plt

    if not scores_by_batch:
        return None
    method = str(method or 'doublet caller').lower()
    method_label = {
        'scrublet': 'Scrublet',
        'scdblfinder': 'scDblFinder',
        'sccomposite': 'scComposite',
        'doubletfinder': 'DoubletFinder',
    }.get(method, method)
    title = title or f'{method_label} doublet score and predicted rate by donor'
    labels = list(scores_by_batch)
    fig, axes = plt.subplots(
        1, 2,
        figsize=(_donor_figure_width(labels, minimum=10.0), 5.5),
        dpi=150,
    )
    values = []
    valid_labels = []
    thresholds = (scrublet_summary or {}).get('threshold_by_batch', {})
    for label in labels:
        current = _finite(scores_by_batch[label])
        if current.size:
            valid_labels.append(str(label))
            values.append(current)
    if not values:
        return None
    axes[0].violinplot(values, positions=np.arange(1, len(values) + 1), showmedians=True,
                       showextrema=False, widths=0.78)
    axes[0].boxplot(values, positions=np.arange(1, len(values) + 1), widths=0.14,
                    showfliers=False, patch_artist=True,
                    boxprops={'facecolor': NATURE_PALETTE[0], 'alpha': 0.55, 'edgecolor': NATURE_PALETTE[0]},
                    medianprops={'color': NATURE_TEXT, 'linewidth': 1.0},
                    whiskerprops={'color': NATURE_PALETTE[0]}, capprops={'color': NATURE_PALETTE[0]})
    has_threshold = False
    for index, label in enumerate(valid_labels, 1):
        threshold = thresholds.get(label, thresholds.get(str(label)))
        try:
            if method == 'scrublet' and threshold is not None and np.isfinite(float(threshold)):
                # A short segment identifies the threshold for this donor
                # without layering numeric labels on neighbouring violins.
                axes[0].hlines(float(threshold), index - 0.34, index + 0.34,
                               color=NATURE_PALETTE[3], linestyle='--', linewidth=1.15)
                has_threshold = True
        except (TypeError, ValueError):
            pass
    axes[0].set_xticks(np.arange(1, len(valid_labels) + 1), valid_labels, rotation=35, ha='right')
    axes[0].set_xlabel('Donor / batch')
    axes[0].set_ylabel(f'{method_label} doublet score')
    axes[0].set_title('Score distribution', loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    if has_threshold:
        from matplotlib.lines import Line2D

        axes[0].legend(
            handles=[Line2D([0], [0], color=NATURE_PALETTE[3], linestyle='--', linewidth=1.15,
                            label='Donor-specific threshold')],
            frameon=False, fontsize=7,
        )

    rates = []
    rate_labels = []
    by_batch = (scrublet_summary or {}).get('by_batch', {})
    for label in valid_labels:
        item = by_batch.get(label, {})
        rate = item.get('detected_doublet_rate', item.get('doublet_rate'))
        if rate is not None:
            rate_labels.append(label)
            rates.append(float(rate) * 100.0)
    if rates:
        bars = axes[1].bar(np.arange(len(rates)), rates, color=NATURE_PALETTE[1], alpha=0.88)
        axes[1].set_xticks(np.arange(len(rates)), rate_labels, rotation=35, ha='right')
        axes[1].set_ylabel('Predicted doublets (%)')
        axes[1].set_xlabel('Donor / batch')
        axes[1].set_title('Predicted doublet rate', loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
        for bar, value in zip(bars, rates):
            axes[1].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(0.02, max(rates) * 0.02),
                         f'{value:.2f}%', ha='center', fontsize=7)
        # A very low detected rate is a review signal, not a caller failure.
        # Scrublet needs its simulated distribution as the decisive evidence;
        # other callers should be checked against their own diagnostics.
        low_rate_review = any(rate < 0.5 for rate in rates)
        if low_rate_review:
            footer_note = (
                'Review: detected doublet rate is unusually low (<0.5%); '
                + (
                    'inspect simulated score distribution and automatic threshold.'
                    if method == 'scrublet' else
                    'inspect the caller-specific scores and classifications.'
                )
            )
        elif len(rates) > 1 and min(rates) > 0 and max(rates) / min(rates) > 5:
            footer_note = 'Review: donor doublet rates differ by >5×.'
        else:
            footer_note = None
        axes[1].set_ylim(top=max(float(max(rates)) * 1.14, 1.0))
    else:
        footer_note = None
        axes[1].text(0.5, 0.5, 'No per-donor rate available', ha='center', va='center', color=NATURE_MUTED)
        axes[1].set_axis_off()
    for ax in axes:
        if ax.axison:
            _style_axis(ax)
    fig.suptitle(title, x=0.06, ha='left', fontsize=12, fontweight='semibold', color=NATURE_TEXT)
    if footer_note:
        # Reserve a figure-level footer so it cannot collide with donor ticks.
        fig.text(0.5, 0.025, footer_note, ha='center', va='bottom', fontsize=7.5,
                 color=NATURE_PALETTE[3])
    layout_rect = (0, 0.10 if footer_note else 0, 1, 0.93)
    fig.tight_layout(pad=1.2, rect=layout_rect)
    fig._native_layout_rect = layout_rect
    return fig



def doublet_obs_vs_sim_figure(scores_by_batch, sim_scores_by_batch, scrublet_summary=None,
                              title='Scrublet doublet score: observed vs simulated'):
    """Per-donor observed vs simulated doublet score distributions + threshold.

    This is the decisive Scrublet QC figure: with threshold=None Scrublet picks
    the automatic threshold from the simulated-doublet score histogram, so the
    threshold can only be judged together with the simulated scores.  Each
    panel shows observed cells, simulated doublets and the donor-specific
    automatic threshold, plus the detected rate and the fraction of simulated
    doublets above the threshold (detectable fraction).
    """
    import matplotlib.pyplot as plt

    labels = [label for label in scores_by_batch if _finite(scores_by_batch[label]).size]
    if not labels:
        return None
    thresholds = (scrublet_summary or {}).get('threshold_by_batch', {})
    by_batch = (scrublet_summary or {}).get('by_batch', {})
    n = len(labels)
    ncols = min(n, 5)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3.6 * ncols + 0.6, 3.4 * nrows + 1.1),
        dpi=150, squeeze=False,
    )
    flat_axes = axes.ravel()
    for index, label in enumerate(labels):
        axis = flat_axes[index]
        observed = _finite(scores_by_batch[label])
        simulated = _finite(sim_scores_by_batch.get(label, []))
        # Normalize both to density so the shapes are comparable even though
        # the simulated set is larger than the observed set.
        axis.hist(observed, bins=42, density=True, color=NATURE_PALETTE[0],
                  alpha=0.55, edgecolor='white', linewidth=0.3,
                  label=f'Observed (n={len(observed):,})')
        if simulated.size:
            axis.hist(simulated, bins=42, density=True, color=NATURE_PALETTE[1],
                      alpha=0.45, edgecolor='white', linewidth=0.3,
                      label=f'Simulated (n={len(simulated):,})')
        threshold = thresholds.get(label, thresholds.get(str(label)))
        try:
            threshold_value = float(threshold)
            if np.isfinite(threshold_value):
                axis.axvline(threshold_value, color=NATURE_PALETTE[3], lw=1.4,
                             ls='--', label=f'Threshold {threshold_value:.3f}')
        except (TypeError, ValueError):
            pass
        item = by_batch.get(label, {}) or {}
        detected = item.get('detected_doublet_rate', item.get('doublet_rate'))
        detectable = item.get('detectable_doublet_fraction')
        sub = f'Detected {detected * 100:.3f}%' if detected is not None else 'Detected n/a'
        if detectable is not None:
            sub += f'  |  Sim above thr {detectable * 100:.1f}%'
        axis.set_title(f'{label}\n{sub}', loc='left', fontsize=8,
                       fontweight='semibold', color=NATURE_TEXT)
        axis.set_xlabel('Doublet score', fontsize=7.5, color=NATURE_TEXT)
        axis.set_ylabel('Density', fontsize=7.5, color=NATURE_TEXT)
        axis.tick_params(labelsize=7, length=3, width=0.7, colors=NATURE_AXIS)
        axis.grid(axis='y', color=NATURE_GRID, linewidth=0.5, alpha=0.7)
        axis.set_axisbelow(True)
        for spine_name, spine in axis.spines.items():
            spine.set_visible(spine_name in ('left', 'bottom'))
            spine.set_color('#98A2B3')
            spine.set_linewidth(0.7)
    for axis in flat_axes[n:]:
        axis.set_visible(False)
    fig.suptitle(title, x=0.06, ha='left', fontsize=12, fontweight='semibold', color=NATURE_TEXT)
    handles, legend_labels = flat_axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, legend_labels, loc='upper center', bbox_to_anchor=(0.5, 0.985),
                   ncol=3, frameon=False, fontsize=8)
    low_rates = [
        (item.get('detected_doublet_rate', item.get('doublet_rate')))
        for item in by_batch.values() if isinstance(item, dict)
    ]
    footer_note = None
    if low_rates and any(rate is not None and float(rate) < 0.005 for rate in low_rates):
        footer_note = (
            'Review: detected doublet rate is unusually low (<0.5%); the automatic '
            'threshold is derived from the simulated distribution - compare it with '
            'the simulated scores above before judging the call.'
        )
    if footer_note:
        fig.text(0.5, 0.015, footer_note, ha='center', va='bottom', fontsize=7.5,
                 color=NATURE_PALETTE[3])
        layout_rect = (0, 0.08, 1, 0.90)
    else:
        layout_rect = (0, 0.02, 1, 0.94)
    fig.tight_layout(pad=1.2, rect=layout_rect)
    fig._native_layout_rect = layout_rect
    return fig


def library_size_by_batch_figure(raw_library_sizes, normalized_library_sizes, obs,
                                 batch_key, target_sum, title='Library size before and after normalization'):
    """Draw donor-aware library-size panels on raw and normalized scales."""
    import matplotlib.pyplot as plt

    if batch_key not in obs.columns:
        return None
    labels = obs[batch_key].astype(str)
    batches = _natural_labels(labels.unique())
    raw_groups = []
    normalized_groups = []
    used = []
    raw_library_sizes = np.asarray(raw_library_sizes, dtype=float)
    normalized_library_sizes = np.asarray(normalized_library_sizes, dtype=float)
    for batch in batches:
        mask = labels.to_numpy() == batch
        raw = raw_library_sizes[mask]
        normalized = normalized_library_sizes[mask]
        raw = raw[np.isfinite(raw)]
        normalized = normalized[np.isfinite(normalized)]
        if len(raw) or len(normalized):
            used.append(batch)
            raw_groups.append(np.log10(np.clip(raw, 0, None) + 1.0))
            normalized_groups.append(normalized)
    if not used:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(max(9.0, 0.8 * len(used) + 5.5), 5.0), dpi=150)
    positions = np.arange(1, len(used) + 1)
    axes[0].violinplot(raw_groups, positions=positions, showmedians=True, showextrema=False, widths=0.8)
    axes[0].set_xticks(positions, used, rotation=35, ha='right')
    axes[0].set_xlabel('Donor / batch')
    axes[0].set_ylabel('log10(raw total counts + 1)')
    axes[0].set_title('Before', loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    axes[1].violinplot(normalized_groups, positions=positions, showmedians=True, showextrema=False, widths=0.8)
    axes[1].axhline(float(target_sum), color=NATURE_PALETTE[3], linestyle='--', linewidth=1.0, label=f'Target {target_sum:g}')
    axes[1].set_xticks(positions, used, rotation=35, ha='right')
    axes[1].set_xlabel('Donor / batch')
    axes[1].set_ylabel('Normalized total counts')
    axes[1].set_title('After (linear scale)', loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    axes[1].legend(frameon=False, fontsize=8)
    for ax in axes:
        _style_axis(ax)
    fig.suptitle(title, x=0.06, ha='left', fontsize=12, fontweight='semibold', color=NATURE_TEXT)
    fig.tight_layout(pad=1.2, rect=(0, 0, 1, 0.93))
    return fig


def normalization_batch_medians_figure(raw_library_sizes, normalized_library_sizes, obs,
                                       batch_key, n_genes_key='n_genes_by_counts',
                                       title='Donor-level normalization diagnostics'):
    """Show library-size harmonization while retaining detected-gene differences."""
    import matplotlib.pyplot as plt

    if batch_key not in obs.columns:
        return None
    labels = obs[batch_key].astype(str)
    batches = _natural_labels(labels.unique())
    raw = np.asarray(raw_library_sizes, dtype=float)
    normalized = np.asarray(normalized_library_sizes, dtype=float)
    raw_medians = []
    normalized_medians = []
    genes_medians = []
    for batch in batches:
        mask = labels.to_numpy() == batch
        raw_medians.append(float(np.nanmedian(raw[mask])))
        normalized_medians.append(float(np.nanmedian(normalized[mask])))
        if n_genes_key in obs.columns:
            genes_medians.append(float(pd.to_numeric(obs.loc[mask, n_genes_key], errors='coerce').median()))
        else:
            genes_medians.append(np.nan)
    fig, axes = plt.subplots(1, 3, figsize=(max(11.0, 0.8 * len(batches) + 7.0), 4.8), dpi=150)
    x = np.arange(len(batches))
    axes[0].bar(x, np.log10(np.clip(raw_medians, 0, None) + 1), color=NATURE_PALETTE[6])
    axes[0].set_ylabel('log10(raw median counts + 1)')
    axes[0].set_title('Raw library size', loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    axes[1].bar(x, normalized_medians, color=NATURE_PALETTE[0])
    axes[1].set_ylabel('Median normalized counts')
    axes[1].set_title('Normalized library size', loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    axes[2].bar(x, genes_medians, color=NATURE_PALETTE[1])
    axes[2].set_ylabel('Median detected genes')
    axes[2].set_title('Complexity retained', loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    for ax in axes:
        ax.set_xticks(x, batches, rotation=35, ha='right')
        ax.set_xlabel('Donor / batch')
        _style_axis(ax)
    fig.suptitle(title, x=0.06, ha='left', fontsize=12, fontweight='semibold', color=NATURE_TEXT)
    fig.tight_layout(pad=1.2, rect=(0, 0, 1, 0.93))
    return fig


def hvg_batch_support_figure(support_counts, n_batches, n_selected, title='Batch support among selected HVGs'):
    import matplotlib.pyplot as plt

    support_counts = {int(key): int(value) for key, value in (support_counts or {}).items()}
    if not support_counts:
        return None
    x_values = list(range(1, max(int(n_batches), max(support_counts)) + 1))
    values = [support_counts.get(value, 0) for value in x_values]
    fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=150)
    bars = ax.bar(x_values, values, color=NATURE_PALETTE[0], alpha=0.88)
    ax.set_xticks(x_values, [f'{value}/{int(n_batches)}' for value in x_values])
    ax.set_xlabel('Number of supporting donors / batches')
    ax.set_ylabel('Selected HVGs')
    ax.set_title(title, loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(0.5, n_selected * 0.01),
                f'{int(value)} ({value / max(n_selected, 1) * 100:.1f}%)', ha='center', fontsize=8)
    ax.text(0.0, -0.22, f'Total selected HVGs = {int(sum(values)):,}', transform=ax.transAxes,
            fontsize=8, color=NATURE_MUTED)
    _style_axis(ax)
    fig.tight_layout(pad=1.1)
    return fig


def hvg_technical_composition_figure(composition, n_selected, title='Technical composition of selected HVGs'):
    """Draw selected-HVG technical flags without clipping small-category labels.

    The counts are flags rather than a mutually exclusive partition: a gene
    can, in principle, match more than one technical category.  Percentages
    therefore use all selected HVGs as their denominator.
    """
    import matplotlib.pyplot as plt

    if not composition:
        return None
    labels = list(composition)
    values = [int(composition[label]) for label in labels]
    fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=150)
    bars = ax.bar(labels, values, color=NATURE_PALETTE[:len(labels)], alpha=0.88)
    ax.set_ylabel('Selected HVGs')
    ax.set_xlabel('Gene category')
    ax.set_title(title, loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    # Technical categories are usually a small fraction of all HVGs.  Using
    # n_selected as an offset placed their labels far above short bars and
    # outside the axes.  Size the offset and top padding from plotted values.
    max_value = max(values)
    label_offset = max(0.3, max_value * 0.025)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + label_offset,
                f'{value} ({value / max(n_selected, 1) * 100:.1f}%)', ha='center', fontsize=8)
    ax.set_ylim(bottom=0, top=max(max_value + label_offset * 3.0, 1.0))
    _style_axis(ax)
    # Keep the interpretation note outside the axes, where it cannot collide
    # with gene-category labels or bar annotations.
    fig.text(0.5, 0.025,
             'Technical flags may overlap; percentages use all selected HVGs.',
             ha='center', va='bottom', fontsize=7.5, color=NATURE_MUTED)
    layout_rect = (0, 0.10, 1, 1)
    fig.tight_layout(pad=1.1, rect=layout_rect)
    fig._native_layout_rect = layout_rect
    return fig


def pca_variance_diagnostic_figure(variance_ratio, selected_n_comps, title='PCA variance explained'):
    import matplotlib.pyplot as plt

    variance = np.asarray(variance_ratio, dtype=float)
    variance = variance[np.isfinite(variance)]
    if variance.size == 0:
        return None
    variance = variance[:50]
    x = np.arange(1, len(variance) + 1)
    cumulative = np.cumsum(variance) * 100.0
    fig, ax = plt.subplots(figsize=(7.8, 5.0), dpi=150)
    ax.bar(x, variance * 100.0, color=NATURE_PALETTE[6], alpha=0.9, label='Per-PC variance')
    ax.set_xlabel('Principal component')
    ax.set_ylabel('Variance explained (%)')
    ax2 = ax.twinx()
    ax2.plot(x, cumulative, color=NATURE_PALETTE[0], marker='o', markersize=2.8, linewidth=1.6,
             label='Cumulative variance')
    ax2.set_ylabel('Cumulative variance (%)', color=NATURE_PALETTE[0])
    ax2.tick_params(axis='y', colors=NATURE_PALETTE[0], labelsize=8)
    if selected_n_comps is not None and int(selected_n_comps) <= len(variance):
        selected = int(selected_n_comps)
        ax.axvline(selected, color=NATURE_PALETTE[3], linestyle='--', linewidth=1.1,
                   label=f'Selected PCs = {selected}')
        ax.text(selected + 0.35, ax.get_ylim()[1] * 0.90, f'{selected} PCs', color=NATURE_PALETTE[3], fontsize=8)
    ax.set_title(title, loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    handles, labels = [], []
    for axis in (ax, ax2):
        h, l = axis.get_legend_handles_labels()
        handles.extend(h); labels.extend(l)
    ax.legend(handles, labels, frameon=False, fontsize=8, loc='upper right')
    _style_axis(ax)
    ax2.spines['right'].set_visible(False)
    fig.tight_layout(pad=1.1)
    return fig


def pca_loading_figure(loadings, gene_labels, variance_ratio=None, n_pcs=5, top_n=8,
                       title='Top PCA loading genes'):
    import matplotlib.pyplot as plt

    matrix = np.asarray(loadings, dtype=float)
    genes = np.asarray([str(value) for value in gene_labels])
    if matrix.ndim != 2 or matrix.shape[0] != len(genes) or matrix.shape[1] < 1:
        return None
    n_pcs = min(int(n_pcs), matrix.shape[1])
    top_n = max(1, int(top_n))
    fig, axes = plt.subplots(n_pcs, 2, figsize=(11.5, max(4.4, 1.65 * n_pcs)), dpi=150, squeeze=False)
    for pc in range(n_pcs):
        values = matrix[:, pc]
        finite = np.isfinite(values)
        if not finite.any():
            for ax in axes[pc]:
                ax.set_visible(False)
            continue
        positive_order = np.argsort(np.where(finite, values, -np.inf))[-top_n:][::-1]
        negative_order = np.argsort(np.where(finite, values, np.inf))[:top_n]
        for ax, order, sign, color in (
            (axes[pc, 0], positive_order, 'positive', NATURE_PALETTE[0]),
            (axes[pc, 1], negative_order, 'negative', NATURE_PALETTE[3]),
        ):
            selected = order[np.isfinite(values[order])]
            selected = selected[:top_n]
            selected = selected[::-1]
            ax.barh(np.arange(len(selected)), values[selected], color=color, alpha=0.86)
            ax.set_yticks(np.arange(len(selected)), genes[selected])
            ax.set_xlabel('Loading')
            suffix = ''
            if variance_ratio is not None and pc < len(variance_ratio):
                suffix = f' ({float(variance_ratio[pc]) * 100:.1f}%)'
            ax.set_title(f'PC{pc + 1}{suffix} — {sign}', loc='left', fontsize=9,
                         fontweight='semibold', color=NATURE_TEXT)
            _style_axis(ax)
            ax.grid(axis='x', color=NATURE_GRID, linewidth=0.45, alpha=0.7)
            ax.spines['left'].set_visible(False)
    fig.suptitle(title, x=0.05, ha='left', fontsize=12, fontweight='semibold', color=NATURE_TEXT)
    fig.tight_layout(pad=1.0, rect=(0, 0, 1, 0.96))
    return fig


def embedding_before_after_figure(before, after, labels, category_label='Donor',
                                  titles=('Before integration', 'After integration'),
                                  x_label='UMAP-1', y_label='UMAP-2'):
    """Draw two embeddings with one shared category-color mapping."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    before = np.asarray(before, dtype=float)[:, :2]
    after = np.asarray(after, dtype=float)[:, :2]
    labels = np.asarray([str(value) for value in labels])
    if before.shape != after.shape or len(labels) != len(before):
        raise ValueError('Before/after embedding dimensions do not match.')
    categories = list(dict.fromkeys(labels.tolist()))
    colors = stable_category_colors(categories)
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.8), dpi=150)
    for ax, coordinates, title in zip(axes, (before, after), titles):
        for category in categories:
            mask = labels == category
            ax.scatter(coordinates[mask, 0], coordinates[mask, 1], s=5, alpha=0.72,
                       color=colors[category], linewidths=0, rasterized=True)
        ax.set_title(title, loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        _style_axis(ax, grid=False)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    handles = [Line2D([], [], linestyle='None', marker='o', markersize=4,
                      markerfacecolor=colors[category], markeredgewidth=0, label=category)
               for category in categories[:20]]
    axes[1].legend(handles=handles, title=category_label, frameon=False, fontsize=7,
                   loc='center left', bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout(pad=1.1, rect=(0, 0, 0.82, 1))
    return fig


def stratified_mixing_figure(rows, group_label='cell type', title='Donor mixing within biological groups'):
    import matplotlib.pyplot as plt

    if not rows:
        return None
    frame = pd.DataFrame(rows)
    if frame.empty or 'group' not in frame.columns:
        return None
    frame = frame.sort_values(['n_cells', 'group'], ascending=[False, True]).head(20)
    x = np.arange(len(frame))
    fig, ax = plt.subplots(figsize=(max(8.0, 0.48 * len(frame) + 4.5), 5.0), dpi=150)
    width = 0.36
    if 'abs_asw_batch' in frame:
        ax.bar(x - width / 2, pd.to_numeric(frame['abs_asw_batch'], errors='coerce'), width,
               color=NATURE_PALETTE[3], label='|Batch ASW| (lower is better)')
    if 'same_donor_neighbor_fraction' in frame:
        ax.bar(x + width / 2, pd.to_numeric(frame['same_donor_neighbor_fraction'], errors='coerce'), width,
               color=NATURE_PALETTE[1], label='Same-donor neighbor fraction (lower generally better)')
    ax.set_xticks(x, frame['group'].astype(str), rotation=40, ha='right')
    ax.set_xlabel(group_label)
    ax.set_ylabel('Metric')
    ax.set_title(title, loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    ax.legend(frameon=False, fontsize=7, loc='upper right')
    ax.text(0.0, -0.23, 'Biology-associated structure can reduce donor mixing; interpret with markers and cell composition.',
            transform=ax.transAxes, fontsize=7.5, color=NATURE_MUTED)
    _style_axis(ax)
    fig.tight_layout(pad=1.1, rect=(0, 0.03, 1, 1))
    return fig


def composition_by_group(obs, group_key, label_key, *, max_labels=30):
    """Return count and within-group fraction tables for composition figures."""
    if group_key not in obs.columns or label_key not in obs.columns:
        return None, None
    work = pd.DataFrame({
        'group': obs[group_key].astype(str),
        'label': obs[label_key].astype(str),
    })
    all_counts = pd.crosstab(work['group'], work['label'])
    if all_counts.empty:
        return all_counts, all_counts
    max_labels = max(2, int(max_labels))
    label_totals = all_counts.sum(axis=0).sort_values(ascending=False)
    if len(label_totals) > max_labels:
        kept_labels = list(label_totals.head(max_labels - 1).index)
        other_label = 'Other'
        if other_label in kept_labels:
            other_label = 'Other (remaining labels)'
        counts = all_counts.loc[:, kept_labels].copy()
        # Keep all cells in the denominator and make the omitted labels
        # visible instead of silently re-normalising the displayed Top-N.
        counts[other_label] = all_counts.drop(columns=kept_labels).sum(axis=1)
    else:
        counts = all_counts.loc[:, list(label_totals.index)].copy()
    fractions = counts.div(all_counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    return counts, fractions


def composition_by_group_figure(fractions, group_label='Donor / sample', label_label='Cell type',
                                title='Cell-type composition by donor'):
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    from matplotlib.ticker import PercentFormatter

    if fractions is None or fractions.empty:
        return None
    groups = [str(value) for value in fractions.index]
    labels = [str(value) for value in fractions.columns]
    # Cell-type names can be long.  A top legend competes with the title and
    # gets clipped or overlaps it in a dense annotation, so reserve a stable
    # right-side column sized from the longest visible label instead.
    plot_width = max(7.2, 0.72 * len(groups) + 4.2)
    longest_legend_label = max([len(label_label)] + [len(label) for label in labels])
    legend_width = max(3.0, 0.075 * longest_legend_label + 0.9)
    figure_width = plot_width + legend_width
    layout_right = plot_width / figure_width
    fig, ax = plt.subplots(
        figsize=(figure_width, max(5.4, 0.20 * len(labels) + 3.8)), dpi=150,
    )
    bottom = np.zeros(len(groups), dtype=float)
    colors = stable_category_colors(labels)
    handles = []
    for index, label in enumerate(labels):
        values = fractions[label].to_numpy(dtype=float)
        bars = ax.bar(
            np.arange(len(groups)), values, bottom=bottom,
            color=colors[label], width=0.72, edgecolor='white', linewidth=0.3,
            label=label,
        )
        handles.append(bars[0])
        # There is no legible way to label every 1--2% segment in a 100%
        # stack.  Show values whenever the segment has room for text; the
        # remaining small segments retain their colour and legend entry.
        red, green, blue = to_rgb(colors[label])
        luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
        text_color = NATURE_TEXT if luminance > 0.62 else 'white'
        for bar_index, value in enumerate(values):
            if value >= 0.03:
                ax.text(
                    bar_index, bottom[bar_index] + value / 2, f'{value:.1%}',
                    ha='center', va='center', fontsize=6.8, fontweight='medium',
                    color=text_color, clip_on=True,
                )
        bottom += values
    ax.set_xticks(np.arange(len(groups)), groups, rotation=35, ha='right')
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=0))
    ax.set_xlabel(group_label)
    ax.set_ylabel('Cell proportion')
    ax.set_title(title, loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    # Reversed legend order mirrors the stack from top to bottom.  Keeping it
    # beside the plotting panel leaves title and labels in independent space.
    ax.legend(handles[::-1], labels[::-1], title=label_label, frameon=False,
              fontsize=7, title_fontsize=8, loc='upper left',
              bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    _style_axis(ax)
    ax.margins(x=0.06)
    layout_rect = (0, 0, layout_right, 1)
    fig.tight_layout(pad=1.1, rect=layout_rect)
    # The shared exporter runs tight_layout once more.  Preserve the legend
    # reservation so the saved PNG/SVG matches the reviewed layout.
    fig._native_layout_rect = layout_rect
    return fig


def composition_heatmap_figure(fractions, group_label='Donor / sample', label_label='Cell type',
                               title='Cell-type fraction by donor'):
    import matplotlib.pyplot as plt

    if fractions is None or fractions.empty:
        return None
    matrix = fractions.to_numpy(dtype=float).T
    fig, ax = plt.subplots(figsize=(max(7.5, 0.55 * len(fractions.index) + 4.8),
                                    max(4.8, 0.36 * len(fractions.columns) + 2.8)), dpi=150)
    image = ax.imshow(matrix, aspect='auto', cmap=nature_continuous_cmap(), vmin=0, vmax=max(0.01, float(np.nanmax(matrix))))
    ax.set_xticks(np.arange(len(fractions.index)), [str(value) for value in fractions.index], rotation=40, ha='right')
    ax.set_yticks(np.arange(len(fractions.columns)), [str(value) for value in fractions.columns])
    ax.set_xlabel(group_label)
    ax.set_ylabel(label_label)
    ax.set_title(title, loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.025)
    colorbar.set_label('Fraction', fontsize=8)
    colorbar.outline.set_visible(False)
    _style_axis(ax, grid=False)
    fig.tight_layout(pad=1.1)
    return fig


def cluster_composition_table(obs, cluster_key, donor_key):
    """Create cluster × donor review rows with entropy and dominance flag."""
    if cluster_key not in obs.columns or donor_key not in obs.columns:
        return pd.DataFrame()
    table = pd.crosstab(obs[cluster_key].astype(str), obs[donor_key].astype(str))
    rows = []
    for cluster in table.index:
        row = table.loc[cluster]
        total = int(row.sum())
        fractions = row / max(total, 1)
        major = str(fractions.idxmax()) if len(fractions) else ''
        major_fraction = float(fractions.max()) if len(fractions) else np.nan
        rows.append({
            'cluster': str(cluster),
            'n_cells': total,
            'major_donor': major,
            'major_donor_fraction': round(major_fraction, 4),
            'donor_entropy': round(_entropy(row.to_numpy(dtype=float)), 4),
            'donor_associated': bool(major_fraction >= 0.90),
        })
    return pd.DataFrame(rows)


def annotation_decision_heatmap(score_frame, cluster_order, decisions=None,
                                final_labels=None, display_label_map=None,
                                title='Candidate marker-programme scores by source cluster'):
    """Plot candidate scores with non-overlapping final annotation labels.

    Rows remain marker-programme scores, while the x-axis reports the final
    reviewed cell-type label for each source cluster.  Keeping those two
    concepts separate makes post-hoc display aliases auditable.
    """
    import matplotlib.pyplot as plt
    import textwrap

    if score_frame is None or score_frame.empty:
        return None
    score_frame = score_frame.copy()
    score_frame = score_frame.loc[:, [column for column in cluster_order if column in score_frame.columns]]
    if score_frame.empty:
        return None
    matrix = score_frame.to_numpy(dtype=float)
    display_label_map = display_label_map or {}
    row_labels = [
        display_label_map.get(str(value).replace('score_', ''), str(value).replace('score_', ''))
        for value in score_frame.index
    ]
    decision_lookup = decisions or {}
    final_labels = final_labels or {}
    tick_labels = []
    max_label_lines = 1
    for cluster in score_frame.columns:
        decision = decision_lookup.get(str(cluster), {})
        label = final_labels.get(
            str(cluster), decision.get('final_label', decision.get('top_candidate', '')),
        )
        label = str(label or 'Unassigned')
        wrapped_label = textwrap.fill(label, width=18, break_long_words=False)
        label_lines = wrapped_label.count('\n') + 1
        max_label_lines = max(max_label_lines, label_lines)
        margin = decision.get('score_margin')
        margin_text = ''
        if margin is not None:
            try:
                if np.isfinite(float(margin)):
                    margin_text = f'\nΔ={float(margin):.2f}'
            except (TypeError, ValueError):
                pass
        tick_labels.append(f'Cluster {cluster}\n{wrapped_label}{margin_text}')

    # Give every source cluster a real text column.  The old fixed-width plot
    # put long final labels as free text above the heatmap, where they all
    # overlapped the title and each other.
    figure_width = max(10.5, 1.22 * matrix.shape[1] + 5.0)
    figure_height = max(5.8, 0.31 * matrix.shape[0] + 3.2 + 0.32 * max_label_lines)
    layout_bottom = min(0.42, 0.07 + 0.045 * (max_label_lines + 2))
    fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=150)
    image = ax.imshow(matrix, aspect='auto', cmap=nature_continuous_cmap(), vmin=0, vmax=max(1.0, float(np.nanmax(matrix))))
    ax.set_xticks(np.arange(matrix.shape[1]), tick_labels, rotation=0, ha='center')
    ax.tick_params(axis='x', labelsize=6.7, pad=4)
    ax.set_yticks(np.arange(matrix.shape[0]), row_labels)
    ax.set_xlabel('Source cluster / final displayed cell type', labelpad=8)
    ax.set_ylabel('Candidate marker programme')
    ax.set_title(title, loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.025)
    colorbar.set_label('Mean marker-programme score', fontsize=8)
    colorbar.outline.set_visible(False)
    _style_axis(ax, grid=False)
    layout_rect = (0, layout_bottom, 1, 1)
    fig.tight_layout(pad=1.1, rect=layout_rect)
    fig._native_layout_rect = layout_rect
    return fig


def state_score_heatmap_figure(obs, score_columns, group_key, title='Cell-state score by donor / group'):
    if group_key not in obs.columns:
        return None
    groups = obs[group_key].astype(str)
    selected = [column for column in score_columns if column in obs.columns]
    if not selected:
        return None
    matrix = pd.DataFrame({column: pd.to_numeric(obs[column], errors='coerce') for column in selected})
    matrix[group_key] = groups.to_numpy()
    matrix = matrix.groupby(group_key, observed=True).mean().T
    if matrix.empty:
        return None
    from modules.native_figures import heatmap_figure
    return heatmap_figure(
        matrix.to_numpy(dtype=float),
        x_labels=matrix.columns.tolist(),
        y_labels=[column.replace('state_score_', '') for column in matrix.index],
        title=title, x_label=group_key, y_label='Cell state',
        colorbar_label='Mean score', vmin=0, vmax=1,
    )


def state_score_by_group_figure(obs, score_columns, group_key, title='Cell-state scores by donor / group'):
    import matplotlib.pyplot as plt

    if group_key not in obs.columns:
        return None
    selected = [column for column in score_columns if column in obs.columns]
    groups = _natural_labels(obs[group_key].astype(str).unique())
    if not selected or not groups:
        return None
    ncols = 2
    nrows = int(np.ceil(len(selected) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(max(8.5, 0.75 * len(groups) + 5.0), 3.0 * nrows), dpi=150, squeeze=False)
    group_values = obs[group_key].astype(str)
    for axis, column in zip(axes.ravel(), selected):
        values = []
        used_groups = []
        for group in groups:
            current = _finite(obs.loc[group_values == group, column])
            if len(current):
                used_groups.append(group)
                values.append(current)
        if not values:
            axis.set_visible(False)
            continue
        axis.violinplot(values, positions=np.arange(1, len(values) + 1), showmedians=True, showextrema=False, widths=0.78)
        axis.set_xticks(np.arange(1, len(used_groups) + 1), used_groups, rotation=35, ha='right')
        axis.set_xlabel(group_key)
        axis.set_ylabel('Score')
        axis.set_ylim(0, 1)
        axis.set_title(column.replace('state_score_', ''), loc='left', fontsize=9, fontweight='semibold', color=NATURE_TEXT)
        _style_axis(axis)
    for axis in axes.ravel()[len(selected):]:
        axis.set_visible(False)
    fig.suptitle(title, x=0.06, ha='left', fontsize=12, fontweight='semibold', color=NATURE_TEXT)
    fig.tight_layout(pad=1.1, rect=(0, 0, 1, 0.94))
    return fig


def condition_composition_figure(obs, sample_key, condition_key, label_key,
                                 title='Cell-type composition by condition'):
    """Plot one 100% stacked cell-type composition bar for each condition.

    The bars are calculated from all cells assigned to each condition, so they
    answer the descriptive question "what fraction of the cells in Control or
    IBD have each annotation?".  ``sample_key`` is retained for call-site
    compatibility but intentionally is not used as a statistical replicate.
    Sample-level composition and inference remain separate downstream steps.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    from matplotlib.ticker import PercentFormatter

    required = {condition_key, label_key}
    if not required.issubset(obs.columns):
        return None
    # Filter genuine missing values before converting to strings.  Converting
    # first would render missing metadata as a literal "nan" bar/legend item.
    frame = obs.loc[obs[condition_key].notna() & obs[label_key].notna(), [
        condition_key, label_key,
    ]].copy()
    frame.columns = ['condition', 'label']
    frame['condition'] = frame['condition'].astype(str).str.strip()
    frame['label'] = frame['label'].astype(str).str.strip()
    frame = frame[frame['condition'].ne('') & frame['label'].ne('')]
    if frame.empty:
        return None

    counts = pd.crosstab(frame['condition'], frame['label'])
    if counts.empty:
        return None
    conditions = _natural_labels(counts.index)
    counts = counts.reindex(conditions)
    # Keep the legend usable in large annotations without silently
    # re-normalising away omitted cells.  The explicit Other segment preserves
    # each condition's 100% denominator.
    max_labels = 20
    label_totals = counts.sum(axis=0).sort_values(ascending=False)
    if len(label_totals) > max_labels:
        kept_labels = label_totals.head(max_labels - 1).index.tolist()
        other = counts.drop(columns=kept_labels).sum(axis=1)
        counts = counts.loc[:, kept_labels].copy()
        counts['Other'] = other
    else:
        counts = counts.loc[:, label_totals.index.tolist()]
    fractions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    labels = [str(value) for value in fractions.columns]
    colors = stable_category_colors(labels)

    # Reserve a dedicated, content-aware legend column.  A fixed 27% of the
    # canvas made the actual bars very narrow when cell-type names were long,
    # and placed the legend too close to the stack in browser previews.
    plot_width = max(7.2, 0.85 * len(conditions) + 4.2)
    longest_legend_label = max([len(label_key)] + [len(label) for label in labels])
    legend_width = max(3.0, 0.075 * longest_legend_label + 0.9)
    figure_width = plot_width + legend_width
    layout_right = plot_width / figure_width
    fig, ax = plt.subplots(
        figsize=(figure_width, max(5.4, 0.25 * len(labels) + 3.8)), dpi=150,
    )
    x = np.arange(len(conditions))
    bottom = np.zeros(len(conditions), dtype=float)
    handles = []
    for label in labels:
        values = fractions[label].to_numpy(dtype=float)
        bars = ax.bar(
            x, values, bottom=bottom, width=0.70,
            color=colors[label], edgecolor='white', linewidth=0.45, label=label,
        )
        handles.append(bars[0])
        red, green, blue = to_rgb(colors[label])
        luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
        text_color = NATURE_TEXT if luminance > 0.62 else 'white'
        for bar_index, value in enumerate(values):
            # Keep labels readable in dense compositions; exact values remain
            # available in the annotation output and composition table.
            if value >= 0.025:
                ax.text(
                    x[bar_index], bottom[bar_index] + value / 2,
                    f'{value:.1%}', ha='center', va='center', fontsize=6.8,
                    color=text_color, fontweight='medium', clip_on=True,
                )
        bottom += values
    ax.set_xticks(x, conditions, rotation=0)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=0))
    ax.set_xlabel(condition_key)
    ax.set_ylabel('Cell proportion')
    ax.set_title(title, loc='left', fontsize=10, fontweight='semibold', color=NATURE_TEXT)
    # Reverse the legend so its visual order agrees with the stack from top
    # to bottom, matching conventional 100% composition charts.
    ax.legend(handles[::-1], labels[::-1], title=label_key, frameon=False,
              fontsize=7, title_fontsize=8, loc='upper left',
              bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)
    _style_axis(ax)
    ax.margins(x=0.12)
    layout_rect = (0, 0, layout_right, 1)
    fig.tight_layout(pad=1.1, rect=layout_rect)
    # ``BaseAnalysis.save_matplotlib_figure`` performs a shared final layout
    # pass.  Preserve the right-side legend reservation for the PNG and SVG.
    fig._native_layout_rect = layout_rect
    return fig
