"""固定图型模板共享的非视觉辅助函数。"""

import math
import re
import textwrap
from collections import Counter
from typing import Iterable, Sequence

import numpy as np


def as_1d(values, length=None, default=''):
    if values is None:
        if length is None:
            return np.asarray([], dtype=object)
        return np.asarray([default] * int(length), dtype=object)
    result = np.asarray(list(values), dtype=object).reshape(-1)
    if length is not None and len(result) != int(length):
        raise ValueError(f'字段长度 {len(result)} 与预期 {length} 不一致')
    return result


def finite_numeric(values, fill=np.nan):
    import pandas as pd

    result = pd.to_numeric(pd.Series(values), errors='coerce').to_numpy(dtype=float)
    if not np.isnan(fill):
        result = np.nan_to_num(result, nan=float(fill), posinf=float(fill), neginf=float(fill))
    return result


def label_tick_indices(count: int, maximum: int) -> np.ndarray:
    if maximum <= 0 or count <= 0:
        return np.asarray([], dtype=int)
    if count <= maximum:
        return np.arange(count, dtype=int)
    stride = max(1, int(math.ceil(count / maximum)))
    indices = np.arange(0, count, stride, dtype=int)
    if indices[-1] != count - 1:
        # Keep the last label without creating an adjacent final pair, which
        # is a common source of overlap at physical single-column size.
        indices[-1] = count - 1
    return indices


def wrap_term(term, width=34, max_chars=82):
    text = strip_term_id(term)
    if len(text) > max_chars:
        text = text[:max(1, max_chars - 1)].rstrip() + '…'
    return textwrap.fill(text, width=max(12, int(width)), break_long_words=False)


def strip_term_id(term):
    """Remove database identifiers from display labels, retaining them in data tables."""
    text = re.sub(r'\s+', ' ', str(term or '')).strip()
    text = re.sub(
        r'\s*\((?:GO|KEGG|REACTOME|WP|WIKIPATHWAYS)\s*:[^)]+\)\s*$',
        '', text, flags=re.IGNORECASE,
    )
    text = re.sub(
        r'\s+(?:GO|KEGG|REACTOME|WP|WIKIPATHWAYS)\s*:[A-Za-z0-9_.-]+\s*$',
        '', text, flags=re.IGNORECASE,
    )
    return text.strip()


def _term_query_key(value):
    """Normalise a user term/GO-ID query without changing stored labels."""
    return re.sub(r'[^a-z0-9]+', ' ', str(value or '').casefold()).strip()


def requested_term_rows(frame, queries, *, raw_column='_term', display_column='_display_term'):
    """Return rows matching target pathway names or database IDs.

    Matching is deterministic and deliberately permissive for user-facing
    controls: exact normalised matches are preferred, followed by a unique
    substring match.  The source index is retained so callers can combine
    requested rows with a ranked Top-N remainder.
    """
    if not queries:
        return frame.iloc[0:0].copy(), [], []
    raw_values = frame[raw_column].map(_term_query_key) if raw_column in frame else frame.index.map(_term_query_key)
    display_values = frame[display_column].map(_term_query_key) if display_column in frame else raw_values
    selected = []
    unmatched = []
    for query in queries:
        key = _term_query_key(query)
        if not key:
            continue
        exact = [index for index, value in zip(frame.index, raw_values) if value == key]
        exact += [index for index, value in zip(frame.index, display_values) if value == key and index not in exact]
        candidates = exact
        if not candidates:
            candidates = [index for index, value in zip(frame.index, raw_values) if key in value]
            candidates += [index for index, value in zip(frame.index, display_values)
                           if key in value and index not in candidates]
        if candidates:
            selected.append(candidates[0])
        else:
            unmatched.append(str(query).strip())
    selected = list(dict.fromkeys(selected))
    return frame.loc[selected].copy() if selected else frame.iloc[0:0].copy(), unmatched, selected


def gene_label_order(genes, selected_genes=(), *, strategy='all', membership=None, maximum=30):
    """Choose directly-labelled membership genes while preserving priorities.

    ``selected_genes`` always take priority.  ``membership`` may be a sequence
    of gene lists and is used by the ``shared`` strategy.  The node/edge list
    remains independent from this label cap, so a warning can explain any
    omitted labels without changing the enrichment result.
    """
    genes = list(dict.fromkeys(str(gene).strip() for gene in genes if str(gene).strip()))
    requested = [str(gene).strip() for gene in (selected_genes or ()) if str(gene).strip()]
    requested = [gene for gene in dict.fromkeys(requested) if gene in genes]
    strategy = str(strategy or 'all').strip().lower()
    if strategy == 'none':
        return [], len(genes)
    if strategy == 'selected':
        ordered = requested
    elif strategy == 'shared' and membership is not None:
        frequency = Counter(gene for values in membership for gene in values)
        shared = sorted(genes, key=lambda gene: (-frequency[gene], gene))
        # Explicit targets are meaningful even when a gene occurs in only one
        # selected pathway; the shared heuristic applies only to the remainder.
        ordered = requested + [gene for gene in shared if gene not in requested and frequency[gene] > 1]
    else:
        ordered = requested + [gene for gene in genes if gene not in requested]
    cap = max(0, int(maximum))
    return ordered[:cap], max(0, len(genes) - min(len(ordered), cap))


def gene_members(value):
    """Parse common enrichment membership delimiters into unique gene labels."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, (list, tuple, set, np.ndarray)):
        values = list(value)
    else:
        values = re.split(r'[;,/|]', str(value))
    return list(dict.fromkeys(
        str(item).strip() for item in values
        if str(item).strip() and str(item).strip().lower() not in {'nan', 'none'}
    ))


def compress_redundant_terms(frame, gene_column, *, threshold=0.85, maximum=None):
    """Keep representative terms while preserving the full source table.

    Rows are assumed to be ordered by scientific priority (usually FDR). A
    later term is hidden only when its hit-gene Jaccard similarity to an
    already-kept term reaches ``threshold``. The function returns the selected
    frame and the number of removed redundant rows.
    """
    if gene_column not in frame.columns:
        return frame.head(maximum) if maximum else frame, 0
    selected = []
    selected_sets = []
    removed = 0
    for index, row in frame.iterrows():
        genes = set(gene_members(row[gene_column]))
        redundant = False
        if genes:
            redundant = any(
                len(genes & other) / len(genes | other) >= float(threshold)
                for other in selected_sets if other
            )
        if redundant:
            removed += 1
            continue
        selected.append(index)
        selected_sets.append(genes)
        if maximum and len(selected) >= int(maximum):
            break
    if not selected:
        return frame.head(maximum) if maximum else frame, 0
    return frame.loc[selected].copy(), removed


def robust_symmetric_limit(matrix, default=2.5, ceiling=3.0):
    values = np.asarray(matrix, dtype=float)
    finite = np.abs(values[np.isfinite(values)])
    if finite.size == 0:
        return float(default)
    value = float(np.nanpercentile(finite, 98.0))
    if not np.isfinite(value) or value <= 0:
        value = float(default)
    return max(1.0, min(float(ceiling), value))


def adjust_labels(texts: Sequence, ax, *, arrow_color='#77808C'):
    """优先使用 adjustText；缺失时保持模板生成的确定性偏移。"""
    if not texts:
        return
    try:
        from adjustText import adjust_text

        adjust_text(
            list(texts), ax=ax,
            expand=(1.08, 1.16), force_text=(0.18, 0.35),
            force_static=(0.08, 0.16), max_move=(18, 28),
            ensure_inside_axes=True,
            arrowprops={'arrowstyle': '-', 'color': arrow_color, 'lw': 0.45},
        )
    except ImportError:
        return


def prune_overlapping_labels(texts: Sequence, ax):
    """Hide only the lowest-priority labels that still collide after layout.

    ``adjustText`` is intentionally used first, but at the final 89 mm/600 dpi
    export a long gene symbol can still touch its neighbour.  A deterministic
    post-layout guard keeps the figure publication-safe: labels are created in
    significance order by the templates, so later labels are the ones removed
    when a collision cannot be resolved without moving a label out of bounds.
    The full gene list remains available in the result table.
    """
    if not texts:
        return []
    try:
        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
    except Exception:
        return []
    visible = [text for text in texts if text.get_visible()]
    hidden = []
    while len(visible) > 1:
        boxes = [text.get_window_extent(renderer=renderer).expanded(1.02, 1.08)
                 for text in visible]
        conflict = None
        for left in range(len(visible)):
            for right in range(left + 1, len(visible)):
                if boxes[left].overlaps(boxes[right]):
                    conflict = right
                    break
            if conflict is not None:
                break
        if conflict is None:
            break
        removed = visible.pop(conflict)
        removed.set_visible(False)
        hidden.append(str(removed.get_text()))
    return hidden


def attach_contract(fig, spec, style, *, semantic_warnings: Iterable[str] = (),
                    semantic_errors: Iterable[str] = (), encodings=None,
                    validation_texts=None):
    profile = style.profile(spec)
    fig._nature_spec = spec.to_dict()
    fig._nature_spec_object = spec
    fig._nature_expected_size_mm = (profile.width_mm, profile.height_mm)
    fig._nature_semantic_warnings = [str(item) for item in semantic_warnings if str(item)]
    fig._nature_semantic_errors = [str(item) for item in semantic_errors if str(item)]
    fig._nature_encodings = dict(encodings or {})
    fig._nature_validation_texts = list(validation_texts or [])
    return fig
