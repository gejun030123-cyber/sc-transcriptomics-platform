"""Shared rules for carrying the final PCA dimensionality downstream."""

DEFAULT_FINAL_N_PCS = 25


def resolve_analysis_n_pcs(adata, requested=None, *, representation_key='X_pca',
                           default=DEFAULT_FINAL_N_PCS):
    """Resolve one explicit PC count for correction, neighbours, and clustering.

    ``dimred`` records the selected PCA dimensionality in ``uns['dimred_pca']``.
    Downstream modules honour that upstream choice and also cap it by their own
    requested value and the actual representation width.  This prevents a
    stale/default request such as 30 or 50 PCs from silently reintroducing
    dimensions that were not part of the final PCA analysis.
    """
    import numpy as np

    requested_n_pcs = int(default if requested is None else requested)
    if requested_n_pcs <= 0:
        raise ValueError('n_pcs 必须是正整数')
    if representation_key not in adata.obsm:
        raise ValueError(f'表示 {representation_key} 不存在')

    representation = np.asarray(adata.obsm[representation_key])
    if representation.ndim != 2 or representation.shape[1] < 2:
        raise ValueError(f'表示 {representation_key} 至少需要 2 个维度')
    available_n_pcs = int(representation.shape[1])

    diagnostics = getattr(adata, 'uns', {}).get('dimred_pca', {}) or {}
    upstream_value = diagnostics.get('selected_n_pcs')
    try:
        upstream_n_pcs = int(upstream_value) if upstream_value is not None else None
    except (TypeError, ValueError):
        upstream_n_pcs = None
    if upstream_n_pcs is not None and upstream_n_pcs <= 0:
        upstream_n_pcs = None

    limits = [requested_n_pcs, available_n_pcs]
    if upstream_n_pcs is not None:
        limits.append(upstream_n_pcs)
    used_n_pcs = min(limits)
    if used_n_pcs < 2:
        raise ValueError('可用 PC 数不足 2，无法构建下游表示')

    if upstream_n_pcs is not None:
        source = 'dimred_selected_n_pcs_capped_by_request_and_representation'
    elif requested_n_pcs > available_n_pcs:
        source = 'representation_width_capped'
    else:
        source = 'requested_n_pcs'
    return int(used_n_pcs), {
        'requested_n_pcs': int(requested_n_pcs),
        'upstream_selected_n_pcs': upstream_n_pcs,
        'available_n_pcs': available_n_pcs,
        'used_n_pcs': int(used_n_pcs),
        'source': source,
    }


def resolution_filename_token(resolution):
    """Make a stable, collision-free filename token such as ``0p6`` or ``1p0``."""
    text = f'{float(resolution):.12f}'.rstrip('0').rstrip('.')
    if '.' not in text:
        text += '.0'
    return text.replace('-', 'm').replace('.', 'p')
