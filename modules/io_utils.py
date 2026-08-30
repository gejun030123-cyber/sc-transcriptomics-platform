import os
import numpy as np
import pandas as pd
from config import Config


QC_OBS_COLUMNS = frozenset({
    'total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'size_factor',
    'total_counts_mt', 'total_counts_ribo', 'total_counts_hb',
    'pct_counts_ribo', 'pct_counts_hb',
    'mito_perc', 'ribo_perc', 'hb_perc', '_auto_group',
})

QC_PERCENTAGE_COLUMNS = (
    ('pct_counts_mt', 'mito_perc'),
    ('pct_counts_ribo', 'ribo_perc'),
    ('pct_counts_hb', 'hb_perc'),
)


def restore_scanpy_qc_percentages(adata):
    """Restore ``pct_counts_*`` columns to Scanpy's 0-100 percentage scale.

    OmicVerse stores its ``*_perc`` aliases as fractions and mirrors those
    fractions into Scanpy-named ``pct_counts_*`` columns. Prefer the explicit
    fraction aliases when present so downstream thresholds and plots keep the
    standard Scanpy percentage semantics. Values outside the fractional 0-1
    range are treated as percentages for compatibility with imported data.
    """
    for pct_column, fraction_column in QC_PERCENTAGE_COLUMNS:
        if fraction_column not in adata.obs.columns:
            continue
        values = pd.to_numeric(adata.obs[fraction_column], errors='coerce')
        finite = values[np.isfinite(values)]
        is_fraction = finite.empty or (
            float(finite.min()) >= 0.0 and float(finite.max()) <= 1.0 + 1e-8
        )
        adata.obs[pct_column] = values * 100.0 if is_fraction else values
    return adata

# These fields describe cells, libraries, or technical processing rather than
# the biological grouping that should normally drive DEG/proportion analyses.
# They remain available in the UI for manual selection and for batch/sample
# specific parameters; they are simply not the first generic ``groupby`` pick.
TECHNICAL_OBS_COLUMNS = frozenset({
    'barcode', 'barcodes', 'cell_barcode', 'cell_id', 'cellid', 'cell_index',
    'sample', 'sample_id', 'library_id', 'orig_ident', 'orig.ident',
    'source_type', 'source_description',
    'batch', 'technical_batch', 'sequencing_batch', 'library_batch',
})

_GROUPING_PRIORITIES = {
    'default': (
        'analysis_group', 'condition', 'treatment', 'treatment_group',
        'group', 'sample_group', 'celltype', 'cell_type', 'annotation',
        'cell_type_annotation', 'leiden', 'cluster', 'seurat_clusters',
    ),
    'bulk': (
        'analysis_group', 'condition', 'treatment', 'treatment_group',
        'group', 'sample_group',
    ),
    'proportion': (
        'celltype', 'cell_type', 'annotation', 'cell_type_annotation',
        'leiden', 'cluster', 'condition', 'treatment', 'group',
    ),
    'cluster': (
        'leiden', 'cluster', 'seurat_clusters', 'celltype', 'cell_type',
        'annotation',
    ),
    'condition': (
        'condition', 'treatment', 'treatment_group', 'group', 'sample_group',
        'analysis_group',
    ),
}

_BATCH_PRIORITY = (
    'batch', 'technical_batch', 'sequencing_batch', 'library_batch',
    'sample', 'sample_id', 'library_id', 'orig.ident', 'orig_ident',
)

_SAMPLE_PRIORITY = (
    'sample_id', 'sample', 'library_id', 'orig.ident', 'orig_ident',
)


def _normalise_obs_column_name(column):
    """Normalise only for matching aliases; preserve the original column name."""
    return str(column).strip().lower().replace('-', '_').replace(' ', '_')


def is_technical_obs_column(column):
    """Return whether an obs field is normally technical/identifier metadata."""
    normalised = _normalise_obs_column_name(column)
    return normalised in {
        _normalise_obs_column_name(name) for name in TECHNICAL_OBS_COLUMNS
    }


def _obs_priority(column, priority_names):
    """Give preferred aliases a stable rank, including leiden_0.8-like names."""
    normalised = _normalise_obs_column_name(column)
    aliases = [_normalise_obs_column_name(name) for name in priority_names]
    try:
        return aliases.index(normalised)
    except ValueError:
        for index, alias in enumerate(aliases):
            if alias in {'leiden', 'cluster'} and normalised.startswith(alias + '_'):
                return index
        return len(aliases) + 100


def obs_grouping_info(adata, column, *, max_categories=50,
                      max_numeric_categories=20, require_multiple=False):
    """Describe whether an ``obs`` column is safe for grouping operations.

    Analysis modules frequently receive a free-form ``cluster_key``/``batch_key``
    or ``groupby`` parameter.  Presence alone is not enough: QC metrics such as
    ``log1p_n_genes_by_counts`` are numeric and often nearly one value per cell,
    which turns a crosstab or per-group loop into hundreds of pseudo-groups.
    Keep the rule in one place so modules use the same semantics.
    """
    result = {
        'column': str(column or ''),
        'valid': False,
        'n_unique': 0,
        'reason': '',
    }
    if not column:
        result['reason'] = '未指定分组列'
        return result
    if column not in getattr(adata, 'obs', pd.DataFrame()).columns:
        result['reason'] = f"列 '{column}' 不存在于 adata.obs"
        return result

    values = adata.obs[column]
    n_unique = int(values.nunique(dropna=False))
    result['n_unique'] = n_unique
    if n_unique < (2 if require_multiple else 1):
        result['reason'] = '分组数量不足'
        return result

    if pd.api.types.is_numeric_dtype(values) and n_unique > max_numeric_categories:
        result['reason'] = (
            f"列 '{column}' 是连续/高基数数值列（{n_unique} 个取值），"
            '不能作为分类分组列'
        )
        return result

    # Reject sample/cell identifiers and other almost-one-value-per-row fields;
    # for small datasets the relative limit is stricter than the absolute cap.
    n_obs = max(int(getattr(adata, 'n_obs', len(values))), 1)
    cardinality_limit = int(max_categories)
    relative_limit = max(20, int(0.2 * n_obs))
    if n_unique > cardinality_limit or (n_unique > relative_limit and n_unique > 20):
        safe_limit = min(cardinality_limit, relative_limit) if relative_limit < cardinality_limit else cardinality_limit
        result['reason'] = (
            f"列 '{column}' 的分组数 {n_unique} 超过安全上限 {safe_limit}"
        )
        return result

    result['valid'] = True
    return result


def batch_cluster_overlap(adata, cluster_key, batch_key, *, dominance_threshold=0.90):
    """Summarise whether clusters are dominated by a single technical batch.

    This is a descriptive warning, not a batch-correction score: a genuinely
    biological population can be present in one batch.  It gives downstream
    review code a bounded, explicit signal instead of requiring visual
    inspection of a composition bar plot alone.
    """
    result = {
        'valid': False,
        'cluster_key': str(cluster_key or ''),
        'batch_key': str(batch_key or ''),
        'n_clusters': 0,
        'n_batches': 0,
        'dominance_threshold': float(dominance_threshold),
        'batch_dominated_cluster_count': 0,
        'batch_dominated_cell_fraction': 0.0,
        'max_batch_fraction': None,
        'weighted_max_batch_fraction': None,
        'warnings': [],
    }
    if not cluster_key or not batch_key:
        result['reason'] = '未指定 cluster 或 batch 列'
        return result
    columns = getattr(adata, 'obs', pd.DataFrame()).columns
    if cluster_key not in columns or batch_key not in columns:
        result['reason'] = 'cluster 或 batch 列不存在'
        return result

    clusters = adata.obs[cluster_key].astype(str)
    batches = adata.obs[batch_key].astype(str)
    n_clusters = int(clusters.nunique(dropna=False))
    n_batches = int(batches.nunique(dropna=False))
    result.update({'n_clusters': n_clusters, 'n_batches': n_batches})
    if n_clusters < 2 or n_batches < 2:
        result['reason'] = '至少需要两个 cluster 和两个 batch 才能评估重合'
        return result

    table = pd.crosstab(clusters, batches)
    sizes = table.sum(axis=1).astype(float)
    fractions = table.div(sizes, axis=0)
    max_fractions = fractions.max(axis=1)
    dominated = max_fractions >= float(dominance_threshold)
    dominated_cells = float(sizes.loc[dominated].sum())
    total_cells = max(float(sizes.sum()), 1.0)
    weights = sizes / total_cells
    weighted_max = float((max_fractions * weights).sum())
    result.update({
        'valid': True,
        'batch_dominated_cluster_count': int(dominated.sum()),
        'batch_dominated_cell_fraction': round(dominated_cells / total_cells, 4),
        'max_batch_fraction': round(float(max_fractions.max()), 4),
        'weighted_max_batch_fraction': round(weighted_max, 4),
    })
    if int(dominated.sum()):
        result['warnings'].append(
            f"{int(dominated.sum())}/{n_clusters} 个 cluster 的单一 batch 占比 ≥ "
            f"{float(dominance_threshold):.0%}；需排查 batch-driven cluster。"
        )
    return result


def rank_obs_grouping_candidates(adata, module_name='', *, purpose='groupby',
                                 include_technical=False):
    """Return safe obs grouping columns in a parameter-aware order.

    The old UI used the physical obs column order, which commonly puts
    ``barcode``/``batch``/``sample`` before ``condition`` or ``leiden``.  This
    helper keeps those fields selectable but separates their automatic role:
    biological group fields for ``groupby``/``group_column``, technical fields
    for ``batch_key``, and cluster fields for ``cluster_key``.
    """
    module = str(module_name or '').strip().lower()
    purpose = str(purpose or 'groupby').strip().lower()
    columns = list(getattr(adata, 'obs', pd.DataFrame()).columns)
    if purpose == 'batch_key':
        priority = _BATCH_PRIORITY
        allow_technical = True
    elif purpose == 'sample_key':
        priority = _SAMPLE_PRIORITY
        allow_technical = True
    elif purpose == 'cluster_key':
        priority = (
            _GROUPING_PRIORITIES['proportion']
            if module == 'cell_communication'
            else _GROUPING_PRIORITIES['cluster']
        )
        allow_technical = False
    elif purpose == 'condition_key' or purpose == 'group_column':
        priority = _GROUPING_PRIORITIES['condition']
        allow_technical = bool(include_technical)
    elif module == 'deg' and any(
        _normalise_obs_column_name(column) in {'celltype', 'cell_type', 'annotation'}
        for column in columns
    ):
        priority = (
            'celltype', 'cell_type', 'annotation', 'cell_type_annotation',
            'leiden', 'cluster', 'condition', 'treatment', 'group',
        )
        allow_technical = bool(include_technical)
    elif module == 'bulk_deg' or module.startswith('bulk_'):
        priority = _GROUPING_PRIORITIES['bulk']
        allow_technical = bool(include_technical)
    elif module == 'proportion':
        priority = _GROUPING_PRIORITIES['proportion']
        allow_technical = bool(include_technical)
    else:
        priority = _GROUPING_PRIORITIES['default']
        allow_technical = bool(include_technical)

    candidates = []
    for index, column in enumerate(columns):
        if column in QC_OBS_COLUMNS:
            continue
        if not allow_technical and is_technical_obs_column(column):
            continue
        info = obs_grouping_info(
            adata, column, max_categories=50,
            max_numeric_categories=20, require_multiple=True,
        )
        if not info['valid']:
            continue
        values = adata.obs[column].astype(str)
        counts = values.value_counts()
        # A grouping with singleton levels is not a useful automatic design.
        if len(counts) < 2 or int(counts.min()) < 2:
            continue
        candidates.append((
            _obs_priority(column, priority), index, column,
        ))

    candidates.sort(key=lambda item: (item[0], item[1]))
    return [item[2] for item in candidates]


def resolve_obs_grouping(adata, requested, fallbacks=(), *, max_categories=50,
                         max_numeric_categories=20, require_multiple=False):
    """Return ``(selected_column, info)`` using the first valid candidate.

    ``selected_column`` is ``None`` when neither the requested column nor a
    fallback is safe.  The returned ``info`` also contains the reason for a
    rejected requested column, allowing callers to surface an actionable warning.
    """
    candidates = []
    for candidate in [requested, *fallbacks]:
        if candidate and candidate not in candidates:
            candidates.append(candidate)

    requested_info = obs_grouping_info(
        adata, requested, max_categories=max_categories,
        max_numeric_categories=max_numeric_categories,
        require_multiple=require_multiple,
    )
    for candidate in candidates:
        info = obs_grouping_info(
            adata, candidate, max_categories=max_categories,
            max_numeric_categories=max_numeric_categories,
            require_multiple=require_multiple,
        )
        if info['valid']:
            info = dict(info)
            info['requested_column'] = str(requested or '')
            info['requested_valid'] = bool(requested_info['valid'])
            info['requested_reason'] = requested_info['reason']
            return candidate, info

    requested_info = dict(requested_info)
    requested_info['requested_column'] = str(requested or '')
    requested_info['requested_valid'] = bool(requested_info['valid'])
    requested_info['requested_reason'] = requested_info['reason']
    return None, requested_info


def infer_expression_measurement(adata, input_path=''):
    """Classify a bulk expression matrix without silently changing its scale."""
    normalization = dict(getattr(adata, 'uns', {}).get('normalization', {}) or {})
    if normalization.get('is_log_transformed'):
        return 'log_transformed'
    matrix = adata.X
    values = matrix.data if hasattr(matrix, 'data') else np.asarray(matrix).ravel()
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    hint = str(input_path).lower()
    if any(token in hint for token in ('fpkm', 'tpm', 'rpkm')):
        return 'continuous_expression'
    if values.size == 0:
        return 'raw_counts'
    integer_fraction = float(np.mean(np.isclose(values, np.round(values))))
    return 'raw_counts' if np.min(values) >= 0 and integer_fraction >= 0.995 else 'continuous_expression'


def infer_sample_group_candidates(sample_names):
    """Infer reusable grouping candidates from names ending in a replicate number.

    Examples:
      Ctr_1 -> Ctr
      Ctr_B_1 -> combined=Ctr_B, factor_1=Ctr, factor_2=B
    """
    import re
    from collections import Counter

    parsed = []
    for raw_name in sample_names:
        name = str(raw_name)
        clean = re.sub(
            r'_(count|FPKM|TPM|fpkm|tpm|Counts|normalized)$', '', name)
        tokens = [token for token in re.split(r'[-_]', clean) if token]
        if len(tokens) < 2:
            return []

        # Most matrices use ``Ctrl_1``/``Ctrl-1``, but exported bulk tables
        # also commonly encode the biological replicate without a separator,
        # for example ``EC_CT_BSA1``. Split only a terminal numeric replicate
        # suffix, leaving the experimental factor itself intact (BSA here).
        replicate_match = re.fullmatch(r'(.+?)(?:rep)?(\d+)', tokens[-1], re.I)
        if re.fullmatch(r'(?:rep)?\d+', tokens[-1], re.I):
            group_tokens = tokens[:-1]
        elif replicate_match and replicate_match.group(1):
            group_tokens = tokens[:-1] + [replicate_match.group(1)]
        else:
            return []
        parsed.append((name, group_tokens))

    if not parsed:
        return []

    raw_candidates = []
    max_factors = max(len(tokens) for _, tokens in parsed)

    # The full pre-replicate name preserves treatment x stratum designs.
    raw_candidates.append({
        'key': 'combined',
        'label': '联合分组（推荐）' if max_factors > 1 else '自动分组（推荐）',
        'mapping': {name: '_'.join(tokens) for name, tokens in parsed},
    })

    for factor_idx in range(max_factors):
        if not all(len(tokens) > factor_idx for _, tokens in parsed):
            continue
        raw_candidates.append({
            'key': f'factor_{factor_idx + 1}',
            'label': f'第 {factor_idx + 1} 因素',
            'mapping': {name: tokens[factor_idx] for name, tokens in parsed},
        })

    candidates = []
    seen_mappings = set()
    for candidate in raw_candidates:
        mapping = candidate['mapping']
        signature = tuple(mapping[name] for name, _ in parsed)
        if signature in seen_mappings:
            continue
        seen_mappings.add(signature)
        counts = Counter(mapping.values())
        if len(counts) < 2 or any(count < 2 for count in counts.values()):
            continue
        candidate['values'] = sorted(counts)
        candidate['group_sizes'] = dict(sorted(counts.items()))
        candidates.append(candidate)

    return candidates


def read_expression_matrix(file_path):
    """读取表达矩阵文件，自动检测格式，只保留数值列"""
    import scanpy as sc

    if file_path.endswith('.h5ad'):
        import scanpy as _sc
        adata = _sc.read_h5ad(file_path)
        import re as _re
        # 过滤掉非样本行（如 Exonic.gene.sizes, Start, End 等）
        non_sample_mask = adata.obs.index.to_series().apply(
            lambda x: bool(_re.match(r'^(Exonic|Start|End|Chr|Strand)', str(x), _re.I))
        )
        if non_sample_mask.any() and not non_sample_mask.all():
            adata = adata[~non_sample_mask].copy()
        # 如果同时有 _count 和 _FPKM 列，只保留 _count
        obs_names = [str(n) for n in adata.obs.index]
        count_idx = [i for i, n in enumerate(obs_names) if '_count' in n]
        fpkm_idx = [i for i, n in enumerate(obs_names) if '_FPKM' in n or '_fpkm' in n]
        if count_idx and fpkm_idx:
            adata = adata[count_idx].copy()
        # 清理样本名：去掉 _count/_FPKM/_TPM 后缀
        clean_names = [_re.sub(r'_(count|FPKM|TPM|fpkm|tpm)$', '', str(n)) for n in adata.obs.index]
        if len(set(clean_names)) == len(clean_names):
            adata.obs.index = pd.Index(clean_names)
        return adata

    df = None

    # 尝试 Excel
    if file_path.endswith(('.xlsx', '.xls')):
        try:
            df = pd.read_excel(file_path, index_col=0)
        except Exception:
            pass

    # 尝试 TSV
    if df is None:
        try:
            df = pd.read_csv(file_path, sep='\t', index_col=0, nrows=5)
            if df.shape[1] > 0:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
            else:
                df = None
        except Exception:
            df = None

    # 尝试自动检测
    if df is None:
        try:
            df = pd.read_csv(file_path, sep=None, engine='python', index_col=0, nrows=5)
            if df.shape[1] > 0:
                df = pd.read_csv(file_path, sep=None, engine='python', index_col=0)
            else:
                df = None
        except Exception:
            df = None

    # 标准 CSV
    if df is None:
        df = pd.read_csv(file_path, index_col=0)

    # 只保留数值列（过滤掉基因注释等非表达数据列）
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if len(numeric_cols) == 0:
        # 尝试强制转换所有列为数值
        df = df.apply(pd.to_numeric, errors='coerce')
        df = df.dropna(axis=1, how='all')
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

    if len(numeric_cols) == 0:
        raise ValueError("无法从输入文件中读取数值列。请确认文件格式正确且包含数值型表达数据。")

    # 智能识别样本列：排除非样本的数值列（如 Exonic.gene.sizes, Start, End 等）
    non_sample_keywords = {'exonic', 'start', 'end', 'size', 'length', 'gc', 'position', 'coord'}
    real_sample_cols = []
    for col in numeric_cols:
        col_lower = col.lower().replace('_', '').replace('.', '').replace('-', '')
        if any(kw in col_lower for kw in non_sample_keywords):
            continue
        # 只保留看起来像样本名的列（包含 _count 或 _FPKM 或 _TPM，或者以数字结尾）
        if any(suffix in col for suffix in ('_count', '_FPKM', '_TPM', '_fpkm', '_tpm')):
            real_sample_cols.append(col)
        elif col[-1].isdigit():
            real_sample_cols.append(col)
        else:
            real_sample_cols.append(col)  # 保留无法判断的列

    if len(real_sample_cols) > 0 and len(real_sample_cols) < len(numeric_cols):
        numeric_cols = real_sample_cols

    # 如果样本列同时包含 count 和 FPKM，优先使用 count
    count_cols = [c for c in numeric_cols if '_count' in c.lower()]
    fpkm_cols = [c for c in numeric_cols if '_fpkm' in c.lower() or '_FPKM' in c]
    if count_cols and fpkm_cols:
        # 优先使用 count 数据
        numeric_cols = count_cols
        print(f"[io_utils] 检测到 count 和 FPKM 数据，使用 count 数据 ({len(count_cols)} 列)")

    # 保留基因名注释列（用于 DEG 结果中显示基因名而非 ID）
    gene_name_col = None
    # 优先级：symbol 类 > name 类 > 其他
    symbol_candidates = ['GeneSymbol', 'gene_symbol', 'gene_name', 'symbol', 'Symbol',
                         'Gene Symbol', 'gene', 'Gene', 'GENE_NAME', 'gene_name_x']
    name_candidates = ['Description', 'description', 'gene_description']
    non_numeric_cols = [c for c in df.columns if c not in numeric_cols]
    for candidate in symbol_candidates:
        if candidate in non_numeric_cols:
            gene_name_col = candidate
            break
    if gene_name_col is None:
        for col in non_numeric_cols:
            col_lower = col.lower()
            if 'symbol' in col_lower:
                gene_name_col = col
                break
    if gene_name_col is None:
        for candidate in name_candidates:
            if candidate in non_numeric_cols:
                gene_name_col = candidate
                break
    if gene_name_col is None:
        for col in non_numeric_cols:
            col_lower = col.lower()
            if 'name' in col_lower:
                gene_name_col = col
                break
    gene_names = df[gene_name_col].tolist() if gene_name_col else None

    if len(numeric_cols) < df.shape[1]:
        df = df[numeric_cols]

    # 替换 NaN 为 0
    nan_count = int(df.isna().sum().sum())
    if nan_count > 0:
        print(f"[io_utils] WARNING: 替换 {nan_count} 个 NaN 值为 0")
    df = df.fillna(0)

    adata = sc.AnnData(X=df.values.T, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))

    # 清理样本名：去掉 _count/_FPKM/_TPM 后缀
    import re as _re
    clean_names = [_re.sub(r'_(count|FPKM|TPM|fpkm|tpm)$', '', str(n)) for n in adata.obs.index]
    if len(set(clean_names)) == len(clean_names):
        adata.obs.index = pd.Index(clean_names)

    # 如果找到基因名列，保存到 var 中
    if gene_names is not None:
        adata.var['gene_name'] = gene_names

    # 将 var_names 从 Ensembl ID 映射为基因名
    adata = remap_var_names(adata)

    return adata


def infer_sc_data_format(input_path):
    """Infer supported single-cell input format from path."""
    if not input_path:
        return 'unknown'

    lower = str(input_path).lower()
    if os.path.isdir(input_path):
        names = set(os.listdir(input_path))
        if (
            ('matrix.mtx' in names or 'matrix.mtx.gz' in names)
            and ('barcodes.tsv' in names or 'barcodes.tsv.gz' in names)
            and (
                'features.tsv' in names or 'features.tsv.gz' in names
                or 'genes.tsv' in names or 'genes.tsv.gz' in names
            )
        ):
            return '10x_mtx'
        if lower.endswith('.zarr'):
            return 'zarr'
        return 'directory'

    if lower.endswith('.h5ad'):
        return 'h5ad'
    if lower.endswith(('.h5', '.hdf5')):
        return '10x_h5'
    if lower.endswith('.loom'):
        return 'loom'
    if lower.endswith('.zarr'):
        return 'zarr'
    if lower.endswith(('.csv', '.txt', '.tsv', '.xlsx', '.xls')):
        return 'expression_matrix'
    if lower.endswith(('.mtx', '.mtx.gz')):
        return '10x_mtx'
    return 'unknown'


def _is_raw_count_matrix(matrix):
    """Return whether a matrix is finite, non-negative and integer-valued."""
    values = np.asarray(matrix.data if hasattr(matrix, 'data') else matrix, dtype=float)
    if values.size == 0:
        return True
    return bool(
        np.isfinite(values).all()
        and (values >= 0).all()
        and np.allclose(values, np.rint(values), rtol=0.0, atol=1e-6)
    )


def _ensure_counts_layer(adata, input_format=None):
    """Preserve counts without inventing them from an arbitrary h5ad ``X``.

    Matrix/10x imports conventionally contain raw counts in ``X`` and may be
    promoted after an explicit integer/non-negative check.  A processed h5ad,
    loom or zarr object has no such contract: its ``X`` may be log-normalized,
    residuals or an embedding matrix.  Leave those objects unmodified and mark
    the missing layer so count-dependent modules can fail with an actionable
    message instead of producing plausible-looking but invalid results.
    """
    if 'counts' in adata.layers:
        adata.uns['counts_layer_missing'] = False
        return adata
    raw_native_formats = {'10x_mtx', '10x_h5', 'expression_matrix'}
    if input_format in raw_native_formats and _is_raw_count_matrix(adata.X):
        adata.layers['counts'] = adata.X.copy()
        adata.uns['counts_layer_inferred_from_x'] = True
        adata.uns['counts_layer_missing'] = False
    else:
        adata.uns['counts_layer_missing'] = True
    return adata


def _standardize_imported_adata(adata, input_format=None, species=None, genome=None):
    """Apply lightweight AnnData normalization needed by downstream modules."""
    if hasattr(adata, 'var_names_make_unique'):
        adata.var_names_make_unique()
    if hasattr(adata, 'obs_names_make_unique'):
        adata.obs_names_make_unique()

    adata = remap_var_names(adata)
    adata = _ensure_counts_layer(adata, input_format=input_format)

    if input_format:
        adata.uns['input_format'] = input_format
    if species:
        adata.uns['species'] = species
    if genome:
        adata.uns['genome'] = genome
    return adata


def read_10x_mtx_compat(mtx_dir, var_names='gene_symbols'):
    """Read a 10x Matrix Market directory with mixed gzip/plain files safely.

    10x downloads are sometimes assembled from files with different compression
    states.  ``scanpy.read_10x_mtx`` expects one consistent filename set, so a
    temporary all-gzip view is created when needed.
    """
    import gzip
    import shutil
    import tempfile
    import scanpy as sc

    if not os.path.isdir(mtx_dir):
        raise FileNotFoundError(f"10x 矩阵目录不存在: {mtx_dir}")

    def _pick(*names):
        return next((name for name in names if os.path.isfile(os.path.join(mtx_dir, name))), None)

    matrix = _pick('matrix.mtx.gz', 'matrix.mtx')
    barcodes = _pick('barcodes.tsv.gz', 'barcodes.tsv')
    features = _pick('features.tsv.gz', 'features.tsv', 'genes.tsv.gz', 'genes.tsv')
    if not matrix or not barcodes or not features:
        raise FileNotFoundError(
            '10x 数据不完整：需要 matrix.mtx(.gz)、barcodes.tsv(.gz) 和 '
            'features.tsv(.gz) 或 genes.tsv(.gz)'
        )

    # Use the source directly only when it already has a complete compressed
    # set; otherwise stage all three selected files to avoid omitting the
    # already-compressed member of a mixed set.
    if (matrix.endswith('.gz') and barcodes.endswith('.gz')
            and features == 'features.tsv.gz'):
        return sc.read_10x_mtx(mtx_dir, var_names=var_names, cache=False)

    # Pass the directory explicitly as a defense in depth: an external caller
    # may import this module after a library has already cached ``/tmp``.
    temp_dir = tempfile.mkdtemp(prefix='10x_gzip_', dir=Config.runtime_tmp_dir())
    try:
        for source_name, target_name in (
            (matrix, 'matrix.mtx.gz'),
            (barcodes, 'barcodes.tsv.gz'),
            # New Scanpy releases unconditionally look for features.tsv.gz.
            # A v2 genes.tsv has the same first two fields and is therefore
            # staged under that filename for a version-independent read.
            (features, 'features.tsv.gz'),
        ):
            source = os.path.join(mtx_dir, source_name)
            target = os.path.join(temp_dir, target_name)
            if source_name.startswith('genes.'):
                # Scanpy 1.11 requires the third feature-type column, whereas
                # 10x v2 genes.tsv has only gene ID and symbol.  Add the
                # conventional type while staging it as a v3-style file.
                opener = gzip.open if source_name.endswith('.gz') else open
                with opener(source, 'rt', encoding='utf-8') as src, gzip.open(
                    target, 'wt', encoding='utf-8'
                ) as dst:
                    for line in src:
                        fields = line.rstrip('\r\n').split('\t')
                        if len(fields) == 2:
                            dst.write(f'{fields[0]}\t{fields[1]}\tGene Expression\n')
                        else:
                            dst.write(line if line.endswith(('\n', '\r')) else line + '\n')
                continue
            if source_name.endswith('.gz'):
                shutil.copy2(source, target)
            else:
                with open(source, 'rb') as src, gzip.open(target, 'wb') as dst:
                    shutil.copyfileobj(src, dst)
        return sc.read_10x_mtx(temp_dir, var_names=var_names, cache=False)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def read_single_cell_data(input_path, input_format='auto', species=None, genome=None):
    """
    Read common single-cell input formats into AnnData.

    Supported formats:
      - h5ad
      - 10x mtx directory or matrix.mtx path
      - 10x h5
      - loom
      - zarr
      - expression matrix csv/tsv/txt/xlsx/xls
    """
    import scanpy as sc

    if not input_path:
        raise ValueError("缺少输入路径")
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"输入路径不存在: {input_path}")

    fmt = infer_sc_data_format(input_path) if input_format in ('', None, 'auto') else input_format
    source_path = input_path

    if fmt == 'h5ad':
        adata = sc.read_h5ad(source_path)
    elif fmt == '10x_mtx':
        mtx_dir = source_path if os.path.isdir(source_path) else os.path.dirname(source_path)
        if not mtx_dir:
            mtx_dir = '.'
        adata = read_10x_mtx_compat(mtx_dir, var_names='gene_symbols')
    elif fmt == '10x_h5':
        adata = sc.read_10x_h5(source_path)
    elif fmt == 'loom':
        adata = sc.read_loom(source_path)
    elif fmt == 'zarr':
        import anndata as ad
        adata = ad.read_zarr(source_path)
    elif fmt == 'expression_matrix':
        adata = read_expression_matrix(source_path)
    else:
        raise ValueError(
            f"不支持的单细胞输入格式: {fmt}。"
            "支持 h5ad、10x mtx、10x h5、loom、zarr、csv/tsv/xlsx 表达矩阵。"
        )

    return _standardize_imported_adata(adata, input_format=fmt, species=species, genome=genome)


def write_single_cell_h5ad(input_path, output_path, input_format='auto', species=None, genome=None):
    """Read a supported single-cell input and write a standardized h5ad."""
    adata = read_single_cell_data(
        input_path,
        input_format=input_format,
        species=species,
        genome=genome,
    )
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    if os.path.abspath(input_path) != os.path.abspath(output_path):
        adata.write_h5ad(output_path)
    return adata


def summarize_adata_import(adata, input_format, output_path):
    """Build a JSON-serializable import summary."""
    total_elements = int(adata.n_obs) * int(adata.n_vars)
    if total_elements > 0:
        if hasattr(adata.X, 'nnz'):
            nonzero = int(adata.X.nnz)
        else:
            nonzero = int(np.count_nonzero(adata.X))
        sparsity = round((1 - nonzero / total_elements) * 100, 1)
    else:
        sparsity = 0.0

    file_size_mb = round(os.path.getsize(output_path) / (1024 * 1024), 1) if os.path.exists(output_path) else 0.0
    return {
        'input_format': input_format,
        'n_cells': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'sparsity': sparsity,
        'file_size_mb': file_size_mb,
        'output_file': os.path.basename(output_path),
        'obs_columns': list(map(str, adata.obs.columns[:20])),
        'var_columns': list(map(str, adata.var.columns[:20])),
        'layers': list(map(str, adata.layers.keys())),
        'counts_layer_status': 'present' if 'counts' in adata.layers else 'missing',
        'counts_layer_inferred_from_x': bool(
            adata.uns.get('counts_layer_inferred_from_x', False)
        ),
    }


def remap_var_names(adata):
    """将 adata.var_names 从 Ensembl ID 映射为基因名，处理重复名。
    原始 ID 保存到 adata.var['gene_id']。
    如果 var_names 已经是基因名（非 Ensembl ID），不做任何改动。"""
    import re

    def _is_ensembl_id(name):
        """检查是否为 Ensembl 基因 ID（如 ENSG00000139618）"""
        return bool(re.match(r'^ENS[A-Z]*G\d{5,}', str(name)))

    # 如果 var_names 已经是基因名（非 Ensembl ID），直接返回，避免重复处理
    if adata.n_vars > 0:
        sample_ids = list(adata.var_names[:min(100, adata.n_vars)])
        ensembl_count = sum(1 for g in sample_ids if _is_ensembl_id(g))
        if ensembl_count < len(sample_ids) * 0.5:
            return adata  # 已经是基因名，无需映射

    gene_name_col = None
    for col in ['gene_name', 'GeneSymbol', 'gene_symbol', 'gene_symbols', 'symbol',
                'Symbol', 'Gene Symbol', 'gene', 'Gene', 'GENE_NAME',
                'feature_name', 'gene_name_x']:
        if col in adata.var.columns:
            gene_name_col = col
            break
    if gene_name_col is None:
        for col in adata.var.columns:
            if 'symbol' in col.lower() or 'name' in col.lower():
                gene_name_col = col
                break

    if gene_name_col is None:
        return adata  # 没有映射信息，保持原样

    gene_names = adata.var[gene_name_col].astype(str).tolist()
    # 过滤无效值
    gene_names = [g if g and g != 'nan' and g.strip() else str(adata.var_names[i])
                  for i, g in enumerate(gene_names)]

    # 保留原始 ID
    if 'gene_id' not in adata.var.columns:
        adata.var['gene_id'] = adata.var_names.tolist()

    # 处理重复基因名：第一次出现保留原名，后续加后缀
    seen = {}
    unique_names = []
    for g in gene_names:
        if g in seen:
            seen[g] += 1
            unique_names.append(f"{g}_{seen[g]}")
        else:
            seen[g] = 0
            unique_names.append(g)

    adata.var_names = pd.Index(unique_names)
    return adata


def convert_10x_to_h5ad(mtx_dir, output_path, species=None, genome=None):
    """
    将 10x Genomics 三文件格式转换为 h5ad。
    自动检测 v2 (genes.tsv) 和 v3 (features.tsv) 格式。

    参数:
        mtx_dir: 包含 barcodes/genes/features/matrix 文件的目录
        output_path: h5ad 输出路径
        species: 可选，物种名（如 "human"、"mouse"）
        genome: 可选，基因组版本（如 "GRCh38"、"mm10"）
    返回:
        anndata.AnnData 对象
    """
    adata = read_10x_mtx_compat(mtx_dir, var_names='gene_symbols')
    adata.var_names_make_unique()

    # 保留 Ensembl ID（read_10x_mtx 在 var_names='gene_symbols' 时
    # 将原始 ID 存为 adata.var 的 gene_ids 列）
    if 'gene_ids' not in adata.var.columns:
        import warnings
        warnings.warn("未能从 10x 数据中提取 Ensembl gene IDs，使用当前 var_names 作为 gene_ids")
        adata.var['gene_ids'] = adata.var.index.tolist()

    # 可选元数据
    if species:
        adata.uns['species'] = species
    if genome:
        adata.uns['genome'] = genome

    # 保存原始计数
    adata.layers['counts'] = adata.X.copy()

    adata.write_h5ad(output_path)
    print(f"[io_utils] 10x 数据转换完成: {mtx_dir} -> {output_path} ({adata.n_obs} cells, {adata.n_vars} genes)")
    return adata
