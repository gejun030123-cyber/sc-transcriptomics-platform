import os
import numpy as np
import pandas as pd


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
        if len(tokens) < 2 or not re.fullmatch(r'(?:rep)?\d+', tokens[-1], re.I):
            return []
        parsed.append((name, tokens[:-1]))

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


def _ensure_counts_layer(adata):
    """Preserve raw/imported matrix in counts layer if absent."""
    if 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()
    return adata


def _standardize_imported_adata(adata, input_format=None, species=None, genome=None):
    """Apply lightweight AnnData normalization needed by downstream modules."""
    if hasattr(adata, 'var_names_make_unique'):
        adata.var_names_make_unique()
    if hasattr(adata, 'obs_names_make_unique'):
        adata.obs_names_make_unique()

    adata = remap_var_names(adata)
    adata = _ensure_counts_layer(adata)

    if input_format:
        adata.uns['input_format'] = input_format
    if species:
        adata.uns['species'] = species
    if genome:
        adata.uns['genome'] = genome
    return adata


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
        adata = sc.read_10x_mtx(mtx_dir, var_names='gene_symbols', cache=True)
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
    if not os.path.isdir(mtx_dir):
        raise FileNotFoundError(f"10x 矩阵目录不存在: {mtx_dir}")
    import scanpy as sc

    # Recent Scanpy versions expect compressed 10x filenames by default.  The
    # web uploader also accepts plain .mtx/.tsv files, so create a temporary
    # gzip view when a ZIP contains the uncompressed form.
    read_dir = mtx_dir
    temp_dir = None
    plain_to_gzip = {
        'matrix.mtx': 'matrix.mtx.gz',
        'barcodes.tsv': 'barcodes.tsv.gz',
        'features.tsv': 'features.tsv.gz',
        'genes.tsv': 'genes.tsv.gz',
    }
    if any(os.path.exists(os.path.join(mtx_dir, plain))
           and not os.path.exists(os.path.join(mtx_dir, compressed))
           for plain, compressed in plain_to_gzip.items()):
        import gzip
        import shutil
        import tempfile

        temp_dir = tempfile.mkdtemp(prefix='10x_gzip_')
        for plain, compressed in plain_to_gzip.items():
            source = os.path.join(mtx_dir, plain)
            if not os.path.isfile(source):
                continue
            target = os.path.join(temp_dir, compressed)
            with open(source, 'rb') as src, gzip.open(target, 'wb') as dst:
                shutil.copyfileobj(src, dst)
        read_dir = temp_dir
    try:
        adata = sc.read_10x_mtx(read_dir, var_names='gene_symbols', cache=True)
    finally:
        if temp_dir:
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)
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
