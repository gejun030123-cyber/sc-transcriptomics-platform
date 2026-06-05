import os
import numpy as np
import pandas as pd


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
        except Exception:
            df = None

    # 尝试自动检测
    if df is None:
        try:
            df = pd.read_csv(file_path, sep=None, engine='python', index_col=0, nrows=5)
            if df.shape[1] > 0:
                df = pd.read_csv(file_path, sep=None, engine='python', index_col=0)
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
