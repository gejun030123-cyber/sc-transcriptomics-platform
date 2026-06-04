import os
import numpy as np
import pandas as pd


def read_expression_matrix(file_path):
    """读取表达矩阵文件，自动检测格式（Excel/CSV/TSV/h5ad）"""
    import scanpy as sc

    if file_path.endswith('.h5ad'):
        return sc.read_h5ad(file_path)

    # 尝试 Excel
    if file_path.endswith(('.xlsx', '.xls')):
        try:
            df = pd.read_excel(file_path, index_col=0)
            return sc.AnnData(X=df.values.T, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))
        except Exception:
            pass  # 可能是伪 Excel 文件（实际为 TSV/CSV）

    # 尝试自动检测分隔符
    try:
        df = pd.read_csv(file_path, sep='\t', index_col=0, nrows=5)
        if df.shape[1] > 0:
            df = pd.read_csv(file_path, sep='\t', index_col=0)
            return sc.AnnData(X=df.values.T, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))
    except Exception:
        pass

    try:
        df = pd.read_csv(file_path, sep=None, engine='python', index_col=0, nrows=5)
        if df.shape[1] > 0:
            df = pd.read_csv(file_path, sep=None, engine='python', index_col=0)
            return sc.AnnData(X=df.values.T, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))
    except Exception:
        pass

    # 最后尝试标准 CSV
    df = pd.read_csv(file_path, index_col=0)
    return sc.AnnData(X=df.values.T, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))
