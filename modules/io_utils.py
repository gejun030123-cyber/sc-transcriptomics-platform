import os
import numpy as np
import pandas as pd


def read_expression_matrix(file_path):
    """读取表达矩阵文件，自动检测格式，只保留数值列"""
    import scanpy as sc

    if file_path.endswith('.h5ad'):
        return sc.read_h5ad(file_path)

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

    if len(numeric_cols) < df.shape[1]:
        df = df[numeric_cols]

    # 替换 NaN 为 0
    df = df.fillna(0)

    return sc.AnnData(X=df.values.T, obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))
