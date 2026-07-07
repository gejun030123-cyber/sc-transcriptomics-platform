# modules/cell_markers.py
"""细胞类型 marker 基因库 — 内置参考 + 用户自定义优先"""

# 内置 marker 库：细胞类型 → {positive, negative}
# 数据来源：文献综合，仅作为默认建议
BUILTIN_MARKERS = {
    't_cell': {
        'display': 'T 细胞',
        'positive': ['CD3D', 'CD3E', 'TRAC'],
        'negative': ['MS4A1', 'LST1'],
    },
    'b_cell': {
        'display': 'B 细胞',
        'positive': ['MS4A1', 'CD79A', 'CD79B'],
        'negative': ['CD3D', 'LST1'],
    },
    'nk': {
        'display': 'NK 细胞',
        'positive': ['NKG7', 'GNLY', 'KLRD1'],
        'negative': ['CD3D', 'MS4A1'],
    },
    'monocyte': {
        'display': '单核细胞',
        'positive': ['LST1', 'S100A8', 'S100A9', 'FCN1'],
        'negative': ['MS4A1', 'CD3D'],
    },
    'macrophage': {
        'display': '巨噬细胞',
        'positive': ['C1QA', 'C1QB', 'CD68', 'APOE'],
        'negative': ['MS4A1', 'CD3D'],
    },
    'dendritic_cell': {
        'display': '树突状细胞',
        'positive': ['FCER1A', 'CLEC10A', 'LILRA4'],
        'negative': ['CD3D', 'MS4A1'],
    },
    'epithelial': {
        'display': '上皮细胞',
        'positive': ['EPCAM', 'KRT8', 'KRT18'],
        'negative': ['PTPRC'],
    },
    'endothelial': {
        'display': '内皮细胞',
        'positive': ['PECAM1', 'VWF', 'KDR'],
        'negative': ['PTPRC'],
    },
    'fibroblast': {
        'display': '成纤维细胞',
        'positive': ['COL1A1', 'COL1A2', 'DCN'],
        'negative': ['PTPRC', 'EPCAM'],
    },
    'microglia': {
        'display': '小胶质细胞',
        'positive': ['P2RY12', 'TMEM119', 'CX3CR1', 'AIF1'],
        'negative': ['S100A8', 'FCGR3A'],
    },
    'astrocyte': {
        'display': '星形胶质细胞',
        'positive': ['GFAP', 'AQP4', 'ALDH1L1'],
        'negative': ['PTPRC'],
    },
    'oligodendrocyte': {
        'display': '少突胶质细胞',
        'positive': ['MBP', 'MOG', 'PLP1'],
        'negative': ['PTPRC'],
    },
}


def get_markers(cell_type, user_positive=None, user_negative=None):
    """
    获取指定细胞类型的 marker 基因列表。
    用户自定义 marker 优先，内置 marker 仅作为回退建议。

    参数：
        cell_type: str, 细胞类型 key（如 'microglia', 't_cell'）
        user_positive: list[str]|None, 用户指定的 positive markers
        user_negative: list[str]|None, 用户指定的 negative markers

    返回：
        dict: {
            'cell_type': str,
            'positive': list[str],
            'negative': list[str],
            'source': 'user' | 'builtin' | 'mixed',
            'warnings': list[str],
        }
    """
    builtin = BUILTIN_MARKERS.get(cell_type.lower(), {})
    builtin_pos = builtin.get('positive', [])
    builtin_neg = builtin.get('negative', [])
    warnings = []

    if user_positive and user_negative:
        return {
            'cell_type': cell_type,
            'positive': [g.upper().strip() for g in user_positive],
            'negative': [g.upper().strip() for g in user_negative],
            'source': 'user',
            'warnings': [],
        }

    if user_positive and not user_negative:
        return {
            'cell_type': cell_type,
            'positive': [g.upper().strip() for g in user_positive],
            'negative': [g.upper().strip() for g in (user_negative or builtin_neg)],
            'source': 'mixed',
            'warnings': [] if user_negative else warnings,
        }

    if not builtin_pos:
        warnings.append(f"未找到 '{cell_type}' 的内置 marker，请提供自定义 marker")
        return {
            'cell_type': cell_type,
            'positive': [],
            'negative': [],
            'source': 'none',
            'warnings': warnings,
        }

    return {
        'cell_type': cell_type,
        'positive': builtin_pos,
        'negative': builtin_neg,
        'source': 'builtin',
        'warnings': warnings,
    }


def list_builtin_cell_types():
    """列出所有内置细胞类型."""
    return {
        key: {'display': info['display'], 'n_positive': len(info['positive']),
              'n_negative': len(info['negative'])}
        for key, info in BUILTIN_MARKERS.items()
    }
