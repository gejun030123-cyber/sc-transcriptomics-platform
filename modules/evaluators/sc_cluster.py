# modules/evaluators/sc_cluster.py
"""单细胞分群签名评分器 — 确定性算法，不依赖 LLM 主观判断"""

import os
import logging
import numpy as np

logger = logging.getLogger(__name__)


def _load_adata(adata_path):
    """安全加载 AnnData（backed 模式优先）。返回 (adata, error)."""
    # 先检查路径合法性，避免无效路径触发 scanpy/numba 初始化
    if not os.path.isfile(adata_path):
        return None, f"文件不存在: {adata_path}"
    if os.path.islink(adata_path):
        return None, "不支持符号链接文件"
    try:
        import scanpy as sc
    except ImportError:
        return None, "scanpy 未安装"
    try:
        # 大文件使用 backed 模式
        adata = sc.read(adata_path, backed='r')
        return adata, None
    except Exception:
        try:
            adata = sc.read(adata_path)
            return adata, None
        except Exception as e:
            return None, f"无法读取 AnnData: {e}"


def _normalize_gene_names(var_names, markers):
    """匹配基因名（大小写不敏感）。返回 (found: list, missing: list)."""
    upper_to_orig = {}
    for g in var_names:
        upper_to_orig[g.upper()] = g

    found = []
    missing = []
    for m in markers:
        if m.upper() in upper_to_orig:
            found.append(upper_to_orig[m.upper()])
        else:
            missing.append(m)
    return found, missing


def score_cluster_signature(adata_path, cluster_key, positive_markers, negative_markers):
    """
    对 AnnData 中的每个 cluster 计算细胞类型签名评分。

    返回：
        dict: {
            'cluster_key': str,
            'positive_markers_found': [...],
            'positive_markers_missing': [...],
            'negative_markers_found': [...],
            'negative_markers_missing': [...],
            'warnings': [...],
            'cluster_scores': [
                {
                    'cluster': str,
                    'n_cells': int,
                    'positive_score': float,
                    'negative_score': float,
                    'specificity_score': float,
                    'size_score': float,
                    'batch_penalty': float,
                    'total_score': float,
                    'detected_positive_markers': [...],
                    'missing_positive_markers': [...],
                }
            ]
        }
    """
    adata, err = _load_adata(adata_path)
    if err:
        return {'error': err, 'cluster_key': cluster_key, 'cluster_scores': [], 'warnings': []}

    warnings = []

    # 检查 cluster_key
    if cluster_key not in adata.obs.columns:
        return {
            'error': f"cluster_key '{cluster_key}' 不在 obs 列中。可用列: {list(adata.obs.columns)}",
            'cluster_key': cluster_key,
            'cluster_scores': [],
            'warnings': [],
        }
    from modules.io_utils import obs_grouping_info
    grouping = obs_grouping_info(
        adata, cluster_key, max_categories=50,
        max_numeric_categories=20, require_multiple=False,
    )
    if not grouping['valid']:
        return {
            'error': f"cluster_key '{cluster_key}' 不是有效的分类聚类列：{grouping['reason']}",
            'cluster_key': cluster_key,
            'cluster_scores': [],
            'warnings': [grouping['reason']],
        }

    # 匹配 marker 基因
    pos_found, pos_missing = _normalize_gene_names(adata.var_names, positive_markers)
    neg_found, neg_missing = _normalize_gene_names(adata.var_names, negative_markers)

    if pos_missing:
        warnings.append(f"未找到 positive markers: {pos_missing}")
    if neg_missing:
        warnings.append(f"未找到 negative markers: {neg_missing}")

    if not pos_found and not neg_found:
        return {
            'cluster_key': cluster_key,
            'positive_markers_found': pos_found,
            'positive_markers_missing': pos_missing,
            'negative_markers_found': neg_found,
            'negative_markers_missing': neg_missing,
            'warnings': warnings,
            'cluster_scores': [],
            'error': '没有匹配到任何 marker 基因',
        }

    clusters = sorted(adata.obs[cluster_key].unique().astype(str))

    # 检查 batch 列
    batch_col = None
    for col in ['batch', 'sample', 'sample_id', 'orig.ident']:
        if col in adata.obs.columns:
            batch_col = col
            break

    total_cells = adata.n_obs
    cluster_scores = []

    for cl in clusters:
        mask = adata.obs[cluster_key].astype(str) == cl
        n_cells = int(mask.sum())

        if n_cells == 0:
            continue

        # 1. positive_score: cluster 内 positive marker 平均表达
        pos_exprs = []
        for g in pos_found:
            try:
                expr = adata[mask, g].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray()
                pos_exprs.append(float(np.mean(expr)))
            except Exception:
                pass

        positive_score = float(np.mean(pos_exprs)) if pos_exprs else 0.0

        # 2. specificity_score: cluster 内外表达差异
        spec_scores = []
        for g in pos_found:
            try:
                inner = adata[mask, g].X
                outer = adata[~mask, g].X
                if hasattr(inner, 'toarray'):
                    inner = inner.toarray()
                    outer = outer.toarray()
                inner_mean = float(np.mean(inner))
                outer_mean = float(np.mean(outer))
                if outer_mean > 0:
                    spec = (inner_mean - outer_mean) / (inner_mean + outer_mean + 1e-8)
                    spec_scores.append(max(0.0, spec))
                else:
                    spec_scores.append(0.5)
            except Exception:
                pass

        specificity_score = float(np.mean(spec_scores)) if spec_scores else 0.0

        # 3. negative_score: negative marker 表达惩罚
        neg_exprs = []
        for g in neg_found:
            try:
                expr = adata[mask, g].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray()
                neg_exprs.append(float(np.mean(expr)))
            except Exception:
                pass

        neg_penalty = float(np.mean(neg_exprs)) if neg_exprs else 0.0

        # 4. size_score: cluster 大小是否合理（太小/太大扣分）
        pct = n_cells / total_cells
        if 0.01 <= pct <= 0.45:
            size_score = 1.0
        elif pct < 0.001:
            size_score = 0.2
        elif pct > 0.60:
            size_score = 0.6
        else:
            size_score = max(0.3, 1.0 - abs(pct - 0.1) * 3)

        # 5. batch_penalty: 是否被单一 batch/sample 主导
        batch_penalty = 0.0
        if batch_col and batch_col in adata.obs.columns:
            try:
                batch_counts = adata[mask].obs[batch_col].value_counts()
                if len(batch_counts) > 1:
                    max_pct = batch_counts.max() / batch_counts.sum()
                    if max_pct > 0.95:
                        batch_penalty = 0.3
                    elif max_pct > 0.80:
                        batch_penalty = 0.15
                    else:
                        batch_penalty = 0.0
                else:
                    batch_penalty = 0.5  # 只有一个 batch
            except Exception:
                pass

        # 综合评分
        total_score = float(
            0.35 * positive_score +
            0.25 * specificity_score -
            0.20 * neg_penalty +
            0.10 * size_score -
            0.10 * batch_penalty
        )
        total_score = max(0.0, min(1.0, total_score))

        cluster_scores.append({
            'cluster': cl,
            'n_cells': n_cells,
            'pct_cells': round(n_cells / total_cells, 4),
            'positive_score': round(positive_score, 4),
            'specificity_score': round(specificity_score, 4),
            'negative_score': round(neg_penalty, 4),
            'size_score': round(size_score, 4),
            'batch_penalty': round(batch_penalty, 4),
            'total_score': round(total_score, 4),
            'detected_positive_markers': pos_found,
            'missing_positive_markers': pos_missing,
        })

    # 按 total_score 降序排列
    cluster_scores.sort(key=lambda x: x['total_score'], reverse=True)

    return {
        'cluster_key': cluster_key,
        'positive_markers_found': pos_found,
        'positive_markers_missing': pos_missing,
        'negative_markers_found': neg_found,
        'negative_markers_missing': neg_missing,
        'warnings': warnings,
        'cluster_scores': cluster_scores,
    }


def score_cell_type_signature(adata_path, cluster_key, target_cell_type,
                               positive_markers=None, negative_markers=None):
    """
    高层包装：获取 marker 库 + 调用评分引擎。

    参数：
        adata_path: str, h5ad 文件路径
        cluster_key: str, 聚类键（如 'leiden'）
        target_cell_type: str, 目标细胞类型 key
        positive_markers: list[str]|None, 用户自定义 positive markers
        negative_markers: list[str]|None, 用户自定义 negative markers

    返回：
        dict: 评分结果 + best_cluster + confidence
    """
    from modules.cell_markers import get_markers

    marker_info = get_markers(target_cell_type, positive_markers, negative_markers)

    result = score_cluster_signature(
        adata_path, cluster_key,
        marker_info['positive'], marker_info['negative']
    )

    if result.get('error'):
        result['target_cell_type'] = target_cell_type
        result['best_cluster'] = None
        result['confidence'] = 'none'
        return result

    scores = result.get('cluster_scores', [])

    # 确定最佳 cluster 和置信度
    if not scores:
        confidence = 'none'
        best_cluster = None
    else:
        top = scores[0]
        best_cluster = top['cluster']

        if top['total_score'] >= 0.75:
            confidence = 'high'
        elif top['total_score'] >= 0.50:
            confidence = 'medium'
        elif top['total_score'] >= 0.25:
            confidence = 'low'
        else:
            confidence = 'very_low'

    # 补充 marker 来源信息
    result['target_cell_type'] = target_cell_type
    result['best_cluster'] = best_cluster
    result['confidence'] = confidence
    result['marker_source'] = marker_info['source']
    result['marker_warnings'] = marker_info.get('warnings', [])

    return result
