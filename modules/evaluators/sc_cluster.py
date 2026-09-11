# modules/evaluators/sc_cluster.py
"""单细胞分群签名评分器 — 确定性算法，不依赖 LLM 主观判断。

评分只使用**有界且与表达量尺度无关**的 marker 证据（簇内相对表达、
簇内检出率、簇内外特异性、负向 marker 比例）；cluster 大小与批次占比
（size_score / batch_penalty）作为独立的质控提示单独报告，不再进入
marker 得分，避免用与 marker 证据无关的量改变“置信度”。

表达量尺度在处理前显式解析：
- X 为非负整数 counts → 先在副本上 normalize_total + log1p；
- X 为非负连续值（log1p 标准化结果）→ 直接使用；
- X 含负值（已 scale / Pearson 残差）→ 改用 layers['counts']，否则明确报错。

任何一步失败都会返回 ``error``，不再用 ``except Exception: pass`` 静默把
失败当成 0 分（那会把“没能读取表达量”伪装成“没有 marker 证据”）。
"""

import os
import logging
import numpy as np

logger = logging.getLogger(__name__)

# 在 log1p 尺度上判定“该 marker 在簇内可检出”的最小平均表达量。
MARKER_DETECTION_CUTOFF = 0.1
# 与 normalize 模块默认值一致的文库大小目标。
NORMALIZATION_TARGET_SUM = 1e4
# 用于判定表达量尺度的抽样上限。
_SCALE_SAMPLE_LIMIT = 200000


def _matrix_sample_values(matrix, max_values=_SCALE_SAMPLE_LIMIT):
    """Return a bounded, deterministic sample of stored matrix values."""
    from scipy import sparse

    if sparse.issparse(matrix):
        values = np.asarray(matrix.data).ravel()
    else:
        values = np.asarray(matrix).ravel()
        if values.size > max_values:
            index = np.linspace(0, values.size - 1, int(max_values), dtype=int)
            values = values[index]
    return np.asarray(values, dtype=float)


def _has_negative_values(matrix):
    values = _matrix_sample_values(matrix)
    finite = values[np.isfinite(values)]
    return bool(finite.size and float(finite.min()) < -1e-6)


def _is_count_matrix(matrix):
    values = _matrix_sample_values(matrix)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return False
    return bool(
        float(finite.min()) >= -1e-6
        and np.allclose(finite, np.rint(finite), atol=1e-6)
    )


def _load_adata(adata_path):
    """安全加载 AnnData 并物化到内存。返回 (adata, error).

    The scorer needs per-cell expression and library sizes for the whole
    matrix, so backed mode only added a second (partially broken) slicing
    path; materialising once keeps every read on the same in-memory object.
    """
    if not os.path.isfile(adata_path):
        return None, f"文件不存在: {adata_path}"
    if os.path.islink(adata_path):
        return None, "不支持符号链接文件"
    try:
        import scanpy as sc
    except ImportError:
        return None, "scanpy 未安装"
    backed = None
    try:
        backed = sc.read(adata_path, backed='r')
        adata = backed.to_memory()
    except Exception as exc:
        logger.debug("backed 读取失败，改用完整读取: %s", exc)
        try:
            adata = sc.read(adata_path)
        except Exception as exc2:
            return None, f"无法读取 AnnData: {exc2}"
    finally:
        if backed is not None:
            try:
                backed.file.close()
            except Exception as exc:  # 句柄泄漏不应掩盖真正的读取错误
                logger.debug("关闭 backed 句柄失败: %s", exc)
    return adata, None


def _scorable_expression(adata):
    """Return (adata, scale_label, error) with marker-comparable expression."""
    import scanpy as sc

    if not _has_negative_values(adata.X):
        if _is_count_matrix(adata.X):
            normalized = adata.copy()
            sc.pp.normalize_total(normalized, target_sum=NORMALIZATION_TARGET_SUM)
            sc.pp.log1p(normalized)
            return normalized, 'log1p_normalized_from_X_counts', None
        return adata, 'assumed_log_normalized_X', None

    counts = adata.layers.get('counts')
    if counts is not None and not _has_negative_values(counts):
        normalized = adata.copy()
        normalized.X = counts.copy() if hasattr(counts, 'copy') else np.array(counts)
        sc.pp.normalize_total(normalized, target_sum=NORMALIZATION_TARGET_SUM)
        sc.pp.log1p(normalized)
        return normalized, 'log1p_normalized_from_layers_counts', None

    return None, '', (
        "marker 评分需要非负表达矩阵；当前 X 含负值（已 scale 或 Pearson 残差），"
        "且没有可用的 layers['counts']。请改用 log1p 标准化后的 h5ad。"
    )


def _normalize_gene_names(var_names, markers):
    """匹配基因名（大小写不敏感）。返回 (found: list, missing: list)."""
    upper_to_orig = {}
    for g in var_names:
        upper_to_orig.setdefault(str(g).upper(), str(g))

    found = []
    missing = []
    for m in markers:
        if str(m).upper() in upper_to_orig:
            found.append(upper_to_orig[str(m).upper()])
        else:
            missing.append(m)
    return found, missing


def _marker_matrix(adata, genes):
    """Return (values, columns, error): an n_cells × n_genes float matrix."""
    if not genes:
        return np.zeros((adata.n_obs, 0), dtype=float), {}, None
    index = adata.var_names
    try:
        positions = {gene: int(index.get_loc(gene)) for gene in genes}
    except Exception as exc:
        # get_loc 在基因名重复时返回数组，__int__ 会失败。
        return None, {}, f"无法定位 marker 基因（基因名可能重复）: {exc}"
    columns = sorted(set(positions.values()))
    column_of = {column: order for order, column in enumerate(columns)}
    try:
        sub = adata.X[:, columns]
    except Exception as exc:
        return None, {}, f"无法读取 marker 表达矩阵: {exc}"
    if hasattr(sub, 'toarray'):
        try:
            sub = sub.toarray()
        except Exception as exc:
            return None, {}, f"无法展开稀疏 marker 矩阵: {exc}"
    values = np.asarray(sub, dtype=float)
    if values.ndim != 2 or values.shape[0] != adata.n_obs:
        return None, {}, "marker 表达矩阵形状与细胞数不一致"
    return values, {gene: column_of[position] for gene, position in positions.items()}, None


def _size_score(pct):
    """Descriptive cluster-size heuristic (reported, not part of marker score)."""
    if 0.01 <= pct <= 0.45:
        return 1.0
    if pct < 0.001:
        return 0.2
    if pct > 0.60:
        return 0.6
    return max(0.3, 1.0 - abs(pct - 0.1) * 3)


def score_cluster_signature(adata_path, cluster_key, positive_markers, negative_markers):
    """
    对 AnnData 中的每个 cluster 计算细胞类型签名评分。

    返回：
        dict: {
            'cluster_key': str,
            'expression_scale': str,
            'positive_markers_found': [...],
            'positive_markers_missing': [...],
            'negative_markers_found': [...],
            'negative_markers_missing': [...],
            'warnings': [...],
            'cluster_scores': [
                {
                    'cluster': str,
                    'n_cells': int,
                    'positive_score': float,        # 簇内 marker 平均表达（展示用）
                    'positive_evidence': float,     # 相对最高簇的 0-1 证据
                    'detection_fraction': float,    # 可检出 positive marker 比例
                    'specificity_score': float,
                    'negative_score': float,
                    'negative_ratio': float,
                    'size_score': float,
                    'batch_penalty': float,
                    'qc_flags': [...],
                    'total_score': float,           # 仅由 marker 证据构成（0-1）
                    'detected_positive_markers': [...],
                    'undetected_positive_markers': [...],
                    'missing_positive_markers': [...],
                }
            ]
        }
    """
    empty = {
        'cluster_key': cluster_key,
        'cluster_scores': [],
        'warnings': [],
    }

    adata, err = _load_adata(adata_path)
    if err:
        return {**empty, 'error': err}

    warnings = []

    # 检查 cluster_key
    if cluster_key not in adata.obs.columns:
        return {
            **empty,
            'error': f"cluster_key '{cluster_key}' 不在 obs 列中。可用列: {list(adata.obs.columns)}",
        }
    from modules.io_utils import obs_grouping_info
    grouping = obs_grouping_info(
        adata, cluster_key, max_categories=50,
        max_numeric_categories=20, require_multiple=False,
    )
    if not grouping['valid']:
        return {
            **empty,
            'error': f"cluster_key '{cluster_key}' 不是有效的分类聚类列：{grouping['reason']}",
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
            **empty,
            'positive_markers_found': pos_found,
            'positive_markers_missing': pos_missing,
            'negative_markers_found': neg_found,
            'negative_markers_missing': neg_missing,
            'warnings': warnings,
            'error': '没有匹配到任何 marker 基因',
        }

    adata, expression_scale, scale_error = _scorable_expression(adata)
    if scale_error:
        return {
            **empty,
            'positive_markers_found': pos_found,
            'positive_markers_missing': pos_missing,
            'negative_markers_found': neg_found,
            'negative_markers_missing': neg_missing,
            'warnings': warnings,
            'error': scale_error,
        }

    pos_values, pos_columns, pos_error = _marker_matrix(adata, pos_found)
    if pos_error:
        return {
            **empty,
            'positive_markers_found': pos_found,
            'positive_markers_missing': pos_missing,
            'negative_markers_found': neg_found,
            'negative_markers_missing': neg_missing,
            'warnings': warnings,
            'error': pos_error,
        }
    neg_values, neg_columns, neg_error = _marker_matrix(adata, neg_found)
    if neg_error:
        return {
            **empty,
            'positive_markers_found': pos_found,
            'positive_markers_missing': pos_missing,
            'negative_markers_found': neg_found,
            'negative_markers_missing': neg_missing,
            'warnings': warnings,
            'error': neg_error,
        }

    clusters = sorted(adata.obs[cluster_key].astype(str).unique())
    labels = adata.obs[cluster_key].astype(str).to_numpy()

    # 检查 batch 列（仅作为质控提示）
    batch_col = None
    for col in ['batch', 'sample', 'sample_id', 'orig.ident']:
        if col in adata.obs.columns:
            batch_col = col
            break
    batch_labels = (
        adata.obs[batch_col].astype(str).to_numpy() if batch_col else None
    )

    total_cells = int(adata.n_obs)
    raw_stats = []

    # ── 第一遍：收集每簇的有界证据分量 ────────────────────────────────
    for cluster in clusters:
        mask = labels == cluster
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue

        pos_means = (
            pos_values[mask].mean(axis=0) if pos_columns else np.zeros(0, dtype=float)
        )
        neg_means = (
            neg_values[mask].mean(axis=0) if neg_columns else np.zeros(0, dtype=float)
        )
        positive_score = float(np.mean(pos_means)) if pos_means.size else 0.0
        negative_score = float(np.mean(neg_means)) if neg_means.size else 0.0

        if pos_means.size:
            detected_mask = pos_means >= MARKER_DETECTION_CUTOFF
            detection_fraction = float(detected_mask.mean())
            detected = [
                gene for gene in pos_found
                if pos_means[pos_columns[gene]] >= MARKER_DETECTION_CUTOFF
            ]
            undetected = [gene for gene in pos_found if gene not in detected]
        else:
            detection_fraction = 0.0
            detected, undetected = [], []

        # 簇内 vs 簇外表达差（每个 marker 先算比值再平均，保持 0-1 有界）
        specificity_scores = []
        if pos_columns and n_cells < total_cells:
            inner = pos_values[mask].mean(axis=0)
            outer = pos_values[~mask].mean(axis=0)
            specificity_scores = [
                max(0.0, float((i - o) / (i + o + 1e-8)))
                for i, o in zip(inner, outer)
            ]
        specificity_score = (
            float(np.mean(specificity_scores)) if specificity_scores else 0.0
        )

        pct = n_cells / total_cells if total_cells else 0.0
        size_score = _size_score(pct)

        qc_flags = []
        batch_penalty = 0.0
        if batch_labels is not None:
            values, counts = np.unique(batch_labels[mask], return_counts=True)
            if values.size <= 1:
                batch_penalty = 0.5
                qc_flags.append('single_batch_only')
            else:
                max_pct = float(counts.max()) / float(counts.sum())
                if max_pct > 0.95:
                    batch_penalty = 0.3
                    qc_flags.append('dominated_by_one_batch')
                elif max_pct > 0.80:
                    batch_penalty = 0.15
                    qc_flags.append('mostly_one_batch')
        if pct < 0.01:
            qc_flags.append('small_cluster')
        if pos_columns and detection_fraction < 0.5:
            qc_flags.append('few_detected_markers')

        raw_stats.append({
            'cluster': cluster,
            'n_cells': n_cells,
            'pct_cells': pct,
            'positive_score': positive_score,
            'negative_score': negative_score,
            'detection_fraction': detection_fraction,
            'specificity_score': specificity_score,
            'size_score': size_score,
            'batch_penalty': batch_penalty,
            'qc_flags': qc_flags,
            'detected_positive_markers': detected,
            'undetected_positive_markers': undetected,
        })

    if not raw_stats:
        return {
            **empty,
            'positive_markers_found': pos_found,
            'positive_markers_missing': pos_missing,
            'negative_markers_found': neg_found,
            'negative_markers_missing': neg_missing,
            'warnings': warnings,
            'error': '没有可用于评分的 cluster',
        }

    # ── 第二遍：相对尺度归一化后合成 marker 得分 ─────────────────────
    # 用簇间最高平均表达做分母，使得分与表达量尺度（counts / log1p）无关。
    positive_scale = max(stat['positive_score'] for stat in raw_stats)

    cluster_scores = []
    for stat in raw_stats:
        positive_evidence = (
            stat['positive_score'] / positive_scale if positive_scale > 0 else 0.0
        )
        negative_ratio = (
            stat['negative_score']
            / (stat['negative_score'] + stat['positive_score'] + 1e-8)
        )
        marker_score = (
            0.40 * positive_evidence
            + 0.25 * stat['detection_fraction']
            + 0.35 * stat['specificity_score']
            - 0.30 * negative_ratio
        )
        total_score = float(max(0.0, min(1.0, marker_score)))
        cluster_scores.append({
            'cluster': stat['cluster'],
            'n_cells': stat['n_cells'],
            'pct_cells': round(stat['pct_cells'], 4),
            'positive_score': round(stat['positive_score'], 4),
            'positive_evidence': round(positive_evidence, 4),
            'detection_fraction': round(stat['detection_fraction'], 4),
            'specificity_score': round(stat['specificity_score'], 4),
            'negative_score': round(stat['negative_score'], 4),
            'negative_ratio': round(negative_ratio, 4),
            'size_score': round(stat['size_score'], 4),
            'batch_penalty': round(stat['batch_penalty'], 4),
            'qc_flags': stat['qc_flags'],
            'total_score': round(total_score, 4),
            'detected_positive_markers': stat['detected_positive_markers'],
            'undetected_positive_markers': stat['undetected_positive_markers'],
            'missing_positive_markers': pos_missing,
        })

    # 按 total_score 降序排列
    cluster_scores.sort(key=lambda item: item['total_score'], reverse=True)

    return {
        'cluster_key': cluster_key,
        'expression_scale': expression_scale,
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

    result['target_cell_type'] = target_cell_type
    result['marker_source'] = marker_info['source']
    result['marker_warnings'] = marker_info.get('warnings', [])

    if result.get('error'):
        result['best_cluster'] = None
        result['confidence'] = 'none'
        return result

    scores = result.get('cluster_scores', [])

    # 置信度只由 marker 证据（total_score）决定；cluster 大小与批次占比
    # 作为 qc_flags 单独报告，不再改变置信度等级。
    if not scores:
        confidence = 'none'
        best_cluster = None
        qc_flags = []
    else:
        top = scores[0]
        best_cluster = top['cluster']
        qc_flags = list(top.get('qc_flags') or [])

        if top['total_score'] >= 0.75:
            confidence = 'high'
        elif top['total_score'] >= 0.50:
            confidence = 'medium'
        elif top['total_score'] >= 0.25:
            confidence = 'low'
        else:
            confidence = 'very_low'

    result['best_cluster'] = best_cluster
    result['confidence'] = confidence
    result['qc_flags'] = qc_flags
    return result
