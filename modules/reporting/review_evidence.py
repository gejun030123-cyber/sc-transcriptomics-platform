"""Review evidence cards for analysis outputs."""

import os


SC_REVIEW_MODULES = {"qc", "clustering", "annotation", "deg"}


def _as_dict(obj):
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    if hasattr(obj, "__dict__"):
        return vars(obj)
    return dict(obj)


def _result_file_names(result_files):
    names = []
    for rf in result_files or []:
        item = _as_dict(rf)
        path = item.get("file_path", "")
        label = item.get("label", "")
        category = item.get("category", "")
        names.append(" ".join([os.path.basename(path), label, category]).lower())
    return names


def _has_file(names, *needles):
    return any(all(needle.lower() in name for needle in needles) for name in names)


def _status_rank(status):
    return {"pass": 0, "review": 1, "warning": 2}.get(status, 1)


def _overall(checks):
    if not checks:
        return "review"
    return max((check["status"] for check in checks), key=_status_rank)


def _check(name, status, message, value=None):
    item = {"name": name, "status": status, "message": message}
    if value is not None:
        item["value"] = value
    return item


def _number(value, default=None):
    try:
        if value is None:
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _qc_evidence(summary, names):
    cells_before = _number(summary.get("cells_before"), 0)
    cells_after = _number(summary.get("cells_after"), 0)
    pct_removed = _number(summary.get("pct_removed"))
    if pct_removed is None and cells_before:
        pct_removed = round((cells_before - cells_after) / max(cells_before, 1) * 100, 1)

    checks = []
    if not cells_before:
        checks.append(_check("细胞保留率", "review", "summary 中缺少 QC 前细胞数，需人工确认过滤强度。"))
    elif pct_removed is not None and pct_removed > 70:
        checks.append(_check("细胞保留率", "warning", "QC 移除比例过高，需确认阈值是否过严。", f"{pct_removed}%"))
    elif pct_removed is not None and pct_removed < 1:
        checks.append(_check("细胞保留率", "review", "几乎没有细胞被过滤，需确认 QC 阈值是否过宽。", f"{pct_removed}%"))
    else:
        checks.append(_check("细胞保留率", "pass", "QC 过滤比例处于常见复核范围。", f"{pct_removed}%"))

    if summary.get("novelty_median") is not None:
        checks.append(_check("复杂度指标", "pass", "已记录 novelty score 中位数。", summary.get("novelty_median")))
    else:
        checks.append(_check("复杂度指标", "review", "未记录 novelty score，中间结果仍需结合 QC 图复核。"))

    if _has_file(names, "qc_filter_summary"):
        checks.append(_check("过滤前后图", "pass", "已生成过滤前后 QC 对比图。"))
    else:
        checks.append(_check("过滤前后图", "review", "未发现过滤前后 QC 对比图。"))

    if _has_file(names, "doublet"):
        checks.append(_check("双细胞证据", "pass", "已生成或记录 doublet 相关证据。"))
    else:
        checks.append(_check("双细胞证据", "review", "未发现 doublet score 图，需结合参数确认。"))

    return {"title": "QC 审批证据", "checks": checks}


def _clustering_evidence(summary, names):
    cluster_keys = [k for k in summary if k.startswith("n_clusters_")]
    cluster_counts = [_number(summary.get(k), 0) for k in cluster_keys]
    primary = summary.get("primary_resolution") or summary.get("best_resolution")
    checks = []

    if not cluster_counts:
        checks.append(_check("Cluster 数量", "review", "summary 中没有 cluster 数量，需确认聚类是否完成。"))
    elif max(cluster_counts) < 2:
        checks.append(_check("Cluster 数量", "warning", "仅得到 1 个 cluster，通常不足以支持后续注释/DEG。", int(max(cluster_counts))))
    elif max(cluster_counts) > 50:
        checks.append(_check("Cluster 数量", "review", "cluster 数量较多，需确认是否过度分裂。", int(max(cluster_counts))))
    else:
        checks.append(_check("Cluster 数量", "pass", "cluster 数量处于可复核范围。", int(max(cluster_counts))))

    if primary:
        checks.append(_check("主分辨率", "pass", "已记录主分辨率，可作为下游默认 cluster。", primary))
    else:
        checks.append(_check("主分辨率", "review", "未记录主分辨率，需确认下游使用哪个 cluster key。"))

    if _has_file(names, "cluster_umap_labeled"):
        checks.append(_check("带标签 UMAP", "pass", "已生成带 cluster 标签的 UMAP。"))
    else:
        checks.append(_check("带标签 UMAP", "review", "未发现带标签 cluster UMAP。"))

    if _has_file(names, "sankey") or _has_file(names, "resolution"):
        checks.append(_check("分辨率流向", "pass", "已有多分辨率拆分证据。"))
    else:
        checks.append(_check("分辨率流向", "review", "未发现多分辨率流向证据。"))

    if _has_file(names, "batch", "composition"):
        checks.append(_check("批次组成", "pass", "已有 cluster 批次组成图，可复核 batch 支配风险。"))
    else:
        checks.append(_check("批次组成", "review", "未发现 cluster 批次组成图，需结合 batch UMAP 复核。"))

    return {"title": "聚类审批证据", "checks": checks}


def _annotation_evidence(summary, names):
    counts = summary.get("celltype_counts") or {}
    total = sum(_number(v, 0) for v in counts.values())
    unknown = sum(_number(v, 0) for k, v in counts.items() if str(k).lower() == "unknown")
    unknown_pct = round(unknown / total * 100, 1) if total else None
    checks = []

    n_celltypes = _number(summary.get("n_celltypes"), 0)
    if n_celltypes <= 1:
        checks.append(_check("细胞类型数", "review", "注释细胞类型过少，需确认 marker 或 cluster 设置。", int(n_celltypes)))
    else:
        checks.append(_check("细胞类型数", "pass", "已得到多个细胞类型注释。", int(n_celltypes)))

    if unknown_pct is None:
        checks.append(_check("Unknown 比例", "review", "summary 中缺少 celltype_counts，无法判断 Unknown 比例。"))
    elif unknown_pct > 30:
        checks.append(_check("Unknown 比例", "warning", "Unknown 比例较高，需调整 marker 或重新分群。", f"{unknown_pct}%"))
    elif unknown_pct > 10:
        checks.append(_check("Unknown 比例", "review", "Unknown 比例中等，建议人工复核低置信 cluster。", f"{unknown_pct}%"))
    else:
        checks.append(_check("Unknown 比例", "pass", "Unknown 比例较低。", f"{unknown_pct}%"))

    confidence = _number(summary.get("mean_confidence"))
    margin = _number(summary.get("mean_score_margin"))
    if confidence is not None:
        if confidence < 0.3:
            checks.append(_check("平均置信度", "warning", "注释平均置信度偏低。", confidence))
        elif confidence < 0.5:
            checks.append(_check("平均置信度", "review", "注释平均置信度中等，建议结合 marker 图确认。", confidence))
        else:
            checks.append(_check("平均置信度", "pass", "注释平均置信度较好。", confidence))
    elif margin is not None:
        if margin < 0.05:
            checks.append(_check("Score margin", "warning", "最高分与次高分差距很小，注释不稳定。", margin))
        elif margin < 0.15:
            checks.append(_check("Score margin", "review", "score margin 中等，建议复核边界 cluster。", margin))
        else:
            checks.append(_check("Score margin", "pass", "score margin 支持当前注释。", margin))
    else:
        checks.append(_check("注释置信度", "review", "未启用或未记录注释置信度。"))

    if _has_file(names, "marker", "heatmap") or _has_file(names, "marker", "dotplot") or _has_file(names, "marker_expression"):
        checks.append(_check("Marker 证据图", "pass", "已生成 marker 表达或评分证据图。"))
    else:
        checks.append(_check("Marker 证据图", "review", "未发现 marker 验证图。"))

    if _has_file(names, "annotation_score_umap"):
        checks.append(_check("置信度 UMAP", "pass", "已生成注释置信度 UMAP。"))
    else:
        checks.append(_check("置信度 UMAP", "review", "未发现注释置信度 UMAP。"))

    return {"title": "注释审批证据", "checks": checks}


def _deg_evidence(summary, names):
    total_deg = _number(summary.get("total_deg_genes"), 0)
    n_groups = _number(summary.get("n_groups"), 0)
    checks = []

    if n_groups < 2:
        checks.append(_check("比较分组", "review", "DEG 分组数量不足，需确认 groupby/cluster 设置。", int(n_groups)))
    else:
        checks.append(_check("比较分组", "pass", "DEG 已覆盖多个分组。", int(n_groups)))

    if total_deg <= 0:
        checks.append(_check("DEG 数量", "warning", "未检测到差异基因，需确认统计方法、阈值或输入层。", int(total_deg)))
    elif total_deg < max(n_groups, 1):
        checks.append(_check("DEG 数量", "review", "差异基因数量较少，需结合阈值和表达图复核。", int(total_deg)))
    else:
        checks.append(_check("DEG 数量", "pass", "检测到差异基因，可作为 marker/注释支持证据。", int(total_deg)))

    if _has_file(names, "deg_marker_heatmap") or _has_file(names, "marker_heatmap"):
        checks.append(_check("Marker 热图", "pass", "已生成 cluster marker 热图。"))
    else:
        checks.append(_check("Marker 热图", "review", "未发现 DEG marker 热图。"))

    if _has_file(names, "volcano"):
        checks.append(_check("火山图", "pass", "已生成 DEG 火山图。"))
    else:
        checks.append(_check("火山图", "review", "未发现 DEG 火山图。"))

    if _has_file(names, "significant_counts") or _has_file(names, "counts_bar"):
        checks.append(_check("显著 DEG 数量图", "pass", "已有显著 DEG 数量概览。"))
    else:
        checks.append(_check("显著 DEG 数量图", "review", "未发现显著 DEG 数量概览图。"))

    return {"title": "DEG 审批证据", "checks": checks}


def build_review_evidence(module_name, summary, result_files=None):
    """Build a review evidence card for key single-cell modules."""
    if module_name not in SC_REVIEW_MODULES:
        return None
    summary = summary or {}
    names = _result_file_names(result_files or [])
    builders = {
        "qc": _qc_evidence,
        "clustering": _clustering_evidence,
        "annotation": _annotation_evidence,
        "deg": _deg_evidence,
    }
    card = builders[module_name](summary, names)
    checks = card.get("checks", [])
    card.update({
        "module_name": module_name,
        "status": _overall(checks),
    })
    return card


def build_result_interpretation(module_name, summary, result_files=None):
    """Create concise, evidence-bound interpretation for every result page."""
    summary = summary or {}
    evidence = build_review_evidence(module_name, summary, result_files)
    display = {
        'qc': '质控', 'normalize': '标准化', 'hvg': '高变基因', 'dimred': '降维',
        'batch_correct': '批次校正', 'clustering': '聚类', 'subcluster': '子簇精细分析',
        'annotation': '细胞注释', 'deg': '差异表达', 'trajectory': '轨迹分析',
        'proportion': '细胞比例', 'cell_communication': '细胞通讯',
    }.get(module_name, module_name)
    if evidence:
        cautions = [check['message'] for check in evidence['checks'] if check['status'] != 'pass']
        passes = [check['message'] for check in evidence['checks'] if check['status'] == 'pass']
        return {
            'title': '结果解读',
            'conclusion': f'已完成{display}；当前结果应结合下方证据图和参数进行生物学解释。',
            'evidence': passes[:2],
            'cautions': cautions[:2],
            'next_step': '优先查看标记为“需复核”或“警告”的项目，再决定是否调整参数或进入下一步。',
        }
    n_cells = summary.get('n_cells') or summary.get('n_samples')
    count_text = f'输出包含 {n_cells} 个分析对象。' if n_cells is not None else '已生成该步骤的结果文件。'
    return {
        'title': '结果解读',
        'conclusion': f'已完成{display}；{count_text}',
        'evidence': [],
        'cautions': ['该模块暂无自动判定阈值，请结合图表、参数和实验设计解读。'],
        'next_step': '确认结果符合预期后再进入下游分析。',
    }
