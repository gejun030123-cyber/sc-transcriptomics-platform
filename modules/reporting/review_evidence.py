"""Review evidence cards for analysis outputs."""

import os


SC_REVIEW_MODULES = {"qc", "batch_correct", "clustering", "annotation", "deg", "proportion"}


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
        # Confidence is an optional review layer.  Its absence should not
        # turn an otherwise valid marker annotation into a red/yellow warning.
        checks.append(_check("注释置信度", "pass", "未启用置信度评估（可选项），未将其作为失败条件。"))

    if _has_file(names, "marker", "heatmap") or _has_file(names, "marker", "dotplot") or _has_file(names, "marker_expression"):
        checks.append(_check("Marker 证据图", "pass", "已生成 marker 表达或评分证据图。"))
    else:
        checks.append(_check("Marker 证据图", "review", "未发现 marker 验证图。"))

    if _has_file(names, "annotation_score_umap"):
        checks.append(_check("置信度 UMAP", "pass", "已生成注释置信度 UMAP。"))
    elif confidence is None and margin is None:
        checks.append(_check("置信度 UMAP", "pass", "未启用置信度评估，因此不生成置信度 UMAP。"))
    else:
        checks.append(_check("置信度 UMAP", "review", "未发现注释置信度 UMAP。"))

    marker_selection = summary.get("marker_selection") or {}
    marker_warnings = list(marker_selection.get("warnings") or [])
    if marker_warnings:
        checks.append(_check(
            "Marker 选择提示", "review",
            "；".join(str(item) for item in marker_warnings[:2]),
        ))

    marker_validation = summary.get("marker_validation") or {}
    display_warning = str(marker_validation.get("display_warning") or "").strip()
    if display_warning:
        checks.append(_check("表达图数据尺度", "review", display_warning))

    maturity = summary.get("maturity_evidence") or {}
    if summary.get("marker_set") == "Organoid":
        time_key = str(maturity.get("time_key") or "未检测到时间元数据")
        state_counts = maturity.get("state_counts") or {}
        checks.append(_check(
            "类器官成熟度证据",
            "review" if not state_counts or "maturity_unresolved" in state_counts else "pass",
            "成熟度来自前体/成熟/增殖表达模块；时间列仅作为独立的描述性轴，不直接定义成熟标签。",
            f"时间列: {time_key}; 状态: {state_counts}",
        ))

    celltypist = summary.get("celltypist_reference") or {}
    if celltypist.get("enabled"):
        comparison_counts = celltypist.get("comparison_counts") or {}
        status_counts = celltypist.get("status_counts") or {}
        conflicts = int(comparison_counts.get("conflict", 0) or 0)
        unassigned = int(status_counts.get("unassigned", 0) or 0)
        low_confidence = int(status_counts.get("low_confidence", 0) or 0)
        multi_label = int(status_counts.get("multi_label", 0) or 0)
        if conflicts or unassigned or low_confidence or multi_label:
            checks.append(_check(
                "CellTypist 参考证据",
                "review",
                "CellTypist 仅作为交叉证据；存在冲突、未分配或低置信度细胞，不能覆盖 Marker 最终标签。",
                {
                    "model": celltypist.get("model_name"),
                    "conflict": conflicts,
                    "unassigned": unassigned,
                    "low_confidence": low_confidence,
                    "multi_label": multi_label,
                },
            ))
        else:
            checks.append(_check(
                "CellTypist 参考证据",
                "pass",
                "本次 CellTypist 参考标签与 Marker 结果未发现未分配或冲突信号；仍需结合组织背景复核。",
                celltypist.get("model_name"),
            ))

    quality = summary.get("quality_evidence") or {}
    suspect_doublet = _number(quality.get("n_suspect_doublet"), 0)
    ambient_signal = _number(quality.get("n_ambient_signal"), 0)
    if suspect_doublet:
        checks.append(_check(
            "Doublet 证据", "warning",
            "发现疑似双细胞/混合谱系信号；仅标记复核，不自动删除。",
            int(suspect_doublet),
        ))
    else:
        checks.append(_check("Doublet 证据", "pass", "未发现超过阈值的双谱系信号。"))
    if ambient_signal:
        checks.append(_check(
            "环境 RNA 证据", "review",
            "发现高普遍性异源 marker 信号；这是启发式提示，需结合 empty droplets 或去污染结果确认。",
            int(ambient_signal),
        ))
    else:
        checks.append(_check("环境 RNA 证据", "pass", "未发现超过阈值的环境 RNA 启发式信号。"))

    annotation_version = str(summary.get("annotation_version") or "").strip()
    if annotation_version:
        checks.append(_check(
            "注释版本",
            "pass",
            "已记录注释版本和备注，可与后续人工修订结果区分。",
            annotation_version,
        ))
    else:
        checks.append(_check("注释版本", "review", "未记录 annotation_version，无法区分不同注释批次。"))

    state_evidence = summary.get("cell_state_evidence") or {}
    state_counts = state_evidence.get("state_counts") or {}
    if state_counts:
        checks.append(_check(
            "细胞状态分离",
            "pass",
            "细胞状态与 cell_type 分开保存；一个细胞可同时带有多个状态标记。",
            state_counts,
        ))
    else:
        checks.append(_check("细胞状态分离", "review", "未生成独立 cell_state 证据。"))

    runtime_warnings = list(summary.get("runtime_warnings") or [])
    if runtime_warnings:
        checks.append(_check(
            "运行提示", "review", "；".join(str(item) for item in runtime_warnings[:2]),
        ))

    reviewed = _number(summary.get("n_clusters_reviewed"), 0)
    need_review = _number(summary.get("n_review_clusters"), 0)
    if _has_file(names, "annotation_cluster_review"):
        status = "review" if need_review else "pass"
        message = (
            "已输出逐簇注释复核表；其中存在需人工确认的 cluster。"
            if need_review else "已输出逐簇注释复核表，当前没有自动标记的高风险 cluster。"
        )
        checks.append(_check("逐簇复核表", status, message, f"{int(need_review)}/{int(reviewed)}"))
    else:
        checks.append(_check("逐簇复核表", "review", "未发现 annotation cluster review CSV。"))

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


def _proportion_evidence(summary, names):
    """Make the statistical unit explicit on cell-composition result pages."""
    checks = []
    sample_ready = bool(summary.get("sample_level_inference_ready"))
    sample_key = str(summary.get("sample_key") or "").strip()
    condition_key = str(summary.get("condition_key") or "").strip()
    n_tests = int(_number(summary.get("sample_level_n_tests"), 0) or 0)
    warnings = list(summary.get("sample_level_warnings") or [])

    if sample_ready and n_tests:
        checks.append(_check(
            "统计单位", "pass",
            "已按独立生物学样本计算细胞类型比例并进行 BH-FDR 校正。",
            f"{sample_key} × {condition_key}",
        ))
    else:
        checks.append(_check(
            "统计单位", "review",
            warnings[0] if warnings else "未形成满足重复要求的样本级比较；细胞计数仅作描述性展示。",
        ))

    if _has_file(names, "sample_level_cell_proportions"):
        checks.append(_check("样本级比例表", "pass", "已输出每个样本 × 细胞类型的比例（含零比例）。"))
    else:
        checks.append(_check("样本级比例表", "review", "未输出样本级比例表；请检查 sample_key/condition_key。"))

    if _has_file(names, "sample_level_proportion_tests"):
        checks.append(_check("FDR 统计表", "pass", "已输出每个细胞类型的样本级检验与 BH-FDR。", n_tests))
    else:
        checks.append(_check("FDR 统计表", "review", "没有可报告的样本级 p 值/FDR。"))

    return {"title": "细胞组成审批证据", "checks": checks}


def _batch_correct_evidence(summary, names):
    """Interpret only comparable, before/after batch-integration metrics.

    A lower batch ASW magnitude and same-batch-neighbour fraction support
    improved mixing; a higher neighbour entropy supports it as well.  These
    are evidence signals, not a substitute for checking cell-type retention or
    a study design in which biological condition and batch are not confounded.
    """
    comparison = summary.get("evaluation_comparison") or {}
    before = comparison.get("before") or {}
    after = comparison.get("after") or {}
    delta = comparison.get("delta") or {}
    sampling = comparison.get("sampling") or {}
    scope = str(comparison.get("comparison_scope") or "")
    checks = []

    if not comparison:
        checks.append(_check(
            "前后评价", "review",
            "未记录共享样本的校正前后指标；请启用“评估整合效果”后重新运行，或仅将当前结果作为探索性展示。",
        ))
        return {"title": "批次校正审批证据", "checks": checks}

    if before.get("error") or after.get("error"):
        checks.append(_check(
            "前后评价", "warning",
            "批次校正前后评价未完整完成，不能据此判断整合改善程度。",
        ))
    else:
        checks.append(_check(
            "前后评价", "pass",
            "已使用同一批次分层抽样细胞计算校正前后指标。",
            f"n={sampling.get('used_size', after.get('n_evaluation_cells', '—'))}",
        ))

    n_eval = _number(sampling.get("used_size", after.get("n_evaluation_cells")))
    if n_eval is None:
        checks.append(_check("评价样本量", "review", "未记录评价细胞数，需确认指标的稳定性。"))
    elif n_eval < 20:
        checks.append(_check("评价样本量", "review", "评价细胞数较少，前后差值可能不稳定。", int(n_eval)))
    else:
        checks.append(_check("评价样本量", "pass", "评价样本量已记录；抽样保留了每个非缺失批次。", int(n_eval)))

    abs_before = _number(before.get("abs_asw_batch"))
    abs_after = _number(after.get("abs_asw_batch"))
    abs_delta = _number(delta.get("abs_asw_batch"))
    if scope == "graph_only":
        checks.append(_check(
            "Batch ASW", "review",
            "BBKNN 主要校正邻居图而非 PCA 表示；Batch ASW 的表示空间差值不应单独作为成败依据。",
            f"{abs_before} → {abs_after}" if abs_before is not None and abs_after is not None else None,
        ))
    elif abs_delta is None:
        checks.append(_check("Batch ASW", "review", "缺少可比较的校正前后 Batch ASW。"))
    elif abs_delta <= -0.02:
        checks.append(_check(
            "Batch ASW", "pass",
            "|Batch ASW| 降低，支持批次在整合表示中的可分性下降。",
            f"{abs_before} → {abs_after}",
        ))
    elif abs_delta >= 0.02:
        checks.append(_check(
            "Batch ASW", "warning",
            "|Batch ASW| 升高，批次在整合表示中反而更易分离，需检查方法与参数。",
            f"{abs_before} → {abs_after}",
        ))
    else:
        checks.append(_check(
            "Batch ASW", "review",
            "|Batch ASW| 前后变化很小；需结合邻居混合和 UMAP 判断是否真的改善。",
            f"{abs_before} → {abs_after}",
        ))

    entropy_delta = _number(delta.get("mean_neighbor_batch_entropy"))
    same_batch_delta = _number(delta.get("mean_neighbor_same_batch_fraction"))
    entropy_before = _number(before.get("mean_neighbor_batch_entropy"))
    entropy_after = _number(after.get("mean_neighbor_batch_entropy"))
    same_before = _number(before.get("mean_neighbor_same_batch_fraction"))
    same_after = _number(after.get("mean_neighbor_same_batch_fraction"))
    if entropy_delta is not None and entropy_delta >= 0.02:
        checks.append(_check(
            "邻居批次混合", "pass",
            "邻居批次熵升高，支持局部批次混合改善。",
            f"entropy {entropy_before} → {entropy_after}",
        ))
    elif same_batch_delta is not None and same_batch_delta <= -0.02:
        checks.append(_check(
            "邻居批次混合", "pass",
            "同批次邻居比例下降，支持局部批次混合改善。",
            f"same-batch {same_before} → {same_after}",
        ))
    elif entropy_delta is not None and entropy_delta <= -0.02:
        checks.append(_check(
            "邻居批次混合", "warning",
            "邻居批次熵下降，局部批次混合可能变差。",
            f"entropy {entropy_before} → {entropy_after}",
        ))
    elif same_batch_delta is not None and same_batch_delta >= 0.02:
        checks.append(_check(
            "邻居批次混合", "warning",
            "同批次邻居比例上升，局部批次混合可能变差。",
            f"same-batch {same_before} → {same_after}",
        ))
    else:
        checks.append(_check(
            "邻居批次混合", "review",
            "缺少或未观察到明确的邻居混合改善；请查看前后 UMAP 与批次组成图。",
        ))

    bio_before = _number(before.get("asw_bio"))
    bio_after = _number(after.get("asw_bio"))
    bio_delta = _number(delta.get("asw_bio"))
    bio_key = after.get("bio_label_key") or before.get("bio_label_key")
    if bio_delta is None:
        checks.append(_check(
            "生物结构保留", "review",
            "未提供可比较的生物标签 ASW；无法自动排查过度校正。",
        ))
    elif bio_delta <= -0.05:
        checks.append(_check(
            "生物结构保留", "warning",
            "生物标签 ASW 明显下降，存在过度校正或生物结构被压平的风险。",
            f"{bio_key}: {bio_before} → {bio_after}",
        ))
    elif bio_delta < -0.02:
        checks.append(_check(
            "生物结构保留", "review",
            "生物标签 ASW 略降，建议复核关键细胞类型是否仍能分开。",
            f"{bio_key}: {bio_before} → {bio_after}",
        ))
    else:
        checks.append(_check(
            "生物结构保留", "pass",
            "生物标签 ASW 未见明显下降，支持主要生物结构得到保留。",
            f"{bio_key}: {bio_before} → {bio_after}",
        ))

    if _has_file(names, "batch_evaluation_pre_post"):
        checks.append(_check("前后指标表", "pass", "已生成可下载的前后指标与 delta 表。"))
    else:
        checks.append(_check("前后指标表", "review", "未发现前后指标表，请确认结果文件是否完整保存。"))

    warnings = comparison.get("warnings") or []
    if warnings:
        checks.append(_check("评价边界", "review", str(warnings[0])))

    return {"title": "批次校正审批证据", "checks": checks}


def build_review_evidence(module_name, summary, result_files=None):
    """Build a review evidence card for key single-cell modules."""
    if module_name not in SC_REVIEW_MODULES:
        return None
    summary = summary or {}
    names = _result_file_names(result_files or [])
    builders = {
        "qc": _qc_evidence,
        "batch_correct": _batch_correct_evidence,
        "clustering": _clustering_evidence,
        "annotation": _annotation_evidence,
        "deg": _deg_evidence,
        "proportion": _proportion_evidence,
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
        'sc_timecourse': '单细胞时序动态', 'proportion': '细胞比例', 'cell_communication': '细胞通讯',
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
    if module_name == 'sc_timecourse':
        inference_available = bool(summary.get('inference_available'))
        composition_inference = bool(summary.get('composition_inference_available'))
        gene_inference = bool(summary.get('gene_inference_available'))
        sample_key = str(summary.get('sample_key') or '').strip()
        n_timepoints = summary.get('n_timepoints')
        n_composition = summary.get('n_composition_significant', 0)
        n_genes = summary.get('n_gene_significant', 0)
        cautions = list(summary.get('warnings') or [])
        if not inference_available:
            cautions.insert(0, '当前没有满足设计与重复要求的样本级推断结果；仅可作描述性探索，不应报告推断性显著性。')
        composition_text = (
            f'细胞组成：{n_composition} 个细胞类型时间关联通过 FDR < 0.05。'
            if composition_inference else
            '细胞组成：未形成满足重复要求的样本级显著性检验，仅展示描述性时间趋势。'
        )
        gene_text = (
            f'基因动态：{n_genes} 个细胞类型特异候选基因关联通过 FDR < 0.05。'
            if gene_inference else
            '基因动态：未形成可报告的样本级显著性检验；若有曲线，仅用于描述性探索。'
        )
        unit_text = (
            f"统计单位为 '{sample_key}' 对应的 sample × celltype，而非单个细胞。"
            if sample_key else
            '未提供可用的生物学样本 ID；时间点内细胞汇总不构成独立统计重复。'
        )
        return {
            'title': '时序结果解读',
            'conclusion': f'已完成 {n_timepoints or "多个"} 时间点的样本级细胞组成与基因动态汇总。',
            'evidence': [
                composition_text,
                gene_text,
                unit_text,
            ],
            'cautions': cautions[:3],
            'next_step': '先核对样本重复、技术 batch 与时间的混杂，再结合时间点 UMAP、组成曲线和伪 bulk 热图选择候选机制；轨迹分析应作为补充验证。',
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
