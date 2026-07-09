"""Branch comparison and acceptance report helpers."""

import json
import os
from datetime import datetime, timezone


def _safe_json_loads(value, default):
    if not value:
        return default
    try:
        return json.loads(value)
    except (TypeError, ValueError, json.JSONDecodeError):
        return default


def _score_summary(scores):
    if not scores:
        return {
            "best_score": None,
            "best_score_id": None,
            "recommendation": "",
            "target_label": "",
            "confidence": "",
            "best_cluster": None,
            "score_detail": {},
        }
    best = scores[0]
    score_detail = _safe_json_loads(best.score_json, {})
    return {
        "best_score": best.total_score,
        "best_score_id": best.id,
        "recommendation": best.recommendation,
        "target_label": best.target_label,
        "confidence": score_detail.get("confidence", ""),
        "best_cluster": score_detail.get("best_cluster"),
        "score_detail": score_detail,
    }


def _flatten_params(params):
    flat = {}
    if not isinstance(params, dict):
        return flat
    for key, value in params.items():
        if isinstance(value, dict):
            for sub_key, sub_val in value.items():
                flat[f"{key}.{sub_key}"] = sub_val
        else:
            flat[key] = value
    return flat


def _param_diff(rows):
    flattened = {row["branch_id"]: _flatten_params(row.get("params", {})) for row in rows}
    all_keys = sorted({key for params in flattened.values() for key in params})
    diff = []
    for key in all_keys:
        values = {branch_id: params.get(key) for branch_id, params in flattened.items()}
        unique_values = {json.dumps(v, sort_keys=True, ensure_ascii=False) for v in values.values()}
        if len(unique_values) > 1:
            diff.append({"param": key, "values": values})
    return diff


def build_branch_comparison(project_id, branches, scores_by_branch):
    """Build a comparable view over branch candidates and their scores."""
    rows = []
    for branch in branches:
        if getattr(branch, "deleted", 0):
            continue
        scores = scores_by_branch.get(branch.id, [])
        score_info = _score_summary(scores)
        try:
            params = json.loads(branch.params_json) if branch.params_json else {}
        except (TypeError, ValueError, json.JSONDecodeError):
            params = {}
        rows.append({
            "branch_id": branch.id,
            "branch_name": branch.branch_name,
            "status": branch.status,
            "accepted": bool(branch.accepted),
            "purpose": branch.purpose,
            "params": params,
            "output_adata_path": branch.output_adata_path,
            "created_at": branch.created_at,
            "finished_at": branch.finished_at,
            **score_info,
        })

    rows.sort(
        key=lambda row: (
            1 if row["accepted"] else 0,
            row["best_score"] if row["best_score"] is not None else -1,
            row["finished_at"] or row["created_at"] or "",
        ),
        reverse=True,
    )
    best = next((row for row in rows if row["status"] == "completed" and row["best_score"] is not None), None)
    accepted = next((row for row in rows if row["accepted"]), None)
    return {
        "project_id": project_id,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "best_branch_id": best["branch_id"] if best else None,
        "accepted_branch_id": accepted["branch_id"] if accepted else None,
        "branches": rows,
        "param_diff": _param_diff(rows),
    }


def _branch_report_paths(project_dir, branch_id):
    report_dir = os.path.join(project_dir, "results", "agent")
    os.makedirs(report_dir, exist_ok=True)
    return (
        os.path.join(report_dir, f"accepted_branch_{branch_id}.json"),
        os.path.join(report_dir, f"accepted_branch_{branch_id}.md"),
    )


def write_acceptance_report(project_dir, branch, scores, comparison=None):
    """Persist JSON and Markdown records explaining an accepted branch."""
    score_info = _score_summary(scores)
    comparison = comparison or {}
    json_path, md_path = _branch_report_paths(project_dir, branch.id)
    payload = {
        "accepted_at": datetime.now(timezone.utc).isoformat(),
        "project_id": branch.project_id,
        "branch": branch.to_dict(),
        "score_summary": score_info,
        "comparison": comparison,
    }
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)

    lines = [
        f"# 已采纳候选分支：{branch.branch_name}",
        "",
        f"- Branch ID：`{branch.id}`",
        f"- 项目 ID：`{branch.project_id}`",
        f"- 采纳时间：{payload['accepted_at']}",
        f"- 输出 h5ad：`{branch.output_adata_path or ''}`",
        f"- 目标标签：{score_info.get('target_label') or '未记录'}",
        f"- 最佳 cluster：{score_info.get('best_cluster') if score_info.get('best_cluster') is not None else '未记录'}",
        f"- 综合评分：{score_info.get('best_score') if score_info.get('best_score') is not None else '未记录'}",
        f"- 推荐理由：{score_info.get('recommendation') or '未记录'}",
        "",
        "## 候选对比",
        "",
    ]
    rows = comparison.get("branches", [])
    if rows:
        lines.append("| 分支 | 状态 | 评分 | 推荐 |")
        lines.append("| --- | --- | --- | --- |")
        for row in rows:
            lines.append(
                f"| {row.get('branch_name')} | {row.get('status')} | "
                f"{row.get('best_score') if row.get('best_score') is not None else ''} | "
                f"{row.get('recommendation') or ''} |"
            )
    else:
        lines.append("暂无可用候选对比。")
    lines.extend(["", "## 参数差异", ""])
    diffs = comparison.get("param_diff", [])
    if diffs:
        for item in diffs:
            lines.append(f"- `{item['param']}`：{json.dumps(item['values'], ensure_ascii=False)}")
    else:
        lines.append("未检测到候选间参数差异，或仅有一个候选。")
    lines.append("")
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return {"json_path": json_path, "markdown_path": md_path}
