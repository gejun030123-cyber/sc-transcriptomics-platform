"""AI 主工作台页面。

The workspace deliberately remains a thin page route: chat execution continues
to use ``routes.chat`` and the existing analysis APIs.  This keeps the first
increment safe while giving the AI experience a dedicated, extensible shell.
"""
import os

from flask import Blueprint, render_template, redirect, url_for, flash

from config import Config
from models import Project
from modules.schemas import MODULE_DISPLAY_MAP, STATUS_MAP


workspace_bp = Blueprint("workspace", __name__)


WORKSPACE_ASSAYS = (
    {
        "key": "scrna",
        "title": "单细胞 RNA",
        "short_title": "scRNA-seq",
        "description": "质控、聚类、注释、差异表达与轨迹分析",
        "status": "available",
        "href": "{pid}/sc-analysis",
        "prompt": "检查当前单细胞项目状态，并告诉我下一步建议",
        "accent": "indigo",
    },
    {
        "key": "bulk_rna",
        "title": "Bulk RNA",
        "short_title": "Bulk RNA-seq",
        "description": "QC、标准化、差异表达、富集与时序分析",
        "status": "available",
        "href": "{pid}/bulk-analysis",
        "prompt": "检查当前 Bulk RNA 数据和实验设计，并推荐下一步分析",
        "accent": "rose",
    },
    {
        "key": "wes",
        "title": "WES 外显子组",
        "short_title": "WES",
        "description": "germline、肿瘤-正常与 tumor-only 的准备、运行追踪和结果",
        "status": "available",
        "href": "{pid}/wes",
        "prompt": "我想进行 WES 分析，请先告诉我需要准备哪些文件、样本表和参考资源",
        "accent": "amber",
    },
    {
        "key": "bulk_atac",
        "title": "Bulk ATAC",
        "short_title": "ATAC-seq",
        "description": "预留：QC、peak、差异可及性与 motif 工作流",
        "status": "planned",
        "href": "",
        "prompt": "我想进行 Bulk ATAC-seq 分析，请先告诉我需要准备哪些 FASTQ、重复和参考资源",
        "accent": "teal",
    },
)


def _workspace_uploads(pid):
    uploads_dir = Config.uploads_dir(pid)
    if not os.path.isdir(uploads_dir):
        return []
    files = []
    for name in sorted(os.listdir(uploads_dir)):
        path = os.path.join(uploads_dir, name)
        if os.path.isfile(path):
            files.append({
                "name": name,
                "size_mb": round(os.path.getsize(path) / (1024 ** 2), 1),
                "extension": os.path.splitext(name)[1].lower(),
            })
    return files


@workspace_bp.route("/<pid>/workspace")
def index(pid):
    project = Project.get_by_id(pid)
    if not project:
        flash("项目未找到", "danger")
        return redirect(url_for("main.index"))

    tasks = project.get_tasks()
    # Keep the page payload intentionally small.  Full task details are loaded
    # through the existing result/status APIs when a user opens a card.
    recent_tasks = []
    for task in tasks[:12]:
        recent_tasks.append({
            "id": task.id,
            "module_name": task.module_name,
            "module_label": MODULE_DISPLAY_MAP.get(task.module_name, task.module_name),
            "status": task.status,
            "status_label": STATUS_MAP.get(task.status, task.status),
            "progress": task.progress,
            "progress_message": task.progress_message or "",
        })

    assays = []
    for assay in WORKSPACE_ASSAYS:
        item = dict(assay)
        item["href"] = assay["href"].format(pid=pid) if assay["href"] else ""
        assays.append(item)

    return render_template(
        "workspace.html",
        project=project,
        all_projects=Project.get_all(),
        uploads=_workspace_uploads(pid),
        recent_tasks=recent_tasks,
        assays=assays,
        active_task_count=sum(t.status in {"pending", "running"} for t in tasks),
    )
