import base64
import json
from types import SimpleNamespace


def test_write_plot_gallery_and_project_report(tmp_path, monkeypatch):
    from modules.reporting import project_report

    project_dir = tmp_path / "project"
    plots_dir = project_dir / "plots"
    plots_dir.mkdir(parents=True)
    plot_path = plots_dir / "qc_overview.png"
    plot_path.write_bytes(base64.b64decode(
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII="
    ))
    svg_path = plots_dir / "qc_overview.svg"
    svg_path.write_text("<svg xmlns='http://www.w3.org/2000/svg' width='1' height='1'></svg>", encoding="utf-8")
    csv_path = project_dir / "results" / "qc.csv"
    csv_path.parent.mkdir()
    csv_path.write_text("metric,value\ncells,45\n", encoding="utf-8")

    project = SimpleNamespace(id="proj-1", name="Demo Project")
    task = SimpleNamespace(
        id="task-1",
        module_name="qc",
        status="completed",
        finished_at="2026-07-08T12:00:00",
        output_adata_path=str(project_dir / "intermediate" / "qc_output.h5ad"),
        result_json=json.dumps({"cells_before": 50, "cells_after": 45}),
    )
    files_by_task = {
        "task-1": [
            SimpleNamespace(
                id="plot-1",
                task_id="task-1",
                project_id="proj-1",
                file_type="png",
                category="qc",
                label="QC overview",
                file_path=str(plot_path),
            ),
            SimpleNamespace(
                id="plot-2",
                task_id="task-1",
                project_id="proj-1",
                file_type="svg",
                category="qc",
                label="QC overview",
                file_path=str(svg_path),
            ),
            SimpleNamespace(
                id="csv-1",
                task_id="task-1",
                project_id="proj-1",
                file_type="csv",
                category="table",
                label="QC metrics",
                file_path=str(csv_path),
            ),
        ]
    }

    info = project_report.ensure_project_report(
        project=project,
        project_dir=str(project_dir),
        tasks=[task],
        files_by_task=files_by_task,
        module_display_map={"qc": "质控"},
    )

    gallery_path = project_dir / "results" / "html_plots" / "plot_gallery.html"
    report_path = project_dir / "results" / "reports" / "project_report.md"
    assert info["gallery_path"] == str(gallery_path)
    assert info["report_path"] == str(report_path)
    assert gallery_path.exists()
    assert report_path.exists()

    gallery = gallery_path.read_text(encoding="utf-8")
    assert "QC overview" in gallery
    assert "Plotly.newPlot" not in gallery
    assert "下载 PNG" in gallery
    assert "下载 SVG" in gallery

    report = report_path.read_text(encoding="utf-8")
    assert "# 项目结果报告：Demo Project" in report
    assert "### 质控" in report
    assert "`cells_after`：45" in report
    assert "审批证据：" in report
    assert "细胞保留率" in report
    assert "plot_gallery.html" in report
