import json
from types import SimpleNamespace


def test_write_plot_gallery_and_project_report(tmp_path, monkeypatch):
    from modules.reporting import project_report

    monkeypatch.setattr(project_report, "_plotly_js", lambda: "window.Plotly={newPlot:function(){}};")

    project_dir = tmp_path / "project"
    plots_dir = project_dir / "plots"
    plots_dir.mkdir(parents=True)
    plot_path = plots_dir / "qc_overview.json"
    plot_path.write_text(
        json.dumps({
            "data": [{"type": "bar", "x": ["before", "after"], "y": [50, 45]}],
            "layout": {"title": "QC overview"},
        }),
        encoding="utf-8",
    )
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
                file_type="plotly_json",
                category="qc",
                label="QC overview",
                file_path=str(plot_path),
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
    assert "Plotly.newPlot" in gallery
    assert "全部导出 PNG" in gallery
    assert "exportAllPlots" in gallery
    assert "导出 SVG" in gallery

    report = report_path.read_text(encoding="utf-8")
    assert "# 项目结果报告：Demo Project" in report
    assert "### 质控" in report
    assert "`cells_after`：45" in report
    assert "审批证据：" in report
    assert "细胞保留率" in report
    assert "plot_gallery.html" in report
