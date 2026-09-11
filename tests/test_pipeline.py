# tests/test_pipeline.py
"""Tests for modules/__init__.py — MODULE_REGISTRY, PIPELINE_ORDER, validate_pipeline_order."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest


def test_qc_passes_omicverse_thresholds_with_supported_argument_name():
    """OmicVerse silently ignores ``thresh`` and expects historical ``tresh``."""
    import inspect
    from modules.qc import QCAnalysis

    source = inspect.getsource(QCAnalysis.run)
    assert "tresh={" in source
    assert "\n            thresh={" not in source


class TestModuleRegistry:
    """测试 MODULE_REGISTRY 注册表。"""

    def test_all_modules_registered(self):
        """所有注册模块都在注册表中。"""
        from modules import MODULE_REGISTRY
        assert len(MODULE_REGISTRY) == 32

    def test_registry_values_are_classes(self):
        """注册表的值都是类（有 run 方法）。"""
        from modules import MODULE_REGISTRY
        for name, cls in MODULE_REGISTRY.items():
            assert hasattr(cls, 'run'), f"{name} 缺少 run 方法"

    def test_expected_modules_exist(self):
        """关键模块存在。"""
        from modules import MODULE_REGISTRY
        expected = ['qc', 'normalize', 'hvg', 'dimred', 'clustering',
                    'functional_state', 'sc_timecourse', 'sc_cell_deg', 'sc_cell_go', 'sc_pseudobulk_deg', 'sc_csv_export', 'neighborhood_da',
                    'sc_batch_import', 'virtual_ko', 'bulk_qc', 'bulk_normalize', 'bulk_deg', 'convert_10x']
        for mod in expected:
            assert mod in MODULE_REGISTRY, f"缺少模块: {mod}"


class TestPipelineOrder:
    """测试 PIPELINE_ORDER。"""

    def test_order_contains_all_non_convert_modules(self):
        """流水线顺序包含除专用数据导入入口外的所有模块。"""
        from modules import PIPELINE_ORDER, MODULE_REGISTRY
        non_convert = {k for k in MODULE_REGISTRY if k not in {'convert_10x', 'sc_batch_import'}}
        order_set = set(PIPELINE_ORDER)
        assert non_convert == order_set, f"差异: {non_convert.symmetric_difference(order_set)}"

    def test_no_duplicates(self):
        """流水线顺序无重复。"""
        from modules import PIPELINE_ORDER
        assert len(PIPELINE_ORDER) == len(set(PIPELINE_ORDER))

    def test_sc_before_bulk(self):
        """单细胞模块在 bulk 模块之前。"""
        from modules import PIPELINE_ORDER
        sc_last = max(PIPELINE_ORDER.index(m) for m in PIPELINE_ORDER if not m.startswith('bulk'))
        bulk_first = min(PIPELINE_ORDER.index(m) for m in PIPELINE_ORDER if m.startswith('bulk'))
        assert sc_last < bulk_first


class TestValidatePipelineOrder:
    """测试 validate_pipeline_order 函数。"""

    def test_full_valid_order(self):
        """完整流水线顺序 → 通过。"""
        from modules import validate_pipeline_order, PIPELINE_ORDER
        is_valid, errors = validate_pipeline_order(PIPELINE_ORDER)
        assert is_valid is True
        assert errors == []

    def test_empty_list(self):
        """空列表 → 通过（无依赖违反）。"""
        from modules import validate_pipeline_order
        is_valid, errors = validate_pipeline_order([])
        assert is_valid is True

    def test_single_module_no_deps(self):
        """单个无依赖模块 → 通过。"""
        from modules import validate_pipeline_order
        is_valid, errors = validate_pipeline_order(['qc'])
        assert is_valid is True

    def test_reverse_order_fails(self):
        """反序（违反依赖）→ 失败。"""
        from modules import validate_pipeline_order, PIPELINE_ORDER
        reversed_order = list(reversed(PIPELINE_ORDER))
        is_valid, errors = validate_pipeline_order(reversed_order)
        assert is_valid is False
        assert len(errors) > 0

    def test_normalize_before_qc_fails(self):
        """normalize 在 qc 之前 → 失败。"""
        from modules import validate_pipeline_order
        is_valid, errors = validate_pipeline_order(['normalize', 'qc'])
        assert is_valid is False
        assert any('normalize' in e for e in errors)

    def test_qc_before_normalize_passes(self):
        """qc 在 normalize 之前 → 通过。"""
        from modules import validate_pipeline_order
        is_valid, errors = validate_pipeline_order(['qc', 'normalize'])
        assert is_valid is True

    def test_partial_subset_valid(self):
        """部分模块子集，顺序正确 → 通过。"""
        from modules import validate_pipeline_order
        is_valid, errors = validate_pipeline_order(['qc', 'normalize', 'hvg', 'dimred'])
        assert is_valid is True

    def test_partial_subset_invalid(self):
        """部分模块子集，顺序错误 → 失败。"""
        from modules import validate_pipeline_order
        is_valid, errors = validate_pipeline_order(['dimred', 'hvg'])
        assert is_valid is False

    def test_bulk_deps(self):
        """bulk 模块依赖测试。"""
        from modules import validate_pipeline_order
        # bulk_deg 依赖 bulk_normalize，bulk_normalize 依赖 bulk_qc
        is_valid, errors = validate_pipeline_order(['bulk_qc', 'bulk_normalize', 'bulk_deg'])
        assert is_valid is True

    def test_bulk_deg_before_normalize_fails(self):
        """bulk_deg 在 bulk_normalize 之前 → 失败。"""
        from modules import validate_pipeline_order
        is_valid, errors = validate_pipeline_order(['bulk_deg', 'bulk_normalize'])
        assert is_valid is False


class TestPipelineTypeValidation:
    """测试流水线类型校验（SC vs Bulk）。"""

    def test_sc_pipeline_valid(self):
        """单细胞最小链路 → 有效。"""
        from modules import validate_pipeline_order
        is_valid, errors = validate_pipeline_order(['qc', 'normalize', 'hvg'])
        assert is_valid is True

    def test_sc_pipeline_reverse_fails(self):
        """单细胞反序 → 失败。"""
        from modules import validate_pipeline_order
        is_valid, errors = validate_pipeline_order(['hvg', 'normalize', 'qc'])
        assert is_valid is False

    def test_bulk_pipeline_valid(self):
        """Bulk 最小链路 → 有效。"""
        from modules import validate_pipeline_order
        is_valid, errors = validate_pipeline_order(['bulk_qc', 'bulk_normalize', 'bulk_pca'])
        assert is_valid is True

    def test_bulk_pipeline_reverse_fails(self):
        """Bulk 反序 → 失败。"""
        from modules import validate_pipeline_order
        is_valid, errors = validate_pipeline_order(['bulk_pca', 'bulk_normalize', 'bulk_qc'])
        assert is_valid is False

    def test_sc_module_in_bulk_pipeline_rejected(self):
        """单细胞模块出现在 bulk pipeline 中 → 被 SC_MODULE_NAMES 拒绝。"""
        from modules import SC_MODULE_NAMES, BULK_MODULE_NAMES
        # qc 是单细胞模块，不属于 BULK_MODULE_NAMES
        assert 'qc' in SC_MODULE_NAMES
        assert 'qc' not in BULK_MODULE_NAMES

    def test_bulk_module_in_sc_pipeline_rejected(self):
        """Bulk 模块出现在单细胞 pipeline 中 → 被 SC_MODULE_NAMES 拒绝。"""
        from modules import SC_MODULE_NAMES, BULK_MODULE_NAMES
        # bulk_qc 是 bulk 模块，不属于 SC_MODULE_NAMES
        assert 'bulk_qc' in BULK_MODULE_NAMES
        assert 'bulk_qc' not in SC_MODULE_NAMES


def test_pipeline_worker_persists_successful_module_summary(test_project, monkeypatch):
    """A successful pipeline module must not fail while saving its summary."""
    import json
    import os
    import modules
    import worker
    from config import Config
    from models import AnalysisTask, PipelineRun

    class SuccessfulModule:
        def __init__(self, project_dir, params, progress_callback):
            self.progress_callback = progress_callback

        def run(self, input_path):
            self.progress_callback(50, "fake module running")
            return {
                "output_adata": input_path,
                "result_files": [],
                "summary": {"status": "ok"},
            }

    monkeypatch.setitem(modules.MODULE_REGISTRY, "fake_pipeline_module", SuccessfulModule)
    monkeypatch.setattr(worker, "register_task_outputs", lambda *args, **kwargs: [])
    monkeypatch.setattr(worker, "_write_pipeline_artifacts", lambda *args, **kwargs: None)
    input_path = os.path.join(Config.project_dir(test_project), "uploads", "input.h5ad")
    run = PipelineRun(
        project_id=test_project, name="worker summary test", analysis_type="sc",
        input_path=input_path, modules_json=json.dumps(["fake_pipeline_module"]),
        params_json=json.dumps({"fake_pipeline_module": {}}),
    )
    run.save()

    worker._run_pipeline_run(
        run.id, test_project, ["fake_pipeline_module"],
        {"fake_pipeline_module": {}}, Config.project_dir(test_project), input_path,
    )

    assert PipelineRun.get_by_id(run.id).status == "completed"
    task_id = json.loads(PipelineRun.get_by_id(run.id).task_ids_json)[0]
    task = AnalysisTask.get_by_id(task_id)
    assert task.status == "completed"
    assert json.loads(task.result_json) == {"status": "ok"}


def test_saved_pipeline_template_runs_from_server_side_definition(test_project, monkeypatch):
    """A one-click template must not depend on browser-copied modules/params."""
    import json
    from app import create_app
    from config import Config
    from models import PipelineRun
    import worker

    monkeypatch.setattr(Config, "PLATFORM_ACCESS_PASSWORD", "")
    input_path = os.path.join(Config.project_dir(test_project), "uploads", "input.h5ad")
    with open(input_path, "wb") as handle:
        handle.write(b"placeholder")

    submitted = []
    monkeypatch.setattr(worker, "submit_pipeline_run", lambda **kwargs: submitted.append(kwargs))
    client = create_app().test_client()
    saved = client.post("/api/presets", json={
        "name": "标准单细胞基础流程",
        "description": "QC 到标准化",
        "scope": "project",
        "project_id": test_project,
        "analysis_type": "sc",
        "pipeline": {"modules": ["qc", "normalize"]},
        # Older saved pipeline templates used this flat parameter shape.
        "params": {"nUMIs": 777, "normalize": {"target_sum": 20000}},
    })
    assert saved.status_code == 200
    template_id = saved.get_json()["id"]

    listed = client.get(f"/api/presets?project_id={test_project}&type=sc")
    assert listed.status_code == 200
    template_summary = next(
        item for item in listed.get_json()["presets"] if item["id"] == template_id
    )
    assert template_summary["has_pipeline"] is True
    assert template_summary["pipeline"]["modules"] == ["qc", "normalize"]
    assert "params" not in template_summary

    launched = client.post(
        f"/api/projects/{test_project}/pipeline-templates/{template_id}/run",
        json={
            "input_path": input_path,
            # These untrusted browser fields must be ignored by the endpoint.
            "modules": ["bulk_qc"],
            "params": {"nUMIs": 1},
        },
    )
    assert launched.status_code == 201
    body = launched.get_json()
    assert body["modules"] == ["qc", "normalize"]
    assert body["template_id"] == template_id
    assert len(submitted) == 1
    assert submitted[0]["modules"] == ["qc", "normalize"]

    run = PipelineRun.get_by_id(body["id"])
    assert run is not None
    saved_params = json.loads(run.params_json)
    assert saved_params["qc"]["nUMIs"] == 777
    assert saved_params["normalize"]["target_sum"] == 20000


def test_pipeline_input_api_lists_safe_project_files(test_project, monkeypatch):
    """The launcher receives absolute, project-owned paths instead of names."""
    from app import create_app
    from config import Config

    monkeypatch.setattr(Config, "PLATFORM_ACCESS_PASSWORD", "")
    uploads = os.path.join(Config.project_dir(test_project), "uploads")
    h5ad_path = os.path.join(uploads, "cells.h5ad")
    csv_path = os.path.join(uploads, "counts.csv")
    text_path = os.path.join(uploads, "notes.md")
    for path in (h5ad_path, csv_path, text_path):
        with open(path, "wb") as handle:
            handle.write(b"placeholder")

    client = create_app().test_client()
    sc_response = client.get(f"/api/projects/{test_project}/pipeline-inputs?type=sc")
    bulk_response = client.get(f"/api/projects/{test_project}/pipeline-inputs?type=bulk")

    assert sc_response.status_code == 200
    assert bulk_response.status_code == 200
    sc_files = sc_response.get_json()["files"]
    bulk_files = bulk_response.get_json()["files"]
    assert [item["path"] for item in sc_files] == [h5ad_path]
    assert {item["path"] for item in bulk_files} == {h5ad_path, csv_path}
    assert all(item["path"].startswith(Config.project_dir(test_project)) for item in bulk_files)


def test_pipeline_launcher_pages_offer_template_creation_and_safe_input_picker(test_project, monkeypatch):
    """Both list pages render the controls needed to populate empty launchers."""
    from app import create_app
    from config import Config

    monkeypatch.setattr(Config, "PLATFORM_ACCESS_PASSWORD", "")
    client = create_app().test_client()
    sc_page = client.get(f"/projects/{test_project}/sc-analysis")
    bulk_page = client.get(f"/projects/{test_project}/bulk-analysis")

    assert sc_page.status_code == 200
    assert bulk_page.status_code == 200
    sc_html = sc_page.get_data(as_text=True)
    bulk_html = bulk_page.get_data(as_text=True)
    assert "新建流程模板" in sc_html
    assert "pipeline-inputs?type=sc" in sc_html
    assert "新建流程模板" in bulk_html
    assert "pipeline-inputs?type=bulk" in bulk_html


def test_pipeline_template_rejects_invalid_step_order_at_save_time(test_project, monkeypatch):
    """The new-template dialog receives a useful error before a run is started."""
    from app import create_app
    from config import Config

    monkeypatch.setattr(Config, "PLATFORM_ACCESS_PASSWORD", "")
    response = create_app().test_client().post("/api/presets", json={
        "name": "错误步骤顺序",
        "scope": "project",
        "project_id": test_project,
        "analysis_type": "sc",
        "pipeline": {"modules": ["normalize", "qc"]},
    })

    assert response.status_code == 400
    assert "前置依赖" in response.get_json()["error"]
