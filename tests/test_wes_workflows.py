import json
import gzip
import io
import os

import pytest

from config import Config
from modules.workflows.nextflow import NextflowExecutor, WorkflowNotConfiguredError
from modules.workflows.preflight import load_manifest_file, validate_manifest
from modules.workflows.registry import list_workflows
from modules.workflows.references import (
    get_capture_kit_profile,
    launch_reference_readiness,
    list_capture_kit_profiles,
    list_reference_assets,
    register_capture_kit_profile,
    register_reference_asset,
)
from modules.workflows.storage import get_manifest, list_manifests, register_manifest
from modules.workflows.assets import list_data_assets, register_data_asset
from modules.workflows.artifacts import collect_workflow_artifacts, list_workflow_artifacts
from modules.workflows.content import inspect_sample_files
from modules.workflows.sarek import render_samplesheet, write_launch_bundle
from modules.workflows.validation import (
    evaluate_metrics,
    parse_benchmark_summary,
    write_evaluation_bed,
)


def _somatic_manifest(project_dir):
    normal_r1 = os.path.join(project_dir, "normal_R1.fastq.gz")
    normal_r2 = os.path.join(project_dir, "normal_R2.fastq.gz")
    tumor_r1 = os.path.join(project_dir, "tumor_R1.fastq.gz")
    tumor_r2 = os.path.join(project_dir, "tumor_R2.fastq.gz")
    capture_bed = os.path.join(project_dir, "capture.bed")
    for path in (normal_r1, normal_r2, tumor_r1, tumor_r2):
        with open(path, "wb") as handle:
            handle.write(b"placeholder")
    with open(capture_bed, "w", encoding="utf-8") as handle:
        handle.write("chr1\t0\t100\n")
    return {
        "reference_bundle_id": "grch38-broad-v1",
        "capture_bed_id": "twist-exome-v2",
        "capture_bed_path": capture_bed,
        "samples": [
            {
                "sample_id": "N01",
                "patient_id": "P01",
                "role": "normal",
                "input_type": "fastq",
                "fastq_1": normal_r1,
                "fastq_2": normal_r2,
                "capture_kit": "twist",
            },
            {
                "sample_id": "T01",
                "patient_id": "P01",
                "role": "tumor",
                "input_type": "fastq",
                "fastq_1": tumor_r1,
                "fastq_2": tumor_r2,
                "matched_normal_id": "N01",
                "capture_kit": "twist",
            },
        ],
    }


def test_wes_registry_is_separate_and_declarative():
    workflows = list_workflows("wes")
    assert {item["key"] for item in workflows} == {
        "wes_germline", "wes_somatic", "wes_annotate_only"
    }
    assert all(item["executor"] == "nextflow" for item in workflows)
    assert all(item["status"] == "validation" for item in workflows)


def test_somatic_manifest_preflight_accepts_matched_pair(test_project):
    manifest = _somatic_manifest(Config.project_dir(test_project))
    result = validate_manifest(
        manifest,
        project_dir=Config.project_dir(test_project),
        workflow_key="wes_somatic",
    )
    assert result.valid is True
    assert result.errors == []
    assert result.checks["sample_count"] == 2
    assert result.checks["roles"] == ["normal", "tumor"]


def test_wes_calling_preflight_requires_capture_bed_path(test_project):
    manifest = _somatic_manifest(Config.project_dir(test_project))
    manifest.pop("capture_bed_path")

    result = validate_manifest(
        manifest,
        project_dir=Config.project_dir(test_project),
        workflow_key="wes_somatic",
    )

    assert result.valid is False
    assert "WES calling workflow 缺少 capture_bed_path" in result.errors


def test_preflight_rejects_missing_pair_and_unknown_normal(test_project):
    project_dir = Config.project_dir(test_project)
    manifest = {
        "reference_bundle_id": "grch38-v1",
        "capture_bed_id": "exome-v1",
        "samples": [{
            "sample_id": "T01",
            "patient_id": "P01",
            "role": "tumor",
            "input_type": "fastq",
            "fastq_1": os.path.join(project_dir, "tumor_R1.fastq.gz"),
            "matched_normal_id": "N99",
        }],
    }
    result = validate_manifest(manifest, project_dir=project_dir,
                               workflow_key="wes_somatic")
    assert result.valid is False
    assert any("fastq_2 缺失" in error for error in result.errors)
    assert any("matched_normal_id 不存在" in error for error in result.errors)


def test_tumor_only_requires_reason_and_explicit_confirmation(test_project):
    project_dir = Config.project_dir(test_project)
    manifest = _somatic_manifest(project_dir)
    manifest["samples"] = [manifest["samples"][1]]
    manifest["samples"][0].pop("matched_normal_id")

    rejected = validate_manifest(
        manifest, project_dir=project_dir, workflow_key="wes_somatic"
    )
    assert rejected.valid is False
    assert any("tumor_only_confirmed=true" in error for error in rejected.errors)
    assert any("tumor_only_reason" in error for error in rejected.errors)

    manifest["tumor_only_confirmed"] = True
    manifest["tumor_only_reason"] = "matched normal 不可用，仅用于受限科研分析"
    accepted = validate_manifest(
        manifest, project_dir=project_dir, workflow_key="wes_somatic"
    )
    assert accepted.valid is True
    assert any("tumor-only" in warning for warning in accepted.warnings)


def test_preflight_rejects_mixed_entry_and_patient_mismatch(test_project):
    project_dir = Config.project_dir(test_project)
    manifest = {
        "reference_bundle_id": "grch38-v1",
        "capture_bed_id": "exome-v1",
        "samples": [
            {
                "sample_id": "N01", "patient_id": "P02", "role": "normal",
                "input_type": "fastq", "fastq_1": "n1", "fastq_2": "n2",
                "capture_kit": "kit",
            },
            {
                "sample_id": "T01", "patient_id": "P01", "role": "tumor",
                "fastq_1": "t1", "fastq_2": "t2", "bam": "t.bam",
                "matched_normal_id": "N01", "capture_kit": "kit",
            },
        ],
    }
    result = validate_manifest(manifest, project_dir=project_dir,
                               workflow_key="wes_somatic", require_files=False)
    assert result.valid is False
    assert any("无法识别或不支持: mixed" in error for error in result.errors)
    assert any("patient_id 不一致" in error for error in result.errors)


def test_manifest_loader_supports_json_and_tsv(tmp_path):
    json_path = tmp_path / "manifest.json"
    json_path.write_text(json.dumps([{"sample_id": "S1"}]), encoding="utf-8")
    assert load_manifest_file(str(json_path))["samples"][0]["sample_id"] == "S1"

    tsv_path = tmp_path / "manifest.tsv"
    tsv_path.write_text("sample_id\trole\nS1\tgermline\n", encoding="utf-8")
    assert load_manifest_file(str(tsv_path))["samples"][0]["role"] == "germline"


def test_manifest_registration_is_idempotent(test_project):
    manifest = _somatic_manifest(Config.project_dir(test_project))
    first = register_manifest(test_project, manifest)
    second = register_manifest(test_project, manifest)
    assert second["id"] == first["id"]
    assert second["version"] == first["version"]
    assert len(list_manifests(test_project)) == 1
    assert get_manifest(first["id"], test_project)["manifest"] == manifest


def test_preflight_api_registers_manifest_without_launching(test_project):
    from app import create_app

    client = create_app().test_client()
    manifest = _somatic_manifest(Config.project_dir(test_project))
    response = client.post(
        f"/api/projects/{test_project}/wes/preflight",
        json={"workflow_key": "wes_somatic", "manifest": manifest},
    )
    assert response.status_code == 201
    body = response.get_json()
    assert body["valid"] is True
    assert body["manifest_id"].startswith("manifest_")
    assert "尚未启动" in body["message"]

    workflows = client.get(f"/api/projects/{test_project}/wes/workflows")
    assert workflows.status_code == 200
    assert len(workflows.get_json()["workflows"]) == 3

    prepared = client.post(
        f"/api/projects/{test_project}/wes/runs/prepare",
        json={"workflow_key": "wes_somatic", "manifest_id": body["manifest_id"]},
    )
    assert prepared.status_code == 201
    prepared_body = prepared.get_json()
    assert prepared_body["run"]["status"] == "prepared"
    assert "尚未启动" in prepared_body["message"]
    assert client.get(f"/api/projects/{test_project}/wes/runs").get_json()["runs"]

    launch_response = client.post(
        f"/api/projects/{test_project}/wes/runs/{prepared_body['run']['id']}/launch",
        json={"confirm": True},
    )
    assert launch_response.status_code == 409
    assert "WES_EXECUTOR_ENABLED" in launch_response.get_json()["error"]

    cancel_response = client.post(
        f"/api/projects/{test_project}/wes/runs/{prepared_body['run']['id']}/cancel",
        json={"confirm": True},
    )
    assert cancel_response.status_code == 202
    assert cancel_response.get_json()["run"]["status"] == "cancelled"


def test_user_capture_bed_upload_is_catalogued_as_test_only(test_project):
    from app import create_app

    client = create_app().test_client()
    response = client.post(
        f"/api/projects/{test_project}/wes/capture-kits/upload",
        data={
            "capture_kit_id": "lab-exome-v1",
            "name": "Lab Exome Panel",
            "version": "v1",
            "assembly": "GRCh38",
            "calling_bed": (io.BytesIO(b"chr1\t0\t100\nchr2\t0\t200\n"), "lab.bed"),
        },
        content_type="multipart/form-data",
    )
    assert response.status_code == 201
    body = response.get_json()
    assert body["profile"]["status"] == "test_only"
    assert body["profile"]["production_allowed"] is False
    asset = body["profile"]["assets"]["calling"]
    assert asset["metadata"]["upload_kind"] == "user_capture_bed"
    assert os.path.isfile(asset["file_path"])
    assert asset["file_path"].startswith(Config.wes_upload_root() + os.sep)
    assert get_capture_kit_profile("lab-exome-v1")["status"] == "test_only"


def test_user_capture_bed_upload_rejects_non_bed_and_keeps_catalog_clean(test_project):
    from app import create_app

    client = create_app().test_client()
    response = client.post(
        f"/api/projects/{test_project}/wes/capture-kits/upload",
        data={
            "capture_kit_id": "lab-exome-invalid",
            "name": "Lab Exome Panel",
            "version": "v1",
            "calling_bed": (io.BytesIO(b"not a bed"), "lab.txt"),
        },
        content_type="multipart/form-data",
    )
    assert response.status_code == 400
    assert "只允许上传" in response.get_json()["error"]
    assert get_capture_kit_profile("lab-exome-invalid") is None


def test_user_wes_fastq_upload_registers_pair_assets(test_project):
    from app import create_app

    client = create_app().test_client()
    response = client.post(
        f"/api/projects/{test_project}/wes/data/upload",
        data={
            "sample_id": "S1",
            "patient_id": "P1",
            "role": "germline",
            "input_type": "fastq",
            "capture_kit_id": "lab-exome-v1",
            "fastq_1": (io.BytesIO(b"@r1\nACGT\n+\n!!!!\n"), "S1_R1.fastq"),
            "fastq_2": (io.BytesIO(b"@r1\nTGCA\n+\n!!!!\n"), "S1_R2.fastq"),
        },
        content_type="multipart/form-data",
    )
    assert response.status_code == 201
    body = response.get_json()
    assert set(body["paths"]) == {"fastq_1", "fastq_2"}
    assets = list_data_assets(test_project, sample_id="S1")
    assert len(assets) == 2
    assert {item["metadata"]["mate"] for item in assets} == {"R1", "R2"}
    assert all(os.path.isfile(item["file_path"]) for item in assets)


def test_user_sra_upload_queues_async_conversion_job(test_project, monkeypatch):
    from app import create_app
    import modules.workflows.sra as sra

    monkeypatch.setattr(sra, "submit_sra_job", lambda job_id: True)
    monkeypatch.setattr(Config, "WES_SRA_MIN_FREE_GB", 1)
    client = create_app().test_client()
    response = client.post(
        f"/api/projects/{test_project}/wes/sra/upload",
        data={
            "sample_id": "SRA1",
            "patient_id": "P1",
            "role": "germline",
            "capture_kit_id": "lab-exome-v1",
            "sra": (io.BytesIO(b"SRA-placeholder"), "SRR123.sra"),
        },
        content_type="multipart/form-data",
    )
    assert response.status_code == 202
    body = response.get_json()
    assert body["job"]["status"] == "pending"
    assert body["job"]["capture_kit_id"] == "lab-exome-v1"
    assert body["source_asset"]["artifact_kind"] == "sra"
    job_response = client.get(
        f"/api/projects/{test_project}/wes/sra/jobs/{body['job']['id']}"
    )
    assert job_response.status_code == 200
    assert job_response.get_json()["job"]["sample_id"] == "SRA1"


def test_sra_conversion_registers_paired_fastq_assets(test_project, monkeypatch, tmp_path):
    from app import create_app
    import modules.workflows.sra as sra

    converter = tmp_path / "fake-fasterq-dump"
    converter.write_text(
        "#!/bin/sh\n"
        "out=''\n"
        "while [ \"$#\" -gt 0 ]; do\n"
        "  if [ \"$1\" = \"-O\" ]; then out=\"$2\"; shift 2; else shift; fi\n"
        "done\n"
        "printf '@r1\\nACGT\\n+\\n!!!!\\n@r2\\nTGCA\\n+\\n!!!!\\n' > \"$out/SRR123_1.fastq\"\n"
        "printf '@r1\\nTGCA\\n+\\n!!!!\\n@r2\\nACGT\\n+\\n!!!!\\n' > \"$out/SRR123_2.fastq\"\n",
        encoding="utf-8",
    )
    converter.chmod(0o750)
    monkeypatch.setattr(Config, "WES_SRA_FASTERQ_BIN", str(converter))
    monkeypatch.setattr(Config, "WES_SRA_MIN_FREE_GB", 1)

    client = create_app().test_client()
    response = client.post(
        f"/api/projects/{test_project}/wes/sra/upload",
        data={
            "sample_id": "SRA2",
            "patient_id": "P2",
            "role": "normal",
            "sra": (io.BytesIO(b"SRA-placeholder"), "SRR123.sra"),
        },
        content_type="multipart/form-data",
    )
    assert response.status_code == 202
    job_id = response.get_json()["job"]["id"]

    job = None
    for _ in range(100):
        job = client.get(f"/api/projects/{test_project}/wes/sra/jobs/{job_id}").get_json()["job"]
        if job["status"] in {"completed", "failed"}:
            break
        import time
        time.sleep(0.05)
    assert job["status"] == "completed", job
    assert job["output_paths"]["fastq_1"].endswith("SRA2_R1.fastq.gz")
    assert job["output_paths"]["fastq_2"].endswith("SRA2_R2.fastq.gz")
    assets = list_data_assets(test_project, sample_id="SRA2")
    assert {item["artifact_kind"] for item in assets} == {"sra", "fastq"}
    assert len([item for item in assets if item["artifact_kind"] == "fastq"]) == 2
    assert all(item["metadata"]["upload_kind"] == "sra_conversion"
               for item in assets if item["artifact_kind"] == "fastq")


def test_nextflow_executor_only_prepares_and_refuses_launch():
    from config import Config
    from modules.workflows.registry import get_workflow

    spec = NextflowExecutor().prepare(
        get_workflow("wes_germline"), "p1", "run1", "manifest1", "/tmp/work", "/tmp/results",
        intervals_path="/tmp/capture.bed",
    )
    assert spec.command[:4] == [Config.WES_NEXTFLOW_BIN, "run", "nf-core/sarek", "-r"]
    assert "-params-file" in spec.command
    assert "--wes" not in spec.command
    assert spec.parameters["wes"] is True
    assert spec.parameters["step"] == "mapping"
    assert spec.parameters["tools"] == "haplotypecaller"
    with pytest.raises(WorkflowNotConfiguredError):
        NextflowExecutor().launch(spec)


def test_nextflow_executor_refuses_calling_without_capture_bed():
    from modules.workflows.registry import get_workflow

    with pytest.raises(WorkflowNotConfiguredError, match="capture BED"):
        NextflowExecutor().prepare(
            get_workflow("wes_germline"), "p1", "run1", "manifest1",
            "/tmp/work", "/tmp/results",
        )


def test_nextflow_executor_starts_alignment_inputs_at_variant_calling():
    from modules.workflows.registry import get_workflow

    for input_type in ("bam", "cram"):
        spec = NextflowExecutor().prepare(
            get_workflow("wes_germline"), "p1", "run1", "manifest1",
            "/tmp/work", "/tmp/results", input_type=input_type,
            intervals_path="/tmp/capture.bed",
        )
        assert spec.parameters["step"] == "variant_calling"
        assert spec.input_type == input_type


def test_nextflow_executor_prepares_vcf_annotation_without_capture_bed():
    from modules.workflows.registry import get_workflow

    spec = NextflowExecutor().prepare(
        get_workflow("wes_annotate_only"), "p1", "run1", "manifest1",
        "/tmp/work", "/tmp/results", input_type="vcf",
    )

    assert spec.parameters["step"] == "annotate"
    assert spec.parameters["tools"] == "vep"
    assert "intervals" not in spec.parameters


def test_sarek_samplesheet_and_launch_bundle_are_reviewable(tmp_path):
    manifest = {
        "reference_bundle_id": "grch38-v1",
        "capture_bed_id": "exome-v1",
        "samples": [{
            "sample_id": "T01", "patient_id": "P01", "role": "tumor",
            "input_type": "fastq", "fastq_1": "/data/T_R1.fastq.gz",
            "fastq_2": "/data/T_R2.fastq.gz", "sex": "F",
        }],
    }
    samplesheet = tmp_path / "samplesheet.csv"
    rendered = render_samplesheet(manifest, "wes_somatic", str(samplesheet))
    assert rendered["row_count"] == 1
    lines = samplesheet.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("patient,sex,status,sample,lane,")
    assert "P01,F,1,T01,lane_1" in lines[1]

    bundle = write_launch_bundle(
        str(tmp_path / "run"),
        workflow={"key": "wes_somatic", "release": "3.10.0"},
        run_id="run1", project_id="p1",
        manifest_record={"id": "manifest1", "checksum": "abc", "manifest": manifest},
        launch={"command": ["nextflow", "run"]}, parameters={"capture": "exome-v1"},
    )
    assert set(bundle["files"]) == {"manifest", "samplesheet", "parameters", "provenance"}
    assert os.path.isfile(bundle["checksums_path"])
    assert len(bundle["checksums"]["samplesheet"]) == 64


def test_local_executor_records_exit_code_with_fixed_binary(test_project, tmp_path, monkeypatch):
    from modules.workflows.registry import get_workflow
    from modules.workflows.runs import create_workflow_run
    from modules.workflows.nextflow import NextflowExecutor
    monkeypatch.setattr(Config, "WES_MIN_FREE_GB", 0)

    manifest = register_manifest(test_project, {
        "reference_bundle_id": "grch38-v1", "capture_bed_id": "exome-v1",
        "samples": [{
            "sample_id": "S1", "patient_id": "P1", "role": "germline",
            "input_type": "fastq", "fastq_1": "r1", "fastq_2": "r2",
        }],
    })
    run_root = tmp_path / "run"
    work_dir = run_root / "work"
    results_dir = run_root / "results"
    launch_dir = run_root / "launch"
    launch_dir.mkdir(parents=True)
    samplesheet = launch_dir / "samplesheet.csv"
    samplesheet.write_text("patient,sample\nP1,S1\n", encoding="utf-8")
    capture_bed = launch_dir / "capture.bed"
    capture_bed.write_text("chr1\t0\t100\n", encoding="utf-8")
    run_id = "wesrun_fixed"
    spec = NextflowExecutor(nextflow_bin="/bin/true", enabled=True).prepare(
        get_workflow("wes_germline"), test_project, run_id, manifest["id"],
        str(work_dir), str(results_dir), samplesheet_path=str(samplesheet),
        intervals_path=str(capture_bed),
    )
    with open(spec.params_file_path, "w", encoding="utf-8") as handle:
        json.dump(spec.parameters, handle)
    create_workflow_run(
        test_project, "wes_germline", manifest["id"], run_id=run_id,
        launch=spec.to_dict(), run_dir=str(run_root),
    )
    executor = NextflowExecutor(
        nextflow_bin="/bin/true", enabled=True, enforce_reference_catalog=False
    )
    runtime = executor.launch(spec, project_id=test_project)
    assert runtime["pid"] > 0
    for _ in range(20):
        run = executor.poll(run_id, test_project, collect_artifacts=False)
        if run["status"] == "completed":
            break
    assert run["status"] == "completed"
    assert run["exit_code"] == 0


def test_asset_and_reference_registration_are_typed_and_checksum_aware(test_project):
    project_dir = Config.project_dir(test_project)
    bam_path = os.path.join(project_dir, "sample.bam")
    ref_path = os.path.join(project_dir, "GRCh38.fa")
    open(bam_path, "wb").write(b"bam-placeholder")
    open(ref_path, "wb").write(b"reference-placeholder")

    asset = register_data_asset(
        test_project, artifact_kind="bam", file_path=bam_path,
        sample_id="S1", source_roots=(project_dir,),
    )
    assert asset["artifact_kind"] == "bam"
    assert len(asset["checksum"]) == 64
    assert list_data_assets(test_project, sample_id="S1")[0]["id"] == asset["id"]

    reference = register_reference_asset(
        assembly="GRCh38", bundle_version="test-v1", asset_type="fasta",
        file_path=ref_path, source_roots=(project_dir,),
    )
    assert reference["assembly"] == "GRCh38"
    assert list_reference_assets(assembly="GRCh38")[0]["id"] == reference["id"]


def test_optional_content_check_is_bounded_and_dependency_aware(test_project):
    path = os.path.join(Config.project_dir(test_project), "S1_R1.fastq.gz")
    with gzip.open(path, "wb") as handle:
        handle.write(b"@r1\nACGT\n+\n!!!!\n")
    summary = inspect_sample_files({
        "sample_id": "S1", "input_type": "fastq", "fastq_1": path,
    })
    assert summary["valid"] is True
    assert summary["checks"][0]["checks"]["gzip"] is True


def test_capture_profiles_are_per_kit_and_test_only_is_not_production(test_project, monkeypatch):
    project_dir = Config.project_dir(test_project)
    bed_path = os.path.join(project_dir, "official-test.bed")
    with open(bed_path, "w", encoding="utf-8") as handle:
        handle.write("chr1\t10\t20\nchr1\t30\t50\n")
    monkeypatch.setenv("WES_SOURCE_ROOTS", project_dir)

    profile = register_capture_kit_profile(
        capture_kit_id="official-tiny-test", name="Official tiny test",
        version="v1", assembly="GRCh38", calling_bed_path=bed_path,
        status="test_only", source_roots=(project_dir,),
    )

    assert profile["capture_kit_id"] == "official-tiny-test"
    assert profile["production_allowed"] is False
    assert profile["assets"]["calling"]["metadata"]["bed_summary"]["target_bases"] == 30
    assert len(list_capture_kit_profiles()) == 1


def test_preflight_api_resolves_registered_capture_profile_path(test_project, monkeypatch):
    from app import create_app

    project_dir = Config.project_dir(test_project)
    monkeypatch.setenv("WES_SOURCE_ROOTS", project_dir)
    bed_path = os.path.join(project_dir, "profile.bed")
    with open(bed_path, "w", encoding="utf-8") as handle:
        handle.write("chr1\t0\t100\n")
    register_capture_kit_profile(
        capture_kit_id="kit-profile-v1", name="Kit profile", version="v1",
        assembly="GRCh38", calling_bed_path=bed_path, status="test_only",
        source_roots=(project_dir,),
    )
    manifest = _somatic_manifest(project_dir)
    manifest["capture_bed_id"] = "kit-profile-v1"
    manifest.pop("capture_bed_path")

    response = create_app().test_client().post(
        f"/api/projects/{test_project}/wes/preflight",
        json={"workflow_key": "wes_somatic", "manifest": manifest},
    )

    assert response.status_code == 201
    body = response.get_json()
    assert body["manifest"]["capture_bed_path"] == bed_path
    assert any("test_only" in warning for warning in body["warnings"])


def test_reference_and_capture_catalog_form_a_hard_launch_gate(test_project, monkeypatch):
    project_dir = Config.project_dir(test_project)
    monkeypatch.setenv("WES_SOURCE_ROOTS", project_dir)
    files = {}
    for asset_type in (
        "fasta", "fai", "dict", "dbsnp", "dbsnp_tbi",
        "germline_resource", "germline_resource_tbi",
    ):
        path = os.path.join(project_dir, asset_type + ".txt")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("chr1\t1000\t0\t50\t51\n" if asset_type == "fai" else asset_type + "\n")
        files[asset_type] = path
        register_reference_asset(
            assembly="GRCh38", bundle_version="bundle-v1", asset_type=asset_type,
            file_path=path, status="validated", source_roots=(project_dir,),
        )
    bed = os.path.join(project_dir, "kit.bed")
    with open(bed, "w", encoding="utf-8") as handle:
        handle.write("chr1\t0\t100\n")
    register_capture_kit_profile(
        capture_kit_id="kit-v1", name="Kit", version="v1", assembly="GRCh38",
        calling_bed_path=bed, status="validated", source_roots=(project_dir,),
    )

    ready = launch_reference_readiness(
        reference_bundle_id="bundle-v1", capture_bed_id="kit-v1",
        capture_bed_path=bed, workflow_key="wes_germline",
    )
    assert ready["valid"] is True

    wrong_path = launch_reference_readiness(
        reference_bundle_id="bundle-v1", capture_bed_id="kit-v1",
        capture_bed_path=files["fasta"], workflow_key="wes_germline",
    )
    assert wrong_path["valid"] is False
    assert any("catalog" in error for error in wrong_path["errors"])


def _completed_run_with_results(test_project, tmp_path, *, include_multiqc=True):
    from modules.workflows.runs import create_workflow_run

    manifest = register_manifest(test_project, {
        "reference_bundle_id": "bundle", "capture_bed_id": "kit",
        "samples": [{"sample_id": "S1", "patient_id": "P1", "role": "germline",
                     "input_type": "fastq", "fastq_1": "r1", "fastq_2": "r2"}],
    })
    run_root = tmp_path / ("run_" + os.urandom(4).hex())
    results = run_root / "results"
    vcf_dir = results / "variant_calling" / "haplotypecaller" / "S1"
    vcf_dir.mkdir(parents=True)
    (vcf_dir / "S1.haplotypecaller.filtered.vcf.gz").write_bytes(b"vcf")
    (vcf_dir / "S1.haplotypecaller.filtered.vcf.gz.tbi").write_bytes(b"index")
    if include_multiqc:
        multiqc = results / "multiqc"
        multiqc.mkdir()
        (multiqc / "multiqc_report.html").write_text("<html>QC</html>", encoding="utf-8")
    run_id = "wesrun_" + os.urandom(6).hex()
    run = create_workflow_run(
        test_project, "wes_germline", manifest["id"], run_id=run_id,
        status="completed", run_dir=str(run_root),
    )
    return run


def test_completed_run_collects_primary_vcf_index_and_multiqc(test_project, tmp_path):
    run = _completed_run_with_results(test_project, tmp_path)
    collection = collect_workflow_artifacts(run["id"], test_project)

    assert collection["valid"] is True
    kinds = {item["artifact_kind"] for item in list_workflow_artifacts(run["id"], test_project)}
    assert {"filtered_vcf", "vcf_index", "multiqc_report"} <= kinds
    assert all(len(item["checksum"]) == 64 for item in collection["artifacts"])


def test_missing_required_result_turns_zero_exit_completion_into_failure(test_project, tmp_path):
    from modules.workflows.runs import get_workflow_run

    run = _completed_run_with_results(test_project, tmp_path, include_multiqc=False)
    collection = collect_workflow_artifacts(run["id"], test_project, strict=True)

    assert collection["valid"] is False
    assert "MultiQC" in "; ".join(collection["errors"])
    assert get_workflow_run(run["id"], test_project)["status"] == "failed"


def test_p3_evaluation_bed_and_threshold_parser(tmp_path):
    truth = tmp_path / "truth.bed"
    capture = tmp_path / "capture.bed"
    callable_bed = tmp_path / "callable.bed"
    truth.write_text("chr1\t0\t100\n", encoding="utf-8")
    capture.write_text("chr1\t20\t90\n", encoding="utf-8")
    callable_bed.write_text("chr1\t40\t120\n", encoding="utf-8")
    output = tmp_path / "evaluation.bed"

    summary = write_evaluation_bed(
        [str(truth), str(capture), str(callable_bed)], str(output)
    )
    assert output.read_text(encoding="utf-8") == "chr1\t40\t90\n"
    assert summary["bases"] == 50

    metrics_csv = tmp_path / "happy.summary.csv"
    metrics_csv.write_text(
        "Type,Filter,METRIC.Precision,METRIC.Recall,METRIC.F1_Score,TRUTH.TP,TRUTH.FN,QUERY.FP\n"
        "SNP,PASS,0.995,0.98,0.987,98,2,1\n"
        "INDEL,PASS,0.985,0.91,0.946,91,9,2\n",
        encoding="utf-8",
    )
    metrics = parse_benchmark_summary(str(metrics_csv))
    assert evaluate_metrics(metrics, mode="germline")["status"] == "validated"
    assert evaluate_metrics(metrics, mode="somatic")["status"] == "not_validated"


def test_wes_dashboard_and_project_scoped_artifact_download(test_project, tmp_path):
    from app import create_app

    run = _completed_run_with_results(test_project, tmp_path)
    collection = collect_workflow_artifacts(run["id"], test_project)
    report = next(item for item in collection["artifacts"]
                  if item["artifact_kind"] == "multiqc_report")
    client = create_app().test_client()

    page = client.get(f"/projects/{test_project}/wes")
    assert page.status_code == 200
    assert "WES 工作台" in page.get_data(as_text=True)
    download = client.get(f"/projects/{test_project}/wes/artifacts/{report['id']}")
    assert download.status_code == 200
    assert b"QC" in download.data


def test_p3_restart_marks_untracked_local_process_interrupted(test_project, tmp_path):
    from modules.workflows.nextflow import NextflowExecutor
    from modules.workflows.runs import create_workflow_run

    manifest = register_manifest(test_project, {
        "reference_bundle_id": "bundle", "capture_bed_id": "kit",
        "samples": [{"sample_id": "S1", "patient_id": "P1", "role": "germline",
                     "input_type": "fastq", "fastq_1": "r1", "fastq_2": "r2"}],
    })
    run = create_workflow_run(
        test_project, "wes_germline", manifest["id"], status="running",
        run_dir=str(tmp_path / "restart-run"),
    )

    recovered = NextflowExecutor.poll(run["id"], test_project)
    assert recovered["status"] == "interrupted"
    assert "平台进程重启" in recovered["error_text"]


def test_p3_disk_gate_blocks_launch_before_process_start(test_project, tmp_path, monkeypatch):
    from modules.workflows.nextflow import NextflowExecutor
    from modules.workflows.registry import get_workflow
    from modules.workflows.runs import create_workflow_run

    manifest = register_manifest(test_project, {
        "reference_bundle_id": "bundle", "capture_bed_id": "kit",
        "samples": [{"sample_id": "S1", "patient_id": "P1", "role": "germline",
                     "input_type": "fastq", "fastq_1": "r1", "fastq_2": "r2"}],
    })
    run_root = tmp_path / "disk-run"
    launch_dir = run_root / "launch"
    launch_dir.mkdir(parents=True)
    samplesheet = launch_dir / "samplesheet.csv"
    samplesheet.write_text("patient,sample\nP1,S1\n", encoding="utf-8")
    capture = launch_dir / "capture.bed"
    capture.write_text("chr1\t0\t10\n", encoding="utf-8")
    executor = NextflowExecutor(
        nextflow_bin="/bin/true", enabled=True, enforce_reference_catalog=False
    )
    spec = executor.prepare(
        get_workflow("wes_germline"), test_project, "wesrun_disk", manifest["id"],
        str(run_root / "work"), str(run_root / "results"),
        samplesheet_path=str(samplesheet), intervals_path=str(capture),
    )
    with open(spec.params_file_path, "w", encoding="utf-8") as handle:
        json.dump(spec.parameters, handle)
    create_workflow_run(
        test_project, "wes_germline", manifest["id"], run_id="wesrun_disk",
        launch=spec.to_dict(), run_dir=str(run_root),
    )
    monkeypatch.setattr(Config, "WES_MIN_FREE_GB", 10 ** 9)

    with pytest.raises(WorkflowNotConfiguredError, match="剩余空间"):
        executor.launch(spec, project_id=test_project)


def test_p3_cancel_terminates_the_process_group(test_project, tmp_path, monkeypatch):
    from modules.workflows.nextflow import NextflowExecutor
    from modules.workflows.registry import get_workflow
    from modules.workflows.runs import create_workflow_run

    monkeypatch.setattr(Config, "WES_MIN_FREE_GB", 0)
    fake_nextflow = tmp_path / "fake-nextflow"
    fake_nextflow.write_text("#!/bin/sh\nsleep 30\n", encoding="utf-8")
    fake_nextflow.chmod(0o755)
    manifest = register_manifest(test_project, {
        "reference_bundle_id": "bundle", "capture_bed_id": "kit",
        "samples": [{"sample_id": "S1", "patient_id": "P1", "role": "germline",
                     "input_type": "fastq", "fastq_1": "r1", "fastq_2": "r2"}],
    })
    run_root = tmp_path / "cancel-run"
    launch_dir = run_root / "launch"
    launch_dir.mkdir(parents=True)
    samplesheet = launch_dir / "samplesheet.csv"
    samplesheet.write_text("patient,sample\nP1,S1\n", encoding="utf-8")
    capture = launch_dir / "capture.bed"
    capture.write_text("chr1\t0\t10\n", encoding="utf-8")
    executor = NextflowExecutor(
        nextflow_bin=str(fake_nextflow), enabled=True, enforce_reference_catalog=False
    )
    run_id = "wesrun_cancel"
    spec = executor.prepare(
        get_workflow("wes_germline"), test_project, run_id, manifest["id"],
        str(run_root / "work"), str(run_root / "results"),
        samplesheet_path=str(samplesheet), intervals_path=str(capture),
    )
    with open(spec.params_file_path, "w", encoding="utf-8") as handle:
        json.dump(spec.parameters, handle)
    create_workflow_run(
        test_project, "wes_germline", manifest["id"], run_id=run_id,
        launch=spec.to_dict(), run_dir=str(run_root),
    )
    executor.launch(spec, project_id=test_project)
    requested = executor.cancel(run_id, test_project)
    assert requested["status"] == "cancel_requested"
    for _ in range(30):
        cancelled = executor.poll(run_id, test_project, collect_artifacts=False)
        if cancelled["status"] == "cancelled":
            break
    assert cancelled["status"] == "cancelled"
