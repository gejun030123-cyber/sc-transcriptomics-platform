"""Regression coverage for sample-aware functional-state analysis."""

import json
import os

import numpy as np
import pandas as pd
import pytest


def _functional_adata(tmp_path):
    anndata = pytest.importorskip("anndata")
    rng = np.random.default_rng(17)
    genes = [
        "PPARA", "HNF4A", "ESRRA", "CPT1A", "CPT2", "SLC25A20", "ACADM",
        "ACADVL", "ACADS", "HADHA", "HADHB", "ECHS1", "ETFDH", "ACAA2",
        "ACOX1", "EHHADH", "ACAA1", "HSD17B4", "ABCD3", "SCP2", "FABP1",
        "FABP2", "SLC27A2", "HMGCS2", "HMGCL", "BDH1", "CXCL8", "CCL20",
        "ICAM1", "NFKBIA", "TNFAIP3", "IL1B", "IL6", "CXCL2", "CXCL3", "CCL2",
        "RELA", "NFKB1", "STAT1", "IRF1", "IRF7", "JUN", "FOS",
    ]
    rows, counts = [], []
    for condition in ("Healthy", "IBD"):
        for replicate in range(2):
            sample = f"{condition}_{replicate + 1}"
            for celltype in ("Colonocyte", "Goblet"):
                for _ in range(5):
                    values = rng.poisson(2, len(genes)).astype(float)
                    if condition == "Healthy":
                        values[3:26] += 7  # metabolic programme
                    else:
                        values[26:36] += 8  # inflammatory programme
                    counts.append(values)
                    rows.append({"sample_id": sample, "condition": condition, "celltype": celltype})
    counts = np.asarray(counts)
    normalized = np.log1p(counts / np.maximum(counts.sum(axis=1, keepdims=True), 1) * 10_000)
    adata = anndata.AnnData(normalized, obs=pd.DataFrame(rows), var=pd.DataFrame(index=genes))
    adata.layers["counts"] = counts
    adata.obsm["X_umap"] = rng.normal(size=(adata.n_obs, 2))
    path = tmp_path / "functional_input.h5ad"
    adata.write_h5ad(path)
    return path


def _run(tmp_path, path, **kwargs):
    from modules.functional_state import FunctionalStateAnalysis

    params = {
        "sample_key": "sample_id",
        "condition_key": "condition",
        "celltype_key": "celltype",
        "reference_condition": "Healthy",
        "comparison_condition": "IBD",
        "min_cells_per_sample_celltype": 3,
        "min_samples_per_group": 2,
        "min_gene_set_genes": 3,
        "min_tf_targets": 2,
    }
    params.update(kwargs)
    return FunctionalStateAnalysis(str(tmp_path), params, lambda *_: None).run(str(path))


def test_functional_state_keeps_gene_pathway_and_tf_layers_separate(tmp_path):
    anndata = pytest.importorskip("anndata")
    from modules.functional_state import PPARA_TARGET_MODULE_NAME
    result = _run(tmp_path, _functional_adata(tmp_path))

    assert os.path.exists(result["output_adata"])
    assert result["summary"]["n_pathway_scores"] > 0
    assert result["summary"]["n_expression_features"] > 0
    assert result["summary"]["n_tf_activities"] == 0
    assert result["summary"]["tf_network_configured"] is False
    assert any("不可用" in warning for warning in result["summary"]["warnings"])

    adata = anndata.read_h5ad(result["output_adata"])
    assert f"fs_pathway_{PPARA_TARGET_MODULE_NAME.replace('/', '_').replace(' ', '_')}_score" in adata.obs
    assert not any(column.startswith("fs_tf_") for column in adata.obs.columns)

    result_root = tmp_path / "results" / "functional_state" / "manual_run"
    stats = pd.read_csv(result_root / "functional_state_statistics.csv")
    assert {"gene_expression", "pathway_score"}.issubset(set(stats["feature_type"]))
    assert set(stats["status"]) == {"sample_level_exploratory"}
    assert set(stats["evidence_tier"]) == {"exploratory_n2"}
    assert stats["fdr_bh"].notna().all()
    coverage = pd.read_csv(result_root / "tf_activity_coverage.csv")
    assert set(coverage["status"]) == {"network_not_configured"}
    manifest = json.loads((result_root / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert manifest["tf_network"]["status"] == "not_configured"
    assert manifest["pseudobulk_validation"]["existing_modules"] == ["sc_pseudobulk_deg", "sc_cell_go"]


def test_functional_state_calculates_network_backed_tf_activity_from_project_file(tmp_path):
    anndata = pytest.importorskip("anndata")
    path = _functional_adata(tmp_path)
    uploads = tmp_path / "uploads"
    uploads.mkdir()
    network = pd.DataFrame({
        "source": ["PPARA"] * 5 + ["RELA"] * 5,
        "target": ["CPT1A", "ACADM", "ACOX1", "FABP1", "HMGCS2", "CXCL8", "CCL20", "ICAM1", "NFKBIA", "TNFAIP3"],
        "weight": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    })
    network.to_csv(uploads / "collectri.tsv", sep="\t", index=False)

    result = _run(
        tmp_path, path, tf_network_path="uploads/collectri.tsv", tf_panel="PPARA,RELA",
        heatmap_features="Inflammatory response score, PPARA activity",
    )

    assert result["summary"]["n_tf_activities"] == 2
    adata = anndata.read_h5ad(result["output_adata"])
    assert {"fs_tf_PPARA_activity", "fs_tf_RELA_activity"}.issubset(adata.obs.columns)
    result_root = tmp_path / "results" / "functional_state" / "manual_run"
    coverage = pd.read_csv(result_root / "tf_activity_coverage.csv")
    assert set(coverage["status"]) == {"low_coverage"}
    assert {"network_targets", "detected_targets", "coverage_fraction", "activity_status", "coverage_level"}.issubset(coverage.columns)
    assert set(coverage["coverage_level"]) == {"low_5_9"}
    assert set(coverage["method"]) == {"decoupler_ulm"}
    correlations = pd.read_csv(result_root / "activity_correlations.csv")
    assert "PPARA activity" in set(correlations["x_feature"])
    assert {"condition", "effect_size_interpretation", "n_unique_samples", "repeated_donor_measures"}.issubset(correlations.columns)
    assert "within_celltype_cell_level_descriptive" in set(correlations["analysis_level"])
    assert {"within_celltype_cell_level_descriptive_by_condition", "sample_x_celltype_by_condition"}.issubset(
        set(correlations["analysis_level"])
    )
    condition_rows = correlations[correlations["analysis_level"].eq("sample_x_celltype_by_condition")]
    assert set(condition_rows["condition"]) == {"Healthy", "IBD"}
    assert set(condition_rows["status"]) == {"exploratory_repeated_measures_condition_stratified"}
    aggregate = correlations[correlations["analysis_level"].eq("sample_x_celltype")].iloc[0]
    assert aggregate["status"] == "exploratory_repeated_measures"
    assert bool(aggregate["repeated_donor_measures"])
    manifest = json.loads((result_root / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert manifest["heatmap_display"]["mode"] == "user_ordered_selection"
    assert manifest["heatmap_display"]["features"] == [
        "Acute myeloid chemokine inflammation score", "PPARA activity",
    ]
    assert manifest["correlation"]["x_feature_requested"] == "PPARA activity"
    assert manifest["pathway_scoring_method"]["requested"] == "scanpy_score_genes"
    assert manifest["tf_coverage_policy"]["minimum_detected_targets_for_scoring"] == 5
    assert result["summary"]["n_heatmap_display_features"] == 2
    result_paths = {os.path.basename(item["file_path"]) for item in result["result_files"]}
    assert "functional_state_score_violin.png" in result_paths
    assert "functional_state_cell_level_score_correlation.png" in result_paths
    assert "functional_state_within_celltype_cell_level_correlation.png" in result_paths
    assert "functional_state_within_celltype_condition_cell_level_correlation.png" in result_paths
    assert "functional_state_condition_sample_celltype_score_correlation.png" in result_paths
    assert "functional_state_regulator_pathway_concordance.png" in result_paths
    assert "functional_state_regulator_pathway_concordance_by_condition.png" in result_paths
    dotplot_svg = (tmp_path / "plots" / "functional_state_gene_expression_dotplot.svg").read_text(encoding="utf-8")
    assert "Regulatory TF genes" in dotplot_svg
    expression_vs_activity = pd.read_csv(result_root / "tf_expression_vs_activity.csv")
    assert "PPARA" in set(expression_vs_activity["tf"])
    assert {"activity_status", "coverage_level", "expression_mean_difference", "activity_mean_difference"}.issubset(expression_vs_activity.columns)
    concordance = pd.read_csv(result_root / "regulator_pathway_concordance.csv")
    assert {"rho", "effect_size_interpretation", "analysis_note"}.issubset(concordance.columns)
    assert set(concordance["status"]) == {"exploratory_repeated_measures"}
    concordance_by_condition = pd.read_csv(result_root / "regulator_pathway_concordance_by_condition.csv")
    assert set(concordance_by_condition["condition"]) == {"Healthy", "IBD"}
    assert set(concordance_by_condition["status"]) == {"exploratory_repeated_measures_condition_stratified"}
    assert manifest["regulator_pathway_concordance"]["condition_stratified_table"] == (
        "regulator_pathway_concordance_by_condition.csv"
    )


def test_functional_state_uses_checksum_verified_managed_collectri_snapshot(tmp_path, monkeypatch):
    from config import Config
    from modules.functional_state import _sha256

    path = _functional_adata(tmp_path)
    data_dir = tmp_path / "platform_data"
    resource_dir = data_dir / "functional_state_resources"
    resource_dir.mkdir(parents=True)
    network = pd.DataFrame({
        "source": ["PPARA"] * 5 + ["RELA"] * 5,
        "target": ["CPT1A", "ACADM", "ACOX1", "FABP1", "HMGCS2", "CXCL8", "CCL20", "ICAM1", "NFKBIA", "TNFAIP3"],
        "weight": [1] * 10,
    })
    snapshot = resource_dir / "collectri_human.tsv"
    network.to_csv(snapshot, sep="\t", index=False)
    (resource_dir / "collectri_human.metadata.json").write_text(json.dumps({
        "resource": "CollecTRI", "organism": "human", "sha256": _sha256(snapshot),
    }), encoding="utf-8")
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "FUNCTIONAL_STATE_RESOURCE_DIR", "")

    result = _run(tmp_path, path, tf_network_source="managed_collectri", tf_panel="PPARA,RELA")

    assert result["summary"]["n_tf_activities"] == 2
    root = tmp_path / "results" / "functional_state" / "manual_run"
    coverage = pd.read_csv(root / "tf_activity_coverage.csv")
    assert set(coverage["method"]) == {"decoupler_ulm"}
    manifest = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert manifest["tf_network"]["source"] == "managed_collectri"
    assert manifest["tf_network"]["sha256"] == _sha256(snapshot)


def test_functional_state_rebuilds_log_expression_from_counts_not_residuals(tmp_path):
    anndata = pytest.importorskip("anndata")
    path = _functional_adata(tmp_path)
    adata = anndata.read_h5ad(path)
    adata.X = np.asarray(adata.X) - 4.0
    adata.uns["normalization"] = {"x_contains": "pearson_residuals"}
    adata.write_h5ad(path)

    result = _run(tmp_path, path)
    manifest = json.loads((tmp_path / "results" / "functional_state" / "manual_run" / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert manifest["expression_source"] == "layers[_functional_log1p]"
    assert manifest["expression_qc"]["decision"] == "reconstructed_from_counts_due_to_negative_or_pearson_residuals"
    output = anndata.read_h5ad(result["output_adata"])
    assert "_functional_log1p" in output.layers
    assert float(np.min(output.X)) < 0  # The input X is preserved, not overwritten.


def test_functional_state_locks_tf_coverage_floor_and_custom_layer_contract(tmp_path):
    path = _functional_adata(tmp_path)
    uploads = tmp_path / "uploads"
    uploads.mkdir()
    pd.DataFrame({
        "source": ["PPARA"] * 4,
        "target": ["CPT1A", "ACADM", "ACOX1", "FABP1"],
        "weight": [1] * 4,
    }).to_csv(uploads / "four_targets.tsv", sep="\t", index=False)

    result = _run(
        tmp_path, path, tf_network_path="uploads/four_targets.tsv", tf_panel="PPARA",
        min_tf_targets=2,  # Historical parameter cannot lower the fixed coverage floor.
    )
    root = tmp_path / "results" / "functional_state" / "manual_run"
    coverage = pd.read_csv(root / "tf_activity_coverage.csv")
    assert result["summary"]["n_tf_activities"] == 0
    assert coverage.loc[0, "status"] == "below_minimum_coverage"
    assert coverage.loc[0, "coverage_level"] == "below_5_not_scored"

    with pytest.raises(ValueError, match="custom_layer"):
        _run(tmp_path, path, expression_source="custom_layer", expression_layer="")


def test_functional_state_refuses_technical_batch_as_sample_without_confirmation(tmp_path):
    anndata = pytest.importorskip("anndata")
    path = _functional_adata(tmp_path)
    adata = anndata.read_h5ad(path)
    adata.obs["batch"] = adata.obs["sample_id"].astype(str)
    adata.write_h5ad(path)

    with pytest.raises(ValueError, match="技术 batch"):
        _run(tmp_path, path, sample_key="batch")


def test_functional_state_registration_and_pipeline_dependency():
    from modules import MODULE_REGISTRY, PIPELINE_DEPS
    from modules.schemas import SC_MODULE_NAMES

    assert "functional_state" in MODULE_REGISTRY
    assert "functional_state" in SC_MODULE_NAMES
    assert PIPELINE_DEPS["functional_state"] == ["annotation"]


def test_functional_state_schema_exposes_locked_defaults_and_dynamic_choices():
    from modules.schemas import PARAM_SCHEMAS

    fields = {field["key"]: field for field in PARAM_SCHEMAS["functional_state"]}
    assert fields["correlation_x"]["default"] == "PPARA activity"
    assert fields["reference_condition"]["type"] == "dynamic_select"
    assert fields["comparison_condition"]["depends_on"] == "condition_key"
    assert fields["pathway_scoring_method"]["options"] == ["scanpy_score_genes"]
    assert "within_celltype_cell_level_descriptive" in fields["correlation_levels"]["options"]
    assert "PPARA-associated lipid-oxidation programme score" in fields["correlation_x"]["options"]
    assert fields["concordance_tf_features"]["default"].startswith("PPARA activity")
    assert fields["include_builtin_pathway_panels"]["default"] is True
    assert "HALLMARK_OXIDATIVE_PHOSPHORYLATION" in fields["managed_gene_set_terms"]["options"]
    assert fields["managed_gene_set_terms_advanced"]["type"] == "textarea"
    assert "min_tf_targets" not in fields
    assert fields["confirm_batch_is_biological_sample"]["show_if"] == {
        "sample_key": ["batch", "technical_batch", "sequencing_batch", "library_batch"],
    }


def test_functional_state_resolves_legacy_builtin_pathway_names_to_precise_labels():
    from modules.functional_state import LEGACY_BUILTIN_PATHWAY_NAMES, _resolve_feature_name

    feature_values = {
        f"{current_name} score": np.zeros(2)
        for current_name in LEGACY_BUILTIN_PATHWAY_NAMES.values()
    }
    for legacy_name, current_name in LEGACY_BUILTIN_PATHWAY_NAMES.items():
        assert _resolve_feature_name(f"{legacy_name} score", feature_values) == f"{current_name} score"


def test_functional_state_reads_only_selected_project_local_gmt_terms(tmp_path):
    path = _functional_adata(tmp_path)
    uploads = tmp_path / "uploads"
    uploads.mkdir()
    (uploads / "hallmark.gmt").write_text(
        "HALLMARK_OXPHOS\tna\tCPT1A\tACADM\tACOX1\n"
        "HALLMARK_UNSELECTED\tna\tCXCL8\tCCL20\tIL1B\n",
        encoding="utf-8",
    )

    _run(
        tmp_path, path, pathway_gmt_path="uploads/hallmark.gmt",
        pathway_gmt_terms="HALLMARK_OXPHOS",
    )
    root = tmp_path / "results" / "functional_state" / "manual_run"
    coverage = pd.read_csv(root / "gene_set_coverage.csv")
    assert "HALLMARK_OXPHOS" in set(coverage["gene_set"])
    assert "HALLMARK_UNSELECTED" not in set(coverage["gene_set"])
    manifest = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert manifest["pathway_gmt"]["selected_terms"] == ["HALLMARK_OXPHOS"]


def test_functional_state_reads_checksum_verified_managed_gene_sets_with_coverage_qc(tmp_path, monkeypatch):
    from config import Config
    from modules.gene_set_registry import _sha256

    data_dir = tmp_path / "platform_data"
    gene_set_dir = data_dir / "functional_state_resources" / "gene_sets" / "hallmark"
    gene_set_dir.mkdir(parents=True)
    snapshot = gene_set_dir / "hallmark.gmt"
    detected = ["CPT1A", "CPT2", "SLC25A20", "ACADM", "ACADVL", "ACADS", "HADHA", "HADHB", "ECHS1", "ETFDH"]
    snapshot.write_text(
        "HALLMARK_TEST_GOOD\tna\t" + "\t".join(detected) + "\n"
        + "HALLMARK_TEST_LOW_COVERAGE\tna\t" + "\t".join(detected + ["ACOX1"] + [f"MISSING{i}" for i in range(20)]) + "\n",
        encoding="utf-8",
    )
    registry = {
        "schema_version": 1,
        "libraries": [{
            "key": "hallmark_human", "resource": "MSigDB Hallmark", "version": "test",
            "file": "hallmark/hallmark.gmt", "sha256": _sha256(snapshot),
            "license": "CC BY 4.0", "license_url": "https://example.test/license", "source_url": "https://example.test/source",
        }],
    }
    (gene_set_dir.parent / "gene_set_registry.json").write_text(json.dumps(registry), encoding="utf-8")
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "FUNCTIONAL_STATE_RESOURCE_DIR", "")

    _run(
        tmp_path, _functional_adata(tmp_path),
        managed_gene_set_terms="HALLMARK_TEST_GOOD,HALLMARK_TEST_LOW_COVERAGE",
        include_builtin_pathway_panels=False,
    )

    root = tmp_path / "results" / "functional_state" / "manual_run"
    coverage = pd.read_csv(root / "gene_set_coverage.csv").set_index("gene_set")
    assert coverage.loc["HALLMARK_TEST_GOOD", "source"] == "managed_hallmark_human"
    assert coverage.loc["HALLMARK_TEST_GOOD", "coverage_level"] == "good_coverage"
    assert coverage.loc["HALLMARK_TEST_LOW_COVERAGE", "status"] == "unreliable_coverage_not_scored"
    manifest = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert manifest["managed_gene_sets"]["selected_terms"] == [
        "HALLMARK_TEST_GOOD", "HALLMARK_TEST_LOW_COVERAGE",
    ]
    assert manifest["pathway_resource_versions"]["managed_hallmark_human"]["version"] == "test"
    assert manifest["builtin_quick_panels"]["included"] is False
    assert "builtin_lipid_inflammation" not in manifest["pathway_resource_versions"]


def test_functional_state_refuses_tampered_managed_gene_set_snapshot(tmp_path, monkeypatch):
    from config import Config
    from modules.gene_set_registry import _sha256, load_selected_managed_gene_sets

    data_dir = tmp_path / "platform_data"
    root = data_dir / "functional_state_resources" / "gene_sets"
    root.mkdir(parents=True)
    snapshot = root / "library.gmt"
    snapshot.write_text("TEST\tna\tCPT1A\tACADM\n", encoding="utf-8")
    registry = {
        "schema_version": 1,
        "libraries": [{
            "key": "test", "file": "library.gmt", "sha256": _sha256(snapshot),
        }],
    }
    (root / "gene_set_registry.json").write_text(json.dumps(registry), encoding="utf-8")
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "FUNCTIONAL_STATE_RESOURCE_DIR", "")

    sets, _ = load_selected_managed_gene_sets("TEST")
    assert sets["TEST"] == ("CPT1A", "ACADM")
    snapshot.write_text("TEST\tna\tCPT1A\tACADM\tACOX1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="校验失败"):
        load_selected_managed_gene_sets("TEST")


def test_functional_state_marks_small_inline_custom_signature_as_exploratory(tmp_path):
    _run(
        tmp_path, _functional_adata(tmp_path), min_gene_set_genes=5,
        custom_gene_sets="Custom lipid handling: CPT1A,FABP1,HMGCS2",
    )

    root = tmp_path / "results" / "functional_state" / "manual_run"
    coverage = pd.read_csv(root / "gene_set_coverage.csv")
    row = coverage.loc[coverage["gene_set"].eq("Custom lipid handling")].iloc[0]
    assert row["n_detected"] == 3
    assert row["signature_kind"] == "exploratory_custom_signature"
    assert row["status"] == "exploratory_custom_signature"
    scores = pd.read_csv(root / "sample_celltype_scores.csv")
    assert "Custom lipid handling score" in set(scores["feature"])


def test_functional_state_uses_distinct_inflammation_and_tnf_nfkB_signatures_with_overlap_qc(tmp_path):
    from modules.functional_state import (
        ACUTE_MYELOID_INFLAMMATION_MODULE_NAME,
        INFLAMMATORY_GENES,
        TNF_NFKB_RESPONSE_GENES,
    )

    assert set(INFLAMMATORY_GENES) != set(TNF_NFKB_RESPONSE_GENES)
    _run(tmp_path, _functional_adata(tmp_path))

    root = tmp_path / "results" / "functional_state" / "manual_run"
    qc = pd.read_csv(root / "pathway_gene_set_overlap_qc.csv")
    row = qc.loc[
        qc["gene_set_a"].eq(ACUTE_MYELOID_INFLAMMATION_MODULE_NAME)
        & qc["gene_set_b"].eq("TNF-NFkB response")
    ].iloc[0]
    assert row["n_shared_detected"] < row["n_detected_a"]
    assert row["n_shared_detected"] < row["n_detected_b"]
    assert row["status"] == "distinct_detected_genes"
    manifest = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert manifest["pathway_resource_versions"]["builtin_lipid_inflammation"]["version"].endswith("_v2")
    assert manifest["pathway_gene_set_overlap_qc"]["table"] == "pathway_gene_set_overlap_qc.csv"


def test_functional_state_ibd_epithelial_focus_uses_frozen_hallmark_terms(tmp_path, monkeypatch):
    """The IBD preset must load a finite, auditable local term set only."""
    from config import Config
    from modules.functional_state import IBD_EPITHELIAL_MANAGED_TERMS, _sha256

    path = _functional_adata(tmp_path)
    resource_dir = tmp_path / "platform_data" / "functional_state_resources" / "gene_sets"
    resource_dir.mkdir(parents=True)
    anndata = pytest.importorskip("anndata")
    genes = list(anndata.read_h5ad(path).var_names.astype(str))
    hallmark = resource_dir / "hallmark.gmt"
    hallmark.write_text("\n".join(
        "\t".join([term, "test", *[genes[(index * 3 + offset) % len(genes)] for offset in range(10)]])
        for index, term in enumerate(IBD_EPITHELIAL_MANAGED_TERMS)
    ) + "\n", encoding="utf-8")
    (resource_dir / "gene_set_registry.json").write_text(json.dumps({
        "schema_version": 1,
        "libraries": [{
            "key": "hallmark_human", "resource": "MSigDB Hallmark", "version": "test",
            "file": "hallmark.gmt", "sha256": _sha256(hallmark),
        }],
    }), encoding="utf-8")
    monkeypatch.setattr(Config, "DATA_DIR", str(tmp_path / "platform_data"))
    monkeypatch.setattr(Config, "FUNCTIONAL_STATE_RESOURCE_DIR", "")

    result = _run(
        tmp_path, path, analysis_focus="ibd_organoid_epithelial",
        include_builtin_pathway_panels=False,
    )

    root = tmp_path / "results" / "functional_state" / "manual_run"
    manifest = json.loads((root / "analysis_manifest.json").read_text(encoding="utf-8"))
    stats = pd.read_csv(root / "pathway_statistics.csv")
    assert result["summary"]["n_pathway_scores"] >= len(IBD_EPITHELIAL_MANAGED_TERMS)
    assert manifest["analysis_focus"] == "ibd_organoid_epithelial"
    assert manifest["ibd_epithelial_focus"]["managed_hallmark_terms"] == list(IBD_EPITHELIAL_MANAGED_TERMS)
    assert set(IBD_EPITHELIAL_MANAGED_TERMS).issubset(
        set(stats["feature"].str.replace(" score", "", regex=False))
    )


def test_functional_state_blocks_duplicate_pathway_scores_after_gene_set_qc(tmp_path):
    with pytest.raises(ValueError, match="基因集重复 QC"):
        _run(
            tmp_path, _functional_adata(tmp_path),
            custom_gene_sets=(
                "Duplicate inflammation: IL1B,IL6,CCL2,CXCL8,CXCL2,CXCL3,CCL20,CXCL1,CXCL5,S100A8,S100A9"
            ),
        )
