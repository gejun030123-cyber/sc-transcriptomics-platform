"""A clone must carry the public single-cell resources needed for offline runs."""

from pathlib import Path


def test_bundled_resources_resolve_without_data_directory(tmp_path, monkeypatch):
    from config import Config
    from modules.functional_state import _read_tf_network
    from modules.gene_set_registry import load_selected_managed_gene_sets
    from modules.sc_cell_go import (
        _local_gene_set_path, _read_local_gene_sets, _resolve_local_gene_set_dir,
    )

    monkeypatch.setattr(Config, "DATA_DIR", str(tmp_path / "empty-data"))
    monkeypatch.setattr(Config, "FUNCTIONAL_STATE_RESOURCE_DIR", "")
    monkeypatch.setattr(Config, "SC_CELL_GO_GENE_SET_DIR", "")

    resource_root = Path(Config.functional_state_resource_dir())
    assert resource_root.name == "functional_state_resources"
    assert (resource_root / "gene_sets" / "gene_set_registry.json").is_file()
    assert _resolve_local_gene_set_dir("", tmp_path / "project") == Config.sc_cell_go_gene_set_dir()

    gene_sets, provenance = load_selected_managed_gene_sets(
        "HALLMARK_INFLAMMATORY_RESPONSE"
    )
    assert len(gene_sets["HALLMARK_INFLAMMATORY_RESPONSE"]) >= 10
    assert provenance["HALLMARK_INFLAMMATORY_RESPONSE"]["library_key"] == "hallmark_human"

    local_path = _local_gene_set_path(
        "GO_Biological_Process_2023", Config.sc_cell_go_gene_set_dir()
    )
    assert len(_read_local_gene_sets(local_path)) >= 1000

    wikipathways_path = _local_gene_set_path(
        "WikiPathway_2021_Human", Config.sc_cell_go_gene_set_dir()
    )
    wikipathways = _read_local_gene_sets(wikipathways_path)
    assert len(wikipathways) >= 800
    assert "TP53" in {gene for genes in wikipathways.values() for gene in genes}

    network = _read_tf_network(resource_root / "collectri_human.tsv")
    assert (network["tf"] == "PPARA").sum() >= 5
    assert (resource_root / "collectri_human.metadata.json").is_file()


def test_sc_cell_go_supplemental_kegg_keeps_bundled_gene_sets(tmp_path, monkeypatch):
    from config import Config
    from modules.sc_cell_go import _local_gene_set_path, _read_local_gene_sets

    data_dir = tmp_path / "data"
    kegg_dir = data_dir / "go_gene_sets"
    kegg_dir.mkdir(parents=True)
    (kegg_dir / "KEGG_2021_Human.gmt").write_text(
        "Example KEGG\thttps://www.kegg.jp/\tTP53\tEGFR\n", encoding="utf-8"
    )
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "FUNCTIONAL_STATE_RESOURCE_DIR", "")
    monkeypatch.setattr(Config, "SC_CELL_GO_GENE_SET_DIR", "")

    configured_dir = Config.sc_cell_go_gene_set_dir()
    assert configured_dir == str(kegg_dir)
    assert _read_local_gene_sets(
        _local_gene_set_path("KEGG_2021_Human", configured_dir)
    )["Example KEGG"] == ["TP53", "EGFR"]
    assert len(_read_local_gene_sets(
        _local_gene_set_path("GO_Biological_Process_2023", configured_dir)
    )) >= 1000


def test_sc_cell_go_keeps_bundled_gmts_when_only_tf_network_is_locally_overridden(tmp_path, monkeypatch):
    from config import Config

    data_dir = tmp_path / "data"
    local_tf_root = data_dir / "functional_state_resources"
    local_tf_root.mkdir(parents=True)
    (local_tf_root / "collectri_human.tsv").write_text(
        "source\ttarget\tweight\nPPARA\tCPT1A\t1\n", encoding="utf-8"
    )
    monkeypatch.setattr(Config, "DATA_DIR", str(data_dir))
    monkeypatch.setattr(Config, "FUNCTIONAL_STATE_RESOURCE_DIR", "")
    monkeypatch.setattr(Config, "SC_CELL_GO_GENE_SET_DIR", "")

    assert Config.functional_state_resource_dir() == str(local_tf_root)
    assert (Path(Config.sc_cell_go_gene_set_dir()) / "gene_set_registry.json").is_file()
