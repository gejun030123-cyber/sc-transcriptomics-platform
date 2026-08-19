"""Dependency health checks for platform capabilities."""

import importlib.util


DEPENDENCY_GROUPS = {
    "core": [
        ("flask", "flask"),
        ("flask-cors", "flask_cors"),
        ("requests", "requests"),
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("scipy", "scipy"),
        ("plotly", "plotly"),
        ("psutil", "psutil"),
        ("anndata", "anndata"),
    ],
    "single_cell": [
        ("scanpy", "scanpy"),
        ("omicverse", "omicverse"),
        ("statsmodels", "statsmodels"),
        ("scikit-learn", "sklearn"),
        ("matplotlib", "matplotlib"),
    ],
    "bulk": [
        ("gseapy", "gseapy"),
        ("pydeseq2", "pydeseq2"),
        ("inmoose", "inmoose"),
        ("patsy", "patsy"),
    ],
    "optional": [
        ("liana", "liana"),
        ("harmonypy", "harmonypy"),
        ("bbknn", "bbknn"),
        ("scanorama", "scanorama"),
        ("scvi-tools", "scvi"),
        ("torch", "torch"),
        ("kneed", "kneed"),
        ("kaleido", "kaleido"),
        ("celltypist", "celltypist"),
        ("openai", "openai"),
        ("anthropic", "anthropic"),
    ],
}

MODULE_DEPENDENCIES = {
    "qc": ["scanpy", "omicverse"],
    "normalize": ["scanpy", "omicverse"],
    "hvg": ["scanpy", "omicverse"],
    "dimred": ["scanpy"],
    "batch_correct": ["scanpy"],
    "clustering": ["scanpy"],
    "subcluster": ["scanpy"],
    "qc_reassess": ["scanpy", "pandas"],
    "annotation": ["scanpy"],
    "deg": ["scanpy"],
    "trajectory": ["scanpy"],
    "sc_timecourse": ["anndata", "pandas", "scipy"],
    "proportion": ["anndata", "pandas", "scipy"],
    "cell_communication": ["scanpy", "liana"],
    "sc_batch_import": ["anndata", "pandas", "scanpy"],
    "sc_cell_deg": ["scanpy", "pandas"],
    "sc_cell_go": ["pandas", "gseapy"],
    "sc_pseudobulk_deg": ["anndata", "pandas", "scipy", "statsmodels"],
    "sc_csv_export": ["anndata", "pandas"],
    "bulk_qc": ["anndata", "pandas"],
    "bulk_normalize": ["numpy", "pandas"],
    "bulk_pca": ["anndata", "pandas", "scikit-learn"],
    "bulk_deg": ["pandas", "scipy", "statsmodels", "omicverse"],
    "bulk_heatmap": ["anndata", "pandas", "scanpy", "scipy"],
    # The current implementation uses OmicVerse and local/Enrichr gene-set
    # files.  gseapy is used by subcluster's optional enrichment path, not by
    # the core Bulk enrichment module.
    "bulk_enrichment": ["omicverse"],
    "bulk_timecourse": ["pandas", "scipy", "patsy"],
    "bulk_deg_integration": ["pandas"],
    "convert_10x": ["scanpy"],
}

MODULE_OPTIONAL_DEPENDENCIES = {
    "dimred": [("kneed", "Kneedle 自动 PC 选择")],
    "batch_correct": [
        ("harmonypy", "Harmony 批次校正"),
        ("bbknn", "BBKNN 批次校正"),
        ("scanorama", "Scanorama 批次校正"),
        ("scvi-tools", "scVI/SysVI 深度整合"),
    ],
    "subcluster": [("gseapy", "子簇通路富集")],
    "annotation": [("celltypist", "CellTypist 参考交叉验证")],
    "bulk_deg": [("inmoose", "DESeq2/edgeR/limma 兼容统计方法")],
    "bulk_enrichment": [("gseapy", "兼容旧版 Enrichr 富集路径")],
    "sc_pseudobulk_deg": [
        ("pydeseq2", "DESeq2 pseudobulk 统计；缺失时回退 Welch log2CPM"),
    ],
}


def _installed(import_name):
    return importlib.util.find_spec(import_name) is not None


def dependency_status():
    packages = {}
    groups = {}
    for group, deps in DEPENDENCY_GROUPS.items():
        entries = []
        for package, import_name in deps:
            ok = _installed(import_name)
            packages[package] = ok
            entries.append({
                "package": package,
                "import_name": import_name,
                "installed": ok,
            })
        groups[group] = entries

    modules = {}
    for module_name, deps in MODULE_DEPENDENCIES.items():
        missing = [dep for dep in deps if not packages.get(dep, _installed(dep.replace("-", "_")))]
        optional_missing = []
        for dep, capability in MODULE_OPTIONAL_DEPENDENCIES.get(module_name, []):
            if not packages.get(dep, _installed(dep.replace("-", "_"))):
                optional_missing.append({"package": dep, "capability": capability})
        modules[module_name] = {
            "available": len(missing) == 0,
            "missing": missing,
            "optional_missing": optional_missing,
        }

    return {
        "groups": groups,
        "modules": modules,
        "summary": {
            "total": sum(len(v) for v in DEPENDENCY_GROUPS.values()),
            "missing": sum(1 for ok in packages.values() if not ok),
        },
    }
