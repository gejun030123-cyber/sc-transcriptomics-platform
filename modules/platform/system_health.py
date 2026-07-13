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
        ("celltypist", "celltypist"),
        ("liana", "liana"),
        ("harmonypy", "harmonypy"),
        ("bbknn", "bbknn"),
        ("scanorama", "scanorama"),
        ("scvi-tools", "scvi"),
        ("torch", "torch"),
        ("kneed", "kneed"),
        ("kaleido", "kaleido"),
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
    "annotation": ["scanpy"],
    "deg": ["scanpy"],
    "cell_communication": ["liana"],
    "bulk_qc": ["anndata", "pandas"],
    "bulk_normalize": ["numpy", "pandas"],
    "bulk_deg": ["scipy", "statsmodels"],
    "bulk_enrichment": ["gseapy"],
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
        modules[module_name] = {
            "available": len(missing) == 0,
            "missing": missing,
        }

    return {
        "groups": groups,
        "modules": modules,
        "summary": {
            "total": sum(len(v) for v in DEPENDENCY_GROUPS.values()),
            "missing": sum(1 for ok in packages.values() if not ok),
        },
    }
