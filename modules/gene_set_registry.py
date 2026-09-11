"""Controlled, versioned gene-set resources for functional-state scoring.

The registry is deliberately read-only at analysis time.  An administrator
uses ``scripts/sync_managed_gene_sets.py`` to fetch and freeze source files
under ``Config.functional_state_resource_dir()/gene_sets``.  A web request can
select terms, but can never supply a server path or a download URL.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


REGISTRY_FILENAME = "gene_set_registry.json"
MAX_SELECTED_TERMS = 60
_SAFE_RELATIVE_PATH = re.compile(r"[A-Za-z0-9_.-]+(?:/[A-Za-z0-9_.-]+)*")


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def managed_gene_set_root():
    """Return the administrator-controlled root for frozen gene-set files."""
    from config import Config

    return (Path(Config.functional_state_resource_dir()) / "gene_sets").resolve()


def _safe_resource_file(root, relative_path, label):
    relative_path = str(relative_path or "").strip()
    if not _SAFE_RELATIVE_PATH.fullmatch(relative_path):
        raise ValueError(f"托管基因集 {label} 的文件路径无效。")
    candidate = root / relative_path
    if candidate.is_symlink():
        raise ValueError(f"托管基因集 {label} 不能是符号链接。")
    resolved = candidate.resolve()
    if root not in resolved.parents or not resolved.is_file():
        raise ValueError(f"托管基因集 {label} 必须位于平台资源目录内。")
    return resolved


def load_managed_registry():
    """Load and minimally validate the frozen, administrator-written registry."""
    root = managed_gene_set_root()
    registry_path = root / REGISTRY_FILENAME
    if registry_path.is_symlink():
        raise ValueError("托管基因集注册表不能是符号链接。")
    if not registry_path.is_file():
        raise FileNotFoundError(
            "托管基因集注册表尚未初始化；请由管理员运行 scripts/sync_managed_gene_sets.py。"
        )
    try:
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError("托管基因集注册表格式无效。") from exc
    if not isinstance(registry, dict) or registry.get("schema_version") != 1:
        raise ValueError("托管基因集注册表版本不受支持。")
    libraries = registry.get("libraries")
    if not isinstance(libraries, list) or not libraries:
        raise ValueError("托管基因集注册表不包含可用数据库。")
    return root, registry


def _requested_terms(value):
    raw_values = [value] if isinstance(value, str) else list(value or [])
    terms = []
    for raw_value in raw_values:
        for item in re.split(r"[,;\n]+", str(raw_value or "")):
            term = item.strip()
            if term and term not in terms:
                terms.append(term)
    if len(terms) > MAX_SELECTED_TERMS:
        raise ValueError(f"一次最多选择 {MAX_SELECTED_TERMS} 个托管通路。")
    return terms


def _library_records(registry):
    records = []
    seen = set()
    for record in registry["libraries"]:
        if not isinstance(record, dict):
            raise ValueError("托管基因集注册表中的数据库记录无效。")
        key = str(record.get("key", "")).strip()
        filename = str(record.get("file", "")).strip()
        expected_sha = str(record.get("sha256", "")).strip().lower()
        if not key or key in seen or not filename or not re.fullmatch(r"[a-f0-9]{64}", expected_sha):
            raise ValueError("托管基因集注册表中的数据库元数据无效。")
        seen.add(key)
        records.append(record)
    return records


def _read_selected_terms(path, targets):
    found = {}
    with Path(path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] not in targets:
                continue
            genes = tuple(dict.fromkeys(
                gene.strip().upper() for gene in fields[2:] if gene.strip()
            ))
            if genes:
                found[fields[0]] = genes
    return found


def load_selected_managed_gene_sets(value):
    """Return selected, checksum-verified term definitions and provenance.

    ``value`` contains term names only.  The term-to-file association comes
    solely from the server-managed registry; users cannot select files.
    """
    selected = _requested_terms(value)
    if not selected:
        return {}, {}

    root, registry = load_managed_registry()
    expected_terms = set(selected)
    found, provenance = {}, {}
    for record in _library_records(registry):
        key = str(record["key"])
        path = _safe_resource_file(root, record["file"], key)
        actual_sha = _sha256(path)
        if actual_sha != str(record["sha256"]).lower():
            raise ValueError(f"托管基因集 {key} 校验失败；请由管理员重新同步资源。")
        library_terms = _read_selected_terms(path, expected_terms)
        for term, genes in library_terms.items():
            if term in found:
                raise ValueError(f"托管基因集 term 名称不唯一: {term}。")
            found[term] = genes
            provenance[term] = {
                "source": f"managed_{key}",
                "library_key": key,
                "resource": str(record.get("resource", key)),
                "version": str(record.get("version", "unknown")),
                "file": str(record["file"]),
                "sha256": actual_sha,
                "license": str(record.get("license", "")),
                "license_url": str(record.get("license_url", "")),
                "source_url": str(record.get("source_url", "")),
            }
    missing = [term for term in selected if term not in found]
    if missing:
        raise ValueError("托管基因集未找到所选 term: " + "、".join(missing[:8]))
    return {term: found[term] for term in selected}, {term: provenance[term] for term in selected}


def search_managed_gene_set_terms(query, limit=60):
    """Search term titles locally without exposing genes or resource paths."""
    query = str(query or "").strip().lower()
    capped_limit = max(1, min(int(limit), 200))
    if not query:
        return []
    root, registry = load_managed_registry()
    matches = []
    for record in _library_records(registry):
        key = str(record["key"])
        path = _safe_resource_file(root, record["file"], key)
        if _sha256(path) != str(record["sha256"]).lower():
            raise ValueError(f"托管基因集 {key} 校验失败；请由管理员重新同步资源。")
        with path.open(encoding="utf-8") as handle:
            for raw_line in handle:
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) < 3 or query not in fields[0].lower():
                    continue
                matches.append({
                    "term": fields[0], "library": key,
                    "resource": str(record.get("resource", key)),
                    "version": str(record.get("version", "unknown")),
                    "n_genes": len([gene for gene in fields[2:] if gene.strip()]),
                })
    matches.sort(key=lambda item: (item["library"], item["term"]))
    return matches[:capped_limit]
