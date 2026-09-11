#!/usr/bin/env python3
"""Fetch and freeze public Human pathway resources for functional-state scoring.

This is an administrator-only maintenance command.  It uses fixed official
URLs, keeps source checksums and licence metadata, and writes atomically into
the configured local resource directory.  It never runs from a web request.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import tempfile
import urllib.request
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


HALLMARK_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/"
    "h.all.v2026.1.Hs.symbols.gmt"
)
REACTOME_URL = "https://reactome.org/download/current/ReactomePathways.gmt.zip"
GO_OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"
GO_GAF_URL = "https://current.geneontology.org/annotations/gaf/HUMAN-uniprot.gaf.gz"


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(url, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "sc-transcriptomics-platform/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        with tempfile.NamedTemporaryFile("wb", dir=destination.parent, delete=False) as handle:
            temporary = Path(handle.name)
            shutil.copyfileobj(response, handle)
    try:
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _atomic_text(destination, content):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=destination.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(content)
    try:
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _count_gmt_sets(path):
    count = 0
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if len(line.rstrip("\n").split("\t")) >= 3:
                count += 1
    if not count:
        raise ValueError(f"GMT 文件为空或格式无效: {path}")
    return count


def _extract_reactome(source_zip, destination):
    with zipfile.ZipFile(source_zip) as archive:
        entries = [entry for entry in archive.namelist() if entry.lower().endswith(".gmt")]
        if len(entries) != 1:
            raise ValueError("Reactome 下载包没有唯一的 GMT 文件。")
        with archive.open(entries[0]) as source:
            payload = source.read()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("wb", dir=destination.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(payload)
    try:
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _reactome_release_label(source_zip):
    """Use the GMT member's published archive date, not the sync date, as version."""
    with zipfile.ZipFile(source_zip) as archive:
        entries = [entry for entry in archive.infolist() if entry.filename.lower().endswith(".gmt")]
        if len(entries) != 1:
            raise ValueError("Reactome 下载包没有唯一的 GMT 文件。")
        return "ReactomePathways.gmt_" + "-".join(f"{item:02d}" for item in entries[0].date_time[:3])


def _read_go_terms(obo_path):
    terms, current = {}, None
    with Path(obo_path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                if current and current.get("id") and current.get("name") and not current.get("obsolete"):
                    terms[current["id"]] = current
                current = {"obsolete": False}
                continue
            if not line or line.startswith("["):
                if current and current.get("id") and current.get("name") and not current.get("obsolete"):
                    terms[current["id"]] = current
                current = None
                continue
            if current is None:
                continue
            if line.startswith("id: GO:"):
                current["id"] = line[4:].strip()
            elif line.startswith("name: "):
                current["name"] = line[6:].strip()
            elif line.startswith("namespace: "):
                current["namespace"] = line[len("namespace: "):].strip()
            elif line == "is_obsolete: true":
                current["obsolete"] = True
    if current and current.get("id") and current.get("name") and not current.get("obsolete"):
        terms[current["id"]] = current
    return terms


def _go_release_label(obo_path, gaf_path):
    """Return the ontology release plus the human annotation generation date."""
    ontology_release, annotation_date = "unknown", "unknown"
    with Path(obo_path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("data-version: releases/"):
                ontology_release = line.rsplit("/", 1)[-1]
                break
    with gzip.open(gaf_path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.startswith("!"):
                break
            if raw_line.startswith("!date-generated:"):
                annotation_date = raw_line.split(":", 1)[1].strip().split(" ", 1)[0]
    if ontology_release == "unknown":
        raise ValueError("GO OBO 文件缺少 data-version 发布标识。")
    return ontology_release, annotation_date


def _build_go_gmts(obo_path, gaf_path, output_dir, release_label):
    terms = _read_go_terms(obo_path)
    annotations = defaultdict(set)
    with gzip.open(gaf_path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            if raw_line.startswith("!"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 12 or "NOT" in fields[3].split("|"):
                continue
            term_id, symbol = fields[4].strip(), fields[2].strip().upper()
            if term_id in terms and symbol:
                annotations[term_id].add(symbol)

    namespaces = {
        "biological_process": ("go_bp_human_" + release_label + ".gmt", "GO Biological Process"),
        "molecular_function": ("go_mf_human_" + release_label + ".gmt", "GO Molecular Function"),
        "cellular_component": ("go_cc_human_" + release_label + ".gmt", "GO Cellular Component"),
    }
    records = []
    for namespace, (filename, resource) in namespaces.items():
        lines = []
        for term_id, genes in annotations.items():
            term = terms[term_id]
            if term.get("namespace") != namespace or not genes:
                continue
            label = f"{term_id}|{term['name']}"
            description = f"{resource}; direct HUMAN-uniprot GAF annotations; {release_label}"
            lines.append("\t".join([label, description, *sorted(genes)]))
        lines.sort()
        destination = output_dir / "go" / filename
        _atomic_text(destination, "\n".join(lines) + "\n")
        records.append({
            "key": {"biological_process": "go_bp_human", "molecular_function": "go_mf_human", "cellular_component": "go_cc_human"}[namespace],
            "resource": resource,
            "version": release_label + "_direct_gaf",
            "file": str(destination.relative_to(output_dir)),
            "sha256": _sha256(destination),
            "n_gene_sets": _count_gmt_sets(destination),
            "source_url": GO_GAF_URL,
            "ontology_source_url": GO_OBO_URL,
            "license": "CC BY 4.0",
            "license_url": "https://geneontology.org/docs/go-citation-policy/",
            "annotation_policy": "direct human GAF annotations only; no ancestor propagation",
        })
    return records


def sync(resource_dir):
    root = Path(resource_dir).expanduser().resolve()
    if root.is_symlink():
        raise ValueError("资源目录不能是符号链接。")
    output_dir = root / "gene_sets"
    source_dir = output_dir / "sources"
    source_dir.mkdir(parents=True, exist_ok=True)

    hallmark = output_dir / "hallmark" / "h.all.v2026.1.Hs.symbols.gmt"
    reactome_zip = source_dir / "ReactomePathways.current.gmt.zip"
    reactome = output_dir / "reactome" / "reactome_human_current.gmt"
    go_obo = source_dir / "go-basic.current.obo"
    go_gaf = source_dir / "HUMAN-uniprot.current.gaf.gz"

    _download(HALLMARK_URL, hallmark)
    _download(REACTOME_URL, reactome_zip)
    _extract_reactome(reactome_zip, reactome)
    _download(GO_OBO_URL, go_obo)
    _download(GO_GAF_URL, go_gaf)

    go_release, go_annotation_date = _go_release_label(go_obo, go_gaf)
    reactome_release = _reactome_release_label(reactome_zip)
    libraries = [
        {
            "key": "hallmark_human",
            "resource": "MSigDB Hallmark",
            "version": "2026.1.Hs",
            "file": str(hallmark.relative_to(output_dir)),
            "sha256": _sha256(hallmark),
            "n_gene_sets": _count_gmt_sets(hallmark),
            "source_url": HALLMARK_URL,
            "license": "CC BY 4.0 (per-gene-set terms retained by MSigDB)",
            "license_url": "https://www.gsea-msigdb.org/gsea/msigdb/index.jsp",
        },
        {
            "key": "reactome_human",
            "resource": "Reactome Pathways",
            "version": reactome_release,
            "file": str(reactome.relative_to(output_dir)),
            "sha256": _sha256(reactome),
            "n_gene_sets": _count_gmt_sets(reactome),
            "source_url": REACTOME_URL,
            "license": "CC0 1.0",
            "license_url": "https://reactome.org/license",
        },
    ]
    libraries.extend(_build_go_gmts(go_obo, go_gaf, output_dir, go_release))
    registry = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "organism": "Human",
        "libraries": libraries,
        "source_artifacts": {
            "reactome_zip": {"file": str(reactome_zip.relative_to(output_dir)), "sha256": _sha256(reactome_zip)},
            "go_basic_obo": {"file": str(go_obo.relative_to(output_dir)), "sha256": _sha256(go_obo)},
            "human_gaf": {"file": str(go_gaf.relative_to(output_dir)), "sha256": _sha256(go_gaf), "date_generated": go_annotation_date},
        },
    }
    _atomic_text(output_dir / "gene_set_registry.json", json.dumps(registry, ensure_ascii=False, indent=2) + "\n")
    return registry


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--resource-dir", default=os.environ.get("FUNCTIONAL_STATE_RESOURCE_DIR", "data/functional_state_resources"),
        help="administrator-controlled functional-state resource directory",
    )
    args = parser.parse_args()
    registry = sync(args.resource_dir)
    print(json.dumps({
        "resource_dir": str(Path(args.resource_dir).resolve()),
        "libraries": [{key: item[key] for key in ("key", "version", "n_gene_sets", "sha256")} for item in registry["libraries"]],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
