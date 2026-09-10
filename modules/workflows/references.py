"""Small, checksum-aware WES reference and capture-kit catalog.

The catalog intentionally reuses ``reference_assets``.  A capture-kit profile
is a group of rows whose metadata share ``capture_kit_id``; this is sufficient
for an internal lab and avoids a second hierarchy of mutable database tables.
"""

import gzip
import hashlib
import json
import os
import uuid
from typing import Any, Dict, Iterable, Mapping, Optional

from config import Config
from database import get_conn


REFERENCE_STATUSES = {"registered", "test_only", "validated", "retired"}
CAPTURE_ASSET_TYPES = {
    "calling": "capture_bed_calling",
    "vendor": "capture_bed_vendor",
    "qc": "capture_bed_qc",
}
_BUNDLE_REQUIREMENTS = {
    "wes_germline": {
        "fasta", "fai", "dict", "dbsnp", "dbsnp_tbi",
        "germline_resource", "germline_resource_tbi",
    },
    "wes_somatic": {
        "fasta", "fai", "dict", "pon", "pon_tbi",
        "germline_resource", "germline_resource_tbi",
    },
    "wes_annotate_only": {"vep_cache_manifest"},
}

# A bundle normally uses the deployment-wide iGenomes configuration.  A
# legacy reference such as hs37d5 needs a separately pinned custom-reference
# contract instead: its complete, catalogued file set is rendered into the
# Sarek params file for that run only.  Keeping this declaration on the FASTA
# asset avoids a mutable second reference registry and leaves existing GRCh38
# runs completely unchanged.
_SAREK_REFERENCE_METADATA_KEY = "sarek_reference"
_SAREK_CUSTOM_PARAM_BY_ASSET_TYPE = {
    "fasta": "fasta",
    "fai": "fasta_fai",
    "dict": "dict",
    "dbsnp": "dbsnp",
    "dbsnp_tbi": "dbsnp_tbi",
    "germline_resource": "germline_resource",
    "germline_resource_tbi": "germline_resource_tbi",
    "pon": "pon",
    "pon_tbi": "pon_tbi",
}


def _checksum(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_reference_path(file_path: str, source_roots: Iterable[str] = ()) -> str:
    if not file_path:
        raise ValueError("缺少 reference asset 文件路径")
    supplied = os.path.abspath(str(file_path))
    if os.path.islink(supplied):
        raise ValueError("reference asset 不允许符号链接")
    resolved = os.path.realpath(supplied)
    roots = [os.path.realpath(root) for root in (source_roots or Config.wes_source_roots())]
    if not roots:
        raise ValueError("未配置 WES_SOURCE_ROOTS，不能登记服务器 reference asset")
    if not any(resolved == root or resolved.startswith(root + os.sep) for root in roots):
        raise ValueError("reference asset 不在管理员配置的 WES_SOURCE_ROOTS 下")
    if not os.path.isfile(resolved):
        raise ValueError("reference asset 文件不存在或不是普通文件")
    return resolved


def register_reference_asset(*, species: str = "human", assembly: str,
                             bundle_version: str, asset_type: str, file_path: str,
                             checksum: str = "", source: str = "", license_note: str = "",
                             metadata: Optional[Dict[str, Any]] = None,
                             status: str = "registered",
                             source_roots: Iterable[str] = ()) -> Dict[str, Any]:
    if status not in REFERENCE_STATUSES:
        raise ValueError(f"未知 reference asset 状态: {status}")
    resolved = _resolve_reference_path(file_path, source_roots)
    actual_checksum = _checksum(resolved)
    if checksum and checksum.lower() != actual_checksum:
        raise ValueError("reference asset 提供的 checksum 与文件内容不一致")
    encoded_metadata = json.dumps(metadata or {}, ensure_ascii=False, sort_keys=True)
    conn = get_conn()
    try:
        existing = conn.execute(
            "SELECT id FROM reference_assets WHERE assembly=? AND bundle_version=? "
            "AND asset_type=? AND file_path=? AND checksum=? AND metadata_json=?",
            (assembly, bundle_version, asset_type, resolved, actual_checksum, encoded_metadata),
        ).fetchone()
    finally:
        conn.close()
    if existing:
        return get_reference_asset(existing["id"])
    asset_id = "ref_" + uuid.uuid4().hex[:16]
    conn = get_conn()
    try:
        conn.execute(
            "INSERT INTO reference_assets "
            "(id, species, assembly, bundle_version, asset_type, file_path, checksum, source, "
            "license_note, metadata_json, status) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (asset_id, species, assembly, bundle_version, asset_type, resolved,
             actual_checksum, source, license_note, encoded_metadata, status),
        )
        conn.commit()
    finally:
        conn.close()
    return get_reference_asset(asset_id)


def _decode(row):
    if not row:
        return None
    data = dict(row)
    try:
        data["metadata"] = json.loads(data.pop("metadata_json") or "{}")
    except (TypeError, ValueError, json.JSONDecodeError):
        data["metadata"] = {}
        data.pop("metadata_json", None)
    return data


def get_reference_asset(asset_id: str):
    conn = get_conn()
    try:
        row = conn.execute("SELECT * FROM reference_assets WHERE id=?", (asset_id,)).fetchone()
    finally:
        conn.close()
    return _decode(row)


def list_reference_assets(*, assembly: str = "", bundle_version: str = "",
                          status: str = "", asset_type: str = ""):
    conn = get_conn()
    try:
        query = "SELECT * FROM reference_assets WHERE 1=1"
        params = []
        if assembly:
            query += " AND assembly=?"
            params.append(assembly)
        if bundle_version:
            query += " AND bundle_version=?"
            params.append(bundle_version)
        if status:
            query += " AND status=?"
            params.append(status)
        if asset_type:
            query += " AND asset_type=?"
            params.append(asset_type)
        query += " ORDER BY assembly, bundle_version, asset_type, id"
        rows = conn.execute(query, tuple(params)).fetchall()
    finally:
        conn.close()
    return [_decode(row) for row in rows]


def set_reference_asset_status(asset_id: str, status: str) -> Dict[str, Any]:
    """Promote/retire one administrator-managed immutable asset."""
    if status not in REFERENCE_STATUSES:
        raise ValueError(f"未知 reference asset 状态: {status}")
    asset = get_reference_asset(asset_id)
    if not asset:
        raise ValueError("reference asset 不存在")
    if status == "validated":
        check = verify_reference_asset(asset)
        if not check["valid"]:
            raise ValueError("reference asset 当前校验失败: " + "; ".join(check["errors"]))
        if check.get("bed") and not check["bed"].get("sorted"):
            raise ValueError("capture BED 未排序，不能提升为 validated")
    conn = get_conn()
    try:
        cursor = conn.execute(
            "UPDATE reference_assets SET status=? WHERE id=?", (status, asset_id)
        )
        updated = cursor.rowcount
        conn.commit()
    finally:
        conn.close()
    if not updated:
        raise ValueError("reference asset 不存在")
    return get_reference_asset(asset_id)


def inspect_bed(path: str, *, max_errors: int = 20) -> Dict[str, Any]:
    """Stream a BED and report structural/contig checks without loading it all."""
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    interval_count = 0
    target_bases = 0
    contigs = set()
    errors = []
    previous = None

    def _contig_sort_key(contig: str):
        """Use natural chromosome order (chr9 before chr10), then stable text order."""
        name = contig[3:] if contig.startswith("chr") else contig
        if name.isdigit():
            return (0, int(name))
        special = {"X": 1000, "Y": 1001, "M": 1002, "MT": 1002}
        if name in special:
            return (0, special[name])
        return (1, name)

    sorted_bed = True
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.strip()
                if not line or line.startswith(("#", "track ", "browser ")):
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    if len(errors) < max_errors:
                        errors.append(f"line {line_number}: 少于 3 列")
                    continue
                contig = fields[0].strip()
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    if len(errors) < max_errors:
                        errors.append(f"line {line_number}: start/end 不是整数")
                    continue
                if not contig or start < 0 or end <= start:
                    if len(errors) < max_errors:
                        errors.append(f"line {line_number}: 区间无效")
                    continue
                current = (_contig_sort_key(contig), start, end)
                if previous and current < previous:
                    sorted_bed = False
                previous = current
                contigs.add(contig)
                interval_count += 1
                target_bases += end - start
    except (OSError, UnicodeError) as exc:
        errors.append(f"BED 读取失败: {exc}")
    styles = {"chr" if name.startswith("chr") else "bare" for name in contigs}
    if len(styles) > 1:
        errors.append("BED 混用了 chr 和非 chr contig 命名")
    if not interval_count:
        errors.append("BED 不包含有效区间")
    return {
        "valid": not errors,
        "interval_count": interval_count,
        "target_bases": target_bases,
        "contig_count": len(contigs),
        "contigs": sorted(contigs),
        "contig_style": next(iter(styles), "unknown") if len(styles) <= 1 else "mixed",
        "sorted": sorted_bed,
        "errors": errors,
    }


def verify_reference_asset(asset_or_id: str | Mapping[str, Any],
                           *, source_roots: Iterable[str] = ()) -> Dict[str, Any]:
    asset = (get_reference_asset(asset_or_id) if isinstance(asset_or_id, str)
             else dict(asset_or_id))
    if not asset:
        return {"valid": False, "errors": ["reference asset 不存在"]}
    errors = []
    try:
        path = _resolve_reference_path(asset.get("file_path", ""), source_roots)
    except ValueError as exc:
        return {"valid": False, "asset_id": asset.get("id", ""), "errors": [str(exc)]}
    actual = _checksum(path)
    if actual != asset.get("checksum"):
        errors.append("文件 checksum 已变化")
    details: Dict[str, Any] = {}
    if asset.get("asset_type") in CAPTURE_ASSET_TYPES.values():
        details["bed"] = inspect_bed(path)
        errors.extend(details["bed"]["errors"])
    return {
        "valid": not errors,
        "asset_id": asset.get("id", ""),
        "status": asset.get("status", ""),
        "path": path,
        "checksum": actual,
        "errors": errors,
        **details,
    }


def register_capture_kit_profile(*, capture_kit_id: str, name: str, version: str,
                                 assembly: str, calling_bed_path: str,
                                 vendor_bed_path: str = "", qc_bed_path: str = "",
                                 status: str = "test_only", source: str = "",
                                 license_note: str = "", padding_bp: int = 0,
                                 metadata: Optional[Dict[str, Any]] = None,
                                 source_roots: Iterable[str] = ()) -> Dict[str, Any]:
    """Register one kit profile; only the calling BED is mandatory."""
    capture_kit_id = str(capture_kit_id or "").strip()
    if not capture_kit_id or not name or not version:
        raise ValueError("capture kit 需要稳定 ID、名称和版本")
    if status not in {"test_only", "validated"}:
        raise ValueError("capture kit 初始状态只能是 test_only 或 validated")
    paths = {
        "calling": calling_bed_path,
        "vendor": vendor_bed_path,
        "qc": qc_bed_path,
    }
    assets = []
    for role, path in paths.items():
        if not path:
            continue
        bed = inspect_bed(_resolve_reference_path(path, source_roots))
        if not bed["valid"]:
            raise ValueError(f"{role} BED 校验失败: {'; '.join(bed['errors'])}")
        item_metadata = {
            **(metadata or {}),
            "capture_kit_id": capture_kit_id,
            "capture_kit_name": name,
            "capture_kit_version": version,
            "bed_role": role,
            "padding_bp": int(padding_bp) if role == "calling" else 0,
            "bed_summary": {key: bed[key] for key in (
                "interval_count", "target_bases", "contig_count", "contig_style", "sorted"
            )},
        }
        assets.append(register_reference_asset(
            assembly=assembly,
            bundle_version=f"capture:{capture_kit_id}:{version}",
            asset_type=CAPTURE_ASSET_TYPES[role], file_path=path,
            source=source, license_note=license_note, metadata=item_metadata,
            status=status, source_roots=source_roots,
        ))
    return get_capture_kit_profile(capture_kit_id)


def list_capture_kit_profiles(*, assembly: str = "", include_retired: bool = False):
    grouped: Dict[str, Dict[str, Any]] = {}
    for asset in list_reference_assets(assembly=assembly):
        metadata = asset.get("metadata") or {}
        kit_id = str(metadata.get("capture_kit_id") or "")
        role = str(metadata.get("bed_role") or "")
        if not kit_id or role not in CAPTURE_ASSET_TYPES:
            continue
        if not include_retired and asset.get("status") == "retired":
            continue
        profile = grouped.setdefault(kit_id, {
            "capture_kit_id": kit_id,
            "name": metadata.get("capture_kit_name", kit_id),
            "version": metadata.get("capture_kit_version", ""),
            "assembly": asset.get("assembly", ""),
            "assets": {},
        })
        profile["assets"][role] = asset
    profiles = []
    for profile in grouped.values():
        calling = profile["assets"].get("calling")
        profile["status"] = calling.get("status") if calling else "incomplete"
        profile["production_allowed"] = bool(calling and calling.get("status") == "validated")
        profiles.append(profile)
    return sorted(profiles, key=lambda item: (item["name"], item["version"], item["capture_kit_id"]))


def get_capture_kit_profile(capture_kit_id: str, *, include_retired: bool = False):
    """Return the active profile by default; retired assets must not be launch inputs."""
    return next((item for item in list_capture_kit_profiles(include_retired=include_retired)
                 if item["capture_kit_id"] == capture_kit_id), None)


def _resolve_capture_calling_asset(capture_bed_id: str):
    direct = get_reference_asset(capture_bed_id)
    if direct and direct.get("asset_type") == CAPTURE_ASSET_TYPES["calling"]:
        return direct
    profile = get_capture_kit_profile(capture_bed_id)
    return (profile.get("assets") or {}).get("calling") if profile else None


def reference_bundle_sarek_parameters(bundle_version: str, workflow_key: str) -> Dict[str, Any]:
    """Return immutable, catalog-derived Sarek reference parameters.

    Most deployed bundles continue to inherit the existing global iGenomes
    settings and therefore return an empty dict.  A custom bundle opts in by
    placing ``{"mode": "custom"}`` under ``sarek_reference`` on its FASTA
    asset.  All paths then come from checksum-tracked catalog assets, never
    from the browser or a manifest.
    """
    assets = list_reference_assets(bundle_version=str(bundle_version or ""), status="validated")
    fasta = next((asset for asset in assets if asset.get("asset_type") == "fasta"), None)
    declaration = dict(((fasta or {}).get("metadata") or {}).get(
        _SAREK_REFERENCE_METADATA_KEY
    ) or {})
    if not declaration:
        return {"configured": False, "parameters": {}, "errors": []}
    if declaration.get("mode") != "custom":
        return {
            "configured": True,
            "parameters": {},
            "errors": ["sarek_reference.mode 仅支持 custom"],
        }

    required = _BUNDLE_REQUIREMENTS.get(workflow_key, set())
    by_type = {asset.get("asset_type"): asset for asset in assets}
    missing = sorted(required - set(by_type))
    errors = []
    if missing:
        errors.append("custom Sarek bundle 缺少 validated assets: " + ", ".join(missing))
    params: Dict[str, Any] = {
        # Sarek documents this combination for a fully local custom
        # reference.  ``None`` is intentionally serialized as JSON null.
        "genome": None,
        "igenomes_ignore": True,
    }
    for asset_type, parameter in _SAREK_CUSTOM_PARAM_BY_ASSET_TYPE.items():
        asset = by_type.get(asset_type)
        if asset and asset_type in required:
            try:
                params[parameter] = Config.validate_wes_reference_path(asset.get("file_path", ""))
            except ValueError as exc:
                errors.append(f"{asset_type} 路径无效: {exc}")
    return {"configured": True, "parameters": params, "errors": errors}


def reference_bundle_readiness(bundle_version: str, workflow_key: str) -> Dict[str, Any]:
    required = _BUNDLE_REQUIREMENTS.get(workflow_key, set())
    assets = list_reference_assets(bundle_version=bundle_version, status="validated")
    present = {asset["asset_type"] for asset in assets}
    missing = sorted(required - present)
    checks = [verify_reference_asset(asset) for asset in assets if asset["asset_type"] in required]
    errors = [error for check in checks for error in check.get("errors", [])]
    if missing:
        errors.append("reference bundle 缺少 validated assets: " + ", ".join(missing))
    assemblies = {asset.get("assembly") for asset in assets if asset.get("asset_type") in required}
    if len(assemblies) > 1:
        errors.append("reference bundle 混用了 assembly")
    by_type = {asset["asset_type"]: asset for asset in assets}
    sarek_reference = reference_bundle_sarek_parameters(bundle_version, workflow_key)
    errors.extend(sarek_reference["errors"])
    configured_base = str(Config.WES_NEXTFLOW_IGENOMES_BASE or "").strip()
    # Only enforce the deployment base when it is itself inside the currently
    # approved WES roots.  Test/admin bundles may deliberately use a temporary
    # source root (for example during a catalog drill); in that case the
    # deployment's production base must not make the isolated bundle invalid.
    configured_roots = Config.wes_source_roots()
    base_is_approved = False
    if configured_base:
        base = os.path.realpath(configured_base)
        base_is_approved = any(
            base == root or base.startswith(root + os.sep)
            for root in configured_roots
        )
    if configured_base and base_is_approved:
        base = os.path.realpath(configured_base)
        for asset_type in required - {"vep_cache_manifest"}:
            asset = by_type.get(asset_type)
            path = os.path.realpath((asset or {}).get("file_path", ""))
            if asset and not (path == base or path.startswith(base + os.sep)):
                errors.append(f"{asset_type} 不在 WES_NEXTFLOW_IGENOMES_BASE 下")
    if workflow_key == "wes_somatic":
        configured = {
            "pon": Config.WES_NEXTFLOW_PON,
            "germline_resource": Config.WES_NEXTFLOW_GERMLINE_RESOURCE,
        }
        for asset_type, configured_path in configured.items():
            asset = by_type.get(asset_type)
            if configured_path and asset and os.path.realpath(configured_path) != os.path.realpath(asset["file_path"]):
                errors.append(f"WES_NEXTFLOW_{asset_type.upper()} 与 catalog 路径不一致")
    if workflow_key == "wes_annotate_only" and Config.WES_NEXTFLOW_VEP_CACHE:
        manifest = by_type.get("vep_cache_manifest")
        if manifest:
            cache_root = str((manifest.get("metadata") or {}).get("cache_root") or "")
            if cache_root and os.path.realpath(cache_root) != os.path.realpath(Config.WES_NEXTFLOW_VEP_CACHE):
                errors.append("WES_NEXTFLOW_VEP_CACHE 与 catalog manifest 不一致")
    return {
        "bundle_version": bundle_version,
        "workflow_key": workflow_key,
        "valid": not errors,
        "required_asset_types": sorted(required),
        "present_asset_types": sorted(present),
        "assembly": next(iter(assemblies), "") if len(assemblies) <= 1 else "mixed",
        "checks": checks,
        "sarek_reference": sarek_reference,
        "errors": errors,
    }


def launch_reference_readiness(*, reference_bundle_id: str, capture_bed_id: str,
                               capture_bed_path: str, workflow_key: str) -> Dict[str, Any]:
    """Hard production gate tying a manifest to immutable catalog rows."""
    bundle = reference_bundle_readiness(reference_bundle_id, workflow_key)
    errors = list(bundle["errors"])
    capture = None
    if workflow_key != "wes_annotate_only":
        capture = _resolve_capture_calling_asset(capture_bed_id)
        if not capture:
            errors.append("capture_bed_id 未对应已登记的 calling BED/profile")
        else:
            capture_check = verify_reference_asset(capture)
            errors.extend(capture_check["errors"])
            if capture.get("status") != "validated":
                errors.append("capture BED 尚未提升为 validated")
            if os.path.realpath(capture.get("file_path", "")) != os.path.realpath(capture_bed_path or ""):
                errors.append("manifest capture_bed_path 与 catalog 登记路径不一致")
            bundle_assembly = bundle.get("assembly")
            if bundle_assembly and capture.get("assembly") != bundle_assembly:
                errors.append("capture BED assembly 与 reference bundle 不一致")
            fai = next((asset for asset in list_reference_assets(
                bundle_version=reference_bundle_id, status="validated", asset_type="fai"
            )), None)
            if fai and capture_check.get("bed"):
                fai_contigs = set()
                try:
                    with open(fai["file_path"], encoding="utf-8") as handle:
                        for raw in handle:
                            fields = raw.rstrip("\n").split("\t")
                            if len(fields) >= 2 and fields[1].isdigit():
                                fai_contigs.add(fields[0])
                except OSError as exc:
                    errors.append(f"FAI contig 检查失败: {exc}")
                bed_contigs = set(capture_check["bed"].get("contigs") or ())
                if not fai_contigs:
                    errors.append("validated FAI 不包含可解析的 contig")
                elif not bed_contigs.issubset(fai_contigs):
                    errors.append("capture BED contig 不属于 reference FAI")
    return {
        "valid": not errors,
        "reference_bundle": bundle,
        "capture_asset": capture,
        "errors": errors,
    }
