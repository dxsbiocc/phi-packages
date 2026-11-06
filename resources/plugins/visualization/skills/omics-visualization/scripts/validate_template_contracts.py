#!/usr/bin/env python3
"""Validate template input contracts against the omics visualization catalog."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_template_gallery import canonical_text, gallery_document
from viz_catalog import CONTRACTS, ROLE_TYPES, YamlSubsetError, load_families, load_yaml_subset, normalize_name

SKILL_ROOT = Path(__file__).resolve().parents[1]
CATALOG_TEXT_KEYS = ("title", "use_when", "avoid_when", "input_shape")
TEMPLATE_ID = re.compile(r"^# Template-ID:\s*(\S+)\s*$", re.MULTILINE)
VARIANT_PREVIEW = re.compile(r"^preview_[A-Za-z0-9_-]+\.png$")


def catalog_records(errors: list[str]) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    try:
        families = load_families()
    except (OSError, YamlSubsetError) as exc:
        errors.append(f"references/plots.yaml: {exc}")
        return records
    for family in families:
        path = family["catalog"]
        if not path:
            errors.append(f"family {family['family']}: no catalog path in plots.yaml")
            continue
        try:
            document = load_yaml_subset(str(SKILL_ROOT / path))
        except (OSError, YamlSubsetError) as exc:
            errors.append(f"{path}: {exc}")
            continue
        templates = document.get("templates") if isinstance(document, dict) else None
        if not isinstance(templates, list):
            errors.append(f"{path}: no templates list")
            continue
        for record in templates:
            if not isinstance(record, dict) or not isinstance(record.get("id"), str):
                errors.append(f"{path}: a template has no id")
                continue
            if record["id"] in records:
                errors.append(f"{record['id']}: listed twice in the catalogs")
            records[record["id"]] = {**record, "family": family["family"], "catalog": path}
    return records


def check_catalog_record(template_id: str, record: dict[str, Any], errors: list[str]) -> None:
    for key in CATALOG_TEXT_KEYS:
        if not isinstance(record.get(key), str) or not record[key].strip():
            errors.append(f"{template_id}: catalog {key} is empty")
    source = record.get("source")
    if not isinstance(source, str) or not (SKILL_ROOT / source).is_file():
        errors.append(f"{template_id}: source does not exist: {source}")
    else:
        declared = TEMPLATE_ID.search((SKILL_ROOT / source).read_text(encoding="utf-8"))
        if declared is None or declared.group(1) != template_id:
            errors.append(f"{template_id}: {source} declares Template-ID {declared.group(1) if declared else 'nothing'}")
    preview = record.get("preview")
    if not isinstance(preview, str) or not (SKILL_ROOT / preview).is_file() or Path(preview).name != "preview.png":
        errors.append(f"{template_id}: preview is not a shipped preview.png: {preview}")
    examples = record.get("examples")
    if not isinstance(examples, list) or not examples:
        errors.append(f"{template_id}: catalog lists no example data")
    else:
        for example in examples:
            if not isinstance(example, str) or not (SKILL_ROOT / example).is_file():
                errors.append(f"{template_id}: example does not exist: {example}")
    variants = record.get("variants")
    if variants is None:
        return
    if not isinstance(variants, list):
        errors.append(f"{template_id}: variants must be a list")
        return
    seen: set[str] = set()
    for index, item in enumerate(variants):
        where = f"{template_id} variant {index + 1}"
        if not isinstance(item, dict):
            errors.append(f"{where}: must be a mapping")
            continue
        image = item.get("image")
        if not isinstance(image, str) or not (SKILL_ROOT / image).is_file():
            errors.append(f"{where}: image does not exist: {image}")
        elif not VARIANT_PREVIEW.match(Path(image).name):
            errors.append(f"{where}: image must be preview_<name>.png: {image}")
        elif image in seen:
            errors.append(f"{where}: duplicate image {image}")
        else:
            seen.add(image)
        example = item.get("example")
        if not isinstance(example, str) or not (SKILL_ROOT / example).is_file():
            errors.append(f"{where}: example does not exist: {example}")
        if not isinstance(item.get("config"), str) or not str(item.get("config")).strip():
            errors.append(f"{where}: config is empty")
        if not isinstance(item.get("shows"), str) or not str(item.get("shows")).strip():
            errors.append(f"{where}: shows is empty")


def validate_contracts(contracts_path: Path) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    contracts = json.loads(contracts_path.read_text(encoding="utf-8"))
    if contracts.get("schema_version") != 2:
        errors.append("template_contracts.json schema_version must be 2")
    role_types = contracts.get("role_types")
    if not isinstance(role_types, dict) or not role_types:
        errors.append("role_types must be a non-empty object")
        role_types = {}
    for role, kind in role_types.items():
        if kind not in ROLE_TYPES:
            errors.append(f"role {role}: type {kind!r} is not one of {', '.join(ROLE_TYPES)}")

    records = catalog_records(errors)
    for template_id, record in records.items():
        check_catalog_record(template_id, record, errors)

    seen: set[str] = set()
    used_roles: set[str] = set()
    templates = contracts.get("templates")
    if not isinstance(templates, list) or not templates:
        errors.append("contracts must contain a non-empty templates list")
        templates = []
    for index, template in enumerate(templates):
        if not isinstance(template, dict) or not isinstance(template.get("id"), str):
            errors.append(f"template at index {index} has no id")
            continue
        template_id = template["id"]
        if template_id in seen:
            errors.append(f"{template_id}: duplicate template id")
        seen.add(template_id)
        if template_id not in records:
            errors.append(f"{template_id}: not found in references/catalog")
        unexpected = sorted(set(template) - {"id", "roles", "risks"})
        if unexpected:
            errors.append(f"{template_id}: unexpected keys {', '.join(unexpected)}")
        roles = template.get("roles")
        if not isinstance(roles, dict):
            errors.append(f"{template_id}: roles must be an object")
            continue
        required = roles.get("required")
        optional = roles.get("optional")
        if not isinstance(required, list) or not required:
            errors.append(f"{template_id}: roles.required must be a non-empty list")
            required = []
        if not isinstance(optional, list):
            errors.append(f"{template_id}: roles.optional must be a list")
            optional = []
        unknown = sorted((set(required) | set(optional)) - set(role_types))
        if unknown:
            errors.append(f"{template_id}: roles without a type in role_types: {', '.join(unknown)}")
        both = sorted(set(required) & set(optional))
        if both:
            errors.append(f"{template_id}: roles both required and optional: {', '.join(both)}")
        used_roles.update(required)
        used_roles.update(optional)

    unused = sorted(set(role_types) - used_roles)
    if unused:
        warnings.append(f"role_types lists roles no template uses: {', '.join(unused)}")

    gallery, gallery_errors = gallery_document(SKILL_ROOT)
    errors.extend(gallery_errors)
    gallery_path = SKILL_ROOT / "references" / "template_gallery.json"
    if not gallery_path.is_file():
        errors.append("references/template_gallery.json is missing; run scripts/build_template_gallery.py")
    elif gallery_path.read_text(encoding="utf-8") != canonical_text(gallery):
        errors.append("references/template_gallery.json is stale; run scripts/build_template_gallery.py")

    names: dict[str, list[str]] = {}
    for template_id, record in records.items():
        for alias in record.get("aliases") or []:
            if isinstance(alias, str):
                names.setdefault(normalize_name(alias), []).append(template_id)
    shared = sorted((alias, sorted(set(ids))) for alias, ids in names.items() if len(set(ids)) > 1)
    for alias, ids in shared:
        warnings.append(f"alias {alias!r} names several templates: {', '.join(ids)}")

    by_family: dict[str, dict[str, Any]] = {}
    for template_id, record in sorted(records.items()):
        row = by_family.setdefault(record["family"], {"family": record["family"], "contracted": 0, "catalog_total": 0, "missing": []})
        row["catalog_total"] += 1
        if template_id in seen:
            row["contracted"] += 1
        else:
            row["missing"].append(template_id)
            errors.append(f"{template_id}: no input contract in template_contracts.json")
    for row in by_family.values():
        row["coverage_ratio"] = round(row["contracted"] / row["catalog_total"], 4) if row["catalog_total"] else 0.0
    contracted = sum(row["contracted"] for row in by_family.values())
    return {
        "ok": not errors,
        "contracts": str(contracts_path),
        "template_count": len(templates),
        "catalog_template_count": len(records),
        "coverage": {
            "contracted": contracted,
            "catalog_total": len(records),
            "coverage_ratio": round(contracted / len(records), 4) if records else 0.0,
            "by_family": [by_family[family] for family in sorted(by_family)],
        },
        "errors": errors,
        "warnings": warnings,
    }


def render_text(payload: dict[str, Any], max_missing_per_family: int) -> str:
    coverage = payload["coverage"]
    lines = [
        f"{'PASS' if payload['ok'] else 'FAIL'}: {payload['contracts']}",
        f"contracts: {coverage['contracted']} / {coverage['catalog_total']} ({coverage['coverage_ratio']:.1%})",
        "coverage by family:",
    ]
    for row in coverage["by_family"]:
        line = f"- {row['family']}: {row['contracted']} / {row['catalog_total']} ({row['coverage_ratio']:.1%})"
        missing = row["missing"][:max_missing_per_family]
        if missing:
            remaining = len(row["missing"]) - len(missing)
            line += f"; missing: {', '.join(missing)}{f', ... +{remaining}' if remaining > 0 else ''}"
        lines.append(line)
    lines.extend(f"warning: {warning}" for warning in payload["warnings"])
    lines.extend(f"error: {error}" for error in payload["errors"])
    return "\n".join(lines)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate template input contracts against the omics visualization catalog.")
    parser.add_argument("--contracts", default=CONTRACTS, help="Path to template_contracts.json.")
    parser.add_argument(
        "--max-missing-per-family",
        type=int,
        default=8,
        help="Maximum missing contract ids to show per family in text output.",
    )
    parser.add_argument(
        "--fail-under",
        type=float,
        default=0.0,
        help="Optional minimum total coverage ratio, between 0 and 1.",
    )
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    payload = validate_contracts(Path(args.contracts).expanduser().resolve())
    if payload["coverage"]["coverage_ratio"] < args.fail_under:
        payload["ok"] = False
        payload["errors"].append(
            f"coverage {payload['coverage']['coverage_ratio']:.1%} is below required {args.fail_under:.1%}"
        )
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print(render_text(payload, args.max_missing_per_family))
    return 0 if payload["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
