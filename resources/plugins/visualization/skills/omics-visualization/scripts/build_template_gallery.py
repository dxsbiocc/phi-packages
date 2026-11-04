#!/usr/bin/env python3
"""Build the plugin-page template gallery from the omics-visualization catalogs.

The gallery lists every catalog template and its preview.png. Layout audit
images under scripts/layouts/ are not plot templates and stay out of the file.
Phi's plugin page reads references/template_gallery.json when the plugin is open.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

SKILL_ROOT = Path(__file__).resolve().parents[1]
GALLERY_PATH = SKILL_ROOT / "references" / "template_gallery.json"
PREVIEW_PATH = re.compile(r"^scripts/[a-z0-9_]+/[a-z0-9_]+/preview\.png$")
TEMPLATE_ID = re.compile(r"^[a-z][a-z0-9-]{0,80}$")


def unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] in {'"', "'"} and value[-1] == value[0]:
        return value[1:-1]
    return value


def parse_families(plots_path: Path) -> list[dict[str, str]]:
    families: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    in_groups = False
    for line in plots_path.read_text(encoding="utf-8").splitlines():
        if line == "groups:":
            in_groups = True
            continue
        if not in_groups:
            continue
        if line and not line.startswith(" "):
            break
        group_match = re.match(r"^  ([a-z0-9_]+):\s*$", line)
        if group_match:
            if current is not None:
                families.append(current)
            current = {"id": group_match.group(1)}
            continue
        if current is None:
            continue
        title_match = re.match(r'^    title:\s*(.+?)\s*$', line)
        if title_match and "title" not in current:
            current["title"] = unquote(title_match.group(1))
        catalog_match = re.match(r'^    catalog:\s*(.+?)\s*$', line)
        if catalog_match and "catalog" not in current:
            current["catalog"] = unquote(catalog_match.group(1))
    if current is not None:
        families.append(current)
    return families


def parse_templates(catalog_path: Path) -> list[dict[str, str]]:
    templates: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    for line in catalog_path.read_text(encoding="utf-8").splitlines():
        id_match = re.match(r"\s*-\s+id:\s*(.+?)\s*$", line)
        if id_match:
            current = {"id": unquote(id_match.group(1))}
            templates.append(current)
            continue
        if current is None:
            continue
        for key in ("title", "preview"):
            attr_match = re.match(rf"\s*{key}:\s*(.+?)\s*$", line)
            if attr_match and key not in current:
                current[key] = unquote(attr_match.group(1))
    return templates


def gallery_document(skill_root: Path = SKILL_ROOT) -> tuple[dict[str, Any], list[str]]:
    errors: list[str] = []
    plots_path = skill_root / "references" / "plots.yaml"
    families = parse_families(plots_path)
    document_families: list[dict[str, Any]] = []
    seen_ids: set[str] = set()

    for family in families:
        family_id = family.get("id", "")
        title = family.get("title", "")
        catalog = family.get("catalog", "")
        if not family_id or not title or not catalog:
            errors.append(f"plots.yaml group {family_id or '<missing>'} needs an id, title, and catalog")
            continue
        expected_catalog = f"references/catalog/{family_id}.yaml"
        if catalog != expected_catalog:
            errors.append(f"{family_id}: catalog path {catalog} != {expected_catalog}")
        catalog_path = skill_root / catalog
        if not catalog_path.is_file():
            errors.append(f"{family_id}: missing catalog {catalog}")
            continue
        templates: list[dict[str, str]] = []
        for record in parse_templates(catalog_path):
            template_id = record.get("id", "")
            template_title = record.get("title", "")
            preview = record.get("preview", "")
            if not TEMPLATE_ID.fullmatch(template_id):
                errors.append(f"{template_id or '<missing>'}: template id is not a gallery id")
                continue
            if template_id in seen_ids:
                errors.append(f"{template_id}: duplicate template id")
                continue
            seen_ids.add(template_id)
            if not template_title:
                errors.append(f"{template_id}: missing title")
                continue
            if not PREVIEW_PATH.fullmatch(preview):
                errors.append(f"{template_id}: preview is not a template preview.png")
                continue
            if not (skill_root / preview).is_file():
                errors.append(f"{template_id}: preview file does not exist: {preview}")
                continue
            templates.append({"id": template_id, "title": template_title, "preview": preview})
        document_families.append({"id": family_id, "title": title, "templates": templates})

    catalog_dir = skill_root / "references" / "catalog"
    if catalog_dir.is_dir():
        listed = {family["id"] for family in families}
        for catalog_path in sorted(catalog_dir.glob("*.yaml")):
            if catalog_path.stem not in listed:
                errors.append(f"{catalog_path.name} is not listed in references/plots.yaml")

    document = {
        "schemaVersion": 1,
        "skill": "omics-visualization",
        "previewBackground": "transparent",
        "families": document_families,
    }
    return document, errors


def canonical_text(document: dict[str, Any]) -> str:
    return json.dumps(document, ensure_ascii=False, indent=2) + "\n"


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write or check references/template_gallery.json.")
    parser.add_argument("--check", action="store_true", help="Exit 1 when the gallery file is missing or stale.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    document, errors = gallery_document()
    for error in errors:
        print(f"error: {error}", file=sys.stderr)
    if errors:
        return 1
    text = canonical_text(document)
    if args.check:
        if not GALLERY_PATH.is_file():
            print("error: references/template_gallery.json is missing", file=sys.stderr)
            return 1
        if GALLERY_PATH.read_text(encoding="utf-8") != text:
            print("error: references/template_gallery.json is stale", file=sys.stderr)
            return 1
        return 0
    GALLERY_PATH.write_text(text, encoding="utf-8")
    count = sum(len(family["templates"]) for family in document["families"])
    print(f"wrote {GALLERY_PATH.relative_to(SKILL_ROOT)} ({count} templates)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
