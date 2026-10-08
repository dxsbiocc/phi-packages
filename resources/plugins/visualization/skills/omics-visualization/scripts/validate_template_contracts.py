#!/usr/bin/env python3
"""Validate omics visualization routing contracts against the template catalog."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

from route_template import DEFAULT_CONTRACTS, ROLE_HINTS

SKILL_ROOT = Path(__file__).resolve().parents[1]
CATALOG_DIR = SKILL_ROOT / "references" / "catalog"
REQUIRED_TEMPLATE_KEYS = {
    "id",
    "family",
    "title",
    "source",
    "preview",
    "preferred_shapes",
    "roles",
    "intent_keywords",
    "use_when",
    "avoid_when",
}
KNOWN_SHAPES = {
    "association_grid",
    "annotation_sidecars",
    "category_grid",
    "category_distribution",
    "category_magnitude",
    "category_two_series",
    "category_uncertainty",
    "classification_enrichment_tree",
    "colored_hierarchy_area",
    "cytoband_table",
    "dendrogram_matrix",
    "embedding_with_tracks",
    "enrichment_sidecar",
    "enrichment_term_grid",
    "enrichment_grouped_terms",
    "enrichment_ranked_terms",
    "enrichment_zoom_sidecars",
    "enrichment_terms",
    "feature_level_testing",
    "feature_distribution",
    "flow_table",
    "genomic_density_windows",
    "genomic_heatmap_table",
    "genomic_interval_table",
    "genomic_locus_table",
    "genomic_annotation_sidecars",
    "genomic_track_table",
    "genomic_coverage_tracks",
    "grouped_category_value",
    "grouped_distribution",
    "group_split_matrix",
    "hierarchy_edges",
    "hierarchy_area",
    "interval_timeline",
    "long_association",
    "long_table",
    "grouped_correlation",
    "grouped_time_series",
    "mean_difference_table",
    "multi_set_overlap",
    "multi_root_hierarchy",
    "mutation_energy_matrix",
    "network_edges",
    "node_link_sidecars",
    "signed_network_edges",
    "nested_composition",
    "nested_genomic_windows",
    "numeric_matrix",
    "oncoprint_events",
    "one_to_many_association",
    "paired_distribution",
    "parallel_sets_table",
    "part_to_whole",
    "protein_lollipop_table",
    "protein_domain_sidecar",
    "radar_wide_table",
    "rank_time_series",
    "ranked_value",
    "set_membership_summary",
    "set_membership_table",
    "single_axis_series",
    "site_score_table",
    "stacked_composition",
    "survival_table",
    "synteny_blocks",
    "synteny_sidecars",
    "svg_marked_observation",
    "time_series",
    "ternary_components",
    "transcript_feature_table",
    "triangular_correlation",
    "two_numeric_observation",
    "two_set_sample_matrix",
    "wide_many_variables",
    "wide_two_series",
    "wide_numeric_table",
}
KNOWN_ROLES = set(ROLE_HINTS) | {
    "cytoband_table",
    "flow_value",
    "genomic_density_windows",
    "genomic_heatmap_table",
    "genomic_interval_table",
    "genomic_locus_table",
    "genomic_track_table",
    "group_split_matrix",
    "mutation_energy_matrix",
    "numeric_matrix",
    "oncoprint_events",
    "protein_lollipop_table",
    "site_score_table",
    "synteny_blocks",
    "ternary_components",
    "transcript_feature_table",
    "two_set_sample_matrix",
}


def unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] in {'"', "'"} and value[-1] == value[0]:
        return value[1:-1]
    return value


def parse_catalog_records(catalog_dir: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    for catalog_path in sorted(catalog_dir.glob("*.yaml")):
        family = catalog_path.stem
        current: dict[str, str] | None = None
        for line in catalog_path.read_text(encoding="utf-8").splitlines():
            id_match = re.match(r"\s*-\s+id:\s*(.+?)\s*$", line)
            if id_match:
                template_id = unquote(id_match.group(1))
                current = {"id": template_id, "family": family, "catalog": str(catalog_path.relative_to(SKILL_ROOT))}
                records[template_id] = current
                continue
            if current is None:
                continue
            attr_match = re.match(r"\s*(source|preview):\s*(.+?)\s*$", line)
            if attr_match:
                current[attr_match.group(1)] = unquote(attr_match.group(2))
    return records


def validate_contracts(contracts_path: Path, skill_root: Path = SKILL_ROOT) -> dict[str, Any]:
    contracts = json.loads(contracts_path.read_text(encoding="utf-8"))
    catalog_records = parse_catalog_records(skill_root / "references" / "catalog")
    errors: list[str] = []
    warnings: list[str] = []
    seen_ids: set[str] = set()

    templates = contracts.get("templates", [])
    if not isinstance(templates, list) or not templates:
        errors.append("contracts must contain a non-empty templates list")
        templates = []

    for index, template in enumerate(templates):
        if not isinstance(template, dict):
            errors.append(f"template at index {index} must be an object")
            continue
        template_id = str(template.get("id", f"<missing-{index}>"))
        missing = sorted(REQUIRED_TEMPLATE_KEYS - set(template))
        if missing:
            errors.append(f"{template_id}: missing keys {', '.join(missing)}")
        if template_id in seen_ids:
            errors.append(f"{template_id}: duplicate template id")
        seen_ids.add(template_id)

        catalog = catalog_records.get(template_id)
        if not catalog:
            errors.append(f"{template_id}: not found in references/catalog")
        else:
            for key in ("family", "source", "preview"):
                if template.get(key) != catalog.get(key):
                    errors.append(f"{template_id}: {key} differs from catalog ({template.get(key)!r} != {catalog.get(key)!r})")

        for key in ("source", "preview"):
            path_value = template.get(key)
            if isinstance(path_value, str) and not (skill_root / path_value).exists():
                errors.append(f"{template_id}: {key} path does not exist: {path_value}")

        shapes = template.get("preferred_shapes", [])
        if not isinstance(shapes, list) or not shapes:
            errors.append(f"{template_id}: preferred_shapes must be a non-empty list")
        else:
            unknown_shapes = sorted(set(shapes) - KNOWN_SHAPES)
            if unknown_shapes:
                errors.append(f"{template_id}: unknown shapes {', '.join(unknown_shapes)}")

        roles = template.get("roles", {})
        if not isinstance(roles, dict):
            errors.append(f"{template_id}: roles must be an object")
            continue
        required = roles.get("required", [])
        optional = roles.get("optional", [])
        if not isinstance(required, list) or not required:
            errors.append(f"{template_id}: roles.required must be a non-empty list")
        if not isinstance(optional, list):
            errors.append(f"{template_id}: roles.optional must be a list")
        unknown_roles = sorted((set(required) | set(optional)) - KNOWN_ROLES)
        if unknown_roles:
            errors.append(f"{template_id}: unknown roles {', '.join(unknown_roles)}")

        keywords = template.get("intent_keywords", [])
        if not isinstance(keywords, list) or not keywords:
            errors.append(f"{template_id}: intent_keywords must be a non-empty list")
        elif len(keywords) < 3:
            warnings.append(f"{template_id}: fewer than three intent keywords")

    missing_paths = [
        template_id
        for template_id, record in catalog_records.items()
        if not (skill_root / record.get("source", "")).exists() or not (skill_root / record.get("preview", "")).exists()
    ]
    if missing_paths:
        warnings.append(f"{len(missing_paths)} catalog templates have missing source or preview paths")

    contract_ids = sorted(seen_ids)
    missing_contracts_by_family: dict[str, list[str]] = {}
    catalog_counts_by_family: dict[str, int] = {}
    contract_counts_by_family: dict[str, int] = {}
    for template_id, record in sorted(catalog_records.items()):
        family = record["family"]
        catalog_counts_by_family[family] = catalog_counts_by_family.get(family, 0) + 1
        if template_id not in seen_ids:
            missing_contracts_by_family.setdefault(family, []).append(template_id)
    for template in templates:
        if isinstance(template, dict) and isinstance(template.get("family"), str):
            family = template["family"]
            contract_counts_by_family[family] = contract_counts_by_family.get(family, 0) + 1

    family_rows = []
    for family in sorted(catalog_counts_by_family):
        contracted = contract_counts_by_family.get(family, 0)
        catalog_total = catalog_counts_by_family[family]
        family_rows.append(
            {
                "family": family,
                "contracted": contracted,
                "catalog_total": catalog_total,
                "coverage_ratio": round(contracted / catalog_total, 4) if catalog_total else 0.0,
                "missing": missing_contracts_by_family.get(family, []),
            }
        )

    coverage = {
        "contracted": len(contract_ids),
        "catalog_total": len(catalog_records),
        "coverage_ratio": round(len(contract_ids) / len(catalog_records), 4) if catalog_records else 0.0,
        "contract_ids": contract_ids,
        "by_family": family_rows,
    }

    return {
        "ok": not errors,
        "contracts": str(contracts_path),
        "template_count": len(templates),
        "catalog_template_count": len(catalog_records),
        "coverage": coverage,
        "errors": errors,
        "warnings": warnings,
    }


def render_text(payload: dict[str, Any], max_missing_per_family: int) -> str:
    status = "PASS" if payload["ok"] else "FAIL"
    coverage = payload["coverage"]
    lines = [
        f"{status}: {payload['contracts']}",
        f"contracts: {coverage['contracted']} / {coverage['catalog_total']} ({coverage['coverage_ratio']:.1%})",
        "coverage by family:",
    ]
    for row in coverage["by_family"]:
        line = f"- {row['family']}: {row['contracted']} / {row['catalog_total']} ({row['coverage_ratio']:.1%})"
        missing = row["missing"][:max_missing_per_family]
        if missing:
            suffix = ""
            remaining = len(row["missing"]) - len(missing)
            if remaining > 0:
                suffix = f", ... +{remaining}"
            line += f"; missing: {', '.join(missing)}{suffix}"
        lines.append(line)
    for warning in payload["warnings"]:
        lines.append(f"warning: {warning}")
    for error in payload["errors"]:
        lines.append(f"error: {error}")
    return "\n".join(lines)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate template routing contracts against omics visualization catalogs.")
    parser.add_argument("--contracts", default=str(DEFAULT_CONTRACTS), help="Path to template_contracts.json.")
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
