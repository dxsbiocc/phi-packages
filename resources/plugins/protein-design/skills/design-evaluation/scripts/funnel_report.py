#!/usr/bin/env python3
"""Summarize protein design campaign funnel counts, metric quantiles, parent
backbone diversity, and top ranked candidates from DESIGN_FILTER outputs.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

KEY_METRICS = (
    "confidence_score",
    "ranking_score",
    "complex_plddt",
    "plddt",
    "ptm",
    "iptm",
    "ligand_iptm",
    "min_pae_interaction",
    "ca_rmsd",
    "backbone_rmsd",
    "interface_contacts",
    "interface_hbonds",
    "ligand_contacts",
    "clash_count",
    "affinity_pred_value",
    "affinity_probability_binary",
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def _to_float(val: Any) -> float | None:
    if val is None:
        return None
    s = str(val).strip()
    if not s or s.lower() in ("na", "nan", "null", "none"):
        return None
    try:
        v = float(s)
        return None if math.isnan(v) or math.isinf(v) else v
    except ValueError:
        return None


def _quantile(sorted_vals: list[float], q: float) -> float:
    n = len(sorted_vals)
    if n == 1:
        return sorted_vals[0]
    pos = q * (n - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_vals[lo]
    return sorted_vals[lo] + (pos - lo) * (sorted_vals[hi] - sorted_vals[lo])


def _summarize_metrics(rows: list[dict[str, str]]) -> dict[str, dict[str, Any]]:
    stats: dict[str, dict[str, Any]] = {}
    if not rows:
        return stats
    cols = list(rows[0].keys())
    candidate_cols = [c for c in KEY_METRICS if c in cols]
    for col in cols:
        if any(col.endswith(f".{k}") for k in KEY_METRICS) and col not in candidate_cols:
            candidate_cols.append(col)

    for col in candidate_cols:
        vals = sorted(v for r in rows if (v := _to_float(r.get(col))) is not None)
        if not vals:
            continue
        stats[col] = {
            "count": len(vals),
            "min": round(vals[0], 4),
            "median": round(_quantile(vals, 0.5), 4),
            "p90": round(_quantile(vals, 0.9), 4),
            "max": round(vals[-1], 4),
        }
    return stats


def _find_parent_col(rows: list[dict[str, str]]) -> str | None:
    if not rows:
        return None
    for col in rows[0].keys():
        if col == "parent_design_id" or col.endswith(".parent_design_id"):
            return col
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate campaign funnel diagnostics JSON.")
    parser.add_argument("--merged-tsv", "--merged_tsv", dest="merged_tsv", required=True, help="Path to merged.tsv.")
    parser.add_argument("--filtered-tsv", "--filtered_tsv", dest="filtered_tsv", default="", help="Optional filtered.tsv.")
    parser.add_argument("--summary-json", "--summary_json", dest="summary_json", default="", help="Optional summary.json.")
    parser.add_argument("--top-k", "--top_k", dest="top_k", type=int, default=10, help="Top candidates preview count.")
    args = parser.parse_args()

    merged_path = Path(args.merged_tsv)
    if not merged_path.is_file():
        raise SystemExit(f"merged_tsv not found: {merged_path}")

    merged_rows = _read_tsv(merged_path)
    passing_rows = [r for r in merged_rows if str(r.get("passed_filter", "")).strip().lower() == "true"]

    if args.filtered_tsv and Path(args.filtered_tsv).is_file():
        selected_rows = _read_tsv(Path(args.filtered_tsv))
    else:
        selected_rows = passing_rows

    summary_meta: dict[str, Any] = {}
    if args.summary_json and Path(args.summary_json).is_file():
        summary_meta = json.loads(Path(args.summary_json).read_text(encoding="utf-8"))

    parent_col = _find_parent_col(merged_rows)
    unique_parents_total = (
        len({r.get(parent_col, "").strip() for r in merged_rows if r.get(parent_col, "").strip()})
        if parent_col
        else None
    )
    unique_parents_selected = (
        len({r.get(parent_col, "").strip() for r in selected_rows if r.get(parent_col, "").strip()})
        if parent_col
        else None
    )

    total = len(merged_rows)
    passed = len(passing_rows)
    selected = len(selected_rows)
    pass_rate = round(passed / total, 4) if total > 0 else 0.0

    preview_cols = ["design_id"]
    if parent_col:
        preview_cols.append(parent_col)
    for col in (selected_rows[0].keys() if selected_rows else []):
        if col in KEY_METRICS or any(col.endswith(f".{k}") for k in KEY_METRICS):
            if col not in preview_cols:
                preview_cols.append(col)

    top_candidates = [
        {k: r.get(k, "") for k in preview_cols}
        for r in selected_rows[: max(0, args.top_k)]
    ]

    report = {
        "total_designs": total,
        "passed_filter": passed,
        "selected_designs": selected,
        "pass_rate": pass_rate,
        "unique_parents_total": unique_parents_total,
        "unique_parents_selected": unique_parents_selected,
        "filter_expr": summary_meta.get("expr", ""),
        "rank_by": summary_meta.get("rank_by", ""),
        "metric_quantiles_all": _summarize_metrics(merged_rows),
        "metric_quantiles_selected": _summarize_metrics(selected_rows),
        "top_candidates": top_candidates,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
