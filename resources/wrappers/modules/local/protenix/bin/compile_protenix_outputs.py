#!/usr/bin/env python3
"""Compile Protenix prediction outputs (*_summary_confidence_sample_*.json and
*_sample_*.cif) into flat best-sample structures and standardized lineage
tables (TSV + JSON).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from pathlib import Path
from typing import Any

SAMPLE_RE = re.compile(r"^(?P<name>.+)_summary_confidence_sample_(?P<idx>\d+)\.json$")


def _maybe_float(val: Any) -> float | None:
    if val is None or isinstance(val, bool):
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile Protenix outputs into standardized tables.")
    parser.add_argument("--out-dir", required=True, help="Protenix --out_dir directory.")
    parser.add_argument("--manifest", default="", help="Optional lineage manifest JSON.")
    parser.add_argument("--step-name", default="protenix", help="Pipeline step name.")
    parser.add_argument("--out-structures-dir", required=True, help="Flat output directory for top predicted CIFs.")
    parser.add_argument("--out-tsv", required=True, help="Output standardized TSV table path.")
    parser.add_argument("--out-json", required=True, help="Output standardized JSON list path.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    manifest: dict[str, dict[str, Any]] = {}
    if args.manifest and Path(args.manifest).is_file():
        manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))

    # Group summary confidence JSON files by job name
    by_job: dict[str, list[tuple[Path, str]]] = {}
    for cf in sorted(out_dir.rglob("*_summary_confidence_sample_*.json")):
        m = SAMPLE_RE.match(cf.name)
        if not m:
            continue
        job_name = m.group("name")
        sample_idx = m.group("idx")
        by_job.setdefault(job_name, []).append((cf, sample_idx))

    if not by_job:
        raise SystemExit(f"No Protenix *_summary_confidence_sample_*.json files found under {out_dir}")

    flat_dir = Path(args.out_structures_dir)
    flat_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for job_name in sorted(by_job.keys()):
        best_cf: Path | None = None
        best_idx = "0"
        best_data: dict[str, Any] = {}
        best_score = float("-inf")

        for cf, s_idx in by_job[job_name]:
            try:
                data = json.loads(cf.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                continue
            r_score = _maybe_float(data.get("ranking_score"))
            if r_score is None:
                r_score = _maybe_float(data.get("iptm")) or _maybe_float(data.get("ptm")) or float("-inf")
            if best_cf is None or r_score > best_score:
                best_cf = cf
                best_idx = s_idx
                best_data = data
                best_score = r_score

        if best_cf is None:
            continue

        info = manifest.get(job_name, {})
        design_id = str(info.get("design_id") or job_name)
        parent_id = str(info.get("parent_design_id") or design_id)

        cand_cif = best_cf.parent / f"{job_name}_sample_{best_idx}.cif"
        if not cand_cif.is_file():
            cifs = sorted(best_cf.parent.glob(f"*_sample_{best_idx}.cif"))
            cand_cif = cifs[0] if cifs else cand_cif

        copied_struct = ""
        if cand_cif.is_file():
            dst = flat_dir / f"{design_id}.cif"
            shutil.copy2(cand_cif, dst)
            copied_struct = str(dst)

        ranking_score = _maybe_float(best_data.get("ranking_score"))
        ptm = _maybe_float(best_data.get("ptm"))
        iptm = _maybe_float(best_data.get("iptm"))
        plddt = _maybe_float(best_data.get("plddt"))
        gpde = _maybe_float(best_data.get("gpde"))
        has_clash = best_data.get("has_clash")
        complex_plddt = round(plddt / 100.0, 4) if (plddt is not None and plddt > 1.0) else plddt

        rows.append(
            {
                "design_id": design_id,
                "parent_design_id": parent_id,
                "step": args.step_name,
                "structure_path": copied_struct,
                "parent_structure_path": info.get("parent_structure_path", ""),
                "confidence_score": ranking_score if ranking_score is not None else iptm,
                "ranking_score": ranking_score,
                "ptm": ptm,
                "iptm": iptm,
                "plddt": plddt,
                "complex_plddt": complex_plddt,
                "gpde": gpde,
                "has_clash": has_clash,
            }
        )

    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})

    out_json = Path(args.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
