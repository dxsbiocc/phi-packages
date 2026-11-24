#!/usr/bin/env python3
"""Compile Foldseek easy-cluster *_clu.tsv outputs into standardized per-design
lineage tables (TSV + JSON) keyed by design_id.

Parses representative -> member pairs and assigns deterministic cluster_id
labels (ordered largest-first, then lexicographically by representative) along
with cluster_representative, is_representative, and cluster_size columns
compatible with DESIGN_FILTER (--group-by cluster_id).
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

STRUCTURE_SUFFIXES = (".cif.gz", ".pdb.gz", ".cif", ".mmcif", ".pdb")


def _clean_design_id(raw: str) -> str:
    name = Path(raw.strip()).name
    lower = name.lower()
    for suf in STRUCTURE_SUFFIXES:
        if lower.endswith(suf):
            return name[: -len(suf)]
    return name


def compile_cluster(clu_path: Path, step_name: str) -> list[dict[str, Any]]:
    if not clu_path.is_file():
        raise SystemExit(f"Foldseek cluster TSV not found: {clu_path}")

    rep_to_members: dict[str, list[str]] = {}
    for line in clu_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        rep = _clean_design_id(parts[0])
        mem = _clean_design_id(parts[1])
        rep_to_members.setdefault(rep, []).append(mem)

    ordered_reps = sorted(rep_to_members.keys(), key=lambda r: (-len(rep_to_members[r]), r))
    rows: list[dict[str, Any]] = []
    for clu_idx, rep in enumerate(ordered_reps, start=1):
        members = rep_to_members[rep]
        clu_id = f"clu_{clu_idx}"
        size = len(members)
        for mem in members:
            rows.append(
                {
                    "design_id": mem,
                    "step": step_name,
                    "cluster_id": clu_id,
                    "cluster_representative": rep,
                    "is_representative": "true" if mem == rep else "false",
                    "cluster_size": size,
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile Foldseek easy-cluster outputs.")
    parser.add_argument("--input-tsv", required=True, help="Raw Foldseek easy-cluster *_clu.tsv.")
    parser.add_argument("--step-name", default="foldseek_easycluster", help="Pipeline step name.")
    parser.add_argument("--out-tsv", required=True, help="Output standardized TSV path.")
    parser.add_argument("--out-json", required=True, help="Output standardized JSON path.")
    args = parser.parse_args()

    rows = compile_cluster(Path(args.input_tsv), args.step_name)
    if not rows:
        raise SystemExit("Foldseek easy-cluster produced zero rows to compile.")

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
