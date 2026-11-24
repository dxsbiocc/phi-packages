#!/usr/bin/env python3
"""Compile Foldseek easy-search or easy-cluster outputs into standardized
per-design lineage tables (TSV + JSON).

Modes:
  - search:  aggregates best hit per design_id (best_target, max_alntmscore,
             max_qtmscore, max_lddt, min_evalue, num_hits, is_novel).
  - cluster: parses *_clu.tsv (representative -> member) into deterministic
             cluster_id, cluster_representative, is_representative, cluster_size.
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


def _maybe_float(val: str) -> float | None:
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _discover_query_ids(query_dir: Path | None) -> list[str]:
    if query_dir is None or not query_dir.exists():
        return []
    ids: list[str] = []
    if query_dir.is_dir():
        for p in sorted(query_dir.rglob("*")):
            if p.is_file() and any(p.name.lower().endswith(s) for s in STRUCTURE_SUFFIXES):
                ids.append(_clean_design_id(p.name))
    elif query_dir.is_file():
        ids.append(_clean_design_id(query_dir.name))
    return sorted(set(ids))


def compile_search(
    hits_path: Path,
    query_ids: list[str],
    novelty_tm_cutoff: float,
    step_name: str,
) -> list[dict[str, Any]]:
    by_design: dict[str, list[dict[str, Any]]] = {qid: [] for qid in query_ids}
    if hits_path.is_file():
        for line in hits_path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue
            qid = _clean_design_id(parts[0])
            by_design.setdefault(qid, []).append(
                {
                    "target": parts[1].strip(),
                    "fident": _maybe_float(parts[2]),
                    "alnlen": _maybe_float(parts[3]),
                    "evalue": _maybe_float(parts[4]),
                    "bits": _maybe_float(parts[5]),
                    "alntmscore": _maybe_float(parts[6]),
                    "qtmscore": _maybe_float(parts[7]),
                    "ttmscore": _maybe_float(parts[8]),
                    "lddt": _maybe_float(parts[9]),
                }
            )

    rows: list[dict[str, Any]] = []
    for design_id in sorted(by_design.keys()):
        hits = by_design[design_id]
        if not hits:
            rows.append(
                {
                    "design_id": design_id,
                    "step": step_name,
                    "mode": "search",
                    "num_hits": 0,
                    "best_target": "",
                    "max_alntmscore": 0.0,
                    "max_qtmscore": 0.0,
                    "max_lddt": 0.0,
                    "min_evalue": None,
                    "is_novel": "true",
                }
            )
            continue

        best = max(hits, key=lambda h: (h["alntmscore"] or 0.0, h["bits"] or 0.0))
        max_tm = max((h["alntmscore"] for h in hits if h["alntmscore"] is not None), default=0.0)
        max_qtm = max((h["qtmscore"] for h in hits if h["qtmscore"] is not None), default=0.0)
        max_lddt = max((h["lddt"] for h in hits if h["lddt"] is not None), default=0.0)
        evalues = [h["evalue"] for h in hits if h["evalue"] is not None]
        min_eval = min(evalues) if evalues else None

        rows.append(
            {
                "design_id": design_id,
                "step": step_name,
                "mode": "search",
                "num_hits": len(hits),
                "best_target": best["target"],
                "max_alntmscore": round(max_tm, 4),
                "max_qtmscore": round(max_qtm, 4),
                "max_lddt": round(max_lddt, 4),
                "min_evalue": min_eval,
                "is_novel": "true" if max_tm < novelty_tm_cutoff else "false",
            }
        )
    return rows


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

    # Sort clusters largest-first for deterministic cluster_id assignment
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
                    "mode": "cluster",
                    "cluster_id": clu_id,
                    "cluster_representative": rep,
                    "is_representative": "true" if mem == rep else "false",
                    "cluster_size": size,
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile Foldseek search or cluster outputs.")
    parser.add_argument("--mode", choices=("search", "cluster"), required=True)
    parser.add_argument("--input-tsv", required=True, help="Raw Foldseek search hits TSV or *_clu.tsv.")
    parser.add_argument("--query-dir", default="", help="Optional staged query directory to include zero-hit designs.")
    parser.add_argument("--novelty-tm-cutoff", type=float, default=0.5, help="TM-score cutoff below which a design is marked is_novel=true.")
    parser.add_argument("--step-name", default="foldseek", help="Pipeline step name.")
    parser.add_argument("--out-tsv", required=True, help="Output standardized TSV path.")
    parser.add_argument("--out-json", required=True, help="Output standardized JSON path.")
    args = parser.parse_args()

    if args.mode == "search":
        qids = _discover_query_ids(Path(args.query_dir)) if args.query_dir else []
        rows = compile_search(Path(args.input_tsv), qids, args.novelty_tm_cutoff, args.step_name)
    else:
        rows = compile_cluster(Path(args.input_tsv), args.step_name)

    if not rows:
        raise SystemExit("Foldseek produced zero rows to compile.")

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
