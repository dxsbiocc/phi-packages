#!/usr/bin/env python3
"""Parse AutoDock Vina docked .pdbqt files into a normalized score_table.tsv
directly compatible with drug-discovery's chem_docking_eval tool.
"""
import argparse
import csv
import json
import re
from pathlib import Path


REMARK_RE = re.compile(
    r"REMARK\s+VINA\s+RESULT:\s*([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)"
)


def parse_pdbqt_poses(pdbqt_path: Path):
    poses = []
    text = pdbqt_path.read_text(encoding="utf-8", errors="replace")
    ligand_id = pdbqt_path.stem
    if ligand_id.endswith("_out"):
        ligand_id = ligand_id[:-4]

    rank = 0
    for line in text.splitlines():
        m = REMARK_RE.search(line)
        if m:
            rank += 1
            affinity = float(m.group(1))
            rmsd_lb = float(m.group(2))
            rmsd_ub = float(m.group(3))
            poses.append(
                {
                    "ligand_id": ligand_id,
                    "pose_rank": rank,
                    "score": round(affinity, 4),
                    "affinity_kcal_mol": round(affinity, 4),
                    "rmsd_lb": round(rmsd_lb, 4),
                    "rmsd_ub": round(rmsd_ub, 4),
                    "pose_file": pdbqt_path.name,
                }
            )
    return poses


def parse_docked_directory(docked_dir: Path):
    rows = []
    for path in sorted(docked_dir.glob("*.pdbqt")):
        rows.extend(parse_pdbqt_poses(path))
    rows.sort(key=lambda r: (r["score"], r["ligand_id"], r["pose_rank"]))
    return rows


def summarize_best_poses(rows):
    best_by_ligand = {}
    for row in rows:
        lid = row["ligand_id"]
        if lid not in best_by_ligand or row["score"] < best_by_ligand[lid]["score"]:
            best_by_ligand[lid] = row
    best_list = sorted(best_by_ligand.values(), key=lambda r: r["score"])
    return {
        "tool": "autodock_vina",
        "n_poses": len(rows),
        "n_ligands": len(best_by_ligand),
        "best_affinity_kcal_mol": best_list[0]["score"] if best_list else None,
        "top_ligands": best_list[:10],
    }


def main():
    parser = argparse.ArgumentParser(description="Parse AutoDock Vina PDBQT outputs.")
    parser.add_argument("--docked_dir", required=True, help="Directory with *_out.pdbqt files.")
    parser.add_argument("--output_tsv", required=True, help="Output normalized score_table.tsv.")
    parser.add_argument("--output_json", default=None, help="Optional summary JSON path.")
    parser.add_argument("--best_only", action="store_true", help="Write only rank-1 best pose per ligand.")
    args = parser.parse_args()

    rows = parse_docked_directory(Path(args.docked_dir))
    if args.best_only:
        seen = {}
        for row in rows:
            lid = row["ligand_id"]
            if lid not in seen or row["score"] < seen[lid]["score"]:
                seen[lid] = row
        rows = sorted(seen.values(), key=lambda r: (r["score"], r["ligand_id"]))

    cols = ["ligand_id", "pose_rank", "score", "affinity_kcal_mol", "rmsd_lb", "rmsd_ub", "pose_file"]
    out_tsv = Path(args.output_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = summarize_best_poses(rows)
    if args.output_json:
        Path(args.output_json).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
