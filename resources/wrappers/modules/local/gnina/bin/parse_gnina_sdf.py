#!/usr/bin/env python3
"""Parse GNINA docked SDF/SDF.GZ outputs (SD tags: minimizedAffinity, CNNscore,
CNNaffinity, CNN_VS) into a normalized score_table.tsv for chem_docking_eval.
"""
import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path


TAG_HEADER_RE = re.compile(r"^>\s*<([^>]+)>\s*$")


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _to_float(val):
    if val is None:
        return None
    try:
        return float(str(val).strip())
    except ValueError:
        return None


def parse_gnina_sdf(sdf_path: Path):
    with _open_text(sdf_path) as fh:
        content = fh.read()

    raw_records = [rec.strip("\r\n") for rec in content.split("$$$$") if rec.strip()]
    parsed = []
    pose_counter = defaultdict(int)

    for idx, rec in enumerate(raw_records, start=1):
        lines = rec.splitlines()
        title = lines[0].strip() if lines else ""
        tags = {}
        cur_tag = None
        cur_vals = []
        for line in lines[1:]:
            m = TAG_HEADER_RE.match(line)
            if m:
                if cur_tag is not None:
                    tags[cur_tag] = "\n".join(cur_vals).strip()
                cur_tag = m.group(1).strip()
                cur_vals = []
            elif cur_tag is not None:
                if line.strip() == "":
                    tags[cur_tag] = "\n".join(cur_vals).strip()
                    cur_tag = None
                    cur_vals = []
                else:
                    cur_vals.append(line.strip())
        if cur_tag is not None:
            tags[cur_tag] = "\n".join(cur_vals).strip()

        ligand_id = (
            tags.get("ligand_id")
            or tags.get("_Name")
            or tags.get("Name")
            or title
            or f"ligand_{idx}"
        ).strip()
        pose_counter[ligand_id] += 1

        min_aff = _to_float(tags.get("minimizedAffinity"))
        cnn_score = _to_float(tags.get("CNNscore"))
        cnn_aff = _to_float(tags.get("CNNaffinity"))
        cnn_vs = _to_float(tags.get("CNN_VS"))
        if cnn_vs is None and cnn_score is not None and cnn_aff is not None:
            cnn_vs = cnn_score * cnn_aff

        if cnn_aff is not None:
            unified_score = -cnn_aff
        elif min_aff is not None:
            unified_score = min_aff
        elif cnn_score is not None:
            unified_score = -cnn_score
        else:
            unified_score = 0.0

        parsed.append(
            {
                "ligand_id": ligand_id,
                "pose_rank": pose_counter[ligand_id],
                "score": round(unified_score, 4),
                "minimized_affinity": round(min_aff, 4) if min_aff is not None else "",
                "cnn_score": round(cnn_score, 4) if cnn_score is not None else "",
                "cnn_affinity": round(cnn_aff, 4) if cnn_aff is not None else "",
                "cnn_vs": round(cnn_vs, 4) if cnn_vs is not None else "",
                "source_file": sdf_path.name,
            }
        )

    parsed.sort(key=lambda r: (r["score"], r["ligand_id"], r["pose_rank"]))
    return parsed


def summarize_gnina(rows):
    best_by_ligand = {}
    for row in rows:
        lid = row["ligand_id"]
        if lid not in best_by_ligand or row["score"] < best_by_ligand[lid]["score"]:
            best_by_ligand[lid] = row
    best_list = sorted(best_by_ligand.values(), key=lambda r: r["score"])
    return {
        "tool": "gnina",
        "n_poses": len(rows),
        "n_ligands": len(best_by_ligand),
        "best_score": best_list[0]["score"] if best_list else None,
        "top_ligands": best_list[:10],
    }


def main():
    parser = argparse.ArgumentParser(description="Parse GNINA SDF/SDF.GZ results.")
    parser.add_argument("--sdf", required=True, help="Input .sdf or .sdf.gz file.")
    parser.add_argument("--output_tsv", required=True, help="Output normalized score_table.tsv.")
    parser.add_argument("--output_json", default=None, help="Optional summary JSON path.")
    parser.add_argument("--best_only", action="store_true", help="Keep only top pose per ligand.")
    args = parser.parse_args()

    rows = parse_gnina_sdf(Path(args.sdf))
    if args.best_only:
        seen = {}
        for row in rows:
            lid = row["ligand_id"]
            if lid not in seen or row["score"] < seen[lid]["score"]:
                seen[lid] = row
        rows = sorted(seen.values(), key=lambda r: (r["score"], r["ligand_id"]))

    cols = [
        "ligand_id",
        "pose_rank",
        "score",
        "minimized_affinity",
        "cnn_score",
        "cnn_affinity",
        "cnn_vs",
        "source_file",
    ]
    out_tsv = Path(args.output_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = summarize_gnina(rows)
    if args.output_json:
        Path(args.output_json).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
