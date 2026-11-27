#!/usr/bin/env python3
"""Parse fpocket *_info.txt and pocket*_atm.pdb files into a ranked pockets.tsv
table plus a top_pocket_box.json for AutoDock Vina / GNINA docking grid sizing.
"""
import argparse
import csv
import json
import re
from pathlib import Path


POCKET_HEADER_RE = re.compile(r"^Pocket\s+(\d+)\s*:")
KV_RE = re.compile(r"^\s*([^:]+?)\s*:\s*([+-]?\d+(?:\.\d+)?)\s*$")


def parse_pocket_pdb_box(pdb_path: Path, padding: float = 6.0, min_box: float = 18.0):
    if not pdb_path.is_file():
        return None
    xs, ys, zs = [], [], []
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
            try:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
            except ValueError:
                continue
    if not xs:
        return None
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)
    sx = max(min_box, (max(xs) - min(xs)) + 2.0 * padding)
    sy = max(min_box, (max(ys) - min(ys)) + 2.0 * padding)
    sz = max(min_box, (max(zs) - min(zs)) + 2.0 * padding)
    return {
        "center_x": round(cx, 3),
        "center_y": round(cy, 3),
        "center_z": round(cz, 3),
        "size_x": round(sx, 3),
        "size_y": round(sy, 3),
        "size_z": round(sz, 3),
    }


def parse_info_txt(info_path: Path):
    pockets = []
    cur = None
    for line in info_path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = POCKET_HEADER_RE.match(line.strip())
        if m:
            if cur is not None:
                pockets.append(cur)
            cur = {"pocket_id": int(m.group(1))}
            continue
        if cur is not None:
            kv = KV_RE.match(line)
            if kv:
                key = kv.group(1).strip().lower()
                val = float(kv.group(2))
                cur[key] = val
    if cur is not None:
        pockets.append(cur)
    return pockets


def parse_fpocket_dir(fpocket_out: Path, padding: float = 6.0):
    info_files = sorted(fpocket_out.glob("*_info.txt"))
    if not info_files:
        return []
    raw_pockets = parse_info_txt(info_files[0])
    pockets_dir = fpocket_out / "pockets"

    rows = []
    for p in raw_pockets:
        pid = p["pocket_id"]
        box = parse_pocket_pdb_box(pockets_dir / f"pocket{pid}_atm.pdb", padding=padding) or {
            "center_x": "",
            "center_y": "",
            "center_z": "",
            "size_x": "",
            "size_y": "",
            "size_z": "",
        }
        rows.append(
            {
                "pocket_id": pid,
                "score": round(p.get("score", 0.0), 4),
                "druggability_score": round(p.get("druggability score", 0.0), 4),
                "n_alpha_spheres": int(round(p.get("number of alpha spheres", 0))),
                "total_sasa": round(p.get("total sasa", 0.0), 3),
                "volume": round(p.get("volume", 0.0), 3),
                "hydrophobicity_score": round(p.get("hydrophobicity score", 0.0), 4),
                "polarity_score": round(p.get("polarity score", 0.0), 4),
                **box,
            }
        )

    rows.sort(key=lambda r: (-float(r["druggability_score"]), -float(r["score"]), int(r["pocket_id"])))
    return rows


def main():
    parser = argparse.ArgumentParser(description="Parse fpocket *_info.txt and pocket PDBs.")
    parser.add_argument("--fpocket_out", required=True, help="fpocket *_out directory.")
    parser.add_argument("--padding", type=float, default=6.0, help="Bounding-box padding per side (Angstrom).")
    parser.add_argument("--output_tsv", required=True, help="Output ranked pockets TSV.")
    parser.add_argument("--output_json", default=None, help="Optional top pocket docking box JSON.")
    args = parser.parse_args()

    rows = parse_fpocket_dir(Path(args.fpocket_out), padding=args.padding)
    cols = [
        "pocket_id",
        "score",
        "druggability_score",
        "n_alpha_spheres",
        "total_sasa",
        "volume",
        "hydrophobicity_score",
        "polarity_score",
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
    ]
    out_tsv = Path(args.output_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    if args.output_json:
        top = rows[0] if rows else None
        payload = {
            "tool": "fpocket",
            "n_pockets": len(rows),
            "top_pocket": top,
            "all_pockets": rows[:10],
        }
        Path(args.output_json).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
