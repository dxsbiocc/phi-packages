#!/usr/bin/env python3
"""Compile BindCraft final_design_stats.csv and Accepted PDBs into standardized
lineage tables (TSV + JSON) keyed strictly by column name.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any


def _maybe_float(val: Any) -> float | None:
    if val is None:
        return None
    s = str(val).strip()
    if not s or s.lower() in ("na", "nan", "none", "null"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile BindCraft outputs into standardized tables.")
    parser.add_argument("--design-path", required=True, help="BindCraft design_path output directory.")
    parser.add_argument("--step-name", default="bindcraft", help="Pipeline step name.")
    parser.add_argument("--out-structures-dir", required=True, help="Flat output directory for accepted PDBs.")
    parser.add_argument("--out-tsv", required=True, help="Output standardized TSV table path.")
    parser.add_argument("--out-json", required=True, help="Output standardized JSON list path.")
    args = parser.parse_args()

    design_dir = Path(args.design_path)
    stats_csv = design_dir / "final_design_stats.csv"
    if not stats_csv.is_file():
        candidates = sorted(design_dir.rglob("final_design_stats.csv"))
        stats_csv = candidates[0] if candidates else stats_csv
    if not stats_csv.is_file():
        raise SystemExit(f"BindCraft final_design_stats.csv not found under {design_dir}")

    accepted_pdbs: dict[str, Path] = {}
    for search_root in (design_dir / "Accepted", design_dir):
        if search_root.is_dir():
            for pdb in sorted(search_root.rglob("*.pdb")):
                accepted_pdbs.setdefault(pdb.stem, pdb)

    flat_dir = Path(args.out_structures_dir)
    flat_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    with stats_csv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for idx, rec in enumerate(reader, start=1):
            design_id = (rec.get("Design") or rec.get("design_id") or f"bindcraft_{idx}").strip()
            parent_id = design_id.rsplit("_mpnn", 1)[0] if "_mpnn" in design_id else design_id
            seq = (rec.get("Sequence") or rec.get("sequence") or "").strip()

            copied_struct = ""
            if design_id in accepted_pdbs:
                dst = flat_dir / f"{design_id}.pdb"
                shutil.copy2(accepted_pdbs[design_id], dst)
                copied_struct = str(dst)

            plddt = _maybe_float(rec.get("Average_pLDDT") or rec.get("pLDDT"))
            complex_plddt = round(plddt / 100.0, 4) if (plddt is not None and plddt > 1.0) else plddt
            iptm = _maybe_float(rec.get("Average_i_pTM") or rec.get("i_pTM"))
            ptm = _maybe_float(rec.get("Average_pTM") or rec.get("pTM"))
            i_pae = _maybe_float(rec.get("Average_i_pAE") or rec.get("i_pAE"))
            dg = _maybe_float(rec.get("Average_dG") or rec.get("dG"))
            sc = _maybe_float(rec.get("Average_ShapeComplementarity") or rec.get("ShapeComplementarity"))
            iface_res = _maybe_float(rec.get("Average_n_InterfaceResidues") or rec.get("n_InterfaceResidues"))
            iface_hb = _maybe_float(rec.get("Average_n_InterfaceHbonds") or rec.get("n_InterfaceHbonds"))
            clashes = _maybe_float(
                rec.get("Average_Relaxed_Clashes")
                or rec.get("Average_Unrelaxed_Clashes")
                or rec.get("Relaxed_Clashes")
            )

            rows.append(
                {
                    "design_id": design_id,
                    "parent_design_id": parent_id,
                    "step": args.step_name,
                    "sequence": seq,
                    "sequence_length": len(seq),
                    "structure_path": copied_struct,
                    "confidence_score": iptm if iptm is not None else complex_plddt,
                    "complex_plddt": complex_plddt,
                    "iptm": iptm,
                    "ptm": ptm,
                    "i_pae_normalized": i_pae,
                    "dG": dg,
                    "shape_complementarity": sc,
                    "interface_residues": iface_res,
                    "interface_hbonds": iface_hb,
                    "clash_count": clashes,
                }
            )

    if not rows:
        raise SystemExit(f"No design rows parsed from {stats_csv}")

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
