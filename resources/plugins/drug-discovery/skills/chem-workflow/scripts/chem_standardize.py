#!/usr/bin/env python3
"""chem_standardize: Standardize small-molecule structures (TSV/CSV/SMI/SDF),
strip salts/solvents, neutralize ionized groups, canonicalize SMILES, generate
structure keys (InChIKey or graph hash), flag invalid SMILES and stereochemistry,
and deduplicate compounds.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    BACKEND,
    emit,
    fail,
    read_compounds,
    setup_matplotlib,
    standardize_smiles,
    write_artifact,
    write_fallback_png,
    write_tsv,
)


def _plot_standardization(
    out_png: Path,
    n_valid: int,
    n_salts: int,
    n_neutralized: int,
    n_dups: int,
    n_invalid: int,
    n_stereo: int,
    mws: list[float],
) -> None:
    plt = setup_matplotlib()
    if plt is None:
        write_fallback_png(out_png, accent_rgb=(15, 118, 110))
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.3))
    labels = ["Valid", "Salts\nStripped", "Neutralized", "Stereo", "Duplicates", "Invalid"]
    counts = [n_valid, n_salts, n_neutralized, n_stereo, n_dups, n_invalid]
    colors = ["#0f766e", "#0284c7", "#0d9488", "#6366f1", "#f59e0b", "#dc2626"]
    bars = ax1.bar(labels, counts, color=colors, width=0.62, edgecolor="#1e293b", linewidth=0.5)
    for b, val in zip(bars, counts):
        ax1.text(b.get_x() + b.get_width() / 2.0, b.get_height() + max(counts, default=1) * 0.02, str(val),
                 ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    ax1.set_ylabel("Compound Count")
    ax1.set_title(f"Structure Standardization QC ({BACKEND})")
    ax1.set_ylim(0, max(1, max(counts)) * 1.18)

    if mws:
        bins = min(18, max(5, len(set(round(m, 0) for m in mws))))
        ax2.hist(mws, bins=bins, color="#0284c7", alpha=0.82, edgecolor="white")
        ax2.axvline(500.0, color="#dc2626", linestyle="--", linewidth=1.2, label="Ro5 MW = 500")
        ax2.legend(frameon=False)
    ax2.set_xlabel("Parent Molecular Weight (Da)")
    ax2.set_ylabel("Compounds")
    ax2.set_title("Standardized Parent MW Distribution")
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--compounds", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--smiles_col", default=None)
    ap.add_argument("--id_col", default=None)
    ap.add_argument("--keep_duplicates", action="store_true", default=False)
    args = ap.parse_args()

    records = read_compounds(args.compounds, smiles_col=args.smiles_col, id_col=args.id_col)
    if not records:
        fail("No compound records found in input.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    valid_rows: list[dict[str, Any]] = []
    invalid_records: list[dict[str, Any]] = []
    n_salts = 0
    n_neutralized = 0
    n_stereo = 0

    for idx, rec in enumerate(records, start=1):
        cid = rec["compound_id"]
        raw_smi = rec["smiles"]
        res = standardize_smiles(raw_smi)
        if not res["is_valid"]:
            invalid_records.append(
                {
                    "row_index": idx,
                    "compound_id": cid,
                    "original_smiles": raw_smi,
                    "error": res.get("error", "Invalid SMILES"),
                }
            )
            continue
        if res["had_salt"]:
            n_salts += 1
        if res["neutralized"]:
            n_neutralized += 1
        if res["has_stereochemistry"]:
            n_stereo += 1
        valid_rows.append(
            {
                "compound_id": cid,
                "original_smiles": raw_smi,
                "canonical_smiles": res["canonical_smiles"],
                "structure_key": res["structure_key"],
                "formula": res["formula"],
                "mw": res["mw"],
                "heavy_atoms": res["heavy_atoms"],
                "had_salt": res["had_salt"],
                "stripped_fragments": ";".join(res["stripped_fragments"]) if res["stripped_fragments"] else "none",
                "neutralized": res["neutralized"],
                "has_stereochemistry": res["has_stereochemistry"],
            }
        )

    if not valid_rows:
        fail(f"All {len(records)} input SMILES failed validation.")

    if invalid_records:
        warnings.append(f"{len(invalid_records)} of {len(records)} input records had invalid SMILES and were excluded.")

    # Deduplication by canonical structure_key
    key_groups: dict[str, list[dict[str, Any]]] = {}
    for row in valid_rows:
        key_groups.setdefault(row["structure_key"], []).append(row)

    n_duplicates = sum(len(grp) - 1 for grp in key_groups.values() if len(grp) > 1)
    if n_duplicates > 0:
        warnings.append(
            f"Detected {n_duplicates} duplicate structure(s) across {len(valid_rows)} valid compounds"
            + (" (retained because keep_duplicates=true)." if args.keep_duplicates else " (collapsed to unique structure keys).")
        )

    out_rows: list[dict[str, Any]] = []
    if args.keep_duplicates:
        seen_keys: set[str] = set()
        for row in valid_rows:
            k = row["structure_key"]
            grp = key_groups[k]
            is_dup = k in seen_keys
            seen_keys.add(k)
            r_copy = dict(row)
            r_copy["is_duplicate"] = is_dup
            r_copy["duplicate_group_size"] = len(grp)
            r_copy["merged_ids"] = ";".join(g["compound_id"] for g in grp)
            out_rows.append(r_copy)
    else:
        seen_keys = set()
        for row in valid_rows:
            k = row["structure_key"]
            if k in seen_keys:
                continue
            seen_keys.add(k)
            grp = key_groups[k]
            r_copy = dict(row)
            r_copy["is_duplicate"] = False
            r_copy["duplicate_group_size"] = len(grp)
            r_copy["merged_ids"] = ";".join(g["compound_id"] for g in grp)
            out_rows.append(r_copy)

    cols = [
        "compound_id",
        "original_smiles",
        "canonical_smiles",
        "structure_key",
        "formula",
        "mw",
        "heavy_atoms",
        "had_salt",
        "stripped_fragments",
        "neutralized",
        "has_stereochemistry",
        "is_duplicate",
        "duplicate_group_size",
        "merged_ids",
    ]
    params = {"keep_duplicates": bool(args.keep_duplicates), "backend": BACKEND}
    tsv_path = out_dir / "standardized_compounds.tsv"
    write_tsv(
        out_rows,
        cols,
        tsv_path,
        tool="chem_standardize",
        title="Standardized Small-Molecule Compound Table",
        inputs=(args.compounds,),
        parameters=params,
        scientific_context={
            "valueSemantics": "salt_stripped_neutralized_canonical_compounds",
            "backend": BACKEND,
            "nValid": len(valid_rows),
            "nOutput": len(out_rows),
        },
    )

    png_path = out_dir / "standardization_summary.png"
    _plot_standardization(
        png_path,
        len(valid_rows),
        n_salts,
        n_neutralized,
        n_duplicates,
        len(invalid_records),
        n_stereo,
        [float(r["mw"]) for r in out_rows],
    )
    write_artifact(
        png_path,
        tool="chem_standardize",
        title="Compound Standardization QC & MW Distribution",
        inputs=(args.compounds,),
        parameters=params,
    )

    emit(
        {
            "status": "ok",
            "tool": "chem_standardize",
            "backend": BACKEND,
            "n_input": len(records),
            "n_valid": len(valid_rows),
            "n_invalid": len(invalid_records),
            "n_salts_stripped": n_salts,
            "n_neutralized": n_neutralized,
            "n_stereochemistry": n_stereo,
            "n_duplicates": n_duplicates,
            "n_output": len(out_rows),
            "keep_duplicates": bool(args.keep_duplicates),
            "invalid_records": invalid_records[:25],
            "outputs": {
                "standardized_compounds": str(tsv_path),
                "standardization_summary": str(png_path),
            },
            "artifacts": [str(tsv_path), str(png_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
