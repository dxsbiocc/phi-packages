#!/usr/bin/env python3
"""chem_sar_cliff: Structure-Activity Relationship (SAR) & Activity Cliff
analysis: normalize potency (nM/uM/pIC50 -> pIC50), compute Ligand Efficiency
(LE) and Lipophilic Ligand Efficiency (LLE), pairwise ECFP4 Tanimoto similarity,
Guha-Van Drie SALI, and Activity Cliffs.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    BACKEND,
    emit,
    fail,
    featurize_smiles,
    read_compounds,
    resolve_column_name,
    setup_matplotlib,
    tanimoto_similarity,
    write_artifact,
    write_fallback_png,
    write_tsv,
)


def _to_pic50(raw_val: float, scale: str) -> float:
    if scale == "pIC50":
        return float(raw_val)
    if raw_val <= 0.0:
        raise ValueError(f"Concentration potency must be > 0 for scale '{scale}'; got {raw_val}")
    if scale == "nM":
        return 9.0 - math.log10(raw_val)
    if scale == "uM":
        return 6.0 - math.log10(raw_val)
    raise ValueError(f"Unsupported potency_scale: {scale}")


def _plot_sar_landscape(
    out_png: Path,
    rows: list[dict[str, Any]],
    pairs: list[dict[str, Any]],
    sim_threshold: float,
) -> None:
    plt = setup_matplotlib()
    if plt is None:
        write_fallback_png(out_png, accent_rgb=(13, 148, 136))
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.4))
    clogps = [float(r["clogp"]) for r in rows]
    pic50s = [float(r["pic50"]) for r in rows]
    les = [float(r["le"]) for r in rows]

    sc = ax1.scatter(clogps, pic50s, c=les, cmap="plasma", s=58, edgecolors="#0f172a", linewidth=0.5)
    cbar = fig.colorbar(sc, ax=ax1)
    cbar.set_label("Ligand Efficiency (LE)")

    x_min, x_max = min(clogps) - 0.6, max(clogps) + 0.6
    xs = [x_min, x_max]
    for lle_iso, ls in ((3.0, ":"), (5.0, "--"), (7.0, "-.")):
        ax1.plot(xs, [x + lle_iso for x in xs], linestyle=ls, color="#64748b", linewidth=1.0, label=f"LLE = {int(lle_iso)}")

    row_by_id = {r["compound_id"]: r for r in rows}
    for cl in [p for p in pairs if p["is_activity_cliff"]][:5]:
        r1 = row_by_id[cl["compound_a"]]
        r2 = row_by_id[cl["compound_b"]]
        ax1.plot([r1["clogp"], r2["clogp"]], [r1["pic50"], r2["pic50"]], color="#dc2626", linewidth=1.2, alpha=0.75)

    ax1.set_xlabel("cLogP")
    ax1.set_ylabel("pIC50")
    ax1.set_title("Lipophilic Ligand Efficiency (LLE) Landscape")
    ax1.legend(frameon=False, loc="best")

    if pairs:
        sims = [float(p["tanimoto_sim"]) for p in pairs]
        dpics = [float(p["delta_pic50"]) for p in pairs]
        is_cliff = [bool(p["is_activity_cliff"]) for p in pairs]
        norm_idx = [i for i, c in enumerate(is_cliff) if not c]
        cliff_idx = [i for i, c in enumerate(is_cliff) if c]
        if norm_idx:
            ax2.scatter([sims[i] for i in norm_idx], [dpics[i] for i in norm_idx],
                        color="#94a3b8", alpha=0.55, s=28, label="Non-cliff pairs")
        if cliff_idx:
            ax2.scatter([sims[i] for i in cliff_idx], [dpics[i] for i in cliff_idx],
                        color="#dc2626", s=54, edgecolors="#7f1d1d", linewidth=0.6, label="Activity Cliffs")
    ax2.axvline(sim_threshold, color="#0f766e", linestyle="--", linewidth=1.1, label=f"Sim >= {sim_threshold:.2f}")
    ax2.axhline(1.0, color="#f59e0b", linestyle="--", linewidth=1.1, label="|dpIC50| >= 1.0")
    ax2.set_xlim(0.0, 1.03)
    ax2.set_xlabel("ECFP4 Tanimoto Similarity")
    ax2.set_ylabel("|Delta pIC50|")
    ax2.set_title("Structure-Activity Landscape (SALI)")
    ax2.legend(frameon=False, loc="best")

    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--compounds", required=True)
    ap.add_argument("--potency_col", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--potency_scale", choices=["nM", "uM", "pIC50"], default="nM")
    ap.add_argument("--smiles_col", default=None)
    ap.add_argument("--id_col", default=None)
    ap.add_argument("--similarity_threshold", type=float, default=0.55)
    args = ap.parse_args()

    if not (0.0 < args.similarity_threshold < 1.0):
        fail(f"similarity_threshold must lie in (0, 1); got {args.similarity_threshold}.")

    records = read_compounds(args.compounds, smiles_col=args.smiles_col, id_col=args.id_col)
    if not records:
        fail("No compound records found in input.")
    cols = list(records[0].keys())
    p_col = resolve_column_name(cols, [args.potency_col], "potency")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    rows: list[dict[str, Any]] = []
    fps: list[int] = []
    n_skipped = 0

    for rec in records:
        raw_p = str(rec.get(p_col, "")).strip()
        try:
            pot_val = float(raw_p)
            if not math.isfinite(pot_val):
                raise ValueError("non-finite")
            pic50 = _to_pic50(pot_val, args.potency_scale)
        except ValueError:
            n_skipped += 1
            continue
        feat = featurize_smiles(rec["smiles"])
        if not feat["is_valid"]:
            n_skipped += 1
            continue
        fp = feat.pop("fp")
        fps.append(fp)
        heavy = max(1, feat["heavy_atoms"])
        le = (1.37 * pic50) / heavy
        lle = pic50 - feat["clogp"]
        rows.append(
            {
                "compound_id": rec["compound_id"],
                "canonical_smiles": feat["canonical_smiles"],
                "raw_potency": pot_val,
                "potency_scale": args.potency_scale,
                "pic50": round(pic50, 4),
                "heavy_atoms": heavy,
                "mw": feat["mw"],
                "clogp": feat["clogp"],
                "tpsa": feat["tpsa"],
                "qed": feat["qed"],
                "le": round(le, 4),
                "lle": round(lle, 4),
                "murcko_scaffold": feat["murcko_scaffold"],
                "n_alerts": feat["n_alerts"],
            }
        )

    if len(rows) < 2:
        fail(f"SAR cliff analysis requires >= 2 valid compounds with potency values (got {len(rows)}).")
    if n_skipped > 0:
        warnings.append(f"Skipped {n_skipped} record(s) with invalid SMILES or invalid potency values.")

    # Pairwise Tanimoto similarity & Guha-Van Drie SALI
    pairs: list[dict[str, Any]] = []
    n = len(rows)
    for i in range(n):
        for j in range(i + 1, n):
            sim = tanimoto_similarity(fps[i], fps[j])
            dp = abs(rows[i]["pic50"] - rows[j]["pic50"])
            sali = dp / (1.0 - sim + 0.001)
            is_cliff = sim >= float(args.similarity_threshold) and dp >= 1.0
            pairs.append(
                {
                    "compound_a": rows[i]["compound_id"],
                    "compound_b": rows[j]["compound_id"],
                    "pic50_a": rows[i]["pic50"],
                    "pic50_b": rows[j]["pic50"],
                    "delta_pic50": round(dp, 4),
                    "tanimoto_sim": round(sim, 4),
                    "sali": round(sali, 3),
                    "is_activity_cliff": is_cliff,
                    "smiles_a": rows[i]["canonical_smiles"],
                    "smiles_b": rows[j]["canonical_smiles"],
                }
            )

    pairs_sorted = sorted(pairs, key=lambda x: (-int(x["is_activity_cliff"]), -x["sali"], -x["delta_pic50"]))
    cliffs = [p for p in pairs_sorted if p["is_activity_cliff"]]

    sar_cols = [
        "compound_id",
        "canonical_smiles",
        "raw_potency",
        "potency_scale",
        "pic50",
        "heavy_atoms",
        "mw",
        "clogp",
        "tpsa",
        "qed",
        "le",
        "lle",
        "murcko_scaffold",
        "n_alerts",
    ]
    cliff_cols = [
        "compound_a",
        "compound_b",
        "pic50_a",
        "pic50_b",
        "delta_pic50",
        "tanimoto_sim",
        "sali",
        "is_activity_cliff",
        "smiles_a",
        "smiles_b",
    ]
    params = {
        "potency_col": p_col,
        "potency_scale": args.potency_scale,
        "similarity_threshold": float(args.similarity_threshold),
        "backend": BACKEND,
    }

    sar_tsv = out_dir / "sar_metrics.tsv"
    write_tsv(
        rows,
        sar_cols,
        sar_tsv,
        tool="chem_sar_cliff",
        title="SAR Potency, Ligand Efficiency (LE) & Lipophilic Efficiency (LLE)",
        inputs=(args.compounds,),
        parameters=params,
        scientific_context={"valueSemantics": "sar_pic50_le_lle_metrics", "backend": BACKEND},
    )

    cliff_tsv = out_dir / "activity_cliffs.tsv"
    write_tsv(
        pairs_sorted,
        cliff_cols,
        cliff_tsv,
        tool="chem_sar_cliff",
        title="Pairwise ECFP4 Similarity & SALI Activity Cliffs",
        inputs=(args.compounds,),
        parameters=params,
    )

    png_path = out_dir / "sar_lle_cliff_landscape.png"
    _plot_sar_landscape(png_path, rows, pairs_sorted, float(args.similarity_threshold))
    write_artifact(
        png_path,
        tool="chem_sar_cliff",
        title="LLE Landscape & SALI Activity Cliff Plot",
        inputs=(args.compounds,),
        parameters=params,
    )

    max_sali = max((p["sali"] for p in pairs_sorted), default=None)
    top_cliffs = [
        {
            "compound_a": c["compound_a"],
            "compound_b": c["compound_b"],
            "delta_pic50": c["delta_pic50"],
            "tanimoto_sim": c["tanimoto_sim"],
            "sali": c["sali"],
        }
        for c in cliffs[:10]
    ]

    emit(
        {
            "status": "ok",
            "tool": "chem_sar_cliff",
            "backend": BACKEND,
            "potency_col": p_col,
            "potency_scale": args.potency_scale,
            "similarity_threshold": float(args.similarity_threshold),
            "n_compounds": len(rows),
            "n_pairs": len(pairs_sorted),
            "n_activity_cliffs": len(cliffs),
            "max_sali": max_sali,
            "mean_pic50": round(sum(r["pic50"] for r in rows) / len(rows), 4),
            "mean_le": round(sum(r["le"] for r in rows) / len(rows), 4),
            "mean_lle": round(sum(r["lle"] for r in rows) / len(rows), 4),
            "top_cliffs": top_cliffs,
            "outputs": {
                "sar_metrics": str(sar_tsv),
                "activity_cliffs": str(cliff_tsv),
                "sar_lle_cliff_landscape": str(png_path),
            },
            "artifacts": [str(sar_tsv), str(cliff_tsv), str(png_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
