#!/usr/bin/env python3
"""chem_medchem_triage: Medicinal chemistry triage with Lipinski Ro5, Veber,
QED, synthetic complexity (1-10), PAINS & Brenk structural alerts, Bemis-Murcko
scaffolds, and Taylor-Butina ECFP4 diversity clustering.
"""
from __future__ import annotations

import argparse
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
    setup_matplotlib,
    taylor_butina_cluster,
    write_artifact,
    write_fallback_png,
    write_tsv,
)


def _plot_chemical_space(out_png: Path, rows: list[dict[str, Any]], tier_counts: dict[str, int]) -> None:
    plt = setup_matplotlib()
    if plt is None:
        write_fallback_png(out_png, accent_rgb=(2, 132, 199))
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
    mws = [float(r["mw"]) for r in rows]
    logps = [float(r["clogp"]) for r in rows]
    qeds = [float(r["qed"]) for r in rows]
    flagged = [r["n_alerts"] > 0 for r in rows]

    clean_idx = [i for i, f in enumerate(flagged) if not f]
    flag_idx = [i for i, f in enumerate(flagged) if f]
    if clean_idx:
        sc = ax1.scatter(
            [mws[i] for i in clean_idx],
            [logps[i] for i in clean_idx],
            c=[qeds[i] for i in clean_idx],
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
            s=52,
            alpha=0.88,
            edgecolors="#0f172a",
            linewidth=0.5,
            label="Clean",
        )
        cbar = fig.colorbar(sc, ax=ax1)
        cbar.set_label("QED")
    if flag_idx:
        ax1.scatter(
            [mws[i] for i in flag_idx],
            [logps[i] for i in flag_idx],
            marker="X",
            s=72,
            color="#dc2626",
            edgecolors="#7f1d1d",
            linewidth=0.6,
            label="PAINS/Brenk Alert",
        )
    ax1.axvline(500.0, color="#64748b", linestyle="--", linewidth=1.0, label="Ro5 MW=500")
    ax1.axhline(5.0, color="#94a3b8", linestyle=":", linewidth=1.0, label="Ro5 cLogP=5")
    ax1.set_xlabel("Molecular Weight (Da)")
    ax1.set_ylabel("cLogP")
    ax1.set_title("Chemical Space & QED Desirability")
    ax1.legend(frameon=False, loc="best")

    tiers = ["lead_like", "drug_like", "alert_flagged", "beyond_ro5"]
    t_vals = [tier_counts.get(t, 0) for t in tiers]
    t_cols = ["#0f766e", "#0284c7", "#dc2626", "#f59e0b"]
    bars = ax2.bar(tiers, t_vals, color=t_cols, width=0.58, edgecolor="#1e293b", linewidth=0.5)
    for b, v in zip(bars, t_vals):
        ax2.text(b.get_x() + b.get_width() / 2.0, b.get_height() + max(t_vals, default=1) * 0.02, str(v),
                 ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    ax2.set_ylabel("Compound Count")
    ax2.set_title("Medicinal Chemistry Triage Tiers")
    ax2.set_ylim(0, max(1, max(t_vals)) * 1.18)

    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--compounds", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--smiles_col", default=None)
    ap.add_argument("--id_col", default=None)
    ap.add_argument("--cluster_cutoff", type=float, default=0.4)
    args = ap.parse_args()

    if not (0.0 < args.cluster_cutoff < 1.0):
        fail(f"cluster_cutoff must lie strictly in (0, 1); got {args.cluster_cutoff}.")

    records = read_compounds(args.compounds, smiles_col=args.smiles_col, id_col=args.id_col)
    if not records:
        fail("No compound records found in input.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    rows: list[dict[str, Any]] = []
    fps: list[int] = []
    n_invalid = 0
    for rec in records:
        feat = featurize_smiles(rec["smiles"])
        if not feat["is_valid"]:
            n_invalid += 1
            continue
        fp = feat.pop("fp")
        fps.append(fp)
        ro5_viol = (
            (1 if feat["mw"] > 500.0 else 0)
            + (1 if feat["clogp"] > 5.0 else 0)
            + (1 if feat["hbd"] > 5 else 0)
            + (1 if feat["hba"] > 10 else 0)
        )
        veber_pass = feat["rotatable_bonds"] <= 10 and feat["tpsa"] <= 140.0
        if feat["n_alerts"] > 0:
            tier = "alert_flagged"
        elif ro5_viol >= 2:
            tier = "beyond_ro5"
        elif feat["mw"] <= 350.0 and feat["clogp"] <= 3.5 and feat["rotatable_bonds"] <= 7 and feat["qed"] >= 0.5:
            tier = "lead_like"
        else:
            tier = "drug_like"

        rows.append(
            {
                "compound_id": rec["compound_id"],
                "canonical_smiles": feat["canonical_smiles"],
                "structure_key": feat["structure_key"],
                "formula": feat["formula"],
                "mw": feat["mw"],
                "heavy_atoms": feat["heavy_atoms"],
                "clogp": feat["clogp"],
                "tpsa": feat["tpsa"],
                "hbd": feat["hbd"],
                "hba": feat["hba"],
                "rotatable_bonds": feat["rotatable_bonds"],
                "fsp3": feat["fsp3"],
                "ring_count": feat["ring_count"],
                "aromatic_ring_count": feat["aromatic_ring_count"],
                "qed": feat["qed"],
                "sa_score": feat["sa_score"],
                "ro5_violations": ro5_viol,
                "veber_pass": veber_pass,
                "pains_alerts": ";".join(feat["pains_alerts"]) if feat["pains_alerts"] else "none",
                "brenk_alerts": ";".join(feat["brenk_alerts"]) if feat["brenk_alerts"] else "none",
                "n_alerts": feat["n_alerts"],
                "murcko_scaffold": feat["murcko_scaffold"],
                "triage_tier": tier,
            }
        )

    if not rows:
        fail("All input compounds had invalid SMILES.")
    if n_invalid > 0:
        warnings.append(f"Skipped {n_invalid} compound(s) with invalid SMILES.")

    cluster_ids, is_centroid = taylor_butina_cluster(fps, distance_cutoff=float(args.cluster_cutoff))
    for r, cid, is_c in zip(rows, cluster_ids, is_centroid):
        r["cluster_id"] = cid
        r["is_cluster_centroid"] = is_c

    # Scaffold summary aggregation
    scaff_groups: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        scaff_groups.setdefault(r["murcko_scaffold"], []).append(r)

    scaffold_rows: list[dict[str, Any]] = []
    for scaff, grp in sorted(scaff_groups.items(), key=lambda kv: (-len(kv[1]), -sum(x["qed"] for x in kv[1]) / len(kv[1]))):
        best_rep = max(grp, key=lambda x: (x["n_alerts"] == 0, x["qed"]))
        scaffold_rows.append(
            {
                "murcko_scaffold": scaff,
                "n_compounds": len(grp),
                "mean_qed": round(sum(x["qed"] for x in grp) / len(grp), 4),
                "mean_mw": round(sum(x["mw"] for x in grp) / len(grp), 2),
                "mean_clogp": round(sum(x["clogp"] for x in grp) / len(grp), 3),
                "n_lead_like": sum(1 for x in grp if x["triage_tier"] == "lead_like"),
                "n_drug_like": sum(1 for x in grp if x["triage_tier"] == "drug_like"),
                "n_alert_flagged": sum(1 for x in grp if x["triage_tier"] == "alert_flagged"),
                "representative_id": best_rep["compound_id"],
            }
        )

    tier_counts = {
        t: sum(1 for r in rows if r["triage_tier"] == t)
        for t in ("lead_like", "drug_like", "alert_flagged", "beyond_ro5")
    }
    n_ro5_pass = sum(1 for r in rows if r["ro5_violations"] <= 1)
    n_veber_pass = sum(1 for r in rows if r["veber_pass"])
    n_alert_flagged = sum(1 for r in rows if r["n_alerts"] > 0)
    if n_alert_flagged > 0:
        warnings.append(f"{n_alert_flagged} compound(s) triggered PAINS/Brenk structural alerts.")

    triage_cols = [
        "compound_id",
        "canonical_smiles",
        "structure_key",
        "formula",
        "mw",
        "heavy_atoms",
        "clogp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "fsp3",
        "ring_count",
        "aromatic_ring_count",
        "qed",
        "sa_score",
        "ro5_violations",
        "veber_pass",
        "pains_alerts",
        "brenk_alerts",
        "n_alerts",
        "murcko_scaffold",
        "cluster_id",
        "is_cluster_centroid",
        "triage_tier",
    ]
    scaff_cols = [
        "murcko_scaffold",
        "n_compounds",
        "mean_qed",
        "mean_mw",
        "mean_clogp",
        "n_lead_like",
        "n_drug_like",
        "n_alert_flagged",
        "representative_id",
    ]
    params = {"cluster_cutoff": float(args.cluster_cutoff), "backend": BACKEND}

    triage_tsv = out_dir / "medchem_triage.tsv"
    write_tsv(
        rows,
        triage_cols,
        triage_tsv,
        tool="chem_medchem_triage",
        title="Medicinal Chemistry Triage & Diversity Clustering Table",
        inputs=(args.compounds,),
        parameters=params,
        scientific_context={"valueSemantics": "medchem_descriptors_alerts_butina_clusters", "backend": BACKEND},
    )

    scaff_tsv = out_dir / "scaffold_summary.tsv"
    write_tsv(
        scaffold_rows,
        scaff_cols,
        scaff_tsv,
        tool="chem_medchem_triage",
        title="Bemis-Murcko Scaffold Summary",
        inputs=(args.compounds,),
        parameters=params,
    )

    png_path = out_dir / "chemical_space_qed.png"
    _plot_chemical_space(png_path, rows, tier_counts)
    write_artifact(
        png_path,
        tool="chem_medchem_triage",
        title="Chemical Space (MW vs cLogP, QED) & Triage Tier Breakdown",
        inputs=(args.compounds,),
        parameters=params,
    )

    clean_sorted = sorted([r for r in rows if r["n_alerts"] == 0], key=lambda x: (-x["qed"], x["sa_score"]))
    top_candidates = [
        {
            "compound_id": r["compound_id"],
            "qed": r["qed"],
            "mw": r["mw"],
            "clogp": r["clogp"],
            "sa_score": r["sa_score"],
            "triage_tier": r["triage_tier"],
            "cluster_id": r["cluster_id"],
        }
        for r in clean_sorted[:5]
    ]

    emit(
        {
            "status": "ok",
            "tool": "chem_medchem_triage",
            "backend": BACKEND,
            "n_compounds": len(records),
            "n_valid": len(rows),
            "cluster_cutoff": float(args.cluster_cutoff),
            "similarity_cutoff": round(1.0 - float(args.cluster_cutoff), 4),
            "n_clusters": max(cluster_ids),
            "n_scaffolds": len(scaffold_rows),
            "tier_counts": tier_counts,
            "n_ro5_pass": n_ro5_pass,
            "n_veber_pass": n_veber_pass,
            "n_alert_flagged": n_alert_flagged,
            "mean_qed": round(sum(r["qed"] for r in rows) / len(rows), 4),
            "mean_sa_score": round(sum(r["sa_score"] for r in rows) / len(rows), 3),
            "top_candidates": top_candidates,
            "outputs": {
                "medchem_triage": str(triage_tsv),
                "scaffold_summary": str(scaff_tsv),
                "chemical_space_qed": str(png_path),
            },
            "artifacts": [str(triage_tsv), str(scaff_tsv), str(png_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
