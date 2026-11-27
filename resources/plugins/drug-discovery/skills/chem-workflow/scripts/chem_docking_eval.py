#!/usr/bin/env python3
"""chem_docking_eval: Virtual-screening docking & co-folding score evaluation:
join ligand SMILES, collapse multi-pose tables to best pose per compound, filter
PAINS/Brenk alerts, compute size-normalized LE_dock and SILE, optionally fuse
auxiliary score via consensus ranking, and select cluster-diverse top hits.
"""
from __future__ import annotations

import argparse
import math
import os
import statistics
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
    read_table_rows,
    resolve_column_name,
    setup_matplotlib,
    tanimoto_similarity,
    write_artifact,
    write_fallback_png,
    write_tsv,
)


def _z_scores(vals: list[float], higher_is_better: bool) -> list[float]:
    if not vals:
        return []
    signed = [v if higher_is_better else -v for v in vals]
    mu = sum(signed) / len(signed)
    sd = statistics.pstdev(signed) if len(signed) > 1 else 0.0
    if sd <= 1e-12:
        return [0.0 for _ in signed]
    return [(x - mu) / sd for x in signed]


def _rank_scores(vals: list[float], higher_is_better: bool) -> list[int]:
    order = sorted(range(len(vals)), key=lambda i: (-vals[i] if higher_is_better else vals[i], i))
    ranks = [0] * len(vals)
    for r, idx in enumerate(order, start=1):
        ranks[idx] = r
    return ranks


def _plot_docking_pareto(out_png: Path, rows: list[dict[str, Any]]) -> None:
    plt = setup_matplotlib()
    if plt is None:
        write_fallback_png(out_png, accent_rgb=(99, 102, 241))
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.4))
    sel = [r for r in rows if r["selected_diverse_hit"]]
    oth_clean = [r for r in rows if r["passes_medchem_filter"] and not r["selected_diverse_hit"]]
    pains = [r for r in rows if not r["passes_medchem_filter"]]

    if oth_clean:
        ax1.scatter([r["heavy_atoms"] for r in oth_clean], [r["fav_score"] for r in oth_clean],
                    color="#94a3b8", s=40, alpha=0.75, label="Clean Non-Selected")
    if pains:
        ax1.scatter([r["heavy_atoms"] for r in pains], [r["fav_score"] for r in pains],
                    color="#dc2626", marker="X", s=64, label=f"PAINS/Brenk Filtered ({len(pains)})")
    if sel:
        ax1.scatter([r["heavy_atoms"] for r in sel], [r["fav_score"] for r in sel],
                    color="#0f766e", s=68, edgecolors="#0f172a", linewidth=0.6, label=f"Selected Diverse Hits ({len(sel)})", zorder=4)

    max_ha = max(r["heavy_atoms"] for r in rows)
    xs = [5, max(10, max_ha + 2)]
    ax1.plot(xs, [0.30 * x for x in xs], color="#0284c7", linestyle="--", linewidth=1.1, label="LE_dock = 0.30")
    ax1.set_xlabel("Heavy Atom Count (N_heavy)")
    ax1.set_ylabel("Favorable Binding Magnitude")
    ax1.set_title("Docking Score vs Molecular Size (LE Pareto)")
    ax1.legend(frameon=False, fontsize=8.0, loc="best")

    top_plot = [r for r in rows if r["passes_medchem_filter"]][:15]
    if top_plot:
        ids = [r["compound_id"] for r in top_plot][::-1]
        siles = [r["sile"] for r in top_plot][::-1]
        cols = ["#0f766e" if r["selected_diverse_hit"] else "#64748b" for r in top_plot][::-1]
        ax2.barh(ids, siles, color=cols, height=0.65, edgecolor="#1e293b", linewidth=0.4)
    ax2.set_xlabel("Size-Independent Ligand Efficiency (SILE)")
    ax2.set_title("Top Clean Virtual-Screening Hits by SILE")

    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--docking_scores", required=True)
    ap.add_argument("--score_col", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--compounds", default=None)
    ap.add_argument("--higher_is_better", action="store_true", default=False)
    ap.add_argument("--aux_score_col", default=None)
    ap.add_argument("--top_n", type=int, default=25)
    args = ap.parse_args()

    if args.top_n < 1:
        fail(f"top_n must be >= 1; got {args.top_n}.")

    cols, raw_rows = read_table_rows(args.docking_scores)
    id_col = resolve_column_name(cols, ["compound_id", "id", "ligand_id", "name", "molecule_id", "mol_id", "title"], "compound ID")
    score_col = resolve_column_name(cols, [args.score_col], "docking score")
    smi_col = resolve_column_name(
        cols,
        ["smiles", "canonical_smiles", "smi", "ligand_smiles", "structure"],
        "SMILES",
        required=False,
    )
    aux_col = resolve_column_name(cols, [args.aux_score_col], "auxiliary score") if args.aux_score_col else None

    smiles_lookup: dict[str, str] = {}
    if args.compounds:
        for rec in read_compounds(args.compounds):
            if rec["smiles"]:
                smiles_lookup[rec["compound_id"]] = rec["smiles"]

    if smi_col is None and not smiles_lookup:
        fail("Docking scores table has no SMILES column and no --compounds table was provided.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    # Group by compound_id and keep best pose per compound
    best_by_id: dict[str, dict[str, Any]] = {}
    pose_counts: dict[str, int] = {}
    n_skipped_score = 0
    for r in raw_rows:
        cid = r.get(id_col, "").strip()
        if not cid:
            continue
        try:
            sc = float(str(r.get(score_col, "")).strip())
            if not math.isfinite(sc):
                raise ValueError
        except ValueError:
            n_skipped_score += 1
            continue
        aux_val: float | None = None
        if aux_col is not None:
            try:
                av = float(str(r.get(aux_col, "")).strip())
                if math.isfinite(av):
                    aux_val = av
            except ValueError:
                aux_val = None
        smi = (r.get(smi_col, "").strip() if smi_col else "") or smiles_lookup.get(cid, "")
        pose_counts[cid] = pose_counts.get(cid, 0) + 1
        prev = best_by_id.get(cid)
        is_better = (
            prev is None
            or (sc > prev["score"] if args.higher_is_better else sc < prev["score"])
        )
        if is_better:
            best_by_id[cid] = {"compound_id": cid, "score": sc, "aux_score": aux_val, "smiles": smi}

    if n_skipped_score > 0:
        warnings.append(f"Skipped {n_skipped_score} row(s) with non-numeric docking scores.")
    if not best_by_id:
        fail("No valid docking score rows found.")

    eval_rows: list[dict[str, Any]] = []
    fps: list[int] = []
    n_missing_struct = 0
    for cid, item in best_by_id.items():
        if not item["smiles"]:
            n_missing_struct += 1
            continue
        feat = featurize_smiles(item["smiles"])
        if not feat["is_valid"]:
            n_missing_struct += 1
            continue
        fp = feat.pop("fp")
        fps.append(fp)
        heavy = max(1, feat["heavy_atoms"])
        fav = item["score"] if args.higher_is_better else -item["score"]
        le_dock = fav / heavy
        sile = fav / (heavy ** 0.3)
        eval_rows.append(
            {
                "compound_id": cid,
                "canonical_smiles": feat["canonical_smiles"],
                "n_poses": pose_counts[cid],
                "score": round(item["score"], 4),
                "aux_score": round(item["aux_score"], 4) if item["aux_score"] is not None else "",
                "fav_score": round(fav, 4),
                "heavy_atoms": heavy,
                "mw": feat["mw"],
                "clogp": feat["clogp"],
                "qed": feat["qed"],
                "le_dock": round(le_dock, 4),
                "sile": round(sile, 4),
                "pains_alerts": ";".join(feat["pains_alerts"]) if feat["pains_alerts"] else "none",
                "brenk_alerts": ";".join(feat["brenk_alerts"]) if feat["brenk_alerts"] else "none",
                "n_alerts": feat["n_alerts"],
                "passes_medchem_filter": feat["n_alerts"] == 0,
                "murcko_scaffold": feat["murcko_scaffold"],
                "_fp": fp,
            }
        )

    if not eval_rows:
        fail("No compounds had both valid docking scores and valid ligand SMILES.")
    if n_missing_struct > 0:
        warnings.append(f"Excluded {n_missing_struct} compound(s) lacking valid SMILES structures.")

    # Compute primary rank & optional consensus scoring
    scores = [float(r["score"]) for r in eval_rows]
    p_ranks = _rank_scores(scores, bool(args.higher_is_better))
    p_z = _z_scores(scores, bool(args.higher_is_better))

    has_aux = aux_col is not None and all(r["aux_score"] != "" for r in eval_rows)
    if aux_col is not None and not has_aux:
        warnings.append(f"Auxiliary score '{aux_col}' missing on some rows; using primary score for ranking.")

    if has_aux:
        aux_vals = [float(r["aux_score"]) for r in eval_rows]
        aux_hib = any(v > 0.0 for v in aux_vals)
        a_ranks = _rank_scores(aux_vals, aux_hib)
        a_z = _z_scores(aux_vals, aux_hib)
        for i, r in enumerate(eval_rows):
            r["primary_rank"] = p_ranks[i]
            r["consensus_z"] = round(0.5 * (p_z[i] + a_z[i]), 4)
            r["consensus_rank_avg"] = round(0.5 * (p_ranks[i] + a_ranks[i]), 2)
        eval_rows.sort(key=lambda r: (-int(r["passes_medchem_filter"]), r["consensus_rank_avg"], -r["consensus_z"]))
    else:
        for i, r in enumerate(eval_rows):
            r["primary_rank"] = p_ranks[i]
            r["consensus_z"] = round(p_z[i], 4)
            r["consensus_rank_avg"] = float(p_ranks[i])
        eval_rows.sort(key=lambda r: (-int(r["passes_medchem_filter"]), r["primary_rank"], -r["sile"]))

    # Diverse hit selection among clean compounds (Tanimoto < 0.65 first, then backfill)
    selected_fps: list[int] = []
    selected_ids: set[str] = set()
    clean_rows = [r for r in eval_rows if r["passes_medchem_filter"]]
    for r in clean_rows:
        if len(selected_ids) >= args.top_n:
            break
        fp = r["_fp"]
        if all(tanimoto_similarity(fp, sfp) < 0.65 for sfp in selected_fps):
            selected_fps.append(fp)
            selected_ids.add(r["compound_id"])
    # Backfill if fewer than top_n distinct clusters exist
    if len(selected_ids) < min(args.top_n, len(clean_rows)):
        for r in clean_rows:
            if len(selected_ids) >= args.top_n:
                break
            if r["compound_id"] not in selected_ids:
                selected_ids.add(r["compound_id"])

    for idx, r in enumerate(eval_rows, start=1):
        r["eval_rank"] = idx
        r["selected_diverse_hit"] = r["compound_id"] in selected_ids
        r.pop("_fp", None)

    n_pains_filtered = sum(1 for r in eval_rows if not r["passes_medchem_filter"])
    if n_pains_filtered > 0:
        warnings.append(f"Flagged and deprioritized {n_pains_filtered} compound(s) with PAINS/Brenk alerts.")

    out_cols = [
        "eval_rank",
        "compound_id",
        "canonical_smiles",
        "n_poses",
        "score",
        "aux_score",
        "primary_rank",
        "consensus_z",
        "consensus_rank_avg",
        "heavy_atoms",
        "mw",
        "clogp",
        "qed",
        "le_dock",
        "sile",
        "pains_alerts",
        "brenk_alerts",
        "n_alerts",
        "passes_medchem_filter",
        "murcko_scaffold",
        "selected_diverse_hit",
    ]
    params = {
        "score_col": score_col,
        "aux_score_col": aux_col,
        "higher_is_better": bool(args.higher_is_better),
        "top_n": int(args.top_n),
        "backend": BACKEND,
    }
    tsv_path = out_dir / "docking_hit_list.tsv"
    write_tsv(
        eval_rows,
        out_cols,
        tsv_path,
        tool="chem_docking_eval",
        title="Virtual Screening Docking Hit List (LE_dock, SILE & PAINS-Filtered)",
        inputs=(args.docking_scores, args.compounds),
        parameters=params,
        scientific_context={"valueSemantics": "docking_le_sile_diverse_hits", "backend": BACKEND},
    )

    png_path = out_dir / "docking_le_pareto.png"
    _plot_docking_pareto(png_path, eval_rows)
    write_artifact(
        png_path,
        tool="chem_docking_eval",
        title="Docking Ligand Efficiency Pareto & Top SILE Hits",
        inputs=(args.docking_scores, args.compounds),
        parameters=params,
    )

    sel_hits = [r for r in eval_rows if r["selected_diverse_hit"]]
    best_sc = max(r["score"] for r in eval_rows) if args.higher_is_better else min(r["score"] for r in eval_rows)

    emit(
        {
            "status": "ok",
            "tool": "chem_docking_eval",
            "backend": BACKEND,
            "score_col": score_col,
            "aux_score_col": aux_col,
            "higher_is_better": bool(args.higher_is_better),
            "n_input_rows": len(raw_rows),
            "n_unique_compounds": len(eval_rows),
            "n_pains_filtered": n_pains_filtered,
            "n_selected_hits": len(sel_hits),
            "top_n": int(args.top_n),
            "best_score": round(best_sc, 4),
            "mean_le_dock": round(sum(r["le_dock"] for r in eval_rows) / len(eval_rows), 4),
            "mean_sile": round(sum(r["sile"] for r in eval_rows) / len(eval_rows), 4),
            "top_hits": [
                {
                    "eval_rank": r["eval_rank"],
                    "compound_id": r["compound_id"],
                    "score": r["score"],
                    "le_dock": r["le_dock"],
                    "sile": r["sile"],
                    "qed": r["qed"],
                    "consensus_z": r["consensus_z"],
                }
                for r in sel_hits[:10]
            ],
            "outputs": {
                "docking_hit_list": str(tsv_path),
                "docking_le_pareto": str(png_path),
            },
            "artifacts": [str(tsv_path), str(png_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
