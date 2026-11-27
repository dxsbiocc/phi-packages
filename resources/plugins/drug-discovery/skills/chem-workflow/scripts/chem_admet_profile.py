#!/usr/bin/env python3
"""chem_admet_profile: Physicochemical ADMET & safety profiling: BOILED-Egg
(HIA/BBB), Wager CNS MPO (0-6), Delaney ESOL solubility (logS), Pfizer 3/75,
GSK 4/400, and hERG basic-amine lipophilicity liability flagging.
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
    setup_matplotlib,
    write_artifact,
    write_fallback_png,
    write_tsv,
)


def _in_ellipse(x: float, y: float, cx: float, cy: float, a: float, b: float, angle_deg: float) -> bool:
    th = math.radians(angle_deg)
    dx, dy = x - cx, y - cy
    xr = dx * math.cos(th) + dy * math.sin(th)
    yr = -dx * math.sin(th) + dy * math.cos(th)
    return (xr / a) ** 2 + (yr / b) ** 2 <= 1.0


def _boiled_egg_class(tpsa: float, clogp: float) -> tuple[bool, bool, str]:
    """Daina & Zoete (2016) BOILED-Egg HIA (white) & BBB (yolk) ellipses."""
    in_bbb = tpsa >= 0.0 and _in_ellipse(tpsa, clogp, 38.117, 3.177, 41.0305, 2.7785, -0.172)
    in_hia = tpsa >= 0.0 and _in_ellipse(tpsa, clogp, 71.051, 2.292, 71.0405, 4.3700, -1.031)
    if in_bbb:
        return True, True, "BBB+"
    if in_hia:
        return True, False, "HIA+ (BBB-)"
    return False, False, "Poorly Absorbed"


def _ramp_down(val: float, low: float, high: float) -> float:
    if val <= low:
        return 1.0
    if val >= high:
        return 0.0
    return (high - val) / (high - low)


def _cns_mpo_score(mw: float, clogp: float, tpsa: float, hbd: int, has_basic: bool, has_acid: bool) -> float:
    """Wager et al. (2010) 6-component CNS Multiparameter Optimization score (0..6)."""
    s_logp = _ramp_down(clogp, 3.0, 5.0)
    clogd_est = clogp - (1.5 if has_basic else 0.0) - (2.0 if has_acid else 0.0)
    s_logd = _ramp_down(clogd_est, 2.0, 4.0)
    s_mw = _ramp_down(mw, 360.0, 500.0)
    if tpsa <= 20.0 or tpsa >= 120.0:
        s_tpsa = 0.0
    elif tpsa < 40.0:
        s_tpsa = (tpsa - 20.0) / 20.0
    elif tpsa <= 90.0:
        s_tpsa = 1.0
    else:
        s_tpsa = (120.0 - tpsa) / 30.0
    s_hbd = _ramp_down(float(hbd), 0.5, 3.5)
    s_pka = 0.55 if has_basic else 1.0
    return round(s_logp + s_logd + s_mw + s_tpsa + s_hbd + s_pka, 3)


def _delaney_esol(mw: float, clogp: float, rot_bonds: int, aromatic_atoms: int, heavy_atoms: int) -> tuple[float, str]:
    ap = aromatic_atoms / max(1, heavy_atoms)
    logs = 0.16 - 0.63 * clogp - 0.0062 * mw + 0.066 * rot_bonds - 0.74 * ap
    if logs > -2.0:
        cls = "Highly soluble"
    elif logs >= -4.0:
        cls = "Soluble"
    elif logs >= -6.0:
        cls = "Moderately soluble"
    else:
        cls = "Poorly soluble"
    return round(logs, 3), cls


def _ellipse_points(cx: float, cy: float, a: float, b: float, angle_deg: float, n_pts: int = 200) -> tuple[list[float], list[float]]:
    th = math.radians(angle_deg)
    xs, ys = [], []
    for k in range(n_pts + 1):
        t = 2.0 * math.pi * k / n_pts
        px = a * math.cos(t)
        py = b * math.sin(t)
        xs.append(cx + px * math.cos(th) - py * math.sin(th))
        ys.append(cy + px * math.sin(th) + py * math.cos(th))
    return xs, ys


def _plot_admet(out_png: Path, rows: list[dict[str, Any]]) -> None:
    plt = setup_matplotlib()
    if plt is None:
        write_fallback_png(out_png, accent_rgb=(217, 119, 6))
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.4))
    ax1.set_facecolor("#e2e8f0")
    wx, wy = _ellipse_points(71.051, 2.292, 71.0405, 4.3700, -1.031)
    yx, yy = _ellipse_points(38.117, 3.177, 41.0305, 2.7785, -0.172)
    ax1.fill(wx, wy, facecolor="#ffffff", edgecolor="#475569", linewidth=1.2, label="HIA (White)")
    ax1.fill(yx, yy, facecolor="#fde047", edgecolor="#ca8a04", linewidth=1.2, label="BBB (Yolk)")

    cls_colors = {"BBB+": "#0f766e", "HIA+ (BBB-)": "#0284c7", "Poorly Absorbed": "#dc2626"}
    for cls, col in cls_colors.items():
        sub = [r for r in rows if r["boiled_egg_class"] == cls]
        if sub:
            ax1.scatter([r["tpsa"] for r in sub], [r["clogp"] for r in sub],
                        color=col, s=52, edgecolors="#0f172a", linewidth=0.5, label=f"{cls} ({len(sub)})", zorder=4)
    ax1.set_xlim(-5, max(165.0, max(r["tpsa"] for r in rows) + 15.0))
    ax1.set_ylim(min(-3.0, min(r["clogp"] for r in rows) - 0.8), max(7.5, max(r["clogp"] for r in rows) + 0.8))
    ax1.set_xlabel("TPSA (A^2)")
    ax1.set_ylabel("cLogP")
    ax1.set_title("BOILED-Egg Permeation Map")
    ax1.legend(frameon=True, facecolor="white", framealpha=0.85, fontsize=8.0, loc="best")

    logs_vals = [r["esol_logs"] for r in rows]
    mpo_vals = [r["cns_mpo"] for r in rows]
    flagged = [r["pfizer_3_75_risk"] or r["gsk_4_400_risk"] or r["herg_risk"] for r in rows]
    clean_i = [i for i, f in enumerate(flagged) if not f]
    risk_i = [i for i, f in enumerate(flagged) if f]
    if clean_i:
        ax2.scatter([logs_vals[i] for i in clean_i], [mpo_vals[i] for i in clean_i],
                    color="#0f766e", s=54, edgecolors="#0f172a", linewidth=0.5, label="Low ADMET Risk")
    if risk_i:
        ax2.scatter([logs_vals[i] for i in risk_i], [mpo_vals[i] for i in risk_i],
                    color="#dc2626", marker="^", s=64, edgecolors="#7f1d1d", linewidth=0.5, label="3/75, 4/400, or hERG Risk")
    ax2.axvline(-4.0, color="#64748b", linestyle="--", linewidth=1.0, label="ESOL Soluble >= -4")
    ax2.axhline(4.0, color="#0284c7", linestyle=":", linewidth=1.1, label="CNS MPO >= 4.0")
    ax2.set_xlabel("Delaney ESOL logS (mol/L)")
    ax2.set_ylabel("Wager CNS MPO (0-6)")
    ax2.set_ylim(-0.2, 6.4)
    ax2.set_title("Solubility vs CNS MPO & Safety Flags")
    ax2.legend(frameon=False, fontsize=8.0, loc="best")

    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--compounds", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--smiles_col", default=None)
    ap.add_argument("--id_col", default=None)
    args = ap.parse_args()

    records = read_compounds(args.compounds, smiles_col=args.smiles_col, id_col=args.id_col)
    if not records:
        fail("No compound records found in input.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    rows: list[dict[str, Any]] = []
    n_invalid = 0
    for rec in records:
        feat = featurize_smiles(rec["smiles"])
        if not feat["is_valid"]:
            n_invalid += 1
            continue
        feat.pop("fp", None)
        hia_pass, bbb_pass, egg_cls = _boiled_egg_class(feat["tpsa"], feat["clogp"])
        cns_mpo = _cns_mpo_score(
            feat["mw"], feat["clogp"], feat["tpsa"], feat["hbd"], feat["has_basic_amine"], feat["has_acidic_group"]
        )
        esol_logs, esol_cls = _delaney_esol(
            feat["mw"], feat["clogp"], feat["rotatable_bonds"], feat["aromatic_atoms"], feat["heavy_atoms"]
        )
        pfizer_risk = feat["clogp"] > 3.0 and feat["tpsa"] < 75.0
        gsk_risk = feat["mw"] > 400.0 and feat["clogp"] > 4.0
        herg_risk = bool(feat["has_basic_amine"] and feat["clogp"] > 3.7)

        rows.append(
            {
                "compound_id": rec["compound_id"],
                "canonical_smiles": feat["canonical_smiles"],
                "mw": feat["mw"],
                "clogp": feat["clogp"],
                "tpsa": feat["tpsa"],
                "hbd": feat["hbd"],
                "hba": feat["hba"],
                "rotatable_bonds": feat["rotatable_bonds"],
                "has_basic_amine": feat["has_basic_amine"],
                "boiled_egg_hia": hia_pass,
                "boiled_egg_bbb": bbb_pass,
                "boiled_egg_class": egg_cls,
                "cns_mpo": cns_mpo,
                "cns_mpo_pass": cns_mpo >= 4.0,
                "esol_logs": esol_logs,
                "esol_solubility_class": esol_cls,
                "pfizer_3_75_risk": pfizer_risk,
                "gsk_4_400_risk": gsk_risk,
                "herg_risk": herg_risk,
            }
        )

    if not rows:
        fail("All input compounds had invalid SMILES.")
    if n_invalid > 0:
        warnings.append(f"Skipped {n_invalid} compound(s) with invalid SMILES.")

    cols = [
        "compound_id",
        "canonical_smiles",
        "mw",
        "clogp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "has_basic_amine",
        "boiled_egg_hia",
        "boiled_egg_bbb",
        "boiled_egg_class",
        "cns_mpo",
        "cns_mpo_pass",
        "esol_logs",
        "esol_solubility_class",
        "pfizer_3_75_risk",
        "gsk_4_400_risk",
        "herg_risk",
    ]
    params = {"backend": BACKEND}

    tsv_path = out_dir / "admet_profile.tsv"
    write_tsv(
        rows,
        cols,
        tsv_path,
        tool="chem_admet_profile",
        title="Physicochemical ADMET, BOILED-Egg, ESOL & Safety Profile",
        inputs=(args.compounds,),
        parameters=params,
        scientific_context={"valueSemantics": "admet_physicochemical_heuristics", "backend": BACKEND},
    )

    png_path = out_dir / "boiled_egg_admet.png"
    _plot_admet(png_path, rows)
    write_artifact(
        png_path,
        tool="chem_admet_profile",
        title="BOILED-Egg Permeation & ESOL/CNS-MPO ADMET Landscape",
        inputs=(args.compounds,),
        parameters=params,
    )

    egg_counts = {
        "bbb_permeant": sum(1 for r in rows if r["boiled_egg_class"] == "BBB+"),
        "hia_absorbed": sum(1 for r in rows if r["boiled_egg_class"] == "HIA+ (BBB-)"),
        "poorly_absorbed": sum(1 for r in rows if r["boiled_egg_class"] == "Poorly Absorbed"),
    }
    esol_counts: dict[str, int] = {}
    for r in rows:
        esol_counts[r["esol_solubility_class"]] = esol_counts.get(r["esol_solubility_class"], 0) + 1

    n_herg = sum(1 for r in rows if r["herg_risk"])
    if n_herg > 0:
        warnings.append(f"{n_herg} compound(s) flagged for potential hERG liability (basic amine + cLogP > 3.7).")

    emit(
        {
            "status": "ok",
            "tool": "chem_admet_profile",
            "backend": BACKEND,
            "n_compounds": len(records),
            "n_valid": len(rows),
            "boiled_egg_counts": egg_counts,
            "n_cns_mpo_pass": sum(1 for r in rows if r["cns_mpo_pass"]),
            "esol_class_counts": esol_counts,
            "n_pfizer_3_75_flagged": sum(1 for r in rows if r["pfizer_3_75_risk"]),
            "n_gsk_4_400_flagged": sum(1 for r in rows if r["gsk_4_400_risk"]),
            "n_herg_flagged": n_herg,
            "mean_esol_logs": round(sum(r["esol_logs"] for r in rows) / len(rows), 3),
            "mean_cns_mpo": round(sum(r["cns_mpo"] for r in rows) / len(rows), 3),
            "outputs": {
                "admet_profile": str(tsv_path),
                "boiled_egg_admet": str(png_path),
            },
            "artifacts": [str(tsv_path), str(png_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
