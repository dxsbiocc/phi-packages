#!/usr/bin/env python3
"""imm_clonotype_tracking: Compute pairwise Morisita-Horn and Jaccard CDR3AA+V
overlap matrices, identify public clonotypes and enriched 4-mer CDR3 motifs,
and run two-sided Fisher's exact tests + BH FDR to call significantly expanded
vs contracted clonotypes between baseline and follow-up samples.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    bh_fdr,
    emit,
    fail,
    harmonize_repertoire_file,
    normalize_gene_base,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)


def _morisita_horn(p_a: dict[str, float], p_b: dict[str, float]) -> float:
    denom = sum(v * v for v in p_a.values()) + sum(v * v for v in p_b.values())
    if denom <= 0:
        return 0.0
    shared = set(p_a.keys()) & set(p_b.keys())
    num = 2.0 * sum(p_a[k] * p_b[k] for k in shared)
    return float(np.clip(num / denom, 0.0, 1.0))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repertoire", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--baseline_sample", default=None)
    ap.add_argument("--followup_sample", default=None)
    ap.add_argument("--locus", default="TRB")
    ap.add_argument("--top_n", type=int, default=20)
    args = ap.parse_args()

    df, detected_fmt = harmonize_repertoire_file(args.repertoire, fmt="auto")
    warnings: list[str] = []

    target_locus = args.locus.strip().upper()
    if target_locus != "ALL":
        df = df[df["locus"] == target_locus].copy()
    if df.empty:
        fail(f"No records found for locus '{args.locus}'.")

    prod = df[(df["productive"]) & (df["cdr3_aa"].astype(str).str.len() > 0) & (df["duplicate_count"] > 0)].copy()
    if prod.empty:
        fail(f"No productive clonotypes with positive counts remain for locus '{args.locus}'.")

    prod["v_gene"] = prod["v_call"].map(normalize_gene_base)
    prod["clonotype_id"] = np.where(
        prod["v_gene"].astype(str).str.len() > 0,
        prod["cdr3_aa"].astype(str) + "|" + prod["v_gene"].astype(str),
        prod["cdr3_aa"].astype(str),
    )

    sample_ids = list(dict.fromkeys(prod["sample_id"].astype(str).tolist()))
    n_samples = len(sample_ids)
    if n_samples < 2:
        fail(f"Clonotype tracking requires >= 2 samples (found {n_samples}: {sample_ids}).")

    agg = (
        prod.groupby(["sample_id", "clonotype_id", "cdr3_aa", "v_gene"], as_index=False, sort=False)["duplicate_count"]
        .sum()
    )

    sample_counts: dict[str, dict[str, float]] = {}
    sample_freqs: dict[str, dict[str, float]] = {}
    sample_totals: dict[str, float] = {}
    clono_meta: dict[str, tuple[str, str]] = {}

    for sid in sample_ids:
        sub = agg[agg["sample_id"] == sid]
        tot = float(sub["duplicate_count"].sum())
        sample_totals[sid] = tot
        cnt_map: dict[str, float] = {}
        freq_map: dict[str, float] = {}
        for _, row in sub.iterrows():
            cid = str(row["clonotype_id"])
            cval = float(row["duplicate_count"])
            cnt_map[cid] = cval
            freq_map[cid] = cval / tot if tot > 0 else 0.0
            clono_meta[cid] = (str(row["cdr3_aa"]), str(row["v_gene"]))
        sample_counts[sid] = cnt_map
        sample_freqs[sid] = freq_map

    # Pairwise Morisita-Horn and Jaccard overlap
    overlap_rows: list[dict[str, Any]] = []
    mh_mat = np.zeros((n_samples, n_samples), dtype=float)
    offdiag_mh: list[float] = []
    offdiag_jac: list[float] = []

    for i, sa in enumerate(sample_ids):
        keys_a = set(sample_counts[sa].keys())
        for j, sb in enumerate(sample_ids):
            keys_b = set(sample_counts[sb].keys())
            inter = len(keys_a & keys_b)
            union = len(keys_a | keys_b)
            jac = float(inter / union) if union > 0 else 0.0
            mh = _morisita_horn(sample_freqs[sa], sample_freqs[sb])
            mh_mat[i, j] = mh
            if i < j:
                offdiag_mh.append(mh)
                offdiag_jac.append(jac)
            overlap_rows.append(
                {
                    "sample_a": sa,
                    "sample_b": sb,
                    "shared_clonotypes": inter,
                    "union_clonotypes": union,
                    "jaccard_index": round(jac, 6),
                    "morisita_horn_index": round(mh, 6),
                }
            )
    overlap_df = pd.DataFrame(overlap_rows)

    # Public / shared clonotypes across >= 2 samples
    all_cids = sorted(clono_meta.keys())
    public_rows: list[dict[str, Any]] = []
    cid_max_freq: dict[str, float] = {}
    cid_n_samples: dict[str, int] = {}
    for cid in all_cids:
        pres_sids = [s for s in sample_ids if cid in sample_counts[s]]
        freqs = [sample_freqs[s].get(cid, 0.0) for s in sample_ids]
        max_f = max(freqs) if freqs else 0.0
        cid_max_freq[cid] = max_f
        cid_n_samples[cid] = len(pres_sids)
        if len(pres_sids) >= 2:
            cdr3_aa, v_gene = clono_meta[cid]
            public_rows.append(
                {
                    "clonotype_id": cid,
                    "cdr3_aa": cdr3_aa,
                    "v_gene": v_gene,
                    "n_samples_shared": len(pres_sids),
                    "samples": ",".join(pres_sids),
                    "mean_frequency": round(float(np.mean([sample_freqs[s][cid] for s in pres_sids])), 6),
                    "max_frequency": round(max_f, 6),
                }
            )
    public_rows.sort(key=lambda r: (-r["n_samples_shared"], -r["mean_frequency"]))
    top_public = public_rows[: max(1, args.top_n)]

    # Enriched 4-mer CDR3 motifs in shared/expanded vs singleton/low clones
    fg_cids = [c for c in all_cids if cid_n_samples[c] >= 2 or cid_max_freq[c] >= 0.01]
    bg_cids = [c for c in all_cids if c not in set(fg_cids)]
    if not fg_cids or not bg_cids:
        # Fallback: top 25% by max frequency vs lower 75%
        sorted_c = sorted(all_cids, key=lambda c: cid_max_freq[c], reverse=True)
        cut = max(1, len(sorted_c) // 4)
        fg_cids = sorted_c[:cut]
        bg_cids = sorted_c[cut:]

    def _kmers_for_cids(cids: list[str], k: int = 4) -> Counter[str]:
        cnt: Counter[str] = Counter()
        for cid in cids:
            seq = clono_meta[cid][0]
            if len(seq) >= k:
                for m in {seq[idx : idx + k] for idx in range(len(seq) - k + 1)}:
                    cnt[m] += 1
        return cnt

    fg_kmers = _kmers_for_cids(fg_cids, 4)
    bg_kmers = _kmers_for_cids(bg_cids, 4)
    motif_rows: list[dict[str, Any]] = []
    if fg_cids and bg_cids:
        n_fg = len(fg_cids)
        n_bg = len(bg_cids)
        for motif, fg_hit in fg_kmers.items():
            if fg_hit < 2:
                continue
            bg_hit = bg_kmers.get(motif, 0)
            table = [[fg_hit, max(0, n_fg - fg_hit)], [bg_hit, max(0, n_bg - bg_hit)]]
            odds, pval = stats.fisher_exact(table, alternative="greater")
            motif_rows.append(
                {
                    "motif": motif,
                    "fg_clones": fg_hit,
                    "bg_clones": bg_hit,
                    "fg_pct": round(100.0 * fg_hit / n_fg, 2),
                    "bg_pct": round(100.0 * bg_hit / n_bg, 2),
                    "odds_ratio": round(float(odds), 3) if math.isfinite(float(odds)) else 999.0,
                    "p_value": float(pval),
                }
            )
        if motif_rows:
            mfdrs = bh_fdr([r["p_value"] for r in motif_rows])
            for r, q in zip(motif_rows, mfdrs):
                r["fdr"] = round(q, 6)
                r["p_value"] = round(r["p_value"], 6)
            motif_rows.sort(key=lambda r: (r["p_value"], -r["fg_clones"]))
    enriched_motifs = motif_rows[: max(1, min(20, args.top_n))]

    # Longitudinal Fisher's exact differential tracking
    base_sid = args.baseline_sample
    foll_sid = args.followup_sample
    if (base_sid is None) ^ (foll_sid is None):
        fail("Provide both --baseline_sample and --followup_sample, or omit both when exactly 2 samples are present.")
    if base_sid is None and foll_sid is None:
        if n_samples == 2:
            base_sid, foll_sid = sample_ids[0], sample_ids[1]
            warnings.append(f"Auto-selected longitudinal pair: baseline='{base_sid}', followup='{foll_sid}'.")
    else:
        if base_sid not in sample_counts:
            fail(f"baseline_sample '{base_sid}' not found in repertoire samples {sample_ids}.")
        if foll_sid not in sample_counts:
            fail(f"followup_sample '{foll_sid}' not found in repertoire samples {sample_ids}.")
        if base_sid == foll_sid:
            fail("baseline_sample and followup_sample must be distinct.")

    diff_df: pd.DataFrame | None = None
    n_expanded: int | None = None
    n_contracted: int | None = None
    top_diff: list[dict[str, Any]] = []

    if base_sid is not None and foll_sid is not None:
        n0 = float(sample_totals[base_sid])
        n1 = float(sample_totals[foll_sid])
        eps = 0.5 / max(n0, n1, 1.0)
        pair_cids = sorted(set(sample_counts[base_sid].keys()) | set(sample_counts[foll_sid].keys()))
        diff_rows: list[dict[str, Any]] = []
        for cid in pair_cids:
            c0 = float(sample_counts[base_sid].get(cid, 0.0))
            c1 = float(sample_counts[foll_sid].get(cid, 0.0))
            f0 = c0 / n0 if n0 > 0 else 0.0
            f1 = c1 / n1 if n1 > 0 else 0.0
            c0_int = int(round(c0))
            c1_int = int(round(c1))
            n0_int = max(c0_int, int(round(n0)))
            n1_int = max(c1_int, int(round(n1)))
            table = [[c1_int, n1_int - c1_int], [c0_int, n0_int - c0_int]]
            odds, pval = stats.fisher_exact(table, alternative="two-sided")
            l2fc = math.log2((f1 + eps) / (f0 + eps))
            cdr3_aa, v_gene = clono_meta[cid]
            diff_rows.append(
                {
                    "clonotype_id": cid,
                    "cdr3_aa": cdr3_aa,
                    "v_gene": v_gene,
                    "baseline_count": round(c0, 2),
                    "followup_count": round(c1, 2),
                    "baseline_freq": round(f0, 6),
                    "followup_freq": round(f1, 6),
                    "freq_diff": round(f1 - f0, 6),
                    "log2_fold_change": round(l2fc, 4),
                    "odds_ratio": round(float(odds), 4) if math.isfinite(float(odds)) else 999.0,
                    "p_value": float(pval),
                }
            )
        fdrs = bh_fdr([r["p_value"] for r in diff_rows])
        for r, q in zip(diff_rows, fdrs):
            r["fdr"] = round(q, 6)
            r["p_value"] = round(r["p_value"], 6)
            if q < 0.05 and r["followup_freq"] > r["baseline_freq"]:
                r["trend"] = "expanded"
            elif q < 0.05 and r["followup_freq"] < r["baseline_freq"]:
                r["trend"] = "contracted"
            else:
                r["trend"] = "stable"
        diff_rows.sort(key=lambda r: (r["fdr"], -abs(r["freq_diff"])))
        diff_df = pd.DataFrame(diff_rows)
        n_expanded = int((diff_df["trend"] == "expanded").sum())
        n_contracted = int((diff_df["trend"] == "contracted").sum())
        top_diff = diff_rows[: max(1, args.top_n)]

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    params = {
        "locus": target_locus,
        "baseline_sample": base_sid,
        "followup_sample": foll_sid,
        "top_n": int(args.top_n),
        "detected_format": detected_fmt,
    }
    sci_ctx = {"locus": target_locus, "n_samples": n_samples, "n_public_clonotypes": len(public_rows)}

    ov_path = write_tsv(
        overlap_df,
        out_dir / "repertoire_overlap_matrix.tsv",
        tool="imm_clonotype_tracking",
        title="Pairwise Morisita-Horn and Jaccard Repertoire Overlap Matrix",
        inputs=[args.repertoire],
        parameters=params,
        scientific_context=sci_ctx,
    )
    artifacts = [ov_path]

    diff_path: str | None = None
    if diff_df is not None:
        diff_path = write_tsv(
            diff_df,
            out_dir / "clonotype_differential_tracking.tsv",
            tool="imm_clonotype_tracking",
            title=f"Longitudinal Clonotype Expansion & Contraction ({base_sid} vs {foll_sid})",
            inputs=[args.repertoire],
            parameters=params,
            scientific_context=sci_ctx,
        )
        artifacts.append(diff_path)

    # 2-panel figure
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.9))

    im = ax1.imshow(mh_mat, cmap="YlGnBu", vmin=0.0, vmax=1.0)
    ax1.set_xticks(np.arange(n_samples))
    ax1.set_yticks(np.arange(n_samples))
    ax1.set_xticklabels(sample_ids, rotation=30, ha="right")
    ax1.set_yticklabels(sample_ids)
    for i in range(n_samples):
        for j in range(n_samples):
            val = mh_mat[i, j]
            ax1.text(j, i, f"{val:.2f}", ha="center", va="center", color="white" if val > 0.55 else "#0f172a", fontsize=8.5)
    ax1.set_title(f"Morisita-Horn Repertoire Overlap ({target_locus})")
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)

    if diff_df is not None and not diff_df.empty:
        eps_plot = 0.5 / max(sample_totals[base_sid], sample_totals[foll_sid], 1.0)
        x_val = np.log10(diff_df["baseline_freq"].to_numpy(dtype=float, copy=True) + eps_plot)
        y_val = np.log10(diff_df["followup_freq"].to_numpy(dtype=float, copy=True) + eps_plot)
        trends = diff_df["trend"].tolist()
        colors = {"expanded": "#dc2626", "contracted": "#2563eb", "stable": "#94a3b8"}
        for tr_label in ("stable", "contracted", "expanded"):
            mask = np.array([t == tr_label for t in trends])
            if np.any(mask):
                ax2.scatter(
                    x_val[mask],
                    y_val[mask],
                    s=36 if tr_label != "stable" else 18,
                    alpha=0.85 if tr_label != "stable" else 0.55,
                    color=colors[tr_label],
                    label=f"{tr_label} ({int(mask.sum())})",
                )
        lims = [min(float(x_val.min()), float(y_val.min())) - 0.2, max(float(x_val.max()), float(y_val.max())) + 0.2]
        ax2.plot(lims, lims, linestyle="--", color="#334155", linewidth=1.0)
        ax2.set_xlim(lims)
        ax2.set_ylim(lims)
        # Annotate top significant clones
        sig_sub = diff_df[diff_df["trend"] != "stable"].head(6)
        for _, r in sig_sub.iterrows():
            xp = math.log10(float(r["baseline_freq"]) + eps_plot)
            yp = math.log10(float(r["followup_freq"]) + eps_plot)
            ax2.annotate(str(r["cdr3_aa"])[:10], (xp, yp), textcoords="offset points", xytext=(4, 4), fontsize=7.5)
        ax2.set_xlabel(f"log10 Frequency ({base_sid})")
        ax2.set_ylabel(f"log10 Frequency ({foll_sid})")
        ax2.set_title("Longitudinal Clonotype Expansion / Contraction")
        ax2.legend(frameon=False, loc="upper left")
        ax2.grid(alpha=0.2, linestyle="--")
    else:
        # Multi-sample top shared clonotypes heatmap
        top_show = public_rows[: min(12, len(public_rows))]
        if top_show:
            mat = np.array([[sample_freqs[s].get(r["clonotype_id"], 0.0) * 100.0 for s in sample_ids] for r in top_show])
            im2 = ax2.imshow(mat, cmap="Reds", aspect="auto")
            ax2.set_xticks(np.arange(n_samples))
            ax2.set_xticklabels(sample_ids, rotation=30, ha="right")
            ax2.set_yticks(np.arange(len(top_show)))
            ax2.set_yticklabels([r["cdr3_aa"] for r in top_show], fontsize=8)
            ax2.set_title("Top Public Clonotypes (% Abundance)")
            fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
        else:
            ax2.text(0.5, 0.5, "No shared clonotypes across samples", ha="center", va="center", transform=ax2.transAxes)
            ax2.set_title("Public Clonotype Abundance")

    fig.tight_layout()
    fig_file = out_dir / "clonotype_tracking_overlap.png"
    fig.savefig(fig_file, bbox_inches="tight")
    plt.close(fig)
    fig_path = write_artifact(
        fig_file,
        tool="imm_clonotype_tracking",
        title="Clonotype Repertoire Overlap & Longitudinal Tracking",
        inputs=[args.repertoire],
        parameters=params,
        scientific_context=sci_ctx,
    )
    artifacts.append(fig_path)

    emit(
        {
            "status": "ok",
            "tool": "imm_clonotype_tracking",
            "locus": target_locus,
            "n_samples": n_samples,
            "baseline_sample": base_sid,
            "followup_sample": foll_sid,
            "n_public_clonotypes": len(public_rows),
            "n_expanded": n_expanded,
            "n_contracted": n_contracted,
            "mean_morisita_horn": round(float(np.mean(offdiag_mh)), 6) if offdiag_mh else 0.0,
            "mean_jaccard": round(float(np.mean(offdiag_jac)), 6) if offdiag_jac else 0.0,
            "top_public_clonotypes": top_public,
            "enriched_4mer_motifs": enriched_motifs,
            "top_differential_clones": top_diff,
            "outputs": {
                "repertoire_overlap_matrix": ov_path,
                "clonotype_differential_tracking": diff_path,
                "clonotype_tracking_overlap": fig_path,
            },
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
