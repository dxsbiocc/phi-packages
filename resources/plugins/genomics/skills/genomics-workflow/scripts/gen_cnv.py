#!/usr/bin/env python3
"""gen_cnv: CBS/CNVkit segment FGA (Fraction Genome Altered) and recurrent gain/loss summary."""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, load_clinical, read_seg


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seg", required=True)
    parser.add_argument("--gain_threshold", type=float, default=0.2)
    parser.add_argument("--loss_threshold", type=float, default=-0.2)
    parser.add_argument("--amp_threshold", type=float, default=0.8)
    parser.add_argument("--del_threshold", type=float, default=-0.8)
    parser.add_argument("--clinical", default=None)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--group_col", default=None)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    seg = read_seg(args.seg)
    seg["Length"] = np.maximum(0.0, seg["End"] - seg["Start"])
    seg = seg[seg["Length"] > 0].copy()
    if seg.empty:
        fail("segment file contains no positive-length segments")

    samples = sorted(seg["Sample"].unique())
    n_samples = len(samples)

    # Per-sample FGA & gain/loss/amp/del burden
    samp_rows = []
    for s, grp in seg.groupby("Sample"):
        tot_bp = float(grp["Length"].sum())
        gain_bp = float(grp.loc[grp["Segment_Mean"] >= args.gain_threshold, "Length"].sum())
        loss_bp = float(grp.loc[grp["Segment_Mean"] <= args.loss_threshold, "Length"].sum())
        amp_bp = float(grp.loc[grp["Segment_Mean"] >= args.amp_threshold, "Length"].sum())
        del_bp = float(grp.loc[grp["Segment_Mean"] <= args.del_threshold, "Length"].sum())
        fga = (gain_bp + loss_bp) / tot_bp if tot_bp > 0 else 0.0
        samp_rows.append(
            {
                "Sample": str(s),
                "n_segments": int(len(grp)),
                "total_mb": round(tot_bp / 1e6, 2),
                "fga": round(fga, 4),
                "gain_fraction": round(gain_bp / tot_bp, 4) if tot_bp else 0.0,
                "loss_fraction": round(loss_bp / tot_bp, 4) if tot_bp else 0.0,
                "amp_fraction": round(amp_bp / tot_bp, 4) if tot_bp else 0.0,
                "del_fraction": round(del_bp / tot_bp, 4) if tot_bp else 0.0,
            }
        )
    fga_df = pd.DataFrame(samp_rows).sort_values("fga", ascending=False)

    # Chromosome-level recurrent gain/loss summary
    chr_rows = []
    for chrom, grp in seg.groupby("Chromosome"):
        n_gain_s = grp.loc[grp["Segment_Mean"] >= args.gain_threshold, "Sample"].nunique()
        n_loss_s = grp.loc[grp["Segment_Mean"] <= args.loss_threshold, "Sample"].nunique()
        n_amp_s = grp.loc[grp["Segment_Mean"] >= args.amp_threshold, "Sample"].nunique()
        n_del_s = grp.loc[grp["Segment_Mean"] <= args.del_threshold, "Sample"].nunique()
        chr_rows.append(
            {
                "Chromosome": str(chrom),
                "gain_samples": int(n_gain_s),
                "loss_samples": int(n_loss_s),
                "amp_samples": int(n_amp_s),
                "del_samples": int(n_del_s),
                "gain_frequency": round(n_gain_s / n_samples, 4),
                "loss_frequency": round(n_loss_s / n_samples, 4),
            }
        )
    chr_df = pd.DataFrame(chr_rows)
    # Natural chromosome ordering
    def chr_key(c: str) -> tuple[int, str]:
        return (0, f"{int(c):02d}") if c.isdigit() else (1, c)
    chr_df = chr_df.sort_values("Chromosome", key=lambda col: col.map(chr_key))

    group_test = None
    if args.clinical and args.group_col:
        clin, scol = load_clinical(args.clinical, args.sample_col)
        if args.group_col not in clin.columns:
            fail(f"group_col '{args.group_col}' not found in clinical table")
        gmap = dict(zip(clin[scol], clin[args.group_col].astype(str)))
        fga_df[args.group_col] = fga_df["Sample"].map(gmap)
        valid_g = fga_df.dropna(subset=[args.group_col])
        groups = sorted(valid_g[args.group_col].unique())
        if len(groups) >= 2:
            arrays = [valid_g.loc[valid_g[args.group_col] == g, "fga"].to_numpy(dtype=float) for g in groups]
            nonempty = [a for a in arrays if len(a) >= 2]
            if len(nonempty) == 2:
                u_stat, p_val = stats.mannwhitneyu(nonempty[0], nonempty[1], alternative="two-sided")
                group_test = {"group_col": args.group_col, "test": "mannwhitneyu", "stat": round(float(u_stat), 4), "pvalue": float(p_val)}
            elif len(nonempty) > 2:
                h_stat, p_val = stats.kruskal(*nonempty)
                group_test = {"group_col": args.group_col, "test": "kruskal", "stat": round(float(h_stat), 4), "pvalue": float(p_val)}

    fga_path = os.path.join(args.output_dir, "sample_fga.tsv")
    chr_path = os.path.join(args.output_dir, "chromosome_cnv_summary.tsv")
    fga_df.to_csv(fga_path, sep="\t", index=False)
    chr_df.to_csv(chr_path, sep="\t", index=False)

    # Plot FGA + chromosome gain/loss frequencies
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    axes[0].bar(fga_df["Sample"], fga_df["fga"] * 100.0, color="#38BDF8")
    axes[0].set_ylabel("Fraction Genome Altered (%)")
    axes[0].set_title("Per-sample FGA")
    axes[0].tick_params(axis="x", rotation=45, labelsize=7)

    x = np.arange(len(chr_df))
    axes[1].bar(x, chr_df["gain_frequency"] * 100.0, color="#EF4444", alpha=0.8, label="Gain")
    axes[1].bar(x, -chr_df["loss_frequency"] * 100.0, color="#3B82F6", alpha=0.8, label="Loss")
    axes[1].axhline(0.0, color="#334155", linewidth=0.9)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(chr_df["Chromosome"], rotation=45, fontsize=8)
    axes[1].set_ylabel("Altered samples (% Gain / -% Loss)")
    axes[1].set_title("Recurrent Chromosome CNV Frequencies")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "cnv-fga-landscape.png")
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)

    emit(
        {
            "seg": os.path.abspath(args.seg),
            "sample_fga_table": os.path.abspath(fga_path),
            "chromosome_summary_table": os.path.abspath(chr_path),
            "n_samples": n_samples,
            "n_segments": int(len(seg)),
            "median_fga": round(float(fga_df["fga"].median()), 4),
            "top_altered_chromosomes": chr_df.sort_values(
                by=["gain_frequency", "loss_frequency"], ascending=False
            ).head(8).to_dict(orient="records"),
            "group_comparison": group_test,
            "figures": [os.path.abspath(fig_path)],
        }
    )


if __name__ == "__main__":
    main()
