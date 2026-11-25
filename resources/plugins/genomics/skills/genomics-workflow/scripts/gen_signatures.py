#!/usr/bin/env python3
"""gen_signatures: SBS96 trinucleotide catalog extraction and NNLS COSMIC signature refitting."""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import nnls

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    SIX_CLASSES,
    TRANSITIONS,
    all_sbs96_channels,
    canonical_sbs96,
    emit,
    fail,
    pyrimidine_sub,
    read_maf,
    read_table,
    write_artifact,
)


def cosine_sim(a: np.ndarray, b: np.ndarray) -> float:
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na == 0.0 or nb == 0.0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", required=True)
    parser.add_argument("--context_col", default=None)
    parser.add_argument("--ref_signatures", default=None)
    parser.add_argument("--min_snvs", type=int, default=20)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_maf(args.maf)
    all_samples = sorted(df["Tumor_Sample_Barcode"].unique())
    snvs = df[df["Variant_Type"] == "SNP"].copy()

    lower_cols = {c.lower(): c for c in df.columns}
    ctx_col = args.context_col or next(
        (lower_cols[k] for k in ("context", "triplet", "trinucleotide_context", "ref_context") if k in lower_cols),
        None,
    )

    figures: list[str] = []
    warnings: list[str] = []

    if ctx_col is None or ctx_col not in df.columns:
        # Honest fallback: 6-class substitution spectrum only
        snvs["sub6"] = [pyrimidine_sub(r, a) for r, a in zip(snvs["Reference_Allele"], snvs["Tumor_Seq_Allele2"])]
        snvs = snvs.dropna(subset=["sub6"])
        counts = {c: int((snvs["sub6"] == c).sum()) for c in SIX_CLASSES}
        n_ti = sum(counts[c] for c in TRANSITIONS)
        n_tv = sum(counts[c] for c in SIX_CLASSES if c not in TRANSITIONS)
        warnings.append(
            "no trinucleotide Context column found in MAF; SBS96 catalog and COSMIC NNLS refitting skipped (6-class Ti/Tv fallback reported)"
        )
        emit(
            {
                "maf": os.path.abspath(args.maf),
                "mode": "fallback_6class",
                "sbs96": {"skipped": "no trinucleotide context column"},
                "refitting": {"skipped": "no trinucleotide context column"},
                "substitution_counts": counts,
                "titv_ratio": round(n_ti / n_tv, 4) if n_tv > 0 else None,
                "figures": [],
                "artifacts": [],
                "warnings": warnings,
            }
        )
        return

    # Build 96-channel SBS catalog (channels x samples)
    channels = all_sbs96_channels()
    snvs["sbs96"] = [
        canonical_sbs96(r, a, ctx)
        for r, a, ctx in zip(snvs["Reference_Allele"], snvs["Tumor_Seq_Allele2"], snvs[ctx_col])
    ]
    snvs_valid = snvs.dropna(subset=["sbs96"])
    ct = pd.crosstab(snvs_valid["sbs96"], snvs_valid["Tumor_Sample_Barcode"])
    catalog = ct.reindex(index=channels, columns=all_samples, fill_value=0)
    catalog_path = os.path.join(args.output_dir, "sbs96_catalog.tsv")
    catalog.reset_index().rename(columns={"index": "Type", "sbs96": "Type"}).to_csv(catalog_path, sep="\t", index=False)

    # Plot cohort aggregate 96-channel spectrum
    cohort_vec = catalog.sum(axis=1).to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(11.5, 3.6))
    palette = ["#38BDF8", "#1E293B", "#EF4444", "#94A3B8", "#22C55E", "#F472B6"]
    bar_colors = [palette[i // 16] for i in range(96)]
    ax.bar(np.arange(96), cohort_vec, color=bar_colors, width=0.8)
    ax.set_xticks([8 + 16 * i for i in range(6)])
    ax.set_xticklabels(SIX_CLASSES, fontsize=9)
    ax.set_ylabel("SNV count")
    ax.set_title(f"Cohort SBS96 Mutational Spectrum (n={int(cohort_vec.sum())} SNVs)")
    fig.tight_layout()
    spec_fig = os.path.join(args.output_dir, "sbs96-spectrum.png")
    fig.savefig(spec_fig, dpi=150)
    plt.close(fig)
    figures.append(spec_fig)

    refit_summary: dict = {"skipped": "no ref_signatures matrix provided (pass COSMIC v3.4 SBS96 TSV)"}
    if args.ref_signatures:
        ref_df = read_table(args.ref_signatures)
        type_col = "Type" if "Type" in ref_df.columns else ref_df.columns[0]
        ref_df = ref_df.set_index(type_col)
        shared = [c for c in channels if c in ref_df.index]
        if len(shared) < 96:
            fail(f"ref_signatures matched only {len(shared)}/96 canonical SBS96 channels")
        W = ref_df.loc[channels].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        sig_names = list(W.columns)
        W_mat = W.to_numpy(dtype=float)

        exp_rows = []
        cos_rows = []
        for s in all_samples:
            y = catalog[s].to_numpy(dtype=float)
            total_snv = int(y.sum())
            if total_snv < args.min_snvs:
                cos_rows.append({"Tumor_Sample_Barcode": s, "n_snvs": total_snv, "cosine_similarity": np.nan, "status": "below_min_snvs"})
                continue
            weights, _ = nnls(W_mat, y)
            recon = W_mat @ weights
            cos = cosine_sim(y, recon)
            rel = weights / max(float(weights.sum()), 1e-12)
            row = {"Tumor_Sample_Barcode": s}
            for idx, sname in enumerate(sig_names):
                row[sname] = round(float(rel[idx]), 5)
            exp_rows.append(row)
            cos_rows.append({"Tumor_Sample_Barcode": s, "n_snvs": total_snv, "cosine_similarity": round(cos, 4), "status": "fitted"})

        exp_df = pd.DataFrame(exp_rows)
        cos_df = pd.DataFrame(cos_rows)
        exp_path = os.path.join(args.output_dir, "signature_exposures.tsv")
        cos_path = os.path.join(args.output_dir, "signature_reconstruction.tsv")
        exp_df.to_csv(exp_path, sep="\t", index=False)
        cos_df.to_csv(cos_path, sep="\t", index=False)

        if len(exp_df) > 0:
            mean_exp = exp_df[sig_names].mean(axis=0).sort_values(ascending=False)
            top_sigs = list(mean_exp.head(8).index)
            fig, ax = plt.subplots(figsize=(max(6.0, len(exp_df) * 0.45), 4.2))
            exp_plot = exp_df.set_index("Tumor_Sample_Barcode")[top_sigs]
            exp_plot.plot(kind="bar", stacked=True, ax=ax, colormap="tab10")
            ax.set_ylabel("Relative signature exposure")
            ax.set_title("Per-sample SBS96 NNLS Signature Exposures")
            ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=8)
            fig.tight_layout()
            exp_fig = os.path.join(args.output_dir, "signature-exposures.png")
            fig.savefig(exp_fig, dpi=150)
            plt.close(fig)
            figures.append(exp_fig)
            median_cos = float(cos_df["cosine_similarity"].dropna().median()) if cos_df["cosine_similarity"].notna().any() else None
            refit_summary = {
                "exposures_table": os.path.abspath(exp_path),
                "reconstruction_table": os.path.abspath(cos_path),
                "n_reference_signatures": len(sig_names),
                "n_samples_fitted": int(len(exp_df)),
                "median_cosine_similarity": round(median_cos, 4) if median_cos is not None else None,
                "top_mean_signatures": {str(k): round(float(v), 4) for k, v in mean_exp.head(8).items()},
            }
        else:
            refit_summary = {"skipped": f"no sample had >= {args.min_snvs} SNVs for NNLS refitting"}

    params_meta = {"min_snvs": args.min_snvs}
    in_files = (args.maf, args.ref_signatures)
    fig_titles = {
        "sbs96-spectrum.png": "Cohort 96-channel SBS mutational spectrum",
        "signature-exposures.png": "Per-sample SBS96 NNLS signature exposures",
    }
    artifacts = [
        write_artifact(f, "gen_signatures", fig_titles.get(os.path.basename(f), os.path.basename(f)), inputs=in_files, parameters=params_meta)
        for f in figures
    ]
    artifacts.append(
        write_artifact(catalog_path, "gen_signatures", "SBS96 mutational catalog (96 channels x samples)", inputs=in_files, parameters=params_meta, shape=(96, len(all_samples) + 1))
    )
    if args.ref_signatures and os.path.exists(os.path.join(args.output_dir, "signature_exposures.tsv")):
        artifacts.append(
            write_artifact(exp_path, "gen_signatures", "Per-sample SBS96 signature exposures", inputs=in_files, parameters=params_meta, shape=exp_df.shape)
        )
        artifacts.append(
            write_artifact(cos_path, "gen_signatures", "SBS96 signature reconstruction cosine similarity", inputs=in_files, parameters=params_meta, shape=cos_df.shape)
        )

    emit(
        {
            "maf": os.path.abspath(args.maf),
            "mode": "sbs96",
            "context_col": ctx_col,
            "sbs96_catalog": os.path.abspath(catalog_path),
            "n_snvs_classified": int(cohort_vec.sum()),
            "n_samples": len(all_samples),
            "refitting": refit_summary,
            "figures": [os.path.abspath(f) for f in figures],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
