#!/usr/bin/env python3
"""imm_repertoire_inspect: Harmonize AIRR, 10x Cell Ranger V(D)J, TRUST4, and
MiXCR repertoire tables into a standardized AIRR-compatible TSV, compute
per-sample QC metrics (productive %, out-of-frame/stop %, chain breakdown,
CDR3 AA length & Kyte-Doolittle hydropathy), and render spectratype QC plots.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    harmonize_repertoire_inputs,
    normalize_gene_base,
    parse_str_list,
    setup_matplotlib,
    write_artifact,
    write_tsv,
)


def _apply_chain_filter(df: pd.DataFrame, chain_filter: str) -> pd.DataFrame:
    if chain_filter == "all":
        return df.copy()
    if chain_filter in {"TCR", "BCR"}:
        return df[df["receptor_type"] == chain_filter].copy()
    return df[df["locus"] == chain_filter].copy()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--format", choices=["auto", "airr", "cellranger_vdj", "trust4", "mixcr"], default="auto")
    ap.add_argument("--chain_filter", choices=["all", "TCR", "BCR", "TRB", "TRA", "IGH"], default="all")
    ap.add_argument("--productive_only", action="store_true", default=False)
    args = ap.parse_args()

    input_paths = parse_str_list(args.inputs)
    if not input_paths:
        fail("At least one input repertoire file is required.")

    raw_df, detected_formats = harmonize_repertoire_inputs(input_paths, fmt=args.format)
    chain_df = _apply_chain_filter(raw_df, args.chain_filter)
    if chain_df.empty:
        fail(f"No records remain after applying chain_filter='{args.chain_filter}'.")

    warnings: list[str] = []
    qc_rows: list[dict[str, Any]] = []

    for sid, sdf in chain_df.groupby("sample_id", sort=False):
        n_rec = int(len(sdf))
        tot_reads = float(sdf["duplicate_count"].sum())
        prod_mask = sdf["productive"].to_numpy(dtype=bool, copy=True)
        oof_mask = sdf["out_of_frame"].to_numpy(dtype=bool, copy=True)
        stop_mask = sdf["stop_codon"].to_numpy(dtype=bool, copy=True)

        n_prod = int(prod_mask.sum())
        prod_pct = 100.0 * n_prod / max(n_rec, 1)
        oof_pct = 100.0 * int(oof_mask.sum()) / max(n_rec, 1)
        stop_pct = 100.0 * int(stop_mask.sum()) / max(n_rec, 1)

        prod_sub = sdf[sdf["productive"]].copy()
        if not prod_sub.empty:
            v_base = prod_sub["v_call"].map(normalize_gene_base)
            clono_keys = prod_sub["cdr3_aa"].astype(str) + "|" + v_base.astype(str)
            unique_clono = int(clono_keys.nunique())
            lens = prod_sub.loc[prod_sub["cdr3_length"] > 0, "cdr3_length"].to_numpy(dtype=float, copy=True)
            hydros = prod_sub["kyte_doolittle_gravy"].dropna().to_numpy(dtype=float, copy=True)
            mean_len = float(np.mean(lens)) if lens.size else None
            sd_len = float(np.std(lens, ddof=1)) if lens.size > 1 else 0.0
            mean_hydro = float(np.mean(hydros)) if hydros.size else None
        else:
            unique_clono = 0
            mean_len = None
            sd_len = None
            mean_hydro = None

        loci_counts = sdf["locus"].value_counts().to_dict()
        tcr_n = int((sdf["receptor_type"] == "TCR").sum())
        bcr_n = int((sdf["receptor_type"] == "BCR").sum())

        qc_flags: list[str] = []
        if prod_pct < 70.0:
            qc_flags.append("low_productive_pct")
            warnings.append(f"Sample '{sid}' has low productive fraction ({prod_pct:.1f}% < 70%).")
        if unique_clono < 10:
            qc_flags.append("low_clonotype_depth")
            warnings.append(f"Sample '{sid}' has fewer than 10 productive clonotypes ({unique_clono}).")
        qc_status = "PASS" if not qc_flags else ",".join(qc_flags)

        qc_rows.append(
            {
                "sample_id": sid,
                "source_format": ",".join(sorted(sdf["source_format"].astype(str).unique())),
                "total_records": n_rec,
                "total_reads_umis": round(tot_reads, 2),
                "unique_clonotypes": unique_clono,
                "productive_records": n_prod,
                "productive_pct": round(prod_pct, 2),
                "out_of_frame_pct": round(oof_pct, 2),
                "stop_codon_pct": round(stop_pct, 2),
                "tcr_pct": round(100.0 * tcr_n / max(n_rec, 1), 2),
                "bcr_pct": round(100.0 * bcr_n / max(n_rec, 1), 2),
                "tra_count": int(loci_counts.get("TRA", 0)),
                "trb_count": int(loci_counts.get("TRB", 0)),
                "igh_count": int(loci_counts.get("IGH", 0)),
                "igk_count": int(loci_counts.get("IGK", 0)),
                "igl_count": int(loci_counts.get("IGL", 0)),
                "mean_cdr3_length": round(mean_len, 2) if mean_len is not None else np.nan,
                "sd_cdr3_length": round(sd_len, 2) if sd_len is not None else np.nan,
                "mean_kd_hydropathy": round(mean_hydro, 3) if mean_hydro is not None else np.nan,
                "qc_status": qc_status,
            }
        )

    out_df = chain_df[chain_df["productive"]].copy() if args.productive_only else chain_df.copy()
    if out_df.empty:
        fail("No productive records remain after --productive_only filtering.")

    # Recompute normalized frequencies per sample on the exported table
    for sid, idx in out_df.groupby("sample_id").groups.items():
        tot = float(out_df.loc[idx, "duplicate_count"].sum())
        out_df.loc[idx, "frequency"] = out_df.loc[idx, "duplicate_count"] / tot if tot > 0 else 0.0

    qc_df = pd.DataFrame(qc_rows)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    params = {
        "format": args.format,
        "chain_filter": args.chain_filter,
        "productive_only": bool(args.productive_only),
    }
    sci_ctx = {"n_samples": int(len(qc_df)), "chain_filter": args.chain_filter}

    harm_path = write_tsv(
        out_df,
        out_dir / "harmonized_repertoire.tsv",
        tool="imm_repertoire_inspect",
        title="Harmonized AIRR-Compatible Immune Repertoire Table",
        inputs=input_paths,
        parameters=params,
        scientific_context=sci_ctx,
    )
    qc_path = write_tsv(
        qc_df,
        out_dir / "repertoire_qc_summary.tsv",
        tool="imm_repertoire_inspect",
        title="Per-Sample Immune Repertoire QC Summary",
        inputs=input_paths,
        parameters=params,
        scientific_context=sci_ctx,
    )

    # 2-panel Spectratype & QC figure
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))

    # Panel A: Productive CDR3 AA length spectratype per sample
    prod_for_plot = chain_df[(chain_df["productive"]) & (chain_df["cdr3_length"] > 0)].copy()
    if not prod_for_plot.empty:
        palette = plt.cm.tab10(np.linspace(0, 1, max(2, len(qc_df))))
        for idx, sid in enumerate(qc_df["sample_id"].tolist()):
            sub_s = prod_for_plot[prod_for_plot["sample_id"] == sid]
            if sub_s.empty:
                continue
            grp = sub_s.groupby("cdr3_length")["duplicate_count"].sum()
            tot = float(grp.sum())
            if tot > 0:
                freqs = 100.0 * grp / tot
                ax1.plot(
                    freqs.index.to_numpy(dtype=int, copy=True),
                    freqs.to_numpy(dtype=float, copy=True),
                    marker="o",
                    markersize=4,
                    linewidth=1.8,
                    color=palette[idx % len(palette)],
                    label=f"{sid} (n={int(sub_s.shape[0])})",
                )
        ax1.set_xlabel("CDR3 Amino Acid Length (aa)")
        ax1.set_ylabel("Repertoire Abundance (%)")
        ax1.set_title("Productive CDR3 Length Spectratype")
        ax1.legend(frameon=False, loc="best")
        ax1.grid(alpha=0.25, linestyle="--")
    else:
        ax1.text(0.5, 0.5, "No productive CDR3 sequences", ha="center", va="center", transform=ax1.transAxes)
        ax1.set_title("Productive CDR3 Length Spectratype")

    # Panel B: Productive / Out-of-frame / Stop codon breakdown per sample
    sids = qc_df["sample_id"].tolist()
    y_pos = np.arange(len(sids))
    prod_vals = qc_df["productive_pct"].to_numpy(dtype=float, copy=True)
    oof_vals = qc_df["out_of_frame_pct"].to_numpy(dtype=float, copy=True)
    stop_vals = qc_df["stop_codon_pct"].to_numpy(dtype=float, copy=True)
    other_vals = np.clip(100.0 - prod_vals - oof_vals - stop_vals, 0.0, 100.0)

    ax2.barh(y_pos, prod_vals, color="#0f766e", label="Productive")
    ax2.barh(y_pos, oof_vals, left=prod_vals, color="#f59e0b", label="Out-of-frame")
    ax2.barh(y_pos, stop_vals, left=prod_vals + oof_vals, color="#dc2626", label="Stop codon")
    ax2.barh(y_pos, other_vals, left=prod_vals + oof_vals + stop_vals, color="#94a3b8", label="Other non-prod")
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(sids)
    ax2.set_xlim(0, 102)
    ax2.axvline(70.0, color="#334155", linestyle="--", linewidth=1.0, alpha=0.7)
    ax2.set_xlabel("Percentage of Records (%)")
    ax2.set_title("Repertoire Productivity & Frame QC")
    ax2.legend(frameon=False, loc="lower right")

    fig.tight_layout()
    fig_file = out_dir / "repertoire_spectratype_qc.png"
    fig.savefig(fig_file, bbox_inches="tight")
    plt.close(fig)
    fig_path = write_artifact(
        fig_file,
        tool="imm_repertoire_inspect",
        title="CDR3 Spectratype & Repertoire Productivity QC",
        inputs=input_paths,
        parameters=params,
        scientific_context=sci_ctx,
    )

    emit(
        {
            "status": "ok",
            "tool": "imm_repertoire_inspect",
            "detected_formats": detected_formats,
            "chain_filter": args.chain_filter,
            "productive_only": bool(args.productive_only),
            "n_samples": int(len(qc_df)),
            "total_reads_umis": float(qc_df["total_reads_umis"].sum()),
            "total_unique_clonotypes": int(qc_df["unique_clonotypes"].sum()),
            "mean_productive_pct": round(float(qc_df["productive_pct"].mean()), 2),
            "samples": qc_df.to_dict(orient="records"),
            "outputs": {
                "harmonized_repertoire": harm_path,
                "repertoire_qc_summary": qc_path,
                "repertoire_spectratype_qc": fig_path,
            },
            "artifacts": [harm_path, qc_path, fig_path],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
