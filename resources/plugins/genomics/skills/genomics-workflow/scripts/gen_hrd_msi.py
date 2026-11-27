#!/usr/bin/env python3
"""gen_hrd_msi: scarHRD/Marquard genomic scar HRD (LOH/TAI/LST) and MSI stratification."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    GENOME_CHR_INFO,
    emit,
    fail,
    read_allele_segments,
    read_maf,
    read_table,
    write_artifact,
)

AUTOSOMES = {str(i) for i in range(1, 23)}
MB = 1_000_000
LOH_MIN_BP = 15 * MB
TAI_MIN_BP = 11 * MB
TAI_TELOMERE_TOL_BP = 2 * MB
LST_SMOOTH_MIN_BP = 3 * MB
LST_SEG_MIN_BP = 10 * MB
LST_GAP_MAX_BP = 3 * MB


def _cn_state(row: pd.Series) -> tuple[int, int]:
    return (int(round(float(row["major_cn"]))), int(round(float(row["minor_cn"]))))


def _spans_whole_chr(start: float, end: float, chr_len: int, cent_start: int, cent_end: int) -> bool:
    length = end - start
    if length >= 0.90 * chr_len:
        return True
    return start <= TAI_TELOMERE_TOL_BP and end >= chr_len - TAI_TELOMERE_TOL_BP and start < cent_start and end > cent_end


def count_hrd_loh(chr_df: pd.DataFrame, chr_len: int, cent_start: int, cent_end: int) -> int:
    """Count LOH regions (minor_cn==0, total_cn>=1) strictly > 15 Mb not spanning entire chromosome."""
    segs = chr_df.sort_values(["start", "end"]).reset_index(drop=True)
    is_loh = (segs["minor_cn"].round() == 0) & (segs["total_cn"].round() >= 1)
    merged: list[tuple[float, float]] = []
    for idx, row in segs.iterrows():
        if not bool(is_loh.iloc[idx]):
            continue
        s, e = float(row["start"]), float(row["end"])
        if merged and s - merged[-1][1] <= 1000.0:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    count = 0
    for s, e in merged:
        if (e - s) > LOH_MIN_BP and not _spans_whole_chr(s, e, chr_len, cent_start, cent_end):
            count += 1
    return count


def count_tai(chr_df: pd.DataFrame, chr_len: int, cent_start: int, cent_end: int) -> int:
    """Count allelic-imbalance segments extending to telomere (>=11 Mb, not crossing centromere or whole-arm)."""
    segs = chr_df.sort_values(["start", "end"]).reset_index(drop=True)
    # Merge adjacent segments with identical allelic state before evaluating telomeric imbalance
    merged: list[tuple[float, float, int, int]] = []
    for _, row in segs.iterrows():
        s, e = float(row["start"]), float(row["end"])
        maj, mino = _cn_state(row)
        if merged and merged[-1][2] == maj and merged[-1][3] == mino and s - merged[-1][1] <= 1000.0:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e), maj, mino)
        else:
            merged.append((s, e, maj, mino))

    count = 0
    for s, e, maj, mino in merged:
        if maj == mino:
            continue
        if (e - s) < TAI_MIN_BP:
            continue
        # Must not cross centromere
        if s < cent_start and e > cent_end:
            continue
        near_p_telo = s <= TAI_TELOMERE_TOL_BP and e <= cent_end
        near_q_telo = e >= (chr_len - TAI_TELOMERE_TOL_BP) and s >= cent_start
        if near_p_telo:
            whole_arm = e >= (cent_start - TAI_TELOMERE_TOL_BP) or (e - s) >= 0.90 * cent_start
            if not whole_arm:
                count += 1
        elif near_q_telo:
            q_len = max(1, chr_len - cent_end)
            whole_arm = s <= (cent_end + TAI_TELOMERE_TOL_BP) or (e - s) >= 0.90 * q_len
            if not whole_arm:
                count += 1
    return count


def _count_lst_on_arm(arm_segs: list[tuple[float, float, int, int]]) -> int:
    """Smooth segments < 3 Mb on a single arm and count adjacent breaks (>=10 Mb, gap <= 3 Mb)."""
    if not arm_segs:
        return 0

    def _merge_identical(items: list[tuple[float, float, int, int]], max_gap: float) -> list[tuple[float, float, int, int]]:
        out: list[tuple[float, float, int, int]] = []
        for s, e, maj, mino in sorted(items, key=lambda x: (x[0], x[1])):
            if out and out[-1][2] == maj and out[-1][3] == mino and (s - out[-1][1]) <= max_gap:
                out[-1] = (out[-1][0], max(out[-1][1], e), maj, mino)
            else:
                out.append((s, e, maj, mino))
        return out

    cur = _merge_identical(arm_segs, max_gap=1000.0)
    # Remove short segments < 3 Mb iteratively (shortest first) and re-merge flanking identical states
    while True:
        short_idx = [i for i, (s, e, _, _) in enumerate(cur) if (e - s) < LST_SMOOTH_MIN_BP]
        if not short_idx:
            break
        drop_i = min(short_idx, key=lambda i: cur[i][1] - cur[i][0])
        cur.pop(drop_i)
        cur = _merge_identical(cur, max_gap=LST_GAP_MAX_BP)

    breaks = 0
    for i in range(len(cur) - 1):
        s1, e1, maj1, mino1 = cur[i]
        s2, e2, maj2, mino2 = cur[i + 1]
        if (maj1, mino1) == (maj2, mino2):
            continue
        if (e1 - s1) >= LST_SEG_MIN_BP and (e2 - s2) >= LST_SEG_MIN_BP and (s2 - e1) <= LST_GAP_MAX_BP:
            breaks += 1
    return breaks


def count_lst(chr_df: pd.DataFrame, cent_start: int, cent_end: int) -> int:
    """Count Large-Scale State Transitions on p and q arms of a chromosome."""
    p_segs: list[tuple[float, float, int, int]] = []
    q_segs: list[tuple[float, float, int, int]] = []
    for _, row in chr_df.iterrows():
        s, e = float(row["start"]), float(row["end"])
        maj, mino = _cn_state(row)
        if s < cent_start:
            pe = min(e, float(cent_start))
            if pe > s:
                p_segs.append((s, pe, maj, mino))
        if e > cent_end:
            qs = max(s, float(cent_end))
            if e > qs:
                q_segs.append((qs, e, maj, mino))
    return _count_lst_on_arm(p_segs) + _count_lst_on_arm(q_segs)


def compute_sample_hrd(seg_df: pd.DataFrame, genome: str) -> dict[str, dict[str, int]]:
    chr_info = GENOME_CHR_INFO[genome]
    results: dict[str, dict[str, int]] = {}
    for sample, sgrp in seg_df.groupby("sample"):
        loh = 0
        tai = 0
        lst = 0
        for chrom, cgrp in sgrp.groupby("chrom"):
            if chrom not in AUTOSOMES or chrom not in chr_info:
                continue
            chr_len, cent_start, cent_end = chr_info[chrom]
            loh += count_hrd_loh(cgrp, chr_len, cent_start, cent_end)
            tai += count_tai(cgrp, chr_len, cent_start, cent_end)
            lst += count_lst(cgrp, cent_start, cent_end)
        results[str(sample)] = {
            "loh_score": loh,
            "tai_score": tai,
            "lst_score": lst,
            "hrd_sum": loh + tai + lst,
        }
    return results


def parse_msi_table(path: str, fallback_sample: str | None = None) -> dict[str, dict[str, Any]]:
    """Parse MSIsensor-pro or Mantis summary TSV/CSV into per-sample MSI metrics."""
    df = read_table(path)
    if df.empty:
        fail(f"MSI table is empty: {path}")
    lower = {c.lower(): c for c in df.columns}
    sample_c = next(
        (lower[k] for k in ("sample", "tumor_sample_barcode", "sample_id", "id") if k in lower),
        None,
    )
    pct_c = next(
        (
            lower[k]
            for k in (
                "%",
                "pct",
                "msi_pct",
                "percent_unstable",
                "unstable_pct",
                "msisensor_pct",
                "msisensor_score",
                "msi_score",
            )
            if k in lower
        ),
        None,
    )
    tot_c = next((lower[k] for k in ("total_number_of_sites", "total_sites", "total") if k in lower), None)
    som_c = next(
        (lower[k] for k in ("number_of_somatic_sites", "somatic_sites", "unstable_sites", "somatic") if k in lower),
        None,
    )
    if pct_c is None and not (tot_c and som_c):
        fail("MSI table must contain '%' (or msi_pct/msi_score) or Total_Number_of_Sites + Number_of_Somatic_Sites")

    out: dict[str, dict[str, Any]] = {}
    for idx, row in df.iterrows():
        if sample_c:
            sname = str(row[sample_c]).strip()
        elif len(df) == 1 and fallback_sample:
            sname = fallback_sample
        else:
            sname = Path(path).stem if len(df) == 1 else f"{Path(path).stem}_{idx + 1}"

        total_sites = int(pd.to_numeric(row[tot_c], errors="coerce")) if tot_c and pd.notna(row[tot_c]) else None
        somatic_sites = int(pd.to_numeric(row[som_c], errors="coerce")) if som_c and pd.notna(row[som_c]) else None
        if pct_c:
            raw_pct = str(row[pct_c]).strip().rstrip("%")
            msi_pct = float(pd.to_numeric(raw_pct, errors="coerce"))
        else:
            msi_pct = (100.0 * somatic_sites / total_sites) if (total_sites and somatic_sites is not None) else np.nan
        if not np.isfinite(msi_pct):
            continue
        out[sname] = {
            "msi_pct": round(float(msi_pct), 4),
            "msi_total_sites": total_sites,
            "msi_somatic_sites": somatic_sites,
        }
    return out


def compute_maf_msi_proxy(path: str) -> dict[str, dict[str, Any]]:
    """Compute per-sample frameshift indel rate and indel-to-SNV ratio from a somatic MAF."""
    df = read_maf(path)
    out: dict[str, dict[str, Any]] = {}
    for sample, grp in df.groupby("Tumor_Sample_Barcode"):
        n_total = int(len(grp))
        is_snv = grp["Variant_Type"] == "SNP"
        is_indel = grp["Variant_Type"].isin({"DEL", "INS"}) | grp["Variant_Classification"].isin(
            {"Frame_Shift_Del", "Frame_Shift_Ins", "Frame_Shift_Indel", "In_Frame_Del", "In_Frame_Ins"}
        )
        is_fs = grp["Variant_Classification"].isin({"Frame_Shift_Del", "Frame_Shift_Ins", "Frame_Shift_Indel"})
        n_snv = int(is_snv.sum())
        n_indel = int(is_indel.sum())
        n_fs = int(is_fs.sum())
        fs_rate = n_fs / n_total if n_total > 0 else 0.0
        indel_snv_ratio = (n_indel / n_snv) if n_snv > 0 else float(n_indel)
        out[str(sample)] = {
            "n_snv": n_snv,
            "n_indel": n_indel,
            "n_frameshift_indel": n_fs,
            "frameshift_indel_rate": round(float(fs_rate), 4),
            "indel_to_snv_ratio": round(float(indel_snv_ratio), 4),
        }
    return out


def classify_biomarker(
    hrd_sum: int | None,
    hrd_cutoff: int,
    msi_pct: float | None,
    msi_pct_cutoff: float,
    maf_metrics: dict[str, Any] | None,
) -> tuple[bool | None, str, str, str, str]:
    """Return (hrd_positive, msi_status, msi_source, biomarker_stratum, therapy_rationale)."""
    hrd_pos = None if hrd_sum is None else bool(hrd_sum >= hrd_cutoff)

    if msi_pct is not None:
        msi_status = "MSI-H" if msi_pct >= msi_pct_cutoff else "MSS"
        msi_source = "msi_table"
    elif maf_metrics is not None:
        is_proxy_msi = (
            maf_metrics["n_indel"] >= 5
            and (maf_metrics["indel_to_snv_ratio"] >= 0.25 or maf_metrics["frameshift_indel_rate"] >= 0.20)
        )
        msi_status = "MSI-H_proxy" if is_proxy_msi else "MSS_proxy"
        msi_source = "maf_indel_proxy"
    else:
        msi_status = "MSI_unknown"
        msi_source = "none"

    if msi_status == "MSI-H":
        if hrd_pos is True:
            return hrd_pos, msi_status, msi_source, "MSI_H_HRD_pos", "ICI_and_PARPi_evaluation"
        return hrd_pos, msi_status, msi_source, "MSI_H", "immunotherapy_candidate"
    if msi_status == "MSI-H_proxy":
        if hrd_pos is True:
            return hrd_pos, msi_status, msi_source, "MSI_H_proxy_HRD_pos", "confirm_MSI_and_PARPi"
        return hrd_pos, msi_status, msi_source, "MSI_H_proxy", "confirm_MSI_immunotherapy"

    if hrd_pos is True:
        stratum = "HRD_pos_MSS" if msi_status in {"MSS", "MSS_proxy"} else "HRD_pos"
        return hrd_pos, msi_status, msi_source, stratum, "PARP_inhibitor_candidate"
    if hrd_pos is False:
        stratum = "HR_proficient_MSS" if msi_status in {"MSS", "MSS_proxy"} else "HR_proficient"
        return hrd_pos, msi_status, msi_source, stratum, "non_HRD_non_MSI"
    if msi_status in {"MSS", "MSS_proxy"}:
        return hrd_pos, msi_status, msi_source, msi_status, "non_MSI_H"
    return hrd_pos, msi_status, msi_source, "undetermined", "insufficient_data"


def render_landscape(
    summary_df: pd.DataFrame,
    out_path: str,
    hrd_cutoff: int,
    msi_pct_cutoff: float,
    has_segments: bool,
    has_msi_table: bool,
    has_maf: bool,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.5))
    plot_df = summary_df.head(40).copy()
    x = np.arange(len(plot_df))

    # Left panel: HRD stacked components (or indel proxy when segments absent)
    if has_segments:
        loh = plot_df["loh_score"].fillna(0).to_numpy(dtype=float)
        tai = plot_df["tai_score"].fillna(0).to_numpy(dtype=float)
        lst = plot_df["lst_score"].fillna(0).to_numpy(dtype=float)
        axes[0].bar(x, loh, color="#0284C7", label="HRD-LOH")
        axes[0].bar(x, tai, bottom=loh, color="#F59E0B", label="TAI")
        axes[0].bar(x, lst, bottom=loh + tai, color="#EF4444", label="LST")
        axes[0].axhline(hrd_cutoff, color="#0F172A", linestyle="--", linewidth=1.2, label=f"Cutoff ({hrd_cutoff})")
        axes[0].set_ylabel("Genomic scar score")
        axes[0].set_title("Per-sample scarHRD Components")
        axes[0].legend(frameon=False, fontsize=8)
    elif has_maf:
        fs = plot_df["frameshift_indel_rate"].fillna(0).to_numpy(dtype=float) * 100.0
        axes[0].bar(x, fs, color="#F97316")
        axes[0].set_ylabel("Frameshift indel rate (%)")
        axes[0].set_title("MAF Frameshift Indel Rate (no CNV segments)")
    else:
        axes[0].text(0.5, 0.5, "No CNV segments provided", ha="center", va="center")
        axes[0].set_axis_off()
    if has_segments or has_maf:
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(plot_df["sample"], rotation=45, ha="right", fontsize=7)

    # Right panel: HRD vs MSI scatter (or MSI barplot if single modality)
    stratum_colors = {
        "HRD_pos_MSS": "#EF4444",
        "HRD_pos": "#EF4444",
        "MSI_H": "#22C55E",
        "MSI_H_HRD_pos": "#A855F7",
        "MSI_H_proxy": "#10B981",
        "MSI_H_proxy_HRD_pos": "#A855F7",
        "HR_proficient_MSS": "#64748B",
        "HR_proficient": "#64748B",
        "MSS": "#64748B",
        "MSS_proxy": "#64748B",
    }
    if has_segments and (has_msi_table or has_maf):
        y_col = "msi_pct" if has_msi_table else "indel_to_snv_ratio"
        y_label = "MSIsensor unstable sites (%)" if has_msi_table else "MAF Indel / SNV ratio"
        for stratum, grp in summary_df.groupby("biomarker_stratum"):
            axes[1].scatter(
                grp["hrd_sum"],
                grp[y_col],
                label=str(stratum),
                color=stratum_colors.get(str(stratum), "#3B82F6"),
                s=52,
                alpha=0.88,
                edgecolors="#0F172A",
                linewidths=0.6,
            )
        axes[1].axvline(hrd_cutoff, color="#0F172A", linestyle="--", linewidth=1.0)
        if has_msi_table:
            axes[1].axhline(msi_pct_cutoff, color="#16A34A", linestyle="--", linewidth=1.0)
        axes[1].set_xlabel("HRD Sum (LOH + TAI + LST)")
        axes[1].set_ylabel(y_label)
        axes[1].set_title("HRD vs MSI Biomarker Stratification")
        axes[1].legend(frameon=False, fontsize=7.5)
    elif has_msi_table:
        vals = plot_df["msi_pct"].fillna(0).to_numpy(dtype=float)
        cols = ["#22C55E" if v >= msi_pct_cutoff else "#64748B" for v in vals]
        axes[1].bar(x, vals, color=cols)
        axes[1].axhline(msi_pct_cutoff, color="#16A34A", linestyle="--", linewidth=1.1, label=f"Cutoff ({msi_pct_cutoff}%)")
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(plot_df["sample"], rotation=45, ha="right", fontsize=7)
        axes[1].set_ylabel("Unstable microsatellites (%)")
        axes[1].set_title("Per-sample MSI Score")
        axes[1].legend(frameon=False, fontsize=8)
    else:
        counts = summary_df["biomarker_stratum"].value_counts()
        axes[1].bar(counts.index.astype(str), counts.values, color="#38BDF8")
        axes[1].set_ylabel("Samples")
        axes[1].set_title("Biomarker Stratum Counts")
        axes[1].tick_params(axis="x", rotation=30, labelsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segments", default=None)
    parser.add_argument("--msi_table", default=None)
    parser.add_argument("--maf", default=None)
    parser.add_argument("--genome", default="hg38", choices=["hg38", "hg19"])
    parser.add_argument("--hrd_cutoff", type=int, default=42)
    parser.add_argument("--msi_pct_cutoff", type=float, default=3.5)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    if not any([args.segments, args.msi_table, args.maf]):
        fail("at least one of --segments, --msi_table, or --maf must be provided")

    os.makedirs(args.output_dir, exist_ok=True)
    warnings: list[str] = []

    hrd_map: dict[str, dict[str, int]] = {}
    if args.segments:
        seg_df = read_allele_segments(args.segments, require_allele_cn=True)
        if seg_df.empty:
            fail("allele-specific segment table contains no valid positive-length segments")
        non_auto = sorted(set(seg_df["chrom"]) - AUTOSOMES)
        if non_auto:
            warnings.append(f"excluded non-autosomal chromosomes from scarHRD counting: {', '.join(non_auto)}")
        hrd_map = compute_sample_hrd(seg_df, args.genome)

    maf_map: dict[str, dict[str, Any]] = {}
    if args.maf:
        maf_map = compute_maf_msi_proxy(args.maf)

    fallback_sample = None
    known_samples = sorted(set(hrd_map) | set(maf_map))
    if len(known_samples) == 1:
        fallback_sample = known_samples[0]

    msi_map: dict[str, dict[str, Any]] = {}
    if args.msi_table:
        msi_map = parse_msi_table(args.msi_table, fallback_sample=fallback_sample)

    all_samples = sorted(set(hrd_map) | set(msi_map) | set(maf_map))
    if not all_samples:
        fail("no samples could be parsed from the provided inputs")

    if not args.segments:
        warnings.append("no allele-specific CNV --segments provided; HRD-LOH/TAI/LST scores omitted")
    if not args.msi_table and args.maf:
        warnings.append("no --msi_table provided; MSI classification uses MAF indel/frameshift surrogate proxy")

    rows: list[dict[str, Any]] = []
    for s in all_samples:
        h = hrd_map.get(s)
        m = msi_map.get(s)
        f = maf_map.get(s)
        hrd_sum = h["hrd_sum"] if h else None
        msi_pct = m["msi_pct"] if m else None
        hrd_pos, msi_status, msi_source, stratum, rationale = classify_biomarker(
            hrd_sum, args.hrd_cutoff, msi_pct, args.msi_pct_cutoff, f
        )
        rows.append(
            {
                "sample": s,
                "loh_score": h["loh_score"] if h else np.nan,
                "tai_score": h["tai_score"] if h else np.nan,
                "lst_score": h["lst_score"] if h else np.nan,
                "hrd_sum": hrd_sum if hrd_sum is not None else np.nan,
                "hrd_positive": hrd_pos,
                "msi_pct": msi_pct if msi_pct is not None else np.nan,
                "msi_status": msi_status,
                "msi_source": msi_source,
                "n_snv": f["n_snv"] if f else np.nan,
                "n_indel": f["n_indel"] if f else np.nan,
                "n_frameshift_indel": f["n_frameshift_indel"] if f else np.nan,
                "frameshift_indel_rate": f["frameshift_indel_rate"] if f else np.nan,
                "indel_to_snv_ratio": f["indel_to_snv_ratio"] if f else np.nan,
                "biomarker_stratum": stratum,
                "therapy_rationale": rationale,
            }
        )

    summary_df = pd.DataFrame(rows)
    table_path = os.path.join(args.output_dir, "hrd_msi_summary.tsv")
    summary_df.to_csv(table_path, sep="\t", index=False)

    fig_path = os.path.join(args.output_dir, "hrd_msi_landscape.png")
    render_landscape(
        summary_df,
        fig_path,
        args.hrd_cutoff,
        args.msi_pct_cutoff,
        bool(args.segments),
        bool(args.msi_table),
        bool(args.maf),
    )

    params_meta = {
        "genome": args.genome,
        "hrd_cutoff": args.hrd_cutoff,
        "msi_pct_cutoff": args.msi_pct_cutoff,
    }
    in_files = (args.segments, args.msi_table, args.maf)
    artifacts = [
        write_artifact(fig_path, "gen_hrd_msi", "Cohort scarHRD and MSI biomarker stratification", inputs=in_files, parameters=params_meta),
        write_artifact(table_path, "gen_hrd_msi", "Per-sample HRD-LOH, TAI, LST, and MSI summary", inputs=in_files, parameters=params_meta, shape=summary_df.shape),
    ]

    hrd_vals = summary_df["hrd_sum"].dropna().to_numpy(dtype=float)
    n_hrd_pos = int((summary_df["hrd_positive"] == True).sum()) if args.segments else None
    n_msi_high = (
        int(summary_df["msi_status"].isin({"MSI-H", "MSI-H_proxy"}).sum())
        if (args.msi_table or args.maf)
        else None
    )
    stratum_counts = {str(k): int(v) for k, v in summary_df["biomarker_stratum"].value_counts().items()}

    # JSON-clean per-sample records (convert NaN to None)
    clean_samples = []
    for rec in summary_df.to_dict(orient="records"):
        clean_samples.append(
            {k: (None if (isinstance(v, float) and not np.isfinite(v)) else v) for k, v in rec.items()}
        )

    emit(
        {
            "genome": args.genome,
            "hrd_cutoff": args.hrd_cutoff,
            "msi_pct_cutoff": args.msi_pct_cutoff,
            "modalities": {
                "segments": bool(args.segments),
                "msi_table": bool(args.msi_table),
                "maf": bool(args.maf),
            },
            "summary_table": os.path.abspath(table_path),
            "n_samples": len(all_samples),
            "n_hrd_positive": n_hrd_pos,
            "hrd_sum_median": round(float(np.median(hrd_vals)), 4) if len(hrd_vals) else None,
            "n_msi_high": n_msi_high,
            "stratum_counts": stratum_counts,
            "samples": clean_samples,
            "figures": [os.path.abspath(fig_path)],
            "artifacts": artifacts,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
