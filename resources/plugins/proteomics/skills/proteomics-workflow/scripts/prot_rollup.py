#!/usr/bin/env python3
"""prot_rollup: aggregate peptide/precursor intensities to protein-level log2 abundances."""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    clean_sample_name,
    contaminant_mask,
    emit,
    fail,
    infer_scale,
    read_table,
)


PEPTIDE_COL_CANDIDATES = [
    "peptide_id",
    "precursor.id",
    "stripped.sequence",
    "peptide sequence",
    "sequence",
    "peptide",
    "modified.sequence",
]
PROTEIN_COL_CANDIDATES = [
    "protein_id",
    "leading razor protein",
    "protein.group",
    "genes",
    "gene",
    "proteins",
    "protein",
    "majority protein ids",
]


def detect_peptide_columns(
    df: pd.DataFrame, peptide_col: str | None, protein_col: str | None
) -> tuple[str, str, list[str]]:
    cols = list(df.columns)
    lower_map = {c.lower(): c for c in cols}

    p_col = peptide_col
    if p_col:
        if p_col not in df.columns:
            fail(f"peptide_col '{p_col}' not found in columns")
    else:
        for cand in PEPTIDE_COL_CANDIDATES:
            if cand in lower_map:
                p_col = lower_map[cand]
                break
        if p_col is None:
            p_col = cols[0]

    prot_col = protein_col
    if prot_col:
        if prot_col not in df.columns:
            fail(f"protein_col '{prot_col}' not found in columns")
    else:
        for cand in PROTEIN_COL_CANDIDATES:
            if cand in lower_map and lower_map[cand] != p_col:
                prot_col = lower_map[cand]
                break
        if prot_col is None:
            others = [c for c in cols if c != p_col]
            if not others:
                fail("cannot identify protein_col in peptide table")
            prot_col = others[0]

    non_sample = {
        p_col.lower(),
        prot_col.lower(),
        "protein.group",
        "protein.ids",
        "protein.names",
        "genes",
        "gene",
        "proteins",
        "protein",
        "leading razor protein",
        "first.protein.description",
        "precursor.id",
        "stripped.sequence",
        "modified.sequence",
        "peptide sequence",
        "sequence",
        "peptide",
        "charge",
        "length",
        "missed cleavages",
        "mass",
        "potential contaminant",
        "reverse",
        "only identified by site",
        "q-value",
        "pep",
        "score",
    }
    lfq_cols = [c for c in cols if c.lower().startswith("lfq intensity ")]
    maxlfq_cols = [c for c in cols if c.lower().endswith(" maxlfq intensity")]
    if lfq_cols:
        sample_cols = lfq_cols
    elif maxlfq_cols:
        sample_cols = maxlfq_cols
    else:
        sample_cols = []
        for c in cols:
            if c.lower() in non_sample:
                continue
            if pd.api.types.is_numeric_dtype(df[c]):
                sample_cols.append(c)
            else:
                coerced = pd.to_numeric(df[c], errors="coerce")
                if coerced.notna().mean() >= 0.5 and coerced.notna().sum() > 0:
                    sample_cols.append(c)
    if len(sample_cols) < 2:
        fail("fewer than 2 quantitative sample columns detected in peptide table")
    return p_col, prot_col, sample_cols


def rrollup_protein(log2_block: np.ndarray) -> np.ndarray:
    """Reference-scaled log2 rollup (Rrollup-style): aligns peptides to best-covered reference."""
    if log2_block.shape[0] == 1:
        return log2_block[0].copy()
    valid_counts = np.sum(np.isfinite(log2_block), axis=1)
    med_ints = np.nanmedian(log2_block, axis=1)
    # Pick reference peptide with highest coverage, breaking ties by highest median intensity
    ref_idx = int(np.lexsort((np.nan_to_num(med_ints, nan=-1e9), valid_counts))[-1])
    ref = log2_block[ref_idx]
    aligned = np.empty_like(log2_block)
    aligned[ref_idx] = ref
    for i in range(log2_block.shape[0]):
        if i == ref_idx:
            continue
        row = log2_block[i]
        both = np.isfinite(row) & np.isfinite(ref)
        if np.any(both):
            shift = float(np.median(ref[both] - row[both]))
        else:
            r_med = float(np.nanmedian(ref)) if np.any(np.isfinite(ref)) else 0.0
            p_med = float(np.nanmedian(row)) if np.any(np.isfinite(row)) else 0.0
            shift = r_med - p_med
        aligned[i] = row + shift
    return np.nanmedian(aligned, axis=0)


def top3_protein(log2_block: np.ndarray) -> np.ndarray:
    """Top-3 (Hi3) mean log2 intensity across the 3 most intense peptides per sample."""
    n_samples = log2_block.shape[1]
    out = np.full(n_samples, np.nan, dtype=float)
    for j in range(n_samples):
        vals = log2_block[:, j]
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        top = np.sort(vals)[-3:]
        out[j] = float(np.mean(top))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--peptide_col", default=None)
    parser.add_argument("--protein_col", default=None)
    parser.add_argument("--shared_peptides", choices=["unique_only", "razor"], default="unique_only")
    parser.add_argument("--method", choices=["rrollup", "top3", "median", "sum"], default="rrollup")
    parser.add_argument("--min_peptides", type=int, default=2)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    df = read_table(args.input)
    p_col, prot_col, sample_cols = detect_peptide_columns(df, args.peptide_col, args.protein_col)

    n_peptides_raw = len(df)
    cmask = contaminant_mask(df, prot_col)
    df = df.loc[~cmask].copy()
    n_contam_removed = n_peptides_raw - len(df)

    # Resolve shared peptides (semicolon-delimited protein lists)
    prot_lists = [
        [p.strip() for p in str(v).split(";") if p.strip() and p.strip().lower() != "nan"]
        for v in df[prot_col]
    ]
    if args.shared_peptides == "unique_only":
        assigned = [plist[0] if len(plist) == 1 else None for plist in prot_lists]
        n_shared_dropped = sum(1 for plist in prot_lists if len(plist) > 1)
    else:
        # Occam's razor: count total peptides per candidate protein, assign shared peptide to most-supported protein
        prot_counts: dict[str, int] = {}
        for plist in prot_lists:
            for p in plist:
                prot_counts[p] = prot_counts.get(p, 0) + 1
        assigned = [
            max(plist, key=lambda p: (prot_counts.get(p, 0), p)) if plist else None
            for plist in prot_lists
        ]
        n_shared_dropped = 0

    keep_mask = np.array([a is not None for a in assigned], dtype=bool)
    df = df.loc[keep_mask].copy()
    assigned_prots = [a for a in assigned if a is not None]
    df["_assigned_protein"] = assigned_prots

    mat = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    mat.columns = [clean_sample_name(c) for c in sample_cols]
    scale_info = infer_scale(mat)
    if scale_info["scale"] == "raw_intensity":
        mat_lin = mat.replace(0, np.nan).where(lambda x: x > 0)
        mat_log2 = np.log2(mat_lin)
    else:
        mat_log2 = mat.copy()
        mat_lin = np.power(2.0, mat_log2)

    mat_log2["_protein"] = assigned_prots
    mat_lin["_protein"] = assigned_prots
    sample_names = list(mat.columns)

    protein_rows = []
    count_rows = []
    n_below_min_pep = 0
    for prot, grp_log2 in mat_log2.groupby("_protein", sort=True):
        n_pep = int(len(grp_log2))
        if n_pep < args.min_peptides:
            n_below_min_pep += 1
            continue
        block = grp_log2[sample_names].to_numpy(dtype=float)
        if args.method == "rrollup":
            vec = rrollup_protein(block)
        elif args.method == "top3":
            vec = top3_protein(block)
        elif args.method == "median":
            vec = np.nanmedian(block, axis=0)
        else:
            lin_block = mat_lin.loc[grp_log2.index, sample_names].to_numpy(dtype=float)
            s = np.nansum(lin_block, axis=0)
            all_na = np.all(~np.isfinite(lin_block), axis=0)
            s[all_na | (s <= 0)] = np.nan
            vec = np.log2(s)

        row_dict = {"protein_id": str(prot)}
        for idx, sname in enumerate(sample_names):
            row_dict[sname] = round(float(vec[idx]), 5) if np.isfinite(vec[idx]) else np.nan
        protein_rows.append(row_dict)
        count_rows.append(
            {
                "protein_id": str(prot),
                "n_peptides": n_pep,
                "n_samples_detected": int(np.sum(np.isfinite(vec))),
                "median_log2_intensity": round(float(np.nanmedian(vec)), 4) if np.any(np.isfinite(vec)) else np.nan,
            }
        )

    if not protein_rows:
        fail(f"no proteins passed min_peptides={args.min_peptides}")

    prot_df = pd.DataFrame(protein_rows)
    counts_df = pd.DataFrame(count_rows).sort_values("n_peptides", ascending=False)

    prot_path = os.path.join(args.output_dir, "protein_matrix.tsv")
    counts_path = os.path.join(args.output_dir, "protein_peptide_counts.tsv")
    prot_df.to_csv(prot_path, sep="\t", index=False)
    counts_df.to_csv(counts_path, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    pep_vals = counts_df["n_peptides"].to_numpy()
    bins = min(25, max(5, int(pep_vals.max() - pep_vals.min() + 1)))
    ax.hist(pep_vals, bins=bins, color="#38BDF8", edgecolor="#0F172A", alpha=0.85)
    ax.set_xlabel("Peptides per protein")
    ax.set_ylabel("Number of proteins")
    ax.set_title(f"Peptide-to-protein rollup ({args.method}, {len(prot_df)} proteins)")
    fig.tight_layout()
    fig_path = os.path.join(args.output_dir, "rollup-peptide-counts.png")
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)

    emit(
        {
            "input": os.path.abspath(args.input),
            "protein_matrix": os.path.abspath(prot_path),
            "peptide_counts_table": os.path.abspath(counts_path),
            "peptide_col": p_col,
            "protein_col": prot_col,
            "input_scale": scale_info["scale"],
            "shared_peptides": args.shared_peptides,
            "method": args.method,
            "min_peptides": args.min_peptides,
            "n_peptides_raw": n_peptides_raw,
            "n_contaminants_removed": n_contam_removed,
            "n_shared_peptides_dropped": n_shared_dropped,
            "n_proteins_below_min_peptides": n_below_min_pep,
            "n_proteins_quantified": int(len(prot_df)),
            "median_peptides_per_protein": round(float(counts_df["n_peptides"].median()), 2),
            "figures": [os.path.abspath(fig_path)],
        }
    )


if __name__ == "__main__":
    main()
