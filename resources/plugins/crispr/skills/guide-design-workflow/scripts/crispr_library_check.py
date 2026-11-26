#!/usr/bin/env python3
"""Pre-synthesis & pre-count audit of pooled CRISPR sgRNA libraries."""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    emit,
    fail,
    find_enzyme_collisions,
    parse_enzyme_list,
    parse_fasta,
    pick_col,
    read_table,
    setup_matplotlib,
    write_artifact,
)


def _load_library_df(path: str | Path, sgrna_col: str, seq_col: str, gene_col: str):
    import pandas as pd

    p = Path(path)
    if not p.is_file():
        fail(f"Library file not found: {path}")
    if p.suffix.lower() in (".fa", ".fasta", ".fna"):
        recs = parse_fasta(p)
        rows = []
        for header, seq in recs:
            parts = header.split("|")
            gid = parts[0]
            gene = parts[1] if len(parts) > 1 else gid.rsplit("_", 1)[0]
            rows.append({"sgRNA": gid, "seq": seq, "gene": gene})
        return pd.DataFrame(rows), "sgRNA", "seq", "gene"

    df = read_table(p).copy()
    cols = list(df.columns)
    sg_c = pick_col(cols, sgrna_col, ["sgrna", "sgRNA", "guide_id", "id", "oligo_id", "name"])
    sq_c = pick_col(cols, seq_col, ["seq", "sequence", "spacer", "sgrna_seq", "guide_seq"])
    gn_c = pick_col(cols, gene_col, ["gene", "gene_id", "target", "target_gene", "symbol"])
    return df, sg_c, sq_c, gn_c


def main() -> None:
    ap = argparse.ArgumentParser(description="Audit pooled CRISPR sgRNA library")
    ap.add_argument("--library", required=True, help="Library TSV/CSV/FASTA")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--sgrna_col", default="")
    ap.add_argument("--seq_col", default="")
    ap.add_argument("--gene_col", default="")
    ap.add_argument("--control_prefix", default="NonTargeting,SafeHarbor,NTC,CTRL")
    ap.add_argument("--avoid_enzymes", action="append", default=None)
    args = ap.parse_args()

    import numpy as np
    import pandas as pd

    df, sg_c, sq_c, gn_c = _load_library_df(args.library, args.sgrna_col, args.seq_col, args.gene_col)
    if df.empty:
        fail("Library table is empty.")

    enzymes = parse_enzyme_list(args.avoid_enzymes)
    ctrl_tokens = [t.strip().lower() for t in args.control_prefix.split(",") if t.strip()]

    df = df.copy()
    df["_id"] = df[sg_c].astype(str).str.strip()
    df["_seq"] = df[sq_c].astype(str).str.strip().str.upper()
    df["_gene"] = df[gn_c].astype(str).str.strip()

    id_counts = df["_id"].value_counts()
    seq_counts = df["_seq"].value_counts()
    mode_len = int(df["_seq"].str.len().mode().iloc[0])

    valid_re = re.compile(r"^[ACGT]+$")
    flagged_rows: list[dict] = []

    n_dup_id = 0
    n_dup_seq = 0
    n_invalid = 0
    n_tttt = 0
    n_enzyme = 0
    n_gc_out = 0
    n_len_out = 0
    gc_vals: list[float] = []

    for _, row in df.iterrows():
        gid = row["_id"]
        seq = row["_seq"]
        gene = row["_gene"]
        issues: list[str] = []

        if id_counts[gid] > 1:
            issues.append("duplicate_id")
            n_dup_id += 1
        if seq_counts[seq] > 1:
            issues.append("duplicate_sequence")
            n_dup_seq += 1
        if not valid_re.fullmatch(seq):
            issues.append("invalid_iupac_base")
            n_invalid += 1
            gc = float("nan")
        else:
            gc = (seq.count("G") + seq.count("C")) / max(1, len(seq))
            gc_vals.append(gc)
            if gc < 0.30 or gc > 0.75:
                issues.append("gc_outlier")
                n_gc_out += 1
        if len(seq) != mode_len:
            issues.append(f"length_mismatch({len(seq)}!={mode_len})")
            n_len_out += 1
        if "TTTT" in seq:
            issues.append("pol3_tttt")
            n_tttt += 1
        enz_hits = find_enzyme_collisions(seq, enzymes)
        if enz_hits:
            issues.append("enzyme_site:" + "+".join(enz_hits))
            n_enzyme += 1

        if issues:
            flagged_rows.append(
                {
                    "sgRNA": gid,
                    "gene": gene,
                    "sequence": seq,
                    "length": len(seq),
                    "gc_content": round(gc, 4) if np.isfinite(gc) else np.nan,
                    "issues": ";".join(issues),
                }
            )

    def _is_ctrl(g: str) -> bool:
        gl = g.lower()
        return any(gl.startswith(tok) or tok in gl for tok in ctrl_tokens)

    df["_is_ctrl"] = df["_gene"].map(_is_ctrl)
    n_controls = int(df["_is_ctrl"].sum())
    ctrl_frac = round(n_controls / max(1, len(df)), 4)

    gene_grp = (
        df.groupby("_gene", as_index=False)
        .agg(
            n_guides=("_id", "count"),
            n_unique_seqs=("_seq", "nunique"),
            is_control=("_is_ctrl", "first"),
        )
        .rename(columns={"_gene": "gene"})
    )
    gene_grp["low_coverage_lt3"] = (~gene_grp["is_control"]) & (gene_grp["n_guides"] < 3)
    n_genes_low = int(gene_grp["low_coverage_lt3"].sum())

    warnings: list[str] = []
    if n_dup_id > 0:
        warnings.append(f"Detected {n_dup_id} row(s) with duplicate sgRNA IDs.")
    if n_dup_seq > 0:
        warnings.append(f"Detected {n_dup_seq} row(s) with duplicate spacer sequences.")
    if n_invalid > 0:
        warnings.append(f"Detected {n_invalid} spacer(s) with non-ACGT bases.")
    if n_tttt > 0:
        warnings.append(f"Detected {n_tttt} spacer(s) containing Pol III TTTT termination motif.")
    if n_enzyme > 0:
        warnings.append(f"Detected {n_enzyme} spacer(s) colliding with cloning restriction sites ({','.join(enzymes)}).")
    if n_controls == 0:
        warnings.append("Zero non-targeting / safe-harbor control sgRNAs detected.")
    elif ctrl_frac < 0.01:
        warnings.append(f"Low control sgRNA representation ({ctrl_frac * 100:.2f}% < 1.0%).")
    if n_genes_low > 0:
        warnings.append(f"{n_genes_low} non-control gene(s) have <3 sgRNAs.")

    if n_dup_id > 0 or n_invalid > 0 or n_dup_seq > 0:
        qc_status = "fail" if (n_dup_id > 0 or n_invalid > 0) else "warn"
    elif warnings:
        qc_status = "warn"
    else:
        qc_status = "pass"

    flagged_df = pd.DataFrame(
        flagged_rows,
        columns=["sgRNA", "gene", "sequence", "length", "gc_content", "issues"],
    )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gene_path = out_dir / "library_qc_per_gene.tsv"
    flag_path = out_dir / "library_flagged_guides.tsv"
    fig_path = out_dir / "library_qc_summary.png"

    gene_grp.to_csv(gene_path, sep="\t", index=False)
    flagged_df.to_csv(flag_path, sep="\t", index=False)

    write_artifact(
        gene_path,
        "crispr_library_check",
        "Pooled Library Per-Gene sgRNA Coverage",
        inputs=[args.library],
        parameters={"control_prefix": args.control_prefix, "avoid_enzymes": enzymes},
        shape=(int(len(gene_grp)), int(gene_grp.shape[1])),
    )
    write_artifact(
        flag_path,
        "crispr_library_check",
        "Flagged sgRNAs in Library Audit",
        inputs=[args.library],
        parameters={"control_prefix": args.control_prefix, "avoid_enzymes": enzymes},
        shape=(int(len(flagged_df)), int(flagged_df.shape[1])),
    )

    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.4), dpi=180)
    non_ctrl_counts = gene_grp[~gene_grp["is_control"]]["n_guides"]
    if not non_ctrl_counts.empty:
        bins = range(1, int(non_ctrl_counts.max()) + 2)
        ax1.hist(non_ctrl_counts, bins=bins, align="left", color="#0d9488", edgecolor="#0f766e", rwidth=0.8)
    ax1.axvline(2.5, color="#dc2626", ls="--", lw=1.0, label="Min 3 guides/gene")
    ax1.set_xlabel("sgRNAs per Target Gene")
    ax1.set_ylabel("Number of Genes")
    ax1.set_title("Library Coverage per Gene")
    ax1.legend(frameon=False, fontsize=8)

    issue_labels = ["Dup ID", "Dup Seq", "Invalid", "TTTT", "Enzyme", "GC Outlier", "Len Skew"]
    issue_vals = [n_dup_id, n_dup_seq, n_invalid, n_tttt, n_enzyme, n_gc_out, n_len_out]
    cols_bar = ["#dc2626" if v > 0 else "#cbd5e1" for v in issue_vals]
    ax2.bar(issue_labels, issue_vals, color=cols_bar)
    ax2.set_ylabel("Flagged Guides")
    ax2.set_title(f"Library QC Issues (Status: {qc_status.upper()})")
    ax2.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "crispr_library_check",
        "Pooled sgRNA Library QC Summary",
        inputs=[args.library],
        parameters={"qc_status": qc_status},
    )

    emit(
        {
            "status": "ok",
            "tool": "crispr_library_check",
            "n_guides": int(len(df)),
            "n_genes": int(gene_grp["gene"].nunique()),
            "n_control_guides": n_controls,
            "control_fraction": ctrl_frac,
            "n_duplicate_ids": n_dup_id,
            "n_duplicate_sequences": n_dup_seq,
            "n_invalid_sequences": n_invalid,
            "n_pol3_tttt": n_tttt,
            "n_enzyme_collisions": n_enzyme,
            "n_gc_outliers": n_gc_out,
            "n_genes_low_coverage": n_genes_low,
            "qc_status": qc_status,
            "outputs": {
                "library_qc_per_gene": str(gene_path),
                "library_flagged_guides": str(flag_path),
                "qc_plot": str(fig_path),
            },
            "artifacts": [str(gene_path), str(flag_path), str(fig_path)],
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
