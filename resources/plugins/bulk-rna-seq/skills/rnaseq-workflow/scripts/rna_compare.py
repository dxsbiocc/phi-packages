#!/usr/bin/env python3
"""rna_compare: integrate multiple DESeq2/edgeR/limma contrast tables into a
unified DEG comparison matrix, direction-concordance summary, pairwise overlap
table, and log2FC correlation plot.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import emit, fail, read_de_table, setup_matplotlib, write_artifact


def parse_labeled_inputs(specs):
    out = []
    for idx, spec in enumerate(specs, start=1):
        if "=" in spec and not Path(spec).is_file():
            label, path = spec.split("=", 1)
            label = label.strip() or f"contrast_{idx}"
        else:
            path = spec
            label = Path(path).stem.replace(".deseq2.results", "").replace(".results", "")
        p = Path(path)
        if not p.is_file():
            fail(f"DE table not found: {p}")
        out.append((label, p))
    if len(out) < 2:
        fail("rna_compare requires at least 2 contrast tables")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True, help="Two or more DE tables, optionally 'label=path'.")
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--log2fc", type=float, default=1.0)
    ap.add_argument("--output_dir", required=True)
    args = ap.parse_args()

    labeled = parse_labeled_inputs(args.inputs)
    tables = {}
    all_genes = set()
    per_contrast = {}

    for label, p in labeled:
        df = read_de_table(p).set_index("gene")
        tables[label] = df
        all_genes.update(df.index.tolist())
        up = set(df.index[(df["padj"] < args.padj) & (df["log2fc"] >= args.log2fc)])
        down = set(df.index[(df["padj"] < args.padj) & (df["log2fc"] <= -args.log2fc)])
        per_contrast[label] = {
            "n_tested": int(len(df)),
            "n_up": int(len(up)),
            "n_down": int(len(down)),
            "n_deg": int(len(up | down)),
            "up_set": up,
            "down_set": down,
        }

    ordered_genes = sorted(all_genes)
    labels = [lbl for lbl, _ in labeled]
    merged = pd.DataFrame(index=ordered_genes)
    merged.index.name = "gene"

    sig_count = np.zeros(len(ordered_genes), dtype=int)
    up_count = np.zeros(len(ordered_genes), dtype=int)
    down_count = np.zeros(len(ordered_genes), dtype=int)

    for lbl in labels:
        df = tables[lbl]
        fc = df["log2fc"].reindex(ordered_genes)
        padj = df["padj"].reindex(ordered_genes)
        merged[f"{lbl}_log2fc"] = fc.round(4)
        merged[f"{lbl}_padj"] = padj
        is_up = (padj < args.padj) & (fc >= args.log2fc)
        is_down = (padj < args.padj) & (fc <= -args.log2fc)
        status = np.where(is_up, "up", np.where(is_down, "down", "ns"))
        merged[f"{lbl}_status"] = status
        up_count += is_up.fillna(False).to_numpy(dtype=int)
        down_count += is_down.fillna(False).to_numpy(dtype=int)
        sig_count += (is_up | is_down).fillna(False).to_numpy(dtype=int)

    merged["n_contrasts_sig"] = sig_count
    merged["n_up"] = up_count
    merged["n_down"] = down_count
    merged["direction_pattern"] = np.where(
        (up_count > 0) & (down_count > 0),
        "discordant",
        np.where(up_count == len(labels), "shared_up",
                 np.where(down_count == len(labels), "shared_down",
                          np.where(sig_count > 1, "partially_shared",
                                   np.where(sig_count == 1, "contrast_specific", "none")))),
    )

    deg_only = merged.loc[merged["n_contrasts_sig"] > 0].sort_values(
        ["n_contrasts_sig", "n_up", "n_down"], ascending=[False, False, False]
    )

    # Pairwise overlap & log2FC Pearson r
    pair_rows = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            a, b = labels[i], labels[j]
            deg_a = per_contrast[a]["up_set"] | per_contrast[a]["down_set"]
            deg_b = per_contrast[b]["up_set"] | per_contrast[b]["down_set"]
            shared_up = len(per_contrast[a]["up_set"] & per_contrast[b]["up_set"])
            shared_down = len(per_contrast[a]["down_set"] & per_contrast[b]["down_set"])
            opp = len((per_contrast[a]["up_set"] & per_contrast[b]["down_set"]) | (per_contrast[a]["down_set"] & per_contrast[b]["up_set"]))
            union_n = len(deg_a | deg_b)
            jaccard = (len(deg_a & deg_b) / union_n) if union_n > 0 else 0.0
            pair_df = merged[[f"{a}_log2fc", f"{b}_log2fc"]].dropna()
            r_val = float(pair_df.iloc[:, 0].corr(pair_df.iloc[:, 1])) if len(pair_df) >= 3 else np.nan
            pair_rows.append(
                {
                    "contrast_a": a,
                    "contrast_b": b,
                    "shared_deg": int(len(deg_a & deg_b)),
                    "shared_up": int(shared_up),
                    "shared_down": int(shared_down),
                    "discordant_deg": int(opp),
                    "jaccard_deg": round(jaccard, 4),
                    "log2fc_pearson_r": round(r_val, 4) if not np.isnan(r_val) else None,
                }
            )

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    raw_inputs = [str(p) for _, p in labeled]
    params = {"contrasts": labels, "padj": args.padj, "log2fc": args.log2fc}

    deg_path = outdir / "multi-contrast-degs.tsv"
    deg_only.to_csv(deg_path, sep="\t")
    pair_df = pd.DataFrame(pair_rows)
    pair_path = outdir / "contrast-pairwise-overlap.tsv"
    pair_df.to_csv(pair_path, sep="\t", index=False)

    # Barplot of up/down DEGs per contrast + scatter for first two contrasts
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    x_pos = np.arange(len(labels))
    ups = [per_contrast[l]["n_up"] for l in labels]
    downs = [per_contrast[l]["n_down"] for l in labels]
    axes[0].bar(x_pos - 0.18, ups, width=0.35, color="#E11D48", label="Up")
    axes[0].bar(x_pos + 0.18, downs, width=0.35, color="#0EA5E9", label="Down")
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    axes[0].set_ylabel("Significant DEGs")
    axes[0].set_title(f"DEG counts (padj<{args.padj}, |log2FC|>={args.log2fc})")
    axes[0].legend(fontsize=8)

    a, b = labels[0], labels[1]
    sc = merged[[f"{a}_log2fc", f"{b}_log2fc", "direction_pattern"]].dropna()
    colors = sc["direction_pattern"].map(
        {
            "shared_up": "#E11D48",
            "shared_down": "#0EA5E9",
            "discordant": "#F59E0B",
            "partially_shared": "#8B5CF6",
            "contrast_specific": "#64748B",
            "none": "#CBD5E1",
        }
    ).fillna("#CBD5E1")
    axes[1].scatter(sc[f"{a}_log2fc"], sc[f"{b}_log2fc"], c=colors, s=18, alpha=0.7)
    axes[1].axhline(0, color="#475569", lw=0.8, ls="--")
    axes[1].axvline(0, color="#475569", lw=0.8, ls="--")
    axes[1].set_xlabel(f"{a} log2FC")
    axes[1].set_ylabel(f"{b} log2FC")
    axes[1].set_title(f"log2FC concordance ({a} vs {b})")
    fig.tight_layout()
    fig_path = outdir / "multi-contrast-comparison.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    artifacts = [
        write_artifact(deg_path, "rna_compare", "Multi-contrast integrated DEG table", inputs=raw_inputs, parameters=params, shape=deg_only.shape),
        write_artifact(pair_path, "rna_compare", "Pairwise contrast overlap and log2FC correlation", inputs=raw_inputs, parameters=params, shape=pair_df.shape),
        write_artifact(fig_path, "rna_compare", "Multi-contrast DEG counts and log2FC concordance", inputs=raw_inputs, parameters=params),
    ]

    summary_clean = {
        lbl: {k: v for k, v in info.items() if not k.endswith("_set")}
        for lbl, info in per_contrast.items()
    }
    emit(
        {
            "tool": "rna_compare",
            "contrasts": labels,
            "per_contrast": summary_clean,
            "n_union_degs": int(len(deg_only)),
            "n_shared_up_all": int((merged["direction_pattern"] == "shared_up").sum()),
            "n_shared_down_all": int((merged["direction_pattern"] == "shared_down").sum()),
            "n_discordant": int((merged["direction_pattern"] == "discordant").sum()),
            "pairwise": pair_rows,
            "outputs": {
                "integrated_degs": str(deg_path),
                "pairwise_overlap": str(pair_path),
                "comparison_figure": str(fig_path),
            },
            "scientificContext": {
                str(deg_path): {
                    "valueSemantics": "de_stats",
                }
            },
            "artifacts": artifacts,
        }
    )


if __name__ == "__main__":
    main()
