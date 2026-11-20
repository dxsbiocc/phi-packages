"""Grade MAG/bin quality from CheckM2 or CheckM tables (MIMAG tiers)."""
import argparse
from pathlib import Path

import pandas as pd

from _common import bin_stem, emit, fail, read_tsv, setup_matplotlib


def load_quality(path):
    df = read_tsv(path)
    cols = set(df.columns)
    if {"Name", "Completeness", "Contamination"} <= cols:
        source = "checkm2"
        name_col = "Name"
    elif {"Bin Id", "Completeness", "Contamination"} <= cols:
        source = "checkm"
        name_col = "Bin Id"
    else:
        fail("unrecognized quality table: expected CheckM2 quality_report.tsv "
             "(Name/Completeness/Contamination) or CheckM qa --tab_table "
             "(Bin Id/Completeness/Contamination)")
    out = pd.DataFrame({
        "bin": df[name_col].astype(str).map(bin_stem),
        "completeness": pd.to_numeric(df["Completeness"], errors="coerce"),
        "contamination": pd.to_numeric(df["Contamination"], errors="coerce"),
    })
    if out[["completeness", "contamination"]].isna().any().any():
        fail("non-numeric completeness/contamination values")
    for src, dst in (("Genome_Size", "genome_size"), ("Contig_N50", "contig_n50"),
                     ("GC_Content", "gc_percent"), ("Total_Contigs", "n_contigs")):
        if src in cols:
            out[dst] = pd.to_numeric(df[src], errors="coerce")
    return source, out


def tier(comp, cont, hq_c, hq_x, mq_c, mq_x):
    if comp >= hq_c and cont <= hq_x:
        return "HQ"
    if comp >= mq_c and cont <= mq_x:
        return "MQ"
    return "LQ"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--hq_completeness", type=float, default=90.0)
    ap.add_argument("--hq_contamination", type=float, default=5.0)
    ap.add_argument("--mq_completeness", type=float, default=50.0)
    ap.add_argument("--mq_contamination", type=float, default=10.0)
    a = ap.parse_args()

    source, q = load_quality(a.input)
    if q.empty:
        fail("quality table has no bins")
    q["tier"] = [
        tier(c, x, a.hq_completeness, a.hq_contamination,
             a.mq_completeness, a.mq_contamination)
        for c, x in zip(q["completeness"], q["contamination"])
    ]

    outdir = Path(a.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    table_path = outdir / "bin-quality.tsv"
    q.to_csv(table_path, sep="\t", index=False)

    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 5.5))
    colors = {"HQ": "#2a9d8f", "MQ": "#e9c46a", "LQ": "#e76f51"}
    for t, sub in q.groupby("tier"):
        ax.scatter(sub["completeness"], sub["contamination"],
                   s=28, alpha=0.8, label=f"{t} (n={len(sub)})", color=colors[t])
    ax.axvline(a.hq_completeness, ls="--", lw=0.8, color="grey")
    ax.axhline(a.hq_contamination, ls="--", lw=0.8, color="grey")
    ax.axvline(a.mq_completeness, ls=":", lw=0.8, color="grey")
    ax.axhline(a.mq_contamination, ls=":", lw=0.8, color="grey")
    ax.set_xlabel("Completeness (%)")
    ax.set_ylabel("Contamination (%)")
    ax.set_title(f"Bin quality ({source})")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig_path = outdir / "bin-quality-scatter.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)

    counts = q["tier"].value_counts().to_dict()
    emit({
        "tool": "mg_bin_qc",
        "source": source,
        "n_bins": int(len(q)),
        "tier_counts": {t: int(counts.get(t, 0)) for t in ("HQ", "MQ", "LQ")},
        "thresholds": {
            "hq": {"completeness": a.hq_completeness, "contamination": a.hq_contamination},
            "mq": {"completeness": a.mq_completeness, "contamination": a.mq_contamination},
        },
        "outputs": {"table": str(table_path), "figure": str(fig_path)},
        "caveat": "MIMAG high-quality also requires 5S/16S/23S rRNA and >=18 tRNAs; "
                  "CheckM/CheckM2 tables do not report these, so HQ here means "
                  "completeness/contamination HQ only.",
    })


if __name__ == "__main__":
    main()