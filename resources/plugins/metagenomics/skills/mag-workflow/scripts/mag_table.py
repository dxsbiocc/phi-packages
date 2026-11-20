"""Join bin quality, GTDB-Tk taxonomy, CoverM abundance and dRep clusters."""
import argparse
from pathlib import Path

import pandas as pd

from _common import bin_stem, emit, fail, read_tsv, setup_matplotlib

RANKS = ["d", "p", "c", "o", "f", "g", "s"]


def load_quality(path):
    df = read_tsv(path)
    cols = set(df.columns)
    if {"bin", "completeness", "contamination"} <= cols:
        df["bin"] = df["bin"].astype(str).map(bin_stem)
        return df
    for name_col in ("Name", "Bin Id"):
        if {name_col, "Completeness", "Contamination"} <= cols:
            out = df.rename(columns={name_col: "bin",
                                     "Completeness": "completeness",
                                     "Contamination": "contamination"})
            out["bin"] = out["bin"].astype(str).map(bin_stem)
            return out
    fail("quality table lacks bin/completeness/contamination columns")


def split_classification(value):
    out = {r: None for r in RANKS}
    if not isinstance(value, str):
        return out
    for part in value.split(";"):
        part = part.strip()
        if len(part) > 3 and part[1:3] == "__" and part[0] in RANKS:
            name = part[3:]
            if name and name not in ("", "uncultured"):
                out[part[0]] = name
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quality", required=True)
    ap.add_argument("--gtdb_summary")
    ap.add_argument("--coverm")
    ap.add_argument("--drep_clusters")
    ap.add_argument("--output_dir", required=True)
    a = ap.parse_args()

    table = load_quality(a.quality)
    joined = ["quality"]

    if a.gtdb_summary:
        g = read_tsv(a.gtdb_summary)
        if not {"user_genome", "classification"} <= set(g.columns):
            fail("GTDB-Tk summary lacks user_genome/classification columns")
        g = g[["user_genome", "classification"]].rename(columns={"user_genome": "bin"})
        g["bin"] = g["bin"].astype(str).map(bin_stem)
        table = table.merge(g, on="bin", how="left")
        ranks = pd.DataFrame([split_classification(v) for v in table["classification"]],
                             index=table.index)
        ranks.columns = [f"{r}__" for r in RANKS]
        table = pd.concat([table, ranks], axis=1)
        joined.append("gtdbtk")

    if a.coverm:
        c = read_tsv(a.coverm)
        if "Genome" not in c.columns:
            fail("CoverM table lacks a Genome column (coverm genome output?)")
        c = c.rename(columns={"Genome": "bin"})
        c["bin"] = c["bin"].astype(str).map(bin_stem)
        table = table.merge(c, on="bin", how="left")
        joined.append("coverm")

    if a.drep_clusters:
        d = read_tsv(a.drep_clusters)
        if "genome" not in d.columns:
            fail("dRep cluster table lacks a genome column (Cdb.csv?)")
        cluster_col = "secondary_cluster" if "secondary_cluster" in d.columns else "cluster"
        if cluster_col not in d.columns:
            fail("dRep table has neither secondary_cluster nor cluster column")
        d = d[["genome", cluster_col]].rename(
            columns={"genome": "bin", cluster_col: "drep_cluster"})
        d["bin"] = d["bin"].astype(str).map(bin_stem)
        table = table.merge(d, on="bin", how="left")
        joined.append("drep")

    outdir = Path(a.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    table_path = outdir / "mag-table.tsv"
    table.to_csv(table_path, sep="\t", index=False)

    figure = None
    if "p__" in table.columns:
        counts = table["p__"].fillna("Unclassified").value_counts()
        plt = setup_matplotlib()
        fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(counts)), 5))
        ax.bar(range(len(counts)), counts.values, color="#457b9d")
        ax.set_xticks(range(len(counts)))
        ax.set_xticklabels(counts.index, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("MAGs")
        ax.set_title("MAGs per phylum (GTDB-Tk)")
        fig.tight_layout()
        figure = outdir / "mags-per-phylum.png"
        fig.savefig(figure, dpi=200)
        plt.close(fig)

    emit({
        "tool": "mg_mag_table",
        "n_bins": int(len(table)),
        "joined": joined,
        "n_classified": int(table["classification"].notna().sum()) if "classification" in table else 0,
        "n_in_drep_clusters": int(table["drep_cluster"].notna().sum()) if "drep_cluster" in table else None,
        "outputs": {"table": str(table_path),
                    **({"figure": str(figure)} if figure else {})},
    })


if __name__ == "__main__":
    main()