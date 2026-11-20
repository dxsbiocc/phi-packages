"""Summarize QUAST (report.tsv) assembly statistics."""
import argparse
from pathlib import Path

from _common import emit, fail, read_tsv

METRICS = [
    "# contigs (>= 0 bp)",
    "# contigs",
    "Total length (>= 0 bp)",
    "Total length",
    "Largest contig",
    "N50",
    "N75",
    "L50",
    "L75",
    "GC (%)",
]
CANON = {
    "# contigs (>= 0 bp)": "contigs",
    "# contigs": "contigs",
    "Total length (>= 0 bp)": "total_length",
    "Total length": "total_length",
    "Largest contig": "largest_contig",
    "N50": "n50",
    "N75": "n75",
    "L50": "l50",
    "L75": "l75",
    "GC (%)": "gc_percent",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    a = ap.parse_args()

    path = Path(a.input)
    if not path.is_file():
        fail(f"input not found: {path}")
    df = read_tsv(path)
    if df.shape[1] < 2 or df.columns[0] != "Assembly":
        fail("not a QUAST report.tsv (first column must be 'Assembly')")

    wanted = {}
    for _, row in df.iterrows():
        metric = str(row.iloc[0]).strip()
        if metric in CANON and CANON[metric] not in wanted:
            wanted[CANON[metric]] = row
    if "n50" not in wanted:
        fail("no N50 row found; is this a QUAST report?")

    assemblies = []
    for col in df.columns[1:]:
        rec = {"assembly": str(col)}
        for key, row in wanted.items():
            try:
                rec[key] = float(row[col])
            except (TypeError, ValueError):
                rec[key] = None
        assemblies.append(rec)

    emit({
        "tool": "mg_assembly_stats",
        "input": str(path),
        "assemblies": assemblies,
        "n_assemblies": len(assemblies),
        "note": "N50/L50 use the >=0 bp (all-contig) rows when present; "
                "QUAST also reports >=1000/5000/... bp variants in the file.",
    })


if __name__ == "__main__":
    main()