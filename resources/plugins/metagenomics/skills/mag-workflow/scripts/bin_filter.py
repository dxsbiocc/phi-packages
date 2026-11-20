"""Copy bins passing completeness/contamination thresholds into a clean set."""
import argparse
import shutil
from pathlib import Path

import pandas as pd

from _common import bin_stem, emit, fail, read_tsv, write_artifact

EXTS = {".fa", ".fasta", ".fna", ".fsa"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bins_dir", required=True)
    ap.add_argument("--quality", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--min_completeness", type=float, default=50.0)
    ap.add_argument("--max_contamination", type=float, default=10.0)
    a = ap.parse_args()

    bins_dir = Path(a.bins_dir)
    if not bins_dir.is_dir():
        fail(f"bins_dir not found: {bins_dir}")
    q = read_tsv(a.quality)
    cols = set(q.columns)
    if {"bin", "completeness", "contamination"} <= cols:
        pass  # mg_bin_qc output
    elif {"Name", "Completeness", "Contamination"} <= cols:
        q = q.rename(columns={"Name": "bin", "Completeness": "completeness",
                              "Contamination": "contamination"})
    elif {"Bin Id", "Completeness", "Contamination"} <= cols:
        q = q.rename(columns={"Bin Id": "bin", "Completeness": "completeness",
                              "Contamination": "contamination"})
    else:
        fail("quality table lacks bin/completeness/contamination columns")
    q["bin"] = q["bin"].astype(str).map(bin_stem)

    fasta = {}
    for p in sorted(bins_dir.iterdir()):
        if p.is_file() and p.suffix.lower() in EXTS:
            fasta[bin_stem(p.name)] = p
    if not fasta:
        fail(f"no .fa/.fasta/.fna files in {bins_dir}")

    missing = [b for b in q["bin"] if b not in fasta]
    if missing:
        fail(f"quality table bins without FASTA in bins_dir: {missing[:10]}"
             + (" ..." if len(missing) > 10 else ""))

    passing = q[(q["completeness"] >= a.min_completeness)
                & (q["contamination"] <= a.max_contamination)]
    outdir = Path(a.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, r in passing.iterrows():
        src = fasta[r["bin"]]
        dst = outdir / src.name
        shutil.copy2(src, dst)
        rows.append({
            "bin": r["bin"], "file": str(dst),
            "completeness": float(r["completeness"]),
            "contamination": float(r["contamination"]),
        })
    manifest = pd.DataFrame(rows, columns=["bin", "file", "completeness", "contamination"])
    manifest_path = outdir / "manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    artifact = write_artifact(manifest_path, "mg_bin_filter", "Filtered bin manifest",
                              skill="mag-workflow", inputs=[a.quality, a.bins_dir],
                              parameters={"min_completeness": a.min_completeness,
                                          "max_contamination": a.max_contamination},
                              shape=manifest.shape)
    emit({
        "tool": "mg_bin_filter",
        "bins_dir": str(bins_dir),
        "n_bins_total": len(fasta),
        "n_passed": len(rows),
        "n_rejected": len(fasta) - len(rows),
        "thresholds": {"min_completeness": a.min_completeness,
                       "max_contamination": a.max_contamination},
        "outputs": {"dir": str(outdir), "manifest": str(manifest_path)},
        "artifacts": [artifact],
    })


if __name__ == "__main__":
    main()