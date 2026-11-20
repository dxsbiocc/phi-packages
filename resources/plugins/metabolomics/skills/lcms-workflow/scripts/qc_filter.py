#!/usr/bin/env python3
"""mb_qc: QC-sample RSD + detection-rate feature filtering."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import (BLANK_PATTERN, emit, fail, qc_samples, read_matrix,
                     setup_matplotlib, write_artifact)

SKILL = "lcms-workflow"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--rsd_threshold", type=float, default=30.0)
    ap.add_argument("--detection_rate", type=float, default=0.8)
    args = ap.parse_args()
    if not 0 <= args.detection_rate <= 1:
        fail("detection_rate must be between 0 and 1")

    meta, X = read_matrix(args.input)
    qc = qc_samples(X.columns)
    if not qc:
        fail("no QC samples found (names containing qc/pool); "
             "QC-based filtering is impossible — report the gap")
    blanks = [c for c in X.columns if BLANK_PATTERN.search(str(c))]
    bio = [c for c in X.columns if c not in qc and c not in blanks]
    if not bio:
        fail("no biological samples besides QC/blank columns")

    Q = X[qc]
    qc_mean = Q.mean(axis=1)
    qc_std = Q.std(axis=1, ddof=1)
    n_qc_det = (Q.notna() & (Q > 0)).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        rsd = np.where(qc_mean > 0, qc_std / qc_mean * 100.0, np.inf)
    rsd = pd.Series(rsd, index=X.index)
    rsd[n_qc_det < 2] = np.inf
    B = X[bio]
    detection = (B.notna() & (B > 0)).sum(axis=1) / len(bio)

    keep = (rsd <= args.rsd_threshold) & (detection >= args.detection_rate)
    removed_rsd = int(((rsd > args.rsd_threshold)).sum())
    removed_det = int(((rsd <= args.rsd_threshold)
                       & (detection < args.detection_rate)).sum())
    if not keep.any():
        fail("no features passed QC filtering; relax thresholds consciously")

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    kept = X[keep]
    matrix_path = out / "filtered-matrix.tsv"
    out_df = pd.concat([meta.loc[kept.index], kept], axis=1)
    out_df.index.name = "feature"
    out_df.to_csv(matrix_path, sep="\t")

    report = pd.DataFrame({
        "rsd_qc_pct": rsd, "detection_rate": detection,
        "n_qc_detected": n_qc_det, "kept": keep,
    })
    report.index.name = "feature"
    report_path = out / "qc-report.tsv"
    report.to_csv(report_path, sep="\t")

    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(6, 4))
    finite = rsd[np.isfinite(rsd)]
    ax.hist(finite.clip(upper=100), bins=40, color="#4878cf", edgecolor="white")
    ax.axvline(args.rsd_threshold, color="#d62728", linestyle="--",
               label=f"threshold {args.rsd_threshold:g}%")
    ax.set_xlabel("QC RSD (%)")
    ax.set_ylabel("features")
    ax.set_title(f"QC RSD distribution ({len(qc)} QC samples)")
    ax.legend()
    fig.tight_layout()
    fig_path = out / "qc-rsd-distribution.png"
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)

    params = {"rsd_threshold": args.rsd_threshold,
              "detection_rate": args.detection_rate}
    artifacts = [
        write_artifact(matrix_path, "mb_qc", "QC-filtered feature matrix",
                       skill=SKILL, inputs=[args.input], parameters=params,
                       shape=kept.shape),
        write_artifact(report_path, "mb_qc", "Per-feature QC report",
                       skill=SKILL, inputs=[args.input], parameters=params,
                       shape=report.shape),
        write_artifact(fig_path, "mb_qc", "QC RSD distribution",
                       skill=SKILL, inputs=[args.input], parameters=params),
    ]
    emit({
        "input": args.input,
        "qc_samples": [str(c) for c in qc],
        "n_qc_samples": len(qc),
        "n_biological_samples": len(bio),
        "blank_samples": [str(c) for c in blanks],
        "rsd_threshold": args.rsd_threshold,
        "detection_rate": args.detection_rate,
        "n_features_before": int(X.shape[0]),
        "n_features_after": int(keep.sum()),
        "removed_by_rsd": removed_rsd,
        "removed_by_detection": removed_det,
        "median_qc_rsd_pct": float(np.median(finite)) if len(finite) else None,
        "output": str(matrix_path),
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()