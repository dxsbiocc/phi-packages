#!/usr/bin/env python3
"""clin_inspect: audit clinical survival cohorts, event polarity, follow-up, and EPV."""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import emit, fail, infer_time_unit, kaplan_meier, parse_event_series, read_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--time_col", required=True)
    parser.add_argument("--event_col", required=True)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--matrix", default=None)
    args = parser.parse_args()

    clin = read_table(args.clinical)
    if args.time_col not in clin.columns:
        fail(f"time_col '{args.time_col}' not found in clinical columns: {list(clin.columns)}")
    if args.event_col not in clin.columns:
        fail(f"event_col '{args.event_col}' not found in clinical columns: {list(clin.columns)}")

    sample_col = args.sample_col or clin.columns[0]
    if sample_col not in clin.columns:
        fail(f"sample_col '{sample_col}' not found in clinical columns: {list(clin.columns)}")

    n_total = int(len(clin))
    raw_time = pd.to_numeric(clin[args.time_col], errors="coerce")
    ev_num, raw_ev_counts = parse_event_series(clin[args.event_col], args.event_col)

    nonpos_mask = raw_time.notna() & (raw_time <= 0)
    n_nonpos = int(nonpos_mask.sum())
    valid_mask = raw_time.notna() & (raw_time > 0) & ev_num.notna()

    t_valid = raw_time.loc[valid_mask].to_numpy(dtype=float)
    e_valid = ev_num.loc[valid_mask].to_numpy(dtype=int)
    n_valid = int(len(t_valid))
    if n_valid == 0:
        fail("0 valid survival records with positive time and non-missing event indicator")

    n_events = int(e_valid.sum())
    n_censored = n_valid - n_events
    censor_frac = round(n_censored / float(n_valid), 4)
    ev_rate = round(n_events / float(n_valid), 4)

    # Reverse Kaplan-Meier median follow-up (flipping event <-> censoring)
    rev_km = kaplan_meier(t_valid, 1 - e_valid)
    med_fup = rev_km["median"]
    unit = infer_time_unit(t_valid, args.time_col)
    epv_10 = int(n_events // 10)

    missingness = {
        str(col): round(float(clin[col].isna().mean()), 4)
        for col in clin.columns
    }

    matrix_overlap = None
    if args.matrix:
        mat = read_table(args.matrix)
        mat_samples = set(str(c).strip() for c in mat.columns[1:])
        clin_samples = set(clin[sample_col].astype(str).str.strip())
        shared = sorted(clin_samples & mat_samples)
        matrix_overlap = {
            "matrix": os.path.abspath(args.matrix),
            "n_features": int(len(mat)),
            "n_matrix_samples": len(mat_samples),
            "n_matched_clinical_samples": len(shared),
        }

    warnings: list[str] = []
    if n_nonpos > 0:
        warnings.append(
            f"{n_nonpos} record(s) have follow-up time <= 0 and will be excluded from KM/Cox models"
        )
    if n_events < 10:
        warnings.append(
            f"only {n_events} events observed (< 10); unpenalized multivariable Cox is underpowered"
        )
    if censor_frac > 0.85:
        warnings.append(
            f"heavy censoring ({censor_frac * 100:.1f}% censored); median survival may not be reached"
        )
    if matrix_overlap and matrix_overlap["n_matched_clinical_samples"] == 0:
        warnings.append(
            f"0 samples overlapped between clinical '{sample_col}' and matrix columns; check sample ID format"
        )

    emit(
        {
            "clinical": os.path.abspath(args.clinical),
            "time_col": args.time_col,
            "event_col": args.event_col,
            "sample_col": sample_col,
            "raw_event_value_counts": raw_ev_counts,
            "n_rows_total": n_total,
            "n_valid_survival": n_valid,
            "n_nonpositive_time_dropped": n_nonpos,
            "n_events": n_events,
            "n_censored": n_censored,
            "event_rate": ev_rate,
            "censoring_fraction": censor_frac,
            "time_summary": {
                "min": round(float(np.min(t_valid)), 4),
                "median": round(float(np.median(t_valid)), 4),
                "q75": round(float(np.percentile(t_valid, 75)), 4),
                "max": round(float(np.max(t_valid)), 4),
                "max_event_time": round(float(np.max(t_valid[e_valid == 1])), 4) if n_events > 0 else None,
            },
            "inferred_time_unit": unit,
            "median_followup_reverse_km": med_fup,
            "epv_max_parameters_at_10": epv_10,
            "missingness": missingness,
            "matrix_overlap": matrix_overlap,
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
