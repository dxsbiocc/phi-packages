#!/usr/bin/env python3
"""mg_function: deterministic operations on HUMAnN-style functional tables."""
import argparse

import pandas as pd

from _common import emit, fail, read_table, write_artifact


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--action", required=True,
                    choices=["normalize", "regroup", "stratified", "unstratified"])
    ap.add_argument("--output", required=True)
    ap.add_argument("--unit", choices=["relative", "cpm"], default="relative")
    ap.add_argument("--mapping")
    args = ap.parse_args()

    df = read_table(args.input)
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.isna().any().any():
        fail("functional table has non-numeric columns")
    rows_in = len(df)
    notes = []

    if args.action == "normalize":
        totals = df.sum(axis=0)
        if (totals == 0).any():
            fail(f"columns sum to zero: {totals[totals == 0].index.tolist()}")
        df = df.div(totals, axis=1)
        if args.unit == "cpm":
            df = df * 1e6
        notes.append("plain column renormalization; for HUMAnN's stratification-aware "
                     "renorm use the humann3 wrapper utilities")
    elif args.action == "regroup":
        if not args.mapping:
            fail("regroup requires --mapping (two-column tsv: feature -> group)")
        mapping = pd.read_csv(args.mapping, sep="\t", header=None, names=["feature", "group"])
        groups = mapping.set_index("feature")["group"]
        assigned = df.index.map(lambda i: groups.get(i.split("|")[0], "UNGROUPED"))
        df = df.groupby(assigned).sum()
        df.index.name = "group"
    elif args.action == "stratified":
        df = df[df.index.str.contains(r"\|", regex=True)]
    elif args.action == "unstratified":
        df = df[~df.index.str.contains(r"\|", regex=True)]

    if df.empty:
        fail(f"action {args.action} left no rows")
    df.to_csv(args.output, sep="\t")
    inputs = [args.input] + ([args.mapping] if args.mapping else [])
    artifact = write_artifact(args.output, "mg_function", f"Functional table ({args.action})",
                              skill="profiling-workflow", inputs=inputs,
                              parameters={"action": args.action, "unit": args.unit},
                              shape=df.shape)
    emit({
        "action": args.action,
        "unit": args.unit if args.action == "normalize" else None,
        "rows_in": rows_in,
        "rows_out": int(len(df)),
        "columns": df.columns.tolist(),
        "notes": notes,
        "output": args.output,
        "artifacts": [artifact],
    })


if __name__ == "__main__":
    main()