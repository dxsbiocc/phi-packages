#!/usr/bin/env python3
"""mb_annotate: deterministic m/z relationships (isotopes, adducts, losses)."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _common import emit, fail, read_matrix, write_artifact

SKILL = "lcms-workflow"
C13_DELTA = 1.003355
ADDUCTS = {
    "pos": {"[M+H]+": 1.007276, "[M+Na]+": 22.989218,
            "[M+K]+": 38.963158, "[M+NH4]+": 18.033823},
    "neg": {"[M-H]-": -1.007276, "[M+Cl]-": 34.969402,
            "[M+FA-H]-": 44.998201},
}
LOSSES = {"H2O": 18.010565, "CO2": 43.989829}
BASE = {"pos": "[M+H]+", "neg": "[M-H]-"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--mode", choices=["pos", "neg"], default="pos")
    ap.add_argument("--ppm", type=float, default=10.0)
    ap.add_argument("--rt_tol", type=float, default=None)
    args = ap.parse_args()

    meta, X = read_matrix(args.input)
    if "mz" not in meta.columns or meta["mz"].isna().all():
        fail("input has no m/z column; annotation needs mz (vendor aliases recognized)")
    has_rt = "rt" in meta.columns and meta["rt"].notna().any()
    rt = meta["rt"].to_numpy(dtype=float) if has_rt else np.zeros(len(meta))
    rt_unit = "min"
    if has_rt and np.nanmedian(rt) > 60:
        rt_unit = "sec"
    if args.rt_tol is not None:
        rt_tol = args.rt_tol
    else:
        rt_tol = 0.2 * (60.0 if rt_unit == "sec" else 1.0)

    df = pd.DataFrame({"mz": meta["mz"].to_numpy(dtype=float), "rt": rt},
                      index=meta.index).dropna(subset=["mz"]).sort_values("mz")
    feats = df.index.to_numpy()
    mz = df["mz"].to_numpy()
    rtv = df["rt"].to_numpy()
    n = len(mz)
    adducts = ADDUCTS[args.mode]
    base_delta = adducts[BASE[args.mode]]
    deltas = []
    for a1, d1 in adducts.items():
        for a2, d2 in adducts.items():
            if d2 > d1:
                deltas.append((d2 - d1, "adduct_pair", f"{a1} / {a2}"))
    deltas.append((C13_DELTA, "isotope_13C", "13C isotope"))
    for name, d in LOSSES.items():
        deltas.append((d, "neutral_loss", f"{name} loss"))
    max_delta = max(d for d, _, _ in deltas)

    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)

    rows = []
    pair_neutral = []
    for i in range(n):
        window = mz[i] + max_delta + args.ppm * (mz[i] + max_delta) / 1e6
        j = i + 1
        while j < n and mz[j] <= window:
            dm = mz[j] - mz[i]
            if has_rt and abs(rtv[j] - rtv[i]) > rt_tol:
                j += 1
                continue
            for delta, kind, label in deltas:
                if abs(dm - delta) <= args.ppm * mz[j] / 1e6:
                    row = {"feature_1": str(feats[i]), "feature_2": str(feats[j]),
                           "relation": kind, "detail": label,
                           "delta_mz": float(dm),
                           "ppm_error": float(abs(dm - delta) / mz[j] * 1e6),
                           "rt_diff": float(abs(rtv[j] - rtv[i])) if has_rt else None}
                    if kind == "adduct_pair":
                        a1 = label.split(" / ")[0]
                        row["neutral_mass"] = float(mz[i] - adducts[a1])
                        pair_neutral.append((i, row["neutral_mass"]))
                    rows.append(row)
                    union(i, j)
            j += 1

    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    neutral_by_root = {}
    for i, nm in pair_neutral:
        neutral_by_root.setdefault(find(i), []).append(nm)
    grows = []
    for gid, (root, members) in enumerate(sorted(groups.items(), key=lambda kv: kv[1][0]), 1):
        if len(members) < 2:
            continue
        neutral = neutral_by_root.get(root)
        if neutral:
            consensus = float(np.median(neutral))
        else:
            consensus = float(np.median(mz[members] - base_delta))
        grows.append({
            "group": f"P{gid:04d}",
            "n_features": len(members),
            "members": ";".join(str(feats[m]) for m in members),
            "neutral_mass": consensus,
            "rt_median": float(np.median(rtv[members])) if has_rt else None,
        })

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    ann = pd.DataFrame(rows, columns=["feature_1", "feature_2", "relation", "detail",
                                      "delta_mz", "ppm_error", "rt_diff", "neutral_mass"])
    ann_path = out / "annotations.tsv"
    ann.to_csv(ann_path, sep="\t", index=False)
    grp = pd.DataFrame(grows, columns=["group", "n_features", "members",
                                       "neutral_mass", "rt_median"])
    grp_path = out / "pseudospectra.tsv"
    grp.to_csv(grp_path, sep="\t", index=False)

    params = {"mode": args.mode, "ppm": args.ppm, "rt_tol": rt_tol}
    artifacts = [
        write_artifact(ann_path, "mb_annotate", "m/z relationship annotations",
                       skill=SKILL, inputs=[args.input], parameters=params,
                       shape=ann.shape),
        write_artifact(grp_path, "mb_annotate", "Pseudospectrum groups",
                       skill=SKILL, inputs=[args.input], parameters=params,
                       shape=grp.shape),
    ]
    counts = ann["relation"].value_counts().to_dict() if len(ann) else {}
    emit({
        "input": args.input,
        "mode": args.mode,
        "ppm": args.ppm,
        "rt_tol": rt_tol,
        "rt_unit": rt_unit,
        "has_rt": bool(has_rt),
        "n_features": n,
        "n_relationships": int(len(ann)),
        "relationships_by_type": {k: int(v) for k, v in counts.items()},
        "n_pseudospectra": int(len(grp)),
        "caveat": "m/z-relationship annotation only (MSI level 3); "
                  "not identification — confirm with MS2 and standards",
        "output_dir": str(out),
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()