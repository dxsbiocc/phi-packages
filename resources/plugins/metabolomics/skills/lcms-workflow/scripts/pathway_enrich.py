#!/usr/bin/env python3
"""mb_pathway: pathway enrichment (ID hypergeometric / mummichog-style m/z)."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from _common import emit, fail, read_tsv, setup_matplotlib, write_artifact

SKILL = "lcms-workflow"
ADDUCTS = {
    "pos": {"[M+H]+": 1.007276, "[M+Na]+": 22.989218,
            "[M+K]+": 38.963158, "[M+NH4]+": 18.033823},
    "neg": {"[M-H]-": -1.007276, "[M+Cl]-": 34.969402,
            "[M+FA-H]-": 44.998201},
}


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def read_list(path):
    ext = Path(path).suffix.lower()
    if ext in {".tsv", ".csv"}:
        df = read_tsv(path)
        return [str(v).strip() for v in df.iloc[:, 0].dropna()]
    return [l.strip() for l in Path(path).read_text().splitlines() if l.strip()]


def read_mapping(path):
    mp = pd.read_csv(path, sep="\t", header=None, names=["pathway", "compound"],
                     dtype=str)
    mp = mp.apply(lambda s: s.str.strip())
    first = mp.iloc[0].str.lower()
    if first["pathway"] in {"pathway", "pathway_id"} and first["compound"] in {
            "compound", "compound_id"}:
        mp = mp.iloc[1:]
    return mp.dropna()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["ids", "masses"])
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--mapping")
    ap.add_argument("--compound_db")
    ap.add_argument("--background")
    ap.add_argument("--ion_mode", choices=["pos", "neg"], default="pos")
    ap.add_argument("--ppm", type=float, default=10.0)
    ap.add_argument("--n_perm", type=int, default=100)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()
    if not args.mapping:
        fail("--mapping (pathway -> compound TSV) is required for both modes")
    mp = read_mapping(args.mapping)
    if mp.empty:
        fail("pathway mapping is empty")
    pathways = {p: set(g["compound"]) for p, g in mp.groupby("pathway")}

    if args.mode == "ids":
        ids = read_list(args.input)
        bg = set(mp["compound"])
        hits = set(ids) & bg
        n_unmatched = len(set(ids) - bg)
        N, n = len(bg), len(hits)
        rows = []
        for pname, members in sorted(pathways.items()):
            K = len(members)
            k = len(hits & members)
            p = float(hypergeom.sf(k - 1, N, K, n)) if n and K else 1.0
            rows.append({"pathway": pname, "pathway_size": K, "overlap": k,
                         "overlap_compounds": ";".join(sorted(hits & members)),
                         "pvalue": p})
        n_input = len(ids)
        n_matched = len(hits)
        background_desc = {"kind": "compounds", "n": N}
        extra = {"n_unmatched_ids": n_unmatched}
    else:
        if not args.compound_db or not args.background:
            fail("mode=masses requires --compound_db and --background")
        db = read_tsv(args.compound_db)
        cols = {str(c).strip().lower(): c for c in db.columns}
        if "compound_id" not in cols or "mass" not in cols:
            fail("compound_db needs compound_id and mass columns")
        db_id = db[cols["compound_id"]].astype(str).to_numpy()
        db_mass = pd.to_numeric(db[cols["mass"]], errors="coerce").to_numpy()
        ok = np.isfinite(db_mass)
        db_id, db_mass = db_id[ok], db_mass[ok]
        order = np.argsort(db_mass)
        db_id, db_mass = db_id[order], db_mass[order]
        sig = np.array([float(v) for v in read_list(args.input)])
        bgm = np.array([float(v) for v in read_list(args.background)])
        if len(sig) > len(bgm):
            fail("significant m/z list is larger than the background list")
        deltas = list(ADDUCTS[args.ion_mode].values())

        def candidates(mz):
            found = set()
            for d in deltas:
                target = mz - d
                tol = args.ppm * mz / 1e6
                lo = np.searchsorted(db_mass, target - tol)
                hi = np.searchsorted(db_mass, target + tol)
                found.update(db_id[lo:hi])
            return found

        bg_sets = [candidates(m) for m in bgm]
        sig_sets = [candidates(m) for m in sig]
        hits = set().union(*sig_sets) if sig_sets else set()
        n_matched = sum(1 for s in sig_sets if s)
        k_obs = {p: len(hits & members) for p, members in pathways.items()}
        rng = np.random.default_rng(args.seed)
        counts = {p: 0 for p in pathways}
        for _ in range(args.n_perm):
            idx = rng.choice(len(bgm), size=len(sig), replace=False)
            perm_hits = set()
            for i in idx:
                perm_hits |= bg_sets[i]
            for p, members in pathways.items():
                if len(perm_hits & members) >= k_obs[p]:
                    counts[p] += 1
        rows = []
        for pname, members in sorted(pathways.items()):
            p = (1 + counts[pname]) / (1 + args.n_perm)
            rows.append({"pathway": pname, "pathway_size": len(members),
                         "overlap": k_obs[pname],
                         "overlap_compounds": ";".join(sorted(hits & members)),
                         "pvalue": float(p)})
        n_input = len(sig)
        background_desc = {"kind": "mz_values", "n": len(bgm)}
        extra = {"n_perm": args.n_perm, "seed": args.seed,
                 "ion_mode": args.ion_mode, "ppm": args.ppm,
                 "n_compound_db": len(db_id)}

    table = pd.DataFrame(rows)
    table["fdr"] = bh_fdr(table["pvalue"].to_numpy())
    table["significant"] = table["fdr"] < 0.05
    table = table.sort_values(["pvalue", "pathway"]).reset_index(drop=True)

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    table_path = out / "pathway-enrichment.tsv"
    table.to_csv(table_path, sep="\t", index=False)

    plt = setup_matplotlib()
    top = table.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, max(3, 0.3 * len(top))))
    colors = ["#d62728" if s else "#9aa0a6" for s in top["significant"]]
    with np.errstate(divide="ignore"):
        vals = -np.log10(top["pvalue"].clip(lower=1e-300))
    ax.barh(top["pathway"], vals, color=colors)
    ax.set_xlabel("-log10 p")
    ax.set_title(f"Pathway enrichment ({args.mode})")
    fig.tight_layout()
    fig_path = out / "pathway-barplot.png"
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)

    params = {"mode": args.mode}
    params.update(extra)
    artifacts = [
        write_artifact(table_path, "mb_pathway", "Pathway enrichment table",
                       skill=SKILL,
                       inputs=[args.input, args.mapping]
                       + ([args.compound_db, args.background] if args.mode == "masses" else []),
                       parameters=params, shape=table.shape),
        write_artifact(fig_path, "mb_pathway", "Top pathways barplot",
                       skill=SKILL,
                       inputs=[args.input, args.mapping]
                       + ([args.compound_db, args.background] if args.mode == "masses" else []),
                       parameters=params),
    ]
    emit({
        "mode": args.mode,
        "input": args.input,
        "n_input": n_input,
        "n_matched_input": n_matched,
        "background": background_desc,
        "n_pathways": len(pathways),
        "n_significant": int(table["significant"].sum()),
        "fdr_threshold": 0.05,
        "top": [{"pathway": r["pathway"], "overlap": int(r["overlap"]),
                 "pvalue": float(r["pvalue"]), "fdr": float(r["fdr"])}
                for _, r in table.head(10).iterrows()],
        **extra,
        "output_dir": str(out),
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()