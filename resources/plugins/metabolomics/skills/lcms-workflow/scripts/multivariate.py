#!/usr/bin/env python3
"""mb_multivariate: PCA (SVD) + PLS-DA (NIPALS) with VIP and permutation Q2."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from _common import emit, fail, read_matrix, read_tsv, setup_matplotlib, write_artifact

SKILL = "lcms-workflow"


def scale_matrix(A, how):
    mean = A.mean(axis=0)
    Ac = A - mean
    std = Ac.std(axis=0, ddof=1)
    if how == "uv":
        s = std
    elif how == "pareto":
        s = np.sqrt(std)
    else:
        s = np.ones_like(std)
    s = np.where(s > 0, s, 1.0)
    return Ac / s, mean, s


def nipals_pls(X, Y, a, tol=1e-10, max_iter=500):
    Xr = X.copy()
    Yr = Y.copy()
    n, p = Xr.shape
    g = Yr.shape[1]
    W = np.zeros((p, a))
    P = np.zeros((p, a))
    Q = np.zeros((g, a))
    T = np.zeros((n, a))
    k = 0
    for k in range(a):
        u = Yr[:, int(np.argmax(Yr.var(axis=0)))].copy()
        t_prev = None
        t = np.zeros(n)
        for _ in range(max_iter):
            w = Xr.T @ u
            nw = np.linalg.norm(w)
            if nw < 1e-14:
                break
            w = w / nw
            t = Xr @ w
            q = Yr.T @ t
            nq = np.linalg.norm(q)
            if nq < 1e-14:
                break
            q = q / nq
            u = Yr @ q
            if t_prev is not None and np.linalg.norm(t - t_prev) <= tol * max(np.linalg.norm(t), 1e-14):
                break
            t_prev = t.copy()
        tt = float(t @ t)
        if tt < 1e-14:
            break
        p_load = Xr.T @ t / tt
        q_load = Yr.T @ t / tt
        Xr = Xr - np.outer(t, p_load)
        Yr = Yr - np.outer(t, q_load)
        W[:, k] = w
        P[:, k] = p_load
        Q[:, k] = q_load
        T[:, k] = t
    else:
        k = a
    return W[:, :k], P[:, :k], Q[:, :k], T[:, :k]


def pls_fit(X, Y, a):
    W, P, Q, T = nipals_pls(X, Y, a)
    k = W.shape[1]
    y0 = Y
    yres = y0 - T @ Q.T
    ss_tot = float((y0 ** 2).sum())
    r2y = 1.0 - float((yres ** 2).sum()) / ss_tot if ss_tot > 0 else 0.0
    B = None
    if k:
        PtW = P.T @ W
        B = W @ np.linalg.solve(PtW, Q.T)
    return {"W": W, "P": P, "Q": Q, "T": T, "B": B, "r2y": r2y, "k": k}


def vip_scores(fit, p):
    W, Q, T = fit["W"], fit["Q"], fit["T"]
    k = W.shape[1]
    if k == 0:
        return np.zeros(p)
    s = (Q ** 2).sum(axis=0) * (T ** 2).sum(axis=0)
    wn = W / np.linalg.norm(W, axis=0)
    total = s.sum()
    if total <= 0:
        return np.zeros(p)
    return np.sqrt(p * ((wn ** 2) @ s) / total)


def one_hot(labels):
    levels = sorted(set(labels))
    Y = np.zeros((len(labels), len(levels)))
    for i, lab in enumerate(labels):
        Y[i, levels.index(lab)] = 1.0
    return Y, levels


def cv_folds(n, seed):
    rng = np.random.default_rng(seed)
    if n <= 50:
        return [[i] for i in range(n)]
    order = rng.permutation(n)
    return [order[j::7].tolist() for j in range(7)]


def q2_cv(A, labels, a, scaling, folds):
    Y, _ = one_hot(labels)
    y_mean_global = Y.mean(axis=0)
    tss = float(((Y - y_mean_global) ** 2).sum())
    press = 0.0
    for fold in folds:
        te = np.array(fold)
        tr = np.array([i for i in range(len(labels)) if i not in set(fold)])
        if len(tr) < 3:
            continue
        Xtr, mean, s = scale_matrix(A[tr], scaling)
        ytr_mean = Y[tr].mean(axis=0)
        fit = pls_fit(Xtr, Y[tr] - ytr_mean, a)
        if fit["B"] is None:
            return None
        Xte = (A[te] - mean) / s
        yhat = Xte @ fit["B"] + ytr_mean
        press += float(((Y[te] - yhat) ** 2).sum())
    if tss <= 0:
        return None
    return 1.0 - press / tss


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--group_col", default="group")
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--scaling", choices=["none", "uv", "pareto"], default="uv")
    ap.add_argument("--n_components", type=int, default=2)
    ap.add_argument("--n_perm", type=int, default=99)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    meta, X = read_matrix(args.input)
    md = read_tsv(args.metadata)
    md = md.set_index(md.columns[0])
    md.index = md.index.astype(str)
    if args.group_col not in md.columns:
        fail(f"group column not in metadata: {args.group_col}")
    missing = [c for c in X.columns if str(c) not in md.index]
    if missing:
        fail(f"samples missing from metadata: {missing}")
    labels = [str(md.loc[str(c), args.group_col]) for c in X.columns]
    levels = sorted(set(labels))
    if len(levels) < 2:
        fail("need at least two groups for PLS-DA")
    A = X.T.to_numpy(dtype=float)
    if np.isnan(A).any():
        fail("matrix has missing values; run mb_preprocess first")
    n, p = A.shape
    a = max(1, min(args.n_components, n - 2, p))
    clamped = a != args.n_components

    As, mean, s = scale_matrix(A, args.scaling)
    U, sv, Vt = np.linalg.svd(As, full_matrices=False)
    pca_scores = U[:, :a] * sv[:a]
    evr = (sv ** 2) / (sv ** 2).sum()
    loadings = Vt[:a].T

    Y, _ = one_hot(labels)
    y_mean = Y.mean(axis=0)
    fit = pls_fit(As, Y - y_mean, a)
    pls_scores = fit["T"]
    vip = vip_scores(fit, p)
    r2y_cum = []
    for k in range(1, fit["k"] + 1):
        sub = pls_fit(As, Y - y_mean, k)
        r2y_cum.append(sub["r2y"])
    folds = cv_folds(n, args.seed)
    q2_cum = [q2_cv(A, labels, k, args.scaling, folds)
              for k in range(1, a + 1)]
    q2_obs = q2_cum[-1] if q2_cum else None

    perm = {"n_perm": 0}
    if args.n_perm > 0 and q2_obs is not None:
        rng = np.random.default_rng(args.seed)
        q2_null = []
        r2_null = []
        for _ in range(args.n_perm):
            yp = labels.copy()
            rng.shuffle(yp)
            Yp, _ = one_hot(yp)
            fp = pls_fit(As, Yp - Yp.mean(axis=0), a)
            r2_null.append(fp["r2y"])
            q2p = q2_cv(A, yp, a, args.scaling, folds)
            if q2p is not None:
                q2_null.append(q2p)
        perm = {
            "n_perm": args.n_perm,
            "q2_null_mean": float(np.mean(q2_null)) if q2_null else None,
            "r2y_null_mean": float(np.mean(r2_null)),
            "q2_p": (1 + sum(1 for v in q2_null if v >= q2_obs)) / (1 + len(q2_null)),
            "r2y_p": (1 + sum(1 for v in r2_null if v >= fit["r2y"])) / (1 + len(r2_null)),
        }

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    samples = [str(c) for c in X.columns]
    group_series = pd.Series(labels, index=samples)

    def scores_df(scores, prefix):
        df = pd.DataFrame(scores, index=samples,
                          columns=[f"{prefix}{i+1}" for i in range(scores.shape[1])])
        df.insert(0, "group", group_series)
        df.index.name = "sample"
        return df

    pca_df = scores_df(pca_scores, "PC")
    pls_df = scores_df(pls_scores, "LV")
    vip_df = pd.DataFrame({"vip": vip}, index=X.index).sort_values("vip", ascending=False)
    vip_df.index.name = "feature"
    load_df = pd.DataFrame(loadings, index=X.index,
                           columns=[f"PC{i+1}" for i in range(loadings.shape[1])])
    load_df.index.name = "feature"

    paths = {}
    for name, df in [("pca-scores.tsv", pca_df), ("plsda-scores.tsv", pls_df),
                     ("vip.tsv", vip_df), ("pca-loadings.tsv", load_df)]:
        df.to_csv(out / name, sep="\t")
        paths[name] = out / name

    plt = setup_matplotlib()
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(levels), 2)))
    figs = []
    for name, df, prefix, title in [
            ("pca-scores.png", pca_df, "PC", "PCA scores"),
            ("plsda-scores.png", pls_df, "LV", "PLS-DA scores")]:
        fig, ax = plt.subplots(figsize=(6, 5))
        for i, lev in enumerate(levels):
            sub = df[df["group"] == lev]
            ax.scatter(sub[f"{prefix}1"], sub[f"{prefix}2"],
                       color=colors[i % len(colors)], label=lev, s=40, edgecolor="white")
        if prefix == "PC":
            ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}%)")
            ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}%)" if len(evr) > 1 else "PC2")
        else:
            ax.set_xlabel("LV1")
            ax.set_ylabel("LV2")
        ax.set_title(title)
        ax.legend()
        ax.axhline(0, color="grey", lw=0.5)
        ax.axvline(0, color="grey", lw=0.5)
        fig.tight_layout()
        fig.savefig(out / name, dpi=150)
        plt.close(fig)
        figs.append(out / name)

    params = {"scaling": args.scaling, "n_components": a,
              "n_perm": args.n_perm, "seed": args.seed,
              "group_col": args.group_col}
    artifacts = [
        write_artifact(paths["pca-scores.tsv"], "mb_multivariate", "PCA scores",
                       skill=SKILL, inputs=[args.input, args.metadata],
                       parameters=params, shape=pca_df.shape),
        write_artifact(paths["plsda-scores.tsv"], "mb_multivariate", "PLS-DA scores",
                       skill=SKILL, inputs=[args.input, args.metadata],
                       parameters=params, shape=pls_df.shape),
        write_artifact(paths["vip.tsv"], "mb_multivariate", "PLS-DA VIP scores",
                       skill=SKILL, inputs=[args.input, args.metadata],
                       parameters=params, shape=vip_df.shape),
        write_artifact(paths["pca-loadings.tsv"], "mb_multivariate", "PCA loadings",
                       skill=SKILL, inputs=[args.input, args.metadata],
                       parameters=params, shape=load_df.shape),
    ] + [
        write_artifact(f, "mb_multivariate", f.stem.replace("-", " ").title(),
                       skill=SKILL, inputs=[args.input, args.metadata],
                       parameters=params)
        for f in figs
    ]
    emit({
        "input": args.input,
        "n_samples": n,
        "n_features": p,
        "groups": {lev: labels.count(lev) for lev in levels},
        "scaling": args.scaling,
        "n_components": a,
        "n_components_clamped": clamped,
        "cv_scheme": "loo" if n <= 50 else "7-fold",
        "pca": {"explained_variance_ratio": [float(v) for v in evr[:a]]},
        "plsda": {
            "r2y_cum": [float(v) for v in r2y_cum],
            "q2_cum": [None if v is None else float(v) for v in q2_cum],
            "permutation": perm,
        },
        "vip_top": [{"feature": str(i), "vip": float(v)}
                    for i, v in vip_df["vip"].head(10).items()],
        "output_dir": str(out),
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()