#!/usr/bin/env python3
"""mg_diversity: alpha/beta diversity, PCoA ordination, PERMANOVA, PERMDISP,
and multi-covariate distance-based PERMANOVA.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from _common import emit, fail, read_matrix, read_metadata, setup_matplotlib, write_artifact


def alpha_indices(counts, values):
    x = np.asarray(counts, dtype=float)
    total = x.sum()
    out = {"observed": None, "chao1": None, "shannon": None, "simpson": None, "invsimpson": None}
    if total <= 0:
        return out
    p = x[x > 0] / total
    out["shannon"] = float(-(p * np.log(p)).sum())
    d = float((p * p).sum())
    out["simpson"] = 1.0 - d
    out["invsimpson"] = 1.0 / d if d > 0 else None
    if values == "counts":
        s_obs = int((x > 0).sum())
        f1 = int((x == 1).sum())
        f2 = int((x == 2).sum())
        out["observed"] = s_obs
        if f2 > 0:
            out["chao1"] = s_obs + f1 * f1 / (2.0 * f2)
        else:
            out["chao1"] = s_obs + f1 * (f1 - 1) / 2.0
    return out


def pcoa(dist):
    n = dist.shape[0]
    d2 = dist * dist
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos = eigvals > 1e-10
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    var = eigvals[pos] / eigvals[pos].sum() if pos.any() else np.array([])
    return coords, var


def permanova(dist, labels, permutations, seed):
    groups = np.unique(labels)
    a, n = len(groups), len(labels)
    if a < 2 or a >= n:
        return None
    d2 = dist * dist
    ss_total = d2.sum() / n
    ss_within = 0.0
    for g in groups:
        idx = np.where(labels == g)[0]
        ng = len(idx)
        if ng > 1:
            ss_within += d2[np.ix_(idx, idx)].sum() / ng
    ss_between = ss_total - ss_within
    if ss_within <= 0:
        return None
    f_obs = (ss_between / (a - 1)) / (ss_within / (n - a))
    r2 = ss_between / ss_total if ss_total > 0 else 0.0
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(permutations):
        perm = rng.permutation(labels)
        ss_w = 0.0
        for g in groups:
            idx = np.where(perm == g)[0]
            ng = len(idx)
            if ng > 1:
                ss_w += d2[np.ix_(idx, idx)].sum() / ng
        ss_b = ss_total - ss_w
        f_perm = (ss_b / (a - 1)) / (ss_w / (n - a)) if ss_w > 0 else np.inf
        if f_perm >= f_obs:
            ge += 1
    return {
        "pseudo_f": float(f_obs),
        "r2": float(r2),
        "p_value": float((1 + ge) / (1 + permutations)),
        "permutations": int(permutations),
        "groups": groups.tolist(),
    }


def permdisp(coords, labels, permutations, seed):
    """Anderson's PERMDISP: test homogeneity of multivariate dispersions
    (distances to group centroids in PCoA space) via seeded permutation F-test.
    """
    groups = np.unique(labels)
    a, n = len(groups), len(labels)
    if a < 2 or a >= n or coords.shape[1] == 0:
        return None

    def _f_stat(lbls):
        dists = np.zeros(n, dtype=float)
        group_means = {}
        for g in groups:
            idx = np.where(lbls == g)[0]
            centroid = coords[idx].mean(axis=0)
            g_dists = np.linalg.norm(coords[idx] - centroid, axis=1)
            dists[idx] = g_dists
            group_means[str(g)] = float(g_dists.mean())
        grand = dists.mean()
        ss_b = sum(len(np.where(lbls == g)[0]) * (group_means[str(g)] - grand) ** 2 for g in groups)
        ss_w = sum(((dists[lbls == g] - group_means[str(g)]) ** 2).sum() for g in groups)
        if ss_w <= 1e-15:
            return None, group_means
        f_val = (ss_b / (a - 1)) / (ss_w / (n - a))
        return float(f_val), group_means

    f_obs, mean_dists = _f_stat(labels)
    if f_obs is None:
        return None
    rng = np.random.default_rng(seed + 17)
    ge = 0
    for _ in range(permutations):
        f_perm, _ = _f_stat(rng.permutation(labels))
        if f_perm is not None and f_perm >= f_obs:
            ge += 1
    return {
        "f_stat": float(f_obs),
        "p_value": float((1 + ge) / (1 + permutations)),
        "permutations": int(permutations),
        "group_mean_distances": mean_dists,
    }


def covariate_permanova(dist, meta, covariates, permutations, seed):
    """Distance-based PERMANOVA (Gower-centered db-RDA pseudo-F) for each
    categorical or numeric covariate in metadata.
    """
    results = []
    for idx_cov, col in enumerate(covariates):
        if col not in meta.columns:
            fail(f"covariate column not in metadata: {col}; columns: {meta.columns.tolist()}")
        series = meta[col]
        valid_mask = series.notna() & (series.astype(str).str.strip() != "")
        valid_idx = np.where(valid_mask.to_numpy())[0]
        n = len(valid_idx)
        if n < 4:
            continue

        sub_dist = dist[np.ix_(valid_idx, valid_idx)]
        d2 = sub_dist * sub_dist
        j = np.eye(n) - np.ones((n, n)) / n
        b = -0.5 * j @ d2 @ j
        ss_total = float(np.trace(b))
        if ss_total <= 1e-12:
            continue

        sub_series = series.iloc[valid_idx]
        num = pd.to_numeric(sub_series, errors="coerce")
        if num.notna().all() and num.nunique() > 2:
            kind = "numeric"
            x_mat = num.to_numpy(dtype=float).reshape(-1, 1)
            x_mat = x_mat - x_mat.mean(axis=0, keepdims=True)
            df_model = 1
        else:
            kind = "categorical"
            cats = sub_series.astype(str)
            if cats.nunique() < 2 or cats.nunique() >= n:
                continue
            dummies = pd.get_dummies(cats, drop_first=True, dtype=float).to_numpy()
            x_mat = dummies - dummies.mean(axis=0, keepdims=True)
            df_model = x_mat.shape[1]

        df_resid = n - 1 - df_model
        if df_resid <= 0:
            continue

        def _calc_f(xm):
            q, _ = np.linalg.qr(xm)
            ss_m = float(np.trace(q.T @ b @ q))
            ss_r = ss_total - ss_m
            if ss_r <= 1e-15:
                return None, 0.0
            return (ss_m / df_model) / (ss_r / df_resid), ss_m / ss_total

        f_obs, r2 = _calc_f(x_mat)
        if f_obs is None:
            continue

        rng = np.random.default_rng(seed + 101 + idx_cov)
        ge = 0
        for _ in range(permutations):
            perm_idx = rng.permutation(n)
            f_perm, _ = _calc_f(x_mat[perm_idx])
            if f_perm is not None and f_perm >= f_obs:
                ge += 1
        results.append(
            {
                "covariate": col,
                "kind": kind,
                "n_samples": int(n),
                "df": int(df_model),
                "r2": float(r2),
                "pseudo_f": float(f_obs),
                "p_value": float((1 + ge) / (1 + permutations)),
                "permutations": int(permutations),
            }
        )
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--values", choices=["counts", "relative"], default="counts")
    ap.add_argument("--metadata")
    ap.add_argument("--group_col", default="group")
    ap.add_argument("--covariates", nargs="+", help="Optional metadata columns for marginal PERMANOVA (e.g. group sex age bmi).")
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--permutations", type=int, default=999)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    matrix = read_matrix(args.input)  # taxa x samples
    samples = matrix.columns.tolist()
    x = matrix.to_numpy().T  # samples x taxa
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = read_metadata(args.metadata, samples) if args.metadata else None
    labels = None
    if meta is not None:
        if args.group_col in meta.columns:
            labels = meta[args.group_col].astype(str).to_numpy()
        elif not args.covariates:
            fail(f"metadata has no column {args.group_col}; columns: {meta.columns.tolist()}")

    alpha = {s: alpha_indices(row, args.values) for s, row in zip(samples, x)}
    alpha_df = pd.DataFrame(alpha).T
    alpha_df.index.name = "sample"
    alpha_path = outdir / "alpha-diversity.tsv"
    alpha_df.to_csv(alpha_path, sep="\t")

    dist = squareform(pdist(x, metric="braycurtis"))
    dist_df = pd.DataFrame(dist, index=samples, columns=samples)
    dist_path = outdir / "bray-curtis.tsv"
    dist_df.to_csv(dist_path, sep="\t")

    coords, var = pcoa(dist)
    axes = [f"PCo{i + 1}" for i in range(coords.shape[1])]
    pcoa_df = pd.DataFrame(coords, index=samples, columns=axes)
    pcoa_path = outdir / "pcoa.tsv"
    pcoa_df.to_csv(pcoa_path, sep="\t")

    perm = permanova(dist, labels, args.permutations, args.seed) if labels is not None else None
    disp = permdisp(coords, labels, args.permutations, args.seed) if labels is not None else None

    cov_list = []
    if args.covariates:
        if meta is None:
            fail("--covariates requires --metadata")
        for item in args.covariates:
            for part in str(item).split(","):
                part = part.strip()
                if part and part not in cov_list:
                    cov_list.append(part)
    cov_results = covariate_permanova(dist, meta, cov_list, args.permutations, args.seed) if cov_list else []

    plt = setup_matplotlib()
    inputs = [args.input] + ([args.metadata] if args.metadata else [])
    params = {
        "values": args.values,
        "permutations": args.permutations,
        "seed": args.seed,
        "covariates": cov_list,
    }
    artifacts = [
        write_artifact(alpha_path, "mg_diversity", "Alpha diversity indices",
                       skill="profiling-workflow", inputs=inputs, parameters=params,
                       shape=alpha_df.shape),
        write_artifact(dist_path, "mg_diversity", "Bray-Curtis distance matrix",
                       skill="profiling-workflow", inputs=inputs, parameters=params,
                       shape=dist_df.shape),
        write_artifact(pcoa_path, "mg_diversity", "PCoA coordinates",
                       skill="profiling-workflow", inputs=inputs, parameters=params,
                       shape=pcoa_df.shape),
    ]

    if cov_results:
        cov_df = pd.DataFrame(cov_results)
        cov_path = outdir / "covariate-permanova.tsv"
        cov_df.to_csv(cov_path, sep="\t", index=False)
        artifacts.append(
            write_artifact(
                cov_path,
                "mg_diversity",
                "Multi-covariate marginal PERMANOVA table",
                skill="profiling-workflow",
                inputs=inputs,
                parameters=params,
                shape=cov_df.shape,
            )
        )

    fig, ax = plt.subplots(figsize=(max(6, len(samples) * 0.5), 4))
    indices = [c for c in alpha_df.columns if alpha_df[c].notna().any()]
    data = [alpha_df[c].dropna().to_numpy() for c in indices]
    ax.boxplot(data, tick_labels=indices)
    ax.set_ylabel("index value")
    ax.set_title(f"Alpha diversity ({args.values})")
    fig.tight_layout()
    p = outdir / "alpha-boxplot.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    artifacts.append(write_artifact(p, "mg_diversity", f"Alpha diversity boxplot ({args.values})",
                                    skill="profiling-workflow", inputs=inputs, parameters=params))

    if coords.shape[1] >= 2:
        fig, ax = plt.subplots(figsize=(6, 5))
        if labels is not None:
            for g in np.unique(labels):
                m = labels == g
                ax.scatter(coords[m, 0], coords[m, 1], label=str(g), s=50, alpha=0.8)
            ax.legend()
        else:
            ax.scatter(coords[:, 0], coords[:, 1], s=50, alpha=0.8)
        for i, s in enumerate(samples):
            ax.annotate(s, (coords[i, 0], coords[i, 1]), fontsize=7)
        ax.set_xlabel(f"PCo1 ({var[0] * 100:.1f}%)")
        ax.set_ylabel(f"PCo2 ({var[1] * 100:.1f}%)")
        title = "PCoA (Bray-Curtis)"
        if perm is not None:
            title += f" | PERMANOVA F={perm['pseudo_f']:.2f}, p={perm['p_value']:.3g}"
        ax.set_title(title)
        fig.tight_layout()
        p = outdir / "pcoa.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        artifacts.append(write_artifact(p, "mg_diversity", "PCoA ordination (Bray-Curtis)",
                                        skill="profiling-workflow", inputs=inputs, parameters=params))

    notes = []
    if args.values == "relative":
        notes.append("observed/chao1 skipped: relative abundances are not counts")
    if perm is not None and disp is not None and perm["p_value"] < 0.05 and disp["p_value"] < 0.05:
        notes.append(
            "PERMDISP is also significant (p < 0.05): group differences in PERMANOVA "
            "may partly reflect unequal within-group dispersion rather than centroid shift alone."
        )
    emit({
        "samples": samples,
        "values": args.values,
        "alpha": alpha,
        "beta_method": "braycurtis",
        "pcoa_variance_explained": var.tolist(),
        "permanova": perm,
        "permdisp": disp,
        "covariate_permanova": cov_results,
        "notes": notes,
        "artifacts": artifacts,
    })


if __name__ == "__main__":
    main()
