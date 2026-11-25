#!/usr/bin/env python3
"""WGCNA-style weighted gene co-expression module detection and trait association."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    bh_fdr,
    emit,
    fail,
    infer_value_kind,
    read_expression_matrix,
    read_metadata,
    setup_matplotlib,
    to_log2_cpm,
    write_artifact,
)


def _bool_arg(val) -> bool:
    return str(val).lower() not in ("false", "0", "no")


def _select_variable_genes(log_expr, top_genes: int):
    import numpy as np

    nonzero_frac = (log_expr > 0).mean(axis=1)
    keep = nonzero_frac >= 0.20
    if int(keep.sum()) < 30:
        keep = log_expr.var(axis=1) > 1e-8
    sub = log_expr.loc[keep]
    if sub.shape[0] < 30:
        fail(f"Too few variable genes ({sub.shape[0]}) for co-expression network construction (need >= 30).")

    arr = sub.to_numpy(dtype=float)
    med = np.median(arr, axis=1, keepdims=True)
    mad = np.median(np.abs(arr - med), axis=1)
    var = np.var(arr, axis=1)
    keep_var = var > 1e-8
    sub = sub.loc[keep_var]
    mad = mad[keep_var]
    var = var[keep_var]
    if sub.shape[0] < 30:
        fail(f"Too few non-constant genes ({sub.shape[0]}) for co-expression network construction.")

    n_keep = min(int(top_genes), sub.shape[0])
    order = np.lexsort((-var, -mad))[:n_keep]
    return sub.iloc[order].copy()


def _evaluate_scale_free(sim_mat, powers: list[int]) -> list[dict]:
    import numpy as np

    rows = []
    n = sim_mat.shape[0]
    for beta in powers:
        adj = np.power(sim_mat, beta)
        np.fill_diagonal(adj, 0.0)
        k = adj.sum(axis=1)
        mean_k = float(np.mean(k))
        median_k = float(np.median(k))
        pos_k = k[k > 1e-12]
        if pos_k.size < 20 or np.max(pos_k) <= np.min(pos_k) * 1.01:
            rows.append(
                {
                    "power": int(beta),
                    "sft_r2": 0.0,
                    "slope": 0.0,
                    "mean_k": round(mean_k, 4),
                    "median_k": round(median_k, 4),
                }
            )
            continue
        counts, edges = np.histogram(pos_k, bins=min(15, max(6, int(np.sqrt(n) // 2))))
        centers = 0.5 * (edges[:-1] + edges[1:])
        valid = (counts > 0) & (centers > 0)
        if int(valid.sum()) < 4:
            r2_signed = 0.0
            slope = 0.0
        else:
            x = np.log10(centers[valid])
            y = np.log10(counts[valid] / float(pos_k.size))
            x_mean = x.mean()
            y_mean = y.mean()
            denom = float(np.sum((x - x_mean) ** 2))
            if denom <= 1e-12:
                slope = 0.0
                r2_signed = 0.0
            else:
                slope = float(np.sum((x - x_mean) * (y - y_mean)) / denom)
                y_hat = y_mean + slope * (x - x_mean)
                ss_tot = float(np.sum((y - y_mean) ** 2))
                ss_res = float(np.sum((y - y_hat) ** 2))
                r2 = max(0.0, 1.0 - ss_res / ss_tot) if ss_tot > 1e-12 else 0.0
                r2_signed = float(-np.sign(slope) * r2) if slope != 0 else 0.0
        rows.append(
            {
                "power": int(beta),
                "sft_r2": round(r2_signed, 4),
                "slope": round(slope, 4),
                "mean_k": round(mean_k, 4),
                "median_k": round(median_k, 4),
            }
        )
    return rows


def _compute_eigengene(expr_sub):
    """Compute oriented 1st PC module eigengene across samples (length N, unit L2 norm)."""
    import numpy as np

    X = expr_sub.to_numpy(dtype=float)  # genes x samples
    means = X.mean(axis=1, keepdims=True)
    stds = X.std(axis=1, ddof=1, keepdims=True)
    stds[stds < 1e-12] = 1.0
    Z = (X - means) / stds  # genes x samples
    # SVD of samples x genes
    u, s, _ = np.linalg.svd(Z.T, full_matrices=False)
    me = u[:, 0]
    avg_expr = Z.mean(axis=0)
    if np.std(me) > 1e-12 and np.std(avg_expr) > 1e-12:
        if np.corrcoef(me, avg_expr)[0, 1] < 0:
            me = -me
    var_explained = float((s[0] ** 2) / max(1e-12, float(np.sum(s ** 2))))
    return me, var_explained


def _extract_numeric_traits(meta, trait_cols: list[str]):
    import numpy as np
    import pandas as pd

    cols = list(trait_cols) if trait_cols else list(meta.columns)
    traits: dict[str, np.ndarray] = {}
    for col in cols:
        if col not in meta.columns:
            fail(f"Trait column '{col}' not found in metadata columns: {list(meta.columns)}")
        s = meta[col]
        num = pd.to_numeric(s, errors="coerce")
        if int(num.notna().sum()) >= max(4, int(0.7 * len(s))) and num.nunique(dropna=True) >= 2:
            traits[col] = num.to_numpy(dtype=float)
            continue
        clean = s.astype(str).str.strip()
        clean = clean.mask(s.isna() | clean.isin(["", "NA", "nan", "None", "null"]))
        levels = [lv for lv in clean.dropna().unique().tolist()]
        if len(levels) == 2:
            pos_level = sorted(levels)[-1]
            arr = np.where(clean.isna(), np.nan, (clean == pos_level).astype(float))
            traits[f"{col}_{pos_level}"] = arr
        elif 3 <= len(levels) <= 5:
            for lv in sorted(levels):
                arr = np.where(clean.isna(), np.nan, (clean == lv).astype(float))
                traits[f"{col}_{lv}"] = arr
    if not traits:
        return pd.DataFrame(index=meta.index)
    return pd.DataFrame(traits, index=meta.index)


def main() -> None:
    ap = argparse.ArgumentParser(description="WGCNA-style weighted gene co-expression network analysis")
    ap.add_argument("--input", required=True, help="Expression matrix TSV/CSV (genes x samples)")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--metadata", default="", help="Optional sample metadata TSV/CSV")
    ap.add_argument("--sample_col", default="", help="Optional sample ID column in metadata")
    ap.add_argument("--trait_cols", action="extend", nargs="+", default=[], help="Optional trait columns in metadata")
    ap.add_argument("--top_genes", type=int, default=1500, help="Top variable genes by MAD to include (default 1500)")
    ap.add_argument("--network_type", choices=["signed", "unsigned"], default="signed")
    ap.add_argument("--power", type=int, default=0, help="Soft-threshold power beta (0 = auto scale-free selection)")
    ap.add_argument("--sft_r2_cutoff", type=float, default=0.80, help="Scale-free R^2 cutoff for auto beta selection")
    ap.add_argument("--min_module_size", type=int, default=25, help="Minimum genes per module (default 25)")
    ap.add_argument(
        "--merge_cut_height",
        type=float,
        default=0.20,
        help="Merge modules whose eigengene correlation exceeds 1 - merge_cut_height (default 0.20 -> r > 0.80)",
    )
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    from scipy import stats

    expr = read_expression_matrix(args.input)
    n_samples = expr.shape[1]
    if n_samples < 6:
        fail(
            f"Co-expression network analysis requires at least 6 samples (recommended >= 12); found {n_samples}."
        )

    warnings: list[str] = []
    if n_samples < 12:
        warnings.append(
            f"Sample size ({n_samples}) is below the WGCNA recommended minimum of 12-15 samples; interpret module correlations cautiously."
        )

    value_kind = infer_value_kind(expr)
    log_expr = to_log2_cpm(expr, value_kind)
    sub_expr = _select_variable_genes(log_expr, max(50, args.top_genes))
    genes = sub_expr.index.tolist()
    samples = sub_expr.columns.tolist()
    n_genes = len(genes)

    effective_min_mod = min(max(5, int(args.min_module_size)), max(5, n_genes // 4))
    if effective_min_mod != args.min_module_size:
        warnings.append(
            f"Adjusted min_module_size from {args.min_module_size} to {effective_min_mod} for {n_genes} variable genes."
        )

    # Pearson correlation across samples
    X = sub_expr.to_numpy(dtype=float)
    corr = np.corrcoef(X)
    corr = np.nan_to_num(corr, nan=0.0)
    np.clip(corr, -1.0, 1.0, out=corr)
    np.fill_diagonal(corr, 1.0)

    if args.network_type == "signed":
        sim = 0.5 + 0.5 * corr
    else:
        sim = np.abs(corr)
    np.clip(sim, 0.0, 1.0, out=sim)

    candidate_powers = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16]
    if args.power > 0 and args.power not in candidate_powers:
        candidate_powers = sorted(set(candidate_powers + [int(args.power)]))
    sft_rows = _evaluate_scale_free(sim, candidate_powers)
    sft_df = pd.DataFrame(sft_rows)

    if args.power > 0:
        selected_power = int(args.power)
        power_selection_mode = "user_specified"
    else:
        passing = [r for r in sft_rows if r["sft_r2"] >= args.sft_r2_cutoff]
        if passing:
            selected_power = int(passing[0]["power"])
            power_selection_mode = "scale_free_r2_threshold"
        else:
            best_row = max(sft_rows, key=lambda r: r["sft_r2"])
            if best_row["sft_r2"] >= 0.55:
                selected_power = int(best_row["power"])
                power_selection_mode = "max_sft_r2_fallback"
            else:
                selected_power = 12 if args.network_type == "signed" else 6
                power_selection_mode = "wgcna_empirical_fallback"
            warnings.append(
                f"No candidate power reached scale-free R^2 >= {args.sft_r2_cutoff:.2f} "
                f"(max R^2 = {best_row['sft_r2']:.3f}); selected power beta = {selected_power} ({power_selection_mode})."
            )

    selected_sft = next((r for r in sft_rows if r["power"] == selected_power), {"sft_r2": 0.0, "mean_k": 0.0})

    # Adjacency and Topological Overlap Matrix (TOM)
    adj = np.power(sim, selected_power)
    np.fill_diagonal(adj, 0.0)
    k_total = adj.sum(axis=1)

    num_tom = adj @ adj.T + adj
    k_min = np.minimum.outer(k_total, k_total)
    denom_tom = k_min + 1.0 - adj
    tom = np.divide(num_tom, np.maximum(denom_tom, 1e-12))
    np.clip(tom, 0.0, 1.0, out=tom)
    np.fill_diagonal(tom, 1.0)
    diss_tom = 1.0 - tom
    np.fill_diagonal(diss_tom, 0.0)
    diss_tom = 0.5 * (diss_tom + diss_tom.T)

    condensed = squareform(diss_tom, checks=False)
    Z_link = linkage(condensed, method="average")

    # Adaptive cut height selection maximizing valid modules >= effective_min_mod (prefer 2..12 modules)
    heights = Z_link[:, 2]
    quantiles = np.linspace(0.35, 0.96, 25)
    candidate_heights = np.unique(np.quantile(heights, quantiles))
    best_labels = None
    best_score = (-1, -1)  # (n_valid_modules_capped, assigned_genes)
    for h in candidate_heights:
        raw_cl = fcluster(Z_link, t=float(h), criterion="distance")
        vc = pd.Series(raw_cl).value_counts()
        valid_cl = vc[vc >= effective_min_mod]
        n_valid = int(len(valid_cl))
        if n_valid == 0 or n_valid > 18:
            continue
        assigned = int(valid_cl.sum())
        capped = n_valid if n_valid <= 12 else max(1, 24 - n_valid)
        score = (capped, assigned)
        if score > best_score:
            best_score = score
            best_labels = raw_cl

    if best_labels is None:
        # Fallback to maxclust=4
        best_labels = fcluster(Z_link, t=min(4, max(2, n_genes // effective_min_mod)), criterion="maxclust")

    # Map small clusters to 0 (unassigned)
    vc = pd.Series(best_labels).value_counts()
    valid_ids = set(vc[vc >= effective_min_mod].index.tolist())
    mod_assign = np.array([int(c) if int(c) in valid_ids else 0 for c in best_labels], dtype=int)

    # If everything fell into 0, force top 2 clusters
    if np.all(mod_assign == 0):
        raw_cl = fcluster(Z_link, t=2, criterion="maxclust")
        mod_assign = raw_cl.astype(int)

    # Iteratively merge modules whose eigengene correlation > 1 - merge_cut_height
    r_merge_thresh = 1.0 - float(args.merge_cut_height)
    while True:
        uniq = sorted([m for m in set(mod_assign.tolist()) if m > 0])
        if len(uniq) <= 1:
            break
        mes = {}
        for m in uniq:
            me, _ = _compute_eigengene(sub_expr.iloc[mod_assign == m])
            mes[m] = me
        me_mat = np.column_stack([mes[m] for m in uniq])
        me_corr = np.corrcoef(me_mat, rowvar=False)
        np.fill_diagonal(me_corr, -1.0)
        max_idx = np.unravel_index(int(np.argmax(me_corr)), me_corr.shape)
        max_r = float(me_corr[max_idx])
        if max_r <= r_merge_thresh:
            break
        m_keep = uniq[min(max_idx[0], max_idx[1])]
        m_drop = uniq[max(max_idx[0], max_idx[1])]
        mod_assign[mod_assign == m_drop] = m_keep

    # Relabel modules M1, M2, ... by descending gene count
    pos_counts = pd.Series(mod_assign[mod_assign > 0]).value_counts()
    ordered_raw_mods = pos_counts.index.tolist()
    mod_name_map = {int(raw_m): f"M{idx + 1}" for idx, raw_m in enumerate(ordered_raw_mods)}
    mod_name_map[0] = "M0_unassigned"
    gene_modules = [mod_name_map[int(m)] for m in mod_assign]

    # Compute final Module Eigengenes, variance explained, kME, and hub genes
    module_names = [f"M{idx + 1}" for idx in range(len(ordered_raw_mods))]
    me_dict: dict[str, np.ndarray] = {}
    me_var_exp: dict[str, float] = {}
    for m_name in module_names:
        mask = np.array([gm == m_name for gm in gene_modules], dtype=bool)
        me, var_exp = _compute_eigengene(sub_expr.iloc[mask])
        me_dict[m_name] = me
        me_var_exp[m_name] = var_exp

    me_df = pd.DataFrame(
        [me_dict[m] for m in module_names],
        index=module_names,
        columns=samples,
    )
    me_df.index.name = "module"

    # Gene-level hub metrics: kTotal, kWithin, and kME against assigned module
    gene_rows = []
    module_summaries = []
    for idx, g in enumerate(genes):
        m_name = gene_modules[idx]
        if m_name == "M0_unassigned":
            k_within = 0.0
            k_me = 0.0
        else:
            same_mod = np.array([gm == m_name for gm in gene_modules], dtype=bool)
            same_mod[idx] = False
            k_within = float(adj[idx, same_mod].sum()) if same_mod.any() else 0.0
            me_vec = me_dict[m_name]
            gx = X[idx, :]
            if np.std(gx) > 1e-12 and np.std(me_vec) > 1e-12:
                k_me = float(np.corrcoef(gx, me_vec)[0, 1])
            else:
                k_me = 0.0
        gene_rows.append(
            {
                "gene": g,
                "module": m_name,
                "kME": round(k_me, 4),
                "kWithin": round(k_within, 4),
                "kTotal": round(float(k_total[idx]), 4),
            }
        )

    gene_df = pd.DataFrame(gene_rows)
    gene_df["is_unassigned"] = gene_df["module"] == "M0_unassigned"
    gene_df = gene_df.sort_values(
        by=["is_unassigned", "module", "kME", "kWithin"],
        ascending=[True, True, False, False],
    ).drop(columns=["is_unassigned"])

    for m_name in module_names:
        sub_g = gene_df[gene_df["module"] == m_name]
        top_hubs = sub_g.head(5)["gene"].tolist()
        module_summaries.append(
            {
                "module": m_name,
                "n_genes": int(len(sub_g)),
                "variance_explained": round(me_var_exp[m_name], 4),
                "top_hub_genes": top_hubs,
            }
        )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sft_path = out_dir / "scale_free_fit.tsv"
    me_path = out_dir / "module_eigengenes.tsv"
    hubs_path = out_dir / "gene_modules_hubs.tsv"
    fig_path = out_dir / "coexpression-modules.png"

    sft_df.to_csv(sft_path, sep="\t", index=False)
    me_df.round(5).to_csv(me_path, sep="\t")
    gene_df.to_csv(hubs_path, sep="\t", index=False)

    write_artifact(
        sft_path,
        "rna_coexpression",
        "WGCNA-Style Scale-Free Topology Fit Across Soft-Threshold Powers",
        inputs=[args.input],
        parameters={"network_type": args.network_type, "selected_power": selected_power},
        shape={"rows": int(len(sft_df)), "cols": int(sft_df.shape[1])},
    )
    write_artifact(
        me_path,
        "rna_coexpression",
        "Co-expression Module Eigengene Matrix",
        inputs=[args.input],
        parameters={"selected_power": selected_power, "n_modules": len(module_names)},
        shape={"rows": int(me_df.shape[0]), "cols": int(me_df.shape[1])},
        scientific_context={
            "dataType": "pathway_scores",
            "valueScale": "normalized",
            "normalization": "module_eigengene_pc1",
            "entityIdType": "module_id",
            "organelle": "nuclear",
        },
    )
    write_artifact(
        hubs_path,
        "rna_coexpression",
        "Gene Co-expression Module Assignments & Intramodular Hub Connectivity (kME)",
        inputs=[args.input],
        parameters={"selected_power": selected_power, "min_module_size": effective_min_mod},
        shape={"rows": int(len(gene_df)), "cols": int(gene_df.shape[1])},
    )

    # Optional module-trait association analysis
    trait_assoc_df = pd.DataFrame()
    trait_assoc_path_str = None
    top_associations = []
    if args.metadata:
        meta = read_metadata(args.metadata, samples=samples, sample_col=args.sample_col or None)
        traits_df = _extract_numeric_traits(meta, args.trait_cols)
        if traits_df.shape[1] == 0:
            warnings.append("Metadata provided, but no numeric or binary/low-cardinality trait columns could be encoded.")
        else:
            assoc_rows = []
            for m_name in module_names:
                me_vec = me_dict[m_name]
                for t_col in traits_df.columns:
                    t_vec = traits_df[t_col].to_numpy(dtype=float)
                    valid = np.isfinite(t_vec) & np.isfinite(me_vec)
                    n_eff = int(valid.sum())
                    if n_eff < 5 or np.std(t_vec[valid]) < 1e-12 or np.std(me_vec[valid]) < 1e-12:
                        r_val, p_val = 0.0, 1.0
                    else:
                        r_res = stats.pearsonr(me_vec[valid], t_vec[valid])
                        r_val = float(r_res.statistic)
                        p_val = float(r_res.pvalue)
                    assoc_rows.append(
                        {
                            "module": m_name,
                            "trait": t_col,
                            "n_samples": n_eff,
                            "cor": round(r_val, 4),
                            "pvalue": p_val,
                        }
                    )
            trait_assoc_df = pd.DataFrame(assoc_rows)
            trait_assoc_df["padj"] = bh_fdr(trait_assoc_df["pvalue"].to_numpy(dtype=float))
            trait_assoc_df = trait_assoc_df.sort_values("pvalue", ascending=True).reset_index(drop=True)
            trait_path = out_dir / "module_trait_associations.tsv"
            trait_assoc_df.to_csv(trait_path, sep="\t", index=False)
            trait_assoc_path_str = str(trait_path)
            write_artifact(
                trait_path,
                "rna_coexpression",
                "Module Eigengene-Trait Pearson Correlations & BH FDR",
                inputs=[args.input, args.metadata],
                parameters={"selected_power": selected_power},
                shape={"rows": int(len(trait_assoc_df)), "cols": int(trait_assoc_df.shape[1])},
            )
            top_associations = [
                {
                    "module": r["module"],
                    "trait": r["trait"],
                    "cor": float(r["cor"]),
                    "pvalue": float(r["pvalue"]),
                    "padj": float(r["padj"]),
                }
                for r in trait_assoc_df.head(8).to_dict(orient="records")
            ]

    # Plot 2-panel summary figure
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2), dpi=180)

    ax1.plot(sft_df["power"], sft_df["sft_r2"], marker="o", color="#2563eb", lw=1.8, label="Signed scale-free $R^2$")
    ax1.axhline(args.sft_r2_cutoff, color="#dc2626", ls="--", lw=1.1, label=f"Cutoff ({args.sft_r2_cutoff:.2f})")
    ax1.axvline(selected_power, color="#059669", ls=":", lw=1.5, label=f"Selected $\\beta={selected_power}$")
    ax1.set_xlabel("Soft-Threshold Power ($\\beta$)")
    ax1.set_ylabel("Scale-Free Topology Fit ($R^2$)")
    ax1.set_ylim(-0.2, 1.05)
    ax1.set_title("Scale-Free Topology Model Fit")
    ax1.legend(loc="lower right", frameon=False, fontsize=8.5)

    if not trait_assoc_df.empty:
        piv_r = trait_assoc_df.pivot(index="module", columns="trait", values="cor").reindex(module_names)
        piv_q = trait_assoc_df.pivot(index="module", columns="trait", values="padj").reindex(module_names)
        mat = piv_r.to_numpy(dtype=float)
        im = ax2.imshow(mat, cmap="RdBu_r", vmin=-1.0, vmax=1.0, aspect="auto")
        ax2.set_xticks(range(piv_r.shape[1]))
        ax2.set_xticklabels(piv_r.columns.tolist(), rotation=30, ha="right", fontsize=8.5)
        ax2.set_yticks(range(piv_r.shape[0]))
        ax2.set_yticklabels(
            [f"{m} (n={int((gene_df['module'] == m).sum())})" for m in piv_r.index.tolist()],
            fontsize=8.5,
        )
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                rv = mat[i, j]
                qv = float(piv_q.iloc[i, j])
                sig = "**" if qv < 0.01 else ("*" if qv < 0.05 else "")
                txt_col = "white" if abs(rv) > 0.55 else "#111827"
                ax2.text(j, i, f"{rv:.2f}{sig}\n({qv:.2g})", ha="center", va="center", fontsize=7.5, color=txt_col)
        fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04, label="Pearson r")
        ax2.set_title("Module-Trait Relationships (* FDR < 0.05)")
    else:
        counts_s = [m["n_genes"] for m in module_summaries]
        labels_s = [f"{m['module']} ({', '.join(m['top_hub_genes'][:2])})" for m in module_summaries]
        y_pos = np.arange(len(module_summaries))
        ax2.barh(y_pos, counts_s, color="#3b82f6", alpha=0.88)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(labels_s, fontsize=8.5)
        ax2.invert_yaxis()
        ax2.set_xlabel("Genes in Module")
        ax2.set_title("Detected Co-expression Modules & Top Hubs")

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "rna_coexpression",
        "Co-expression Scale-Free Fit & Module Summary Plot",
        inputs=[args.input] + ([args.metadata] if args.metadata else []),
        parameters={"selected_power": selected_power, "n_modules": len(module_names)},
    )

    n_unassigned = int((gene_df["module"] == "M0_unassigned").sum())
    emit(
        {
            "status": "ok",
            "tool": "rna_coexpression",
            "inferred_value_kind": value_kind,
            "n_samples": int(n_samples),
            "n_genes_input": int(expr.shape[0]),
            "n_genes_network": int(n_genes),
            "n_genes_assigned": int(n_genes - n_unassigned),
            "n_genes_unassigned": n_unassigned,
            "network_type": args.network_type,
            "selected_power": int(selected_power),
            "power_selection_mode": power_selection_mode,
            "scale_free_r2": float(selected_sft.get("sft_r2", 0.0)),
            "n_modules": int(len(module_names)),
            "modules": module_summaries,
            "top_trait_associations": top_associations,
            "outputs": {
                "scale_free_fit": str(sft_path),
                "module_eigengenes": str(me_path),
                "gene_modules_hubs": str(hubs_path),
                "module_trait_associations": trait_assoc_path_str,
                "summary_plot": str(fig_path),
            },
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
