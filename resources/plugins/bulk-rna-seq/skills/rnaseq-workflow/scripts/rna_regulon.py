#!/usr/bin/env python3
"""rna_regulon: decoupler-style ULM/MLM transcription factor and pathway regulon activity inference."""
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
    setup_matplotlib,
    to_log2_cpm,
    write_artifact,
)

# Curated 14-pathway / key-TF human & mouse canonical signed signaling footprint resource
# (uppercase gene symbols for case-insensitive HGNC/MGI matching; +1.0 activation, -1.0 repression)
BUILTIN_REGULON_NETWORK: dict[str, list[tuple[str, float]]] = {
    "NFkB": [
        ("NFKBIA", 1.0), ("TNFAIP3", 1.0), ("CCL2", 1.0), ("CXCL8", 1.0),
        ("IL6", 1.0), ("BCL2A1", 1.0), ("BIRC3", 1.0), ("ICAM1", 1.0),
        ("VCAM1", 1.0), ("PTGS2", 1.0), ("RELB", 1.0), ("NFKB2", 1.0),
        ("TRAF1", 1.0), ("CCL20", 1.0), ("CXCL1", 1.0), ("CXCL2", 1.0),
        ("SOD2", 1.0), ("IRF1", 1.0),
    ],
    "JAK_STAT": [
        ("STAT1", 1.0), ("STAT2", 1.0), ("IRF1", 1.0), ("IRF9", 1.0),
        ("SOCS1", 1.0), ("SOCS3", 1.0), ("GBP1", 1.0), ("GBP2", 1.0),
        ("CXCL9", 1.0), ("CXCL10", 1.0), ("ISG15", 1.0), ("MX1", 1.0),
        ("OAS1", 1.0), ("IFIT1", 1.0), ("IFIT3", 1.0), ("IDO1", 1.0),
        ("CISH", 1.0), ("PIM1", 1.0),
    ],
    "TGFb": [
        ("SERPINE1", 1.0), ("SMAD7", 1.0), ("CTGF", 1.0), ("CCN2", 1.0),
        ("COL1A1", 1.0), ("COL1A2", 1.0), ("COL3A1", 1.0), ("FN1", 1.0),
        ("ACTA2", 1.0), ("SNAI1", 1.0), ("SNAI2", 1.0), ("CDKN2B", 1.0),
        ("TGFBI", 1.0), ("LTBP2", 1.0), ("PMEPA1", 1.0), ("SKIL", 1.0),
        ("BAMBI", 1.0), ("ID1", -1.0),
    ],
    "TNFa": [
        ("TNF", 1.0), ("TNFAIP3", 1.0), ("NFKBIA", 1.0), ("CCL2", 1.0),
        ("CXCL10", 1.0), ("IL1B", 1.0), ("IL6", 1.0), ("JUNB", 1.0),
        ("FOS", 1.0), ("IER3", 1.0), ("ZC3H12A", 1.0), ("ATF3", 1.0),
        ("GADD45B", 1.0), ("BIRC3", 1.0), ("ICAM1", 1.0), ("MMP9", 1.0),
        ("PLAU", 1.0),
    ],
    "Hypoxia": [
        ("VEGFA", 1.0), ("CA9", 1.0), ("SLC2A1", 1.0), ("LDHA", 1.0),
        ("PDK1", 1.0), ("BNIP3", 1.0), ("BNIP3L", 1.0), ("HK2", 1.0),
        ("ENO1", 1.0), ("ALDOA", 1.0), ("PGK1", 1.0), ("PFKFB3", 1.0),
        ("NDRG1", 1.0), ("LOX", 1.0), ("EGLN3", 1.0), ("ADM", 1.0),
        ("ANGPTL4", 1.0), ("MXI1", 1.0),
    ],
    "p53": [
        ("CDKN1A", 1.0), ("MDM2", 1.0), ("BAX", 1.0), ("BBC3", 1.0),
        ("GADD45A", 1.0), ("FAS", 1.0), ("TP53I3", 1.0), ("DDB2", 1.0),
        ("RRM2B", 1.0), ("SESN1", 1.0), ("SESN2", 1.0), ("ZMAT3", 1.0),
        ("CCNG1", 1.0), ("BTG2", 1.0), ("TNFRSF10B", 1.0), ("BIRC5", -1.0),
        ("CCNB1", -1.0), ("CDK1", -1.0),
    ],
    "MAPK": [
        ("DUSP4", 1.0), ("DUSP6", 1.0), ("SPRY2", 1.0), ("SPRY4", 1.0),
        ("ETV4", 1.0), ("ETV5", 1.0), ("FOS", 1.0), ("FOSL1", 1.0),
        ("JUN", 1.0), ("EGR1", 1.0), ("IER3", 1.0), ("MYC", 1.0),
        ("CCND1", 1.0), ("PHLDA1", 1.0), ("TRIB1", 1.0), ("MAFF", 1.0),
        ("BHLHE40", 1.0),
    ],
    "PI3K": [
        ("AKT1", 1.0), ("PIK3CA", 1.0), ("CCND1", 1.0), ("SLC2A1", 1.0),
        ("HK2", 1.0), ("SREBF1", 1.0), ("FASN", 1.0), ("ACLY", 1.0),
        ("VEGFA", 1.0), ("MYC", 1.0), ("RPS6KB1", 1.0), ("EIF4EBP1", 1.0),
        ("FOXO3", -1.0), ("CDKN1B", -1.0), ("BCL2L11", -1.0), ("TXNIP", -1.0),
        ("DDIT4", -1.0),
    ],
    "WNT": [
        ("AXIN2", 1.0), ("LEF1", 1.0), ("TCF7", 1.0), ("MYC", 1.0),
        ("CCND1", 1.0), ("LGR5", 1.0), ("NKD1", 1.0), ("NOTUM", 1.0),
        ("ASCL2", 1.0), ("SP5", 1.0), ("RNF43", 1.0), ("ZNRF3", 1.0),
        ("DKK1", 1.0), ("MMP7", 1.0), ("BIRC5", 1.0), ("CDX2", 1.0),
        ("SFRP1", -1.0),
    ],
    "Estrogen": [
        ("ESR1", 1.0), ("PGR", 1.0), ("GREB1", 1.0), ("TFF1", 1.0),
        ("XBP1", 1.0), ("FOXA1", 1.0), ("GATA3", 1.0), ("CCND1", 1.0),
        ("MYB", 1.0), ("CA12", 1.0), ("STC2", 1.0), ("CXCL12", 1.0),
        ("IGFBP4", 1.0), ("RET", 1.0), ("BCL2", 1.0), ("NRIP1", 1.0),
        ("AREG", 1.0),
    ],
    "Androgen": [
        ("KLK3", 1.0), ("KLK2", 1.0), ("TMPRSS2", 1.0), ("NKX3-1", 1.0),
        ("FKBP5", 1.0), ("STEAP4", 1.0), ("SLC45A3", 1.0), ("ACSL3", 1.0),
        ("ELOVL5", 1.0), ("FASN", 1.0), ("CAMKK2", 1.0), ("CENPN", 1.0),
        ("ZBTB16", 1.0), ("PMEPA1", 1.0), ("NDRG1", 1.0), ("DDC", -1.0),
    ],
    "VEGF": [
        ("VEGFA", 1.0), ("KDR", 1.0), ("FLT1", 1.0), ("DLL4", 1.0),
        ("ANGPT2", 1.0), ("ESM1", 1.0), ("APLN", 1.0), ("PGF", 1.0),
        ("NRP1", 1.0), ("PECAM1", 1.0), ("CDH5", 1.0), ("TIE1", 1.0),
        ("TEK", 1.0), ("ROBO4", 1.0), ("PLVAP", 1.0), ("RGS5", 1.0),
    ],
    "Trail": [
        ("TNFSF10", 1.0), ("TNFRSF10A", 1.0), ("TNFRSF10B", 1.0), ("CASP8", 1.0),
        ("CASP3", 1.0), ("CASP10", 1.0), ("FADD", 1.0), ("BID", 1.0),
        ("BAX", 1.0), ("BAK1", 1.0), ("PMAIP1", 1.0), ("BBC3", 1.0),
        ("CFLAR", -1.0), ("XIAP", -1.0), ("BIRC2", -1.0), ("MCL1", -1.0),
    ],
    "EGFR": [
        ("EGFR", 1.0), ("EREG", 1.0), ("AREG", 1.0), ("TGFA", 1.0),
        ("HBEGF", 1.0), ("DUSP6", 1.0), ("DUSP4", 1.0), ("SPRY2", 1.0),
        ("FOSL1", 1.0), ("MYC", 1.0), ("CCND1", 1.0), ("MMP1", 1.0),
        ("PLAU", 1.0), ("PHLDA1", 1.0), ("ERRFI1", 1.0), ("ETV5", 1.0),
        ("EPHA2", 1.0),
    ],
}


def _load_network(path: str):
    import pandas as pd

    if not path:
        rows = [
            {"source": src, "target": tgt.upper(), "weight": float(w)}
            for src, pairs in BUILTIN_REGULON_NETWORK.items()
            for tgt, w in pairs
        ]
        return pd.DataFrame(rows), "builtin_canonical_14_signaling_regulons"

    p = Path(path)
    if not p.is_file():
        fail(f"Regulon network file not found: {p}")
    sep = "," if p.name.lower().endswith(".csv") else "\t"
    raw = pd.read_csv(p, sep=sep)
    if raw.empty:
        fail(f"Regulon network file is empty: {p}")

    cols_lower = {c.lower().strip(): c for c in raw.columns}
    src_col = next(
        (cols_lower[c] for c in ("source", "tf", "pathway", "regulon", "regulator") if c in cols_lower),
        None,
    )
    tgt_col = next(
        (cols_lower[c] for c in ("target", "gene", "gene_symbol", "genesymbol", "symbol") if c in cols_lower),
        None,
    )
    w_col = next(
        (cols_lower[c] for c in ("mor", "weight", "mode", "sign", "score") if c in cols_lower),
        None,
    )
    if not src_col or not tgt_col:
        fail(
            f"Regulon network {p.name} requires 'source' and 'target' columns; found {raw.columns.tolist()}"
        )

    net = pd.DataFrame()
    net["source"] = raw[src_col].astype(str).str.strip()
    net["target"] = raw[tgt_col].astype(str).str.strip().str.upper()
    if w_col:
        net["weight"] = pd.to_numeric(raw[w_col], errors="coerce").fillna(1.0)
    else:
        net["weight"] = 1.0
    net = net[(net["source"] != "") & (net["target"] != "") & (net["weight"] != 0.0)]
    net = net.groupby(["source", "target"], as_index=False)["weight"].mean()
    if net.empty:
        fail(f"No valid (source, target, weight) edges parsed from {p}")
    return net, f"custom:{p.name}"


def _load_input_profiles(path: str, stat_col_req: str):
    """Load expression matrix (genes x samples) or DE result table (gene + stat/log2FC).
    Returns (input_mode, df_profiles_upper, warnings).
    """
    import numpy as np
    import pandas as pd

    p = Path(path)
    if not p.is_file():
        fail(f"Input file not found: {p}")
    sep = "," if p.name.lower().endswith(".csv") else "\t"
    raw = pd.read_csv(p, sep=sep)
    if raw.empty or raw.shape[1] < 2:
        fail(f"Input file must have a gene column and at least one numeric column: {p}")

    cols_lower = {c.lower().strip(): c for c in raw.columns}
    de_markers = {
        "log2foldchange", "log2fc", "logfc", "lfc",
        "pvalue", "pval", "p.value", "p_val",
        "padj", "fdr", "adj.p.val", "qvalue", "q_value",
        "stat", "t", "wald", "lr",
    }
    is_de = bool(stat_col_req) or any(k in cols_lower for k in de_markers)
    warnings: list[str] = []

    if is_de:
        gene_col = next(
            (cols_lower[c] for c in ("gene", "gene_id", "geneid", "symbol", "gene_name", "id", "unnamed: 0") if c in cols_lower),
            raw.columns[0],
        )
        chosen_col = None
        if stat_col_req:
            if stat_col_req in raw.columns:
                chosen_col = stat_col_req
            elif stat_col_req.lower().strip() in cols_lower:
                chosen_col = cols_lower[stat_col_req.lower().strip()]
            else:
                fail(f"Requested stat_col '{stat_col_req}' not found in columns: {raw.columns.tolist()}")
        else:
            for cand in ("stat", "t", "wald", "log2foldchange", "log2fc", "logfc", "lfc"):
                if cand in cols_lower:
                    chosen_col = cols_lower[cand]
                    break
        if chosen_col is None:
            p_col = next((cols_lower[c] for c in ("pvalue", "pval", "p.value", "p_val", "padj", "fdr") if c in cols_lower), None)
            fc_col = next((cols_lower[c] for c in ("log2foldchange", "log2fc", "logfc", "lfc") if c in cols_lower), None)
            if p_col and fc_col:
                pvals = np.clip(pd.to_numeric(raw[p_col], errors="coerce").fillna(1.0).to_numpy(dtype=float, copy=True), 1e-300, 1.0)
                fcs = pd.to_numeric(raw[fc_col], errors="coerce").fillna(0.0).to_numpy(dtype=float, copy=True)
                stat_vals = np.sign(fcs) * (-np.log10(pvals))
                chosen_label = "signed_neglog10p"
            else:
                fail(f"Could not detect statistic or log2FC column in DE table: {raw.columns.tolist()}")
        else:
            stat_vals = pd.to_numeric(raw[chosen_col], errors="coerce").to_numpy(dtype=float, copy=True)
            chosen_label = chosen_col

        genes_up = raw[gene_col].astype(str).str.strip().str.upper()
        de_df = pd.DataFrame({"gene": genes_up, chosen_label: stat_vals})
        de_df = de_df.dropna(subset=[chosen_label])
        de_df = de_df[de_df["gene"] != ""]
        if de_df.index.size == 0:
            fail("No finite numeric values found in DE ranking statistic column.")
        if de_df["gene"].duplicated().any():
            de_df["_abs"] = de_df[chosen_label].abs()
            de_df = de_df.sort_values("_abs", ascending=False).drop_duplicates(subset=["gene"], keep="first").drop(columns=["_abs"])
        prof = de_df.set_index("gene")[[chosen_label]]
        return "de", prof, warnings

    # Matrix mode
    idx_col = raw.columns[0]
    mat = raw.set_index(idx_col)
    mat.index = mat.index.astype(str).str.strip().str.upper()
    mat.columns = mat.columns.astype(str).str.strip()
    num = mat.apply(pd.to_numeric, errors="coerce")
    if num.isna().any().any():
        bad = num.columns[num.isna().any()].tolist()
        fail(f"Non-numeric values in expression matrix columns: {bad}")
    if num.index.has_duplicates:
        num["_mean"] = num.mean(axis=1)
        num = num.sort_values("_mean", ascending=False)
        num = num[~num.index.duplicated(keep="first")].drop(columns=["_mean"])

    vkind = infer_value_kind(num)
    log_df = to_log2_cpm(num, vkind)
    X = log_df.to_numpy(dtype=float, copy=True)
    if log_df.shape[1] >= 2:
        mu = X.mean(axis=1, keepdims=True)
        sd = X.std(axis=1, ddof=1, keepdims=True)
        sd = np.where(sd > 1e-12, sd, 1.0)
        Z = (X - mu) / sd
    else:
        Z = X - X.mean(axis=0, keepdims=True)
        warnings.append("Single-sample expression profile provided; skipped cohort gene-wise z-scoring.")
    prof = pd.DataFrame(Z, index=log_df.index, columns=log_df.columns)
    return "matrix", prof, warnings


def _run_ulm(X_mat, Y_mat):
    """Univariate Linear Model: fit each regulon column independently across G shared genes."""
    import numpy as np
    from scipy import stats

    G, P = X_mat.shape
    df = G - 2
    if df < 2:
        fail(f"Need at least 4 shared genes for ULM regression (found {G}).")

    Y_c = Y_mat - Y_mat.mean(axis=0, keepdims=True)
    t_mat = np.zeros((P, Y_mat.shape[1]), dtype=float)
    p_mat = np.ones((P, Y_mat.shape[1]), dtype=float)

    for p_idx in range(P):
        x = X_mat[:, p_idx]
        x_c = x - x.mean()
        ssx = float(np.dot(x_c, x_c))
        if ssx <= 1e-12:
            continue
        beta = (x_c @ Y_c) / ssx  # shape (S,)
        resid = Y_c - np.outer(x_c, beta)
        rss = np.sum(resid ** 2, axis=0)
        se = np.sqrt(np.maximum(rss / float(df) / ssx, 1e-24))
        t_vals = np.where(rss <= 1e-24, 0.0, beta / se)
        p_vals = 2.0 * stats.t.sf(np.abs(t_vals), df=df)
        t_mat[p_idx, :] = t_vals
        p_mat[p_idx, :] = np.clip(np.where(np.isfinite(p_vals), p_vals, 1.0), 0.0, 1.0)

    return t_mat, p_mat


def _run_mlm(X_mat, Y_mat):
    """Ridge-stabilized Multivariate Linear Model across all P regulons simultaneously."""
    import numpy as np
    from scipy import stats

    G, P = X_mat.shape
    df = max(2, G - P - 1)
    X_c = X_mat - X_mat.mean(axis=0, keepdims=True)
    Y_c = Y_mat - Y_mat.mean(axis=0, keepdims=True)

    XtX = X_c.T @ X_c
    trace_avg = float(np.trace(XtX)) / max(1, P)
    cond = float(np.linalg.cond(XtX)) if P > 1 else 1.0
    ridge_scale = 1e-3 if cond > 1e6 or G <= P + 2 else 1e-5
    lam = ridge_scale * max(trace_avg, 1e-6)
    A = XtX + lam * np.eye(P, dtype=float)
    A_inv = np.linalg.inv(A)

    B_hat = A_inv @ (X_c.T @ Y_c)  # (P, S)
    resid = Y_c - X_c @ B_hat
    rss = np.sum(resid ** 2, axis=0)  # (S,)
    sigma2 = rss / float(df)

    cov_diag = np.maximum(np.diag(A_inv @ XtX @ A_inv), 1e-24)  # (P,)
    se_mat = np.sqrt(np.outer(cov_diag, sigma2))
    t_mat = np.where(se_mat > 1e-24, B_hat / se_mat, 0.0)
    p_mat = 2.0 * stats.t.sf(np.abs(t_mat), df=df)
    p_mat = np.clip(np.where(np.isfinite(p_mat), p_mat, 1.0), 0.0, 1.0)
    return t_mat, p_mat


def main() -> None:
    ap = argparse.ArgumentParser(description="ULM/MLM transcription factor and pathway regulon activity inference")
    ap.add_argument("--input", required=True, help="Expression matrix (genes x samples) or DE table TSV/CSV")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--network", default="", help="Optional signed regulon TSV/CSV (source, target, optional mor/weight)")
    ap.add_argument("--method", choices=["ulm", "mlm"], default="ulm", help="Linear model method (default ulm)")
    ap.add_argument("--min_targets", type=int, default=5, help="Minimum matched target genes per regulon (default 5)")
    ap.add_argument("--stat_col", default="", help="Optional DE statistic column name")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd

    net_df, net_source = _load_network(args.network)
    input_mode, prof_df, warnings = _load_input_profiles(args.input, args.stat_col)

    sample_genes = [str(g) for g in prof_df.index[: min(200, len(prof_df.index))]]
    ens_frac = sum(1 for g in sample_genes if g.startswith(("ENSG0", "ENSMUSG0", "ENSRNOG0"))) / max(1, len(sample_genes))
    if ens_frac > 0.5 and not args.network:
        fail(
            "Input appears to use Ensembl IDs; convert to HGNC/MGI gene symbols with rna_id_map first "
            "or pass a matching --network."
        )

    prof_genes = set(prof_df.index.tolist())
    net_sub = net_df[net_df["target"].isin(prof_genes)].copy()
    counts = net_sub.groupby("source")["target"].nunique()
    min_t = max(3, int(args.min_targets))
    valid_sources = sorted(counts[counts >= min_t].index.tolist())
    if not valid_sources:
        fail(
            f"No regulons retained >= {min_t} overlapping targets with input genes "
            f"(max overlap={int(counts.max()) if not counts.empty else 0})."
        )

    net_sub = net_sub[net_sub["source"].isin(valid_sources)].copy()
    piv = net_sub.pivot_table(index="target", columns="source", values="weight", aggfunc="mean", fill_value=0.0)
    # Retain columns with non-zero variance across shared targets
    shared_genes = sorted(piv.index.tolist())
    piv = piv.loc[shared_genes]
    keep_srcs = [s for s in piv.columns if float(np.std(piv[s].to_numpy(dtype=float, copy=True))) > 1e-12]
    if not keep_srcs:
        fail("All retained regulons have zero target-weight variance over shared genes.")
    piv = piv[keep_srcs]

    X_mat = piv.to_numpy(dtype=float, copy=True)
    Y_mat = prof_df.loc[shared_genes].to_numpy(dtype=float, copy=True)
    samples = prof_df.columns.tolist()

    if args.method == "mlm" and len(shared_genes) <= len(keep_srcs) + 2:
        warnings.append(
            f"Shared gene count ({len(shared_genes)}) is close to regulon count ({len(keep_srcs)}); "
            "MLM relies on ridge stabilization."
        )

    if args.method == "ulm":
        t_mat, p_mat = _run_ulm(X_mat, Y_mat)
    else:
        t_mat, p_mat = _run_mlm(X_mat, Y_mat)

    fdr_mat = np.column_stack([bh_fdr(p_mat[:, j]) for j in range(p_mat.shape[1])])

    act_df = pd.DataFrame(t_mat, index=keep_srcs, columns=samples)
    act_df.index.name = "source"
    pval_df = pd.DataFrame(p_mat, index=keep_srcs, columns=samples)
    pval_df.index.name = "source"
    fdr_df = pd.DataFrame(fdr_mat, index=keep_srcs, columns=samples)
    fdr_df.index.name = "source"

    # Rank top regulons by mean |t| across samples (or |t| in DE mode)
    summary_rows = []
    for i, src in enumerate(keep_srcs):
        scores_i = t_mat[i, :]
        pvals_i = p_mat[i, :]
        fdrs_i = fdr_mat[i, :]
        best_j = int(np.argmax(np.abs(scores_i)))
        mean_score = float(np.mean(scores_i))
        n_sig_samples = int(np.sum(fdrs_i < 0.05))
        summary_rows.append(
            {
                "source": src,
                "n_targets": int(counts.loc[src]),
                "score": round(float(scores_i[best_j]) if len(samples) == 1 else mean_score, 4),
                "max_abs_score": round(float(np.max(np.abs(scores_i))), 4),
                "pvalue": float(pvals_i[best_j]),
                "padj": float(fdrs_i[best_j]),
                "n_samples_fdr05": n_sig_samples,
                "direction": "activated" if (scores_i[best_j] if len(samples) == 1 else mean_score) >= 0 else "repressed",
            }
        )
    summary_rows.sort(key=lambda r: (r["padj"], -r["max_abs_score"]))

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    act_path = out_dir / "regulon_activities.tsv"
    pval_path = out_dir / "regulon_pvalues.tsv"
    fdr_path = out_dir / "regulon_fdr.tsv"
    fig_path = out_dir / "regulon_activity_heatmap.png"

    act_df.round(5).to_csv(act_path, sep="\t")
    pval_df.to_csv(pval_path, sep="\t")
    fdr_df.to_csv(fdr_path, sep="\t")

    inputs_list = [args.input] + ([args.network] if args.network else [])
    params = {"method": args.method, "min_targets": min_t, "network_source": net_source}
    write_artifact(
        act_path,
        "rna_regulon",
        f"Regulon Activity t-Statistics ({args.method.upper()})",
        inputs=inputs_list,
        parameters=params,
        shape=act_df.shape,
        scientific_context={
            "dataType": "pathway_scores",
            "valueScale": "zscore",
            "normalization": f"regulon_{args.method}_tstat",
            "entityIdType": "regulon_id",
            "organelle": "nuclear",
        },
    )
    write_artifact(
        pval_path,
        "rna_regulon",
        f"Regulon Activity Two-Sided p-Values ({args.method.upper()})",
        inputs=inputs_list,
        parameters=params,
        shape=pval_df.shape,
    )
    write_artifact(
        fdr_path,
        "rna_regulon",
        f"Regulon Activity Benjamini-Hochberg FDR ({args.method.upper()})",
        inputs=inputs_list,
        parameters=params,
        shape=fdr_df.shape,
    )

    # Render heatmap (matrix mode >=2 samples) or diverging barplot (DE / single-sample)
    plt = setup_matplotlib()
    top_srcs = [r["source"] for r in summary_rows[: min(30, len(summary_rows))]]
    sub_act = act_df.loc[top_srcs]
    sub_fdr = fdr_df.loc[top_srcs]

    if len(samples) >= 2:
        fig_w = max(7.0, min(14.0, 3.8 + 0.35 * len(samples)))
        fig_h = max(4.5, min(12.0, 2.0 + 0.32 * len(top_srcs)))
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=180)
        mat = sub_act.to_numpy(dtype=float, copy=True)
        vmax = max(2.0, float(np.nanpercentile(np.abs(mat), 98)))
        im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_yticks(range(len(top_srcs)))
        ax.set_yticklabels(top_srcs, fontsize=8.5)
        if len(samples) <= 30:
            ax.set_xticks(range(len(samples)))
            ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
        else:
            ax.set_xticks([])
            ax.set_xlabel(f"Samples (n={len(samples)})")
        if len(samples) <= 18 and len(top_srcs) <= 25:
            q_mat = sub_fdr.to_numpy(dtype=float, copy=True)
            for i in range(mat.shape[0]):
                for j in range(mat.shape[1]):
                    if q_mat[i, j] < 0.05:
                        txt_c = "white" if abs(mat[i, j]) > 0.55 * vmax else "#111827"
                        ax.text(j, i, "*" if q_mat[i, j] >= 0.01 else "**", ha="center", va="center", fontsize=8, color=txt_c)
        fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03, label=f"Regulon Activity (t-stat, {args.method.upper()})")
        ax.set_title(f"Regulon Activity Heatmap ({args.method.upper()}, * FDR < 0.05)")
    else:
        col_name = samples[0]
        df_plot = pd.DataFrame(
            {"source": top_srcs, "score": sub_act[col_name].to_numpy(dtype=float, copy=True), "padj": sub_fdr[col_name].to_numpy(dtype=float, copy=True)}
        ).sort_values("score", ascending=True)
        fig_h = max(4.2, min(11.0, 1.8 + 0.34 * len(df_plot)))
        fig, ax = plt.subplots(figsize=(8.0, fig_h), dpi=180)
        colors = ["#dc2626" if s >= 0 else "#2563eb" for s in df_plot["score"]]
        y_pos = np.arange(len(df_plot))
        ax.barh(y_pos, df_plot["score"], color=colors, alpha=0.88)
        ax.axvline(0.0, color="#374151", lw=1.0)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(
            [f"{r.source} ({'*' if r.padj < 0.05 else 'ns'}, q={r.padj:.2g})" for r in df_plot.itertuples()],
            fontsize=8.5,
        )
        ax.set_xlabel(f"Regulon Activity t-Statistic ({args.method.upper()})")
        ax.set_title(f"Regulon Activity Profile ({col_name})")

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "rna_regulon",
        f"Regulon Activity Heatmap ({args.method.upper()})",
        inputs=inputs_list,
        parameters=params,
    )

    n_sig_fdr = int(sum(1 for r in summary_rows if r["padj"] < 0.05))
    emit(
        {
            "status": "ok",
            "tool": "rna_regulon",
            "input_mode": input_mode,
            "method": args.method,
            "network_source": net_source,
            "n_samples": int(len(samples)),
            "n_genes_matched": int(len(shared_genes)),
            "n_regulons_tested": int(len(keep_srcs)),
            "n_significant_fdr": n_sig_fdr,
            "top_regulons": summary_rows[:14],
            "outputs": {
                "activities_table": str(act_path),
                "pvalues_table": str(pval_path),
                "fdr_table": str(fdr_path),
                "heatmap_plot": str(fig_path),
            },
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
