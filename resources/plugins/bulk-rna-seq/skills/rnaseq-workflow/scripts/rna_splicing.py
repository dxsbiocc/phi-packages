#!/usr/bin/env python3
"""rna_splicing: harmonize rMATS, SUPPA2, and generic alternative splicing results with reading-frame impact."""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import bh_fdr, emit, fail, setup_matplotlib, write_artifact

SUPPA_CLASS_MAP = {
    "SE": "SE",
    "A5": "A5SS",
    "A5SS": "A5SS",
    "A3": "A3SS",
    "A3SS": "A3SS",
    "MX": "MXE",
    "MXE": "MXE",
    "RI": "RI",
    "AF": "AF",
    "AL": "AL",
}


def _expand_inputs(raw_specs: list[str]) -> list[Path]:
    paths: list[Path] = []
    for spec in raw_specs:
        for token in str(spec).split(","):
            t = token.strip()
            if not t:
                continue
            p = Path(t)
            if not p.is_file():
                fail(f"Splicing input file not found: {p}")
            paths.append(p)
    if not paths:
        fail("At least one splicing_inputs file is required.")
    return paths


def _infer_rmats_class(path: Path, cols_lower: dict[str, str], df) -> str:
    fname = path.name.upper()
    for code in ("A5SS", "A3SS", "MXE", "SE", "RI"):
        if re.search(rf"(^|[^A-Z]){code}([^A-Z]|$)", fname):
            return code
    if "exonstart_0base" in cols_lower and "exonend" in cols_lower:
        return "SE"
    if "1stexonstart_0base" in cols_lower or "2ndexonstart_0base" in cols_lower:
        return "MXE"
    if "riexonstart_0base" in cols_lower or "riexonend" in cols_lower:
        return "RI"
    if "longexonstart_0base" in cols_lower and "shortes" in cols_lower:
        # Distinguish A5SS vs A3SS by strand and shared splice site on row 0
        row0 = df.iloc[0]
        strand = str(row0[cols_lower["strand"]]).strip() if "strand" in cols_lower else "+"
        ls = row0[cols_lower["longexonstart_0base"]]
        le = row0[cols_lower["longexonend"]]
        ss = row0[cols_lower["shortes"]]
        if ls == ss:
            return "A5SS" if strand == "+" else "A3SS"
        return "A3SS" if strand == "+" else "A5SS"
    return "OTHER"


def _mean_psi_list(series):
    import numpy as np

    vals = []
    for raw in series:
        parts = [p.strip() for p in str(raw).split(",") if p.strip() not in ("", "NA", "nan", "None")]
        nums = []
        for p in parts:
            try:
                nums.append(float(p))
            except ValueError:
                continue
        vals.append(float(np.mean(nums)) if nums else float("nan"))
    return vals


def _parse_suppa_event(ev_id: str) -> tuple[str, str, float]:
    """Parse SUPPA2 event_id e.g. ENSG00000141510.18;SE:chr17:100-150:220-300:+"""
    s = str(ev_id).strip()
    if ";" not in s:
        return s, "OTHER", float("nan")
    gene_part, coord_part = s.split(";", 1)
    gene = gene_part.strip()
    toks = coord_part.split(":")
    raw_cls = toks[0].strip().upper() if toks else "OTHER"
    ev_cls = SUPPA_CLASS_MAP.get(raw_cls, raw_cls or "OTHER")
    exon_len = float("nan")
    if ev_cls == "SE" and len(toks) >= 4:
        try:
            ex_start = int(toks[2].split("-")[1])
            ex_end = int(toks[3].split("-")[0])
            if ex_end >= ex_start:
                exon_len = float(ex_end - ex_start + 1)
        except Exception:
            pass
    elif ev_cls == "RI" and len(toks) >= 4:
        try:
            parts_ri = toks[3].split("-")
            if len(parts_ri) == 2:
                i_start, i_end = int(parts_ri[0]), int(parts_ri[1])
                if i_end > i_start:
                    exon_len = float(i_end - i_start - 1)
        except Exception:
            pass
    return gene, ev_cls, exon_len


def _parse_single_splicing_file(path: Path):
    import numpy as np
    import pandas as pd

    sep = "," if path.name.lower().endswith(".csv") else "\t"
    raw = pd.read_csv(path, sep=sep)
    if raw.empty:
        fail(f"Splicing table is empty: {path}")

    cols_lower = {c.lower().strip(): c for c in raw.columns}

    # 1. Detect rMATS format
    if "incleveldifference" in cols_lower and ("genesymbol" in cols_lower or "geneid" in cols_lower):
        ev_cls = _infer_rmats_class(path, cols_lower, raw)
        gene_col = cols_lower.get("genesymbol") or cols_lower.get("geneid")
        id_col = cols_lower.get("id")
        genes = raw[gene_col].astype(str).str.replace('"', "", regex=False).str.strip()
        if "geneid" in cols_lower:
            gene_ids = raw[cols_lower["geneid"]].astype(str).str.replace('"', "", regex=False).str.strip()
            genes = genes.where(genes != "", gene_ids)

        dpsi = pd.to_numeric(raw[cols_lower["incleveldifference"]], errors="coerce")
        pvals = pd.to_numeric(raw[cols_lower["pvalue"]], errors="coerce") if "pvalue" in cols_lower else pd.Series(np.nan, index=raw.index)
        fdrs = pd.to_numeric(raw[cols_lower["fdr"]], errors="coerce") if "fdr" in cols_lower else pd.Series(np.nan, index=raw.index)

        exon_len = pd.Series(np.nan, index=raw.index, dtype=float)
        if ev_cls == "SE" and "exonstart_0base" in cols_lower and "exonend" in cols_lower:
            exon_len = pd.to_numeric(raw[cols_lower["exonend"]], errors="coerce") - pd.to_numeric(raw[cols_lower["exonstart_0base"]], errors="coerce")
        elif ev_cls in ("A5SS", "A3SS") and "longexonstart_0base" in cols_lower and "shortes" in cols_lower:
            l_len = pd.to_numeric(raw[cols_lower["longexonend"]], errors="coerce") - pd.to_numeric(raw[cols_lower["longexonstart_0base"]], errors="coerce")
            s_len = pd.to_numeric(raw[cols_lower["shortee"]], errors="coerce") - pd.to_numeric(raw[cols_lower["shortes"]], errors="coerce")
            exon_len = (l_len - s_len).abs()
        elif ev_cls == "RI":
            if "downstreames" in cols_lower and "upstreamee" in cols_lower:
                exon_len = pd.to_numeric(raw[cols_lower["downstreames"]], errors="coerce") - pd.to_numeric(raw[cols_lower["upstreamee"]], errors="coerce")
            elif "riexonstart_0base" in cols_lower and "riexonend" in cols_lower:
                exon_len = pd.to_numeric(raw[cols_lower["riexonend"]], errors="coerce") - pd.to_numeric(raw[cols_lower["riexonstart_0base"]], errors="coerce")
        elif ev_cls == "MXE" and "1stexonstart_0base" in cols_lower and "1stexonend" in cols_lower:
            l1 = pd.to_numeric(raw[cols_lower["1stexonend"]], errors="coerce") - pd.to_numeric(raw[cols_lower["1stexonstart_0base"]], errors="coerce")
            if "2ndexonstart_0base" in cols_lower and "2ndexonend" in cols_lower:
                l2 = pd.to_numeric(raw[cols_lower["2ndexonend"]], errors="coerce") - pd.to_numeric(raw[cols_lower["2ndexonstart_0base"]], errors="coerce")
                diff = (l1 - l2).abs()
                exon_len = diff.where(diff > 0, l1)
            else:
                exon_len = l1

        ev_ids = [
            f"{ev_cls}:{g}:{raw.iloc[i][id_col] if id_col else i}"
            for i, g in enumerate(genes)
        ]
        psi1 = _mean_psi_list(raw[cols_lower["inclevel1"]]) if "inclevel1" in cols_lower else [float("nan")] * len(raw)
        psi2 = _mean_psi_list(raw[cols_lower["inclevel2"]]) if "inclevel2" in cols_lower else [float("nan")] * len(raw)

        out = pd.DataFrame(
            {
                "event_id": ev_ids,
                "gene": genes.to_numpy(copy=True),
                "event_class": ev_cls,
                "psi_group1": psi1,
                "psi_group2": psi2,
                "dpsi": dpsi.to_numpy(dtype=float, copy=True),
                "pvalue": pvals.to_numpy(dtype=float, copy=True),
                "fdr": fdrs.to_numpy(dtype=float, copy=True),
                "exon_len_bp": exon_len.to_numpy(dtype=float, copy=True),
                "source_format": "rMATS",
                "source_file": path.name,
            }
        )
        return out

    # 2. Detect SUPPA2 .dpsi or generic differential splicing format
    dpsi_col = next(
        (cols_lower[c] for c in ("dpsi", "delta_psi", "deltapsi", "psi_diff", "incleveldifference") if c in cols_lower),
        None,
    )
    if dpsi_col is None:
        # SUPPA2 sometimes prefixes dpsi column like "<cond1>-<cond2>_dPSI"
        dpsi_col = next((c for c in raw.columns if c.lower().endswith("_dpsi") or "dpsi" in c.lower()), None)
    p_col = next(
        (cols_lower[c] for c in ("pvalue", "p-val", "pval", "p_val", "p.value", "p") if c in cols_lower),
        None,
    )
    if p_col is None:
        p_col = next((c for c in raw.columns if "p-val" in c.lower() or c.lower().endswith("_pval")), None)
    fdr_col = next((cols_lower[c] for c in ("fdr", "padj", "qvalue", "q_value", "adj.p.val") if c in cols_lower), None)

    if dpsi_col is None:
        fail(f"Could not find dPSI / IncLevelDifference column in splicing input: {path}")

    ev_col = next((cols_lower[c] for c in ("event_id", "event", "id", "unnamed: 0") if c in cols_lower), raw.columns[0])
    gene_col = next((cols_lower[c] for c in ("gene", "genesymbol", "gene_name", "symbol", "gene_id", "geneid") if c in cols_lower), None)
    cls_col = next((cols_lower[c] for c in ("event_class", "event_type", "type", "class") if c in cols_lower), None)
    len_col = next((cols_lower[c] for c in ("exon_len_bp", "exon_length", "exon_len", "length") if c in cols_lower), None)

    raw_evs = raw[ev_col].astype(str).str.strip().tolist()
    is_suppa = any(";" in ev and ":" in ev for ev in raw_evs[: min(20, len(raw_evs))])

    parsed_genes = []
    parsed_classes = []
    parsed_lens = []
    for idx, ev in enumerate(raw_evs):
        if is_suppa:
            g_s, c_s, l_s = _parse_suppa_event(ev)
        else:
            g_s, c_s, l_s = ev, "OTHER", float("nan")
        if gene_col:
            gv = str(raw.iloc[idx][gene_col]).strip()
            if gv and gv.lower() != "nan":
                g_s = gv
        if cls_col:
            cv = str(raw.iloc[idx][cls_col]).strip().upper()
            c_s = SUPPA_CLASS_MAP.get(cv, cv or c_s)
        if len_col:
            lv = pd.to_numeric(raw.iloc[idx][len_col], errors="coerce")
            if math.isfinite(lv):
                l_s = float(lv)
        parsed_genes.append(g_s)
        parsed_classes.append(c_s)
        parsed_lens.append(l_s)

    # Generic exonStart_0base / exonEnd fallback
    if not len_col and "exonstart_0base" in cols_lower and "exonend" in cols_lower:
        ex_lens = pd.to_numeric(raw[cols_lower["exonend"]], errors="coerce") - pd.to_numeric(raw[cols_lower["exonstart_0base"]], errors="coerce")
        parsed_lens = ex_lens.to_numpy(dtype=float, copy=True).tolist()

    dpsi_vals = pd.to_numeric(raw[dpsi_col], errors="coerce").to_numpy(dtype=float, copy=True)
    pvals = pd.to_numeric(raw[p_col], errors="coerce").to_numpy(dtype=float, copy=True) if p_col else np.full(len(raw), np.nan)
    if fdr_col:
        fdrs = pd.to_numeric(raw[fdr_col], errors="coerce").to_numpy(dtype=float, copy=True)
    else:
        valid_p = np.where(np.isfinite(pvals), np.clip(pvals, 0.0, 1.0), 1.0)
        fdrs = bh_fdr(valid_p)

    return pd.DataFrame(
        {
            "event_id": raw_evs,
            "gene": parsed_genes,
            "event_class": parsed_classes,
            "psi_group1": np.nan,
            "psi_group2": np.nan,
            "dpsi": dpsi_vals,
            "pvalue": pvals,
            "fdr": fdrs,
            "exon_len_bp": parsed_lens,
            "source_format": "SUPPA2" if is_suppa else "generic",
            "source_file": path.name,
        }
    )


def main() -> None:
    ap = argparse.ArgumentParser(description="Harmonize and annotate differential alternative splicing events")
    ap.add_argument("--splicing_inputs", nargs="+", required=True, help="One or more rMATS/SUPPA2/generic splicing tables")
    ap.add_argument("--output_dir", required=True, help="Output directory")
    ap.add_argument("--dpsi_cutoff", type=float, default=0.10, help="Minimum |dPSI| cutoff (default 0.10)")
    ap.add_argument("--fdr_cutoff", type=float, default=0.05, help="Maximum FDR cutoff (default 0.05)")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd

    input_paths = _expand_inputs(args.splicing_inputs)
    frames = [_parse_single_splicing_file(p) for p in input_paths]
    all_df = pd.concat(frames, ignore_index=True)
    all_df = all_df.dropna(subset=["dpsi"]).copy()
    if all_df.empty:
        fail("No valid splicing events with finite dPSI found across input files.")

    # Fill missing FDR via BH on p-value if any file lacked FDR
    missing_fdr = ~np.isfinite(all_df["fdr"].to_numpy(dtype=float, copy=True))
    if missing_fdr.any():
        p_fill = np.where(np.isfinite(all_df.loc[missing_fdr, "pvalue"]), all_df.loc[missing_fdr, "pvalue"], 1.0)
        all_df.loc[missing_fdr, "fdr"] = bh_fdr(p_fill)

    ex_len = all_df["exon_len_bp"].to_numpy(dtype=float, copy=True)
    coding_impact = []
    for lv in ex_len:
        if math.isfinite(lv) and lv > 0:
            coding_impact.append("frameshift_candidate" if int(round(lv)) % 3 != 0 else "in_frame_modular")
        else:
            coding_impact.append("unknown")
    all_df["coding_impact"] = coding_impact
    all_df["abs_dpsi"] = all_df["dpsi"].abs()
    all_df["direction"] = np.where(all_df["dpsi"] >= 0, "inclusion", "skipping_or_exclusion")
    all_df["significant"] = (all_df["abs_dpsi"] >= float(args.dpsi_cutoff)) & (all_df["fdr"] <= float(args.fdr_cutoff))

    sig_df = all_df[all_df["significant"]].sort_values(["fdr", "abs_dpsi"], ascending=[True, False]).reset_index(drop=True)
    warnings: list[str] = []
    if sig_df.empty:
        warnings.append(
            f"Zero splicing events passed |dPSI| >= {args.dpsi_cutoff} and FDR <= {args.fdr_cutoff}."
        )

    # Event class summary across all parsed events
    class_rows = []
    for ev_cls, grp in all_df.groupby("event_class", sort=True):
        sig_g = grp[grp["significant"]]
        class_rows.append(
            {
                "event_class": str(ev_cls),
                "n_total": int(len(grp)),
                "n_significant": int(len(sig_g)),
                "n_inclusion": int((sig_g["dpsi"] > 0).sum()),
                "n_exclusion": int((sig_g["dpsi"] < 0).sum()),
                "n_frameshift_candidate": int((sig_g["coding_impact"] == "frameshift_candidate").sum()),
                "n_in_frame_modular": int((sig_g["coding_impact"] == "in_frame_modular").sum()),
                "median_abs_dpsi_sig": round(float(sig_g["abs_dpsi"].median()), 4) if len(sig_g) > 0 else 0.0,
            }
        )
    class_df = pd.DataFrame(class_rows).sort_values(["n_significant", "n_total"], ascending=[False, False]).reset_index(drop=True)

    # Multi-hit genes
    if not sig_df.empty:
        gene_rows = []
        for gene, g_df in sig_df.groupby("gene", sort=False):
            classes = sorted(set(g_df["event_class"].tolist()))
            gene_rows.append(
                {
                    "gene": str(gene),
                    "n_sig_events": int(len(g_df)),
                    "event_classes": ",".join(classes),
                    "n_event_classes": int(len(classes)),
                    "n_frameshift_candidates": int((g_df["coding_impact"] == "frameshift_candidate").sum()),
                    "max_abs_dpsi": round(float(g_df["abs_dpsi"].max()), 4),
                    "min_fdr": float(g_df["fdr"].min()),
                }
            )
        gene_df = pd.DataFrame(gene_rows).sort_values(
            ["n_sig_events", "max_abs_dpsi", "min_fdr"], ascending=[False, False, True]
        ).reset_index(drop=True)
    else:
        gene_df = pd.DataFrame(
            columns=["gene", "n_sig_events", "event_classes", "n_event_classes", "n_frameshift_candidates", "max_abs_dpsi", "min_fdr"]
        )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sig_path = out_dir / "splicing_significant_events.tsv"
    cls_path = out_dir / "splicing_event_class_summary.tsv"
    gene_path = out_dir / "splicing_multi_hit_genes.tsv"
    fig_path = out_dir / "splicing_dpsi_volcano.png"

    sig_cols = [
        "event_id", "gene", "event_class", "dpsi", "abs_dpsi", "direction",
        "pvalue", "fdr", "exon_len_bp", "coding_impact",
        "psi_group1", "psi_group2", "source_format", "source_file",
    ]
    sig_df[sig_cols].to_csv(sig_path, sep="\t", index=False)
    class_df.to_csv(cls_path, sep="\t", index=False)
    gene_df.to_csv(gene_path, sep="\t", index=False)

    str_inputs = [str(p) for p in input_paths]
    params = {"dpsi_cutoff": float(args.dpsi_cutoff), "fdr_cutoff": float(args.fdr_cutoff)}
    write_artifact(
        sig_path,
        "rna_splicing",
        "Significant Differential Alternative Splicing Events & Frame Impact",
        inputs=str_inputs,
        parameters=params,
        shape=(int(len(sig_df)), int(len(sig_cols))),
        scientific_context={
            "dataType": "de_stats",
            "effectField": "dpsi",
            "fdrField": "fdr",
        },
    )
    write_artifact(
        cls_path,
        "rna_splicing",
        "Alternative Splicing Event-Class Summary",
        inputs=str_inputs,
        parameters=params,
        shape=class_df.shape,
    )
    write_artifact(
        gene_path,
        "rna_splicing",
        "Multi-Hit Differentially Spliced Genes",
        inputs=str_inputs,
        parameters=params,
        shape=gene_df.shape,
    )

    # 2-panel volcano + event class & frame-shift breakdown
    plt = setup_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 5.0), dpi=180, gridspec_kw={"width_ratios": [1.3, 1.0]})

    dpsi_arr = all_df["dpsi"].to_numpy(dtype=float, copy=True)
    fdr_arr = np.clip(all_df["fdr"].to_numpy(dtype=float, copy=True), 1e-300, 1.0)
    neglogq = -np.log10(fdr_arr)
    sig_mask = all_df["significant"].to_numpy(dtype=bool, copy=True)

    ax1.scatter(dpsi_arr[~sig_mask], neglogq[~sig_mask], s=14, color="#9ca3af", alpha=0.55, label="NS")
    palette = {"SE": "#2563eb", "A5SS": "#f59e0b", "A3SS": "#10b981", "MXE": "#8b5cf6", "RI": "#ef4444", "AF": "#06b6d4", "AL": "#ec4899"}
    for ev_cls in class_df["event_class"].tolist():
        m = sig_mask & (all_df["event_class"] == ev_cls).to_numpy(dtype=bool, copy=True)
        if m.any():
            ax1.scatter(dpsi_arr[m], neglogq[m], s=26, color=palette.get(ev_cls, "#111827"), alpha=0.88, label=f"{ev_cls} (n={int(m.sum())})")
    ax1.axhline(-math.log10(max(1e-300, float(args.fdr_cutoff))), color="#6b7280", ls="--", lw=1.0)
    ax1.axvline(float(args.dpsi_cutoff), color="#6b7280", ls=":", lw=1.0)
    ax1.axvline(-float(args.dpsi_cutoff), color="#6b7280", ls=":", lw=1.0)
    for _, r in sig_df.head(6).iterrows():
        yv = -math.log10(max(1e-300, float(r["fdr"])))
        ax1.annotate(
            f"{r['gene']} ({r['event_class']})",
            (float(r["dpsi"]), yv),
            fontsize=7.3,
            xytext=(4, 3),
            textcoords="offset points",
        )
    ax1.set_xlabel("Delta PSI (dPSI / IncLevelDifference)")
    ax1.set_ylabel("-log10(FDR)")
    ax1.set_title("Differential Splicing Volcano")
    ax1.legend(loc="best", frameon=False, fontsize=8)

    y_pos = np.arange(len(class_df))
    fs_vals = class_df["n_frameshift_candidate"].to_numpy(dtype=float, copy=True)
    inf_vals = class_df["n_in_frame_modular"].to_numpy(dtype=float, copy=True)
    unk_vals = np.maximum(0.0, class_df["n_significant"].to_numpy(dtype=float, copy=True) - fs_vals - inf_vals)
    ax2.barh(y_pos, fs_vals, color="#dc2626", alpha=0.88, label="Frameshift (len % 3 != 0)")
    ax2.barh(y_pos, inf_vals, left=fs_vals, color="#2563eb", alpha=0.88, label="In-frame modular (len % 3 == 0)")
    if np.any(unk_vals > 0):
        ax2.barh(y_pos, unk_vals, left=fs_vals + inf_vals, color="#9ca3af", alpha=0.85, label="Unknown length")
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(class_df["event_class"].tolist(), fontsize=9)
    ax2.invert_yaxis()
    ax2.set_xlabel("Significant Events")
    ax2.set_title("Significant Events by Class & Reading-Frame Impact")
    ax2.legend(loc="lower right", frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    write_artifact(
        fig_path,
        "rna_splicing",
        "Differential Splicing Volcano & Reading-Frame Impact Breakdown",
        inputs=str_inputs,
        parameters=params,
    )

    n_fs = int((sig_df["coding_impact"] == "frameshift_candidate").sum()) if not sig_df.empty else 0
    n_inf = int((sig_df["coding_impact"] == "in_frame_modular").sum()) if not sig_df.empty else 0
    n_multi = int((gene_df["n_sig_events"] >= 2).sum()) if not gene_df.empty else 0

    emit(
        {
            "status": "ok",
            "tool": "rna_splicing",
            "n_input_files": int(len(input_paths)),
            "n_events_total": int(len(all_df)),
            "n_events_significant": int(len(sig_df)),
            "n_genes_significant": int(len(gene_df)),
            "n_multi_hit_genes": n_multi,
            "n_frameshift_candidates": n_fs,
            "n_in_frame_modular": n_inf,
            "event_class_summary": class_df.to_dict(orient="records"),
            "top_events": sig_df[sig_cols].head(10).to_dict(orient="records"),
            "top_multi_hit_genes": gene_df.head(10).to_dict(orient="records"),
            "outputs": {
                "significant_events": str(sig_path),
                "event_class_summary": str(cls_path),
                "multi_hit_genes": str(gene_path),
                "volcano_plot": str(fig_path),
            },
            "warnings": warnings,
        }
    )


if __name__ == "__main__":
    main()
