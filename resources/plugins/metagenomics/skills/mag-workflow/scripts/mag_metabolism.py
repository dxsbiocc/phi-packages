#!/usr/bin/env python3
"""mg_mag_metabolism: summarize MAG metabolic module completeness (from KO
matrices / eggNOG-mapper or gapseq pathway tables) and genome-scale metabolic
model (GEM) complexity from CarveMe / gapseq SBML (.xml / .xml.gz) files.
"""
import argparse
import gzip
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

from _common import bin_stem, emit, fail, setup_matplotlib, write_artifact

# Curated core microbial KEGG modules (KO sets for deterministic completeness scoring)
DEFAULT_MODULES = {
    "M00001_Glycolysis": ["K00844", "K01810", "K00850", "K01623", "K01803", "K00134", "K00927", "K01834", "K01689", "K00873"],
    "M00004_Pentose_phosphate": ["K00036", "K01057", "K00033", "K01807", "K01783", "K00615", "K00616"],
    "M00009_TCA_cycle": ["K01647", "K01681", "K00031", "K00164", "K01902", "K00239", "K01676", "K00024"],
    "M00086_Beta_oxidation": ["K00249", "K01692", "K00022", "K00632"],
    "M00175_Nitrogen_fixation": ["K02588", "K02586", "K02591"],
    "M00529_Denitrification": ["K00370", "K00368", "K04561", "K00376"],
    "M00530_DNRA": ["K00370", "K00362", "K03385"],
    "M00596_Sulfate_reduction": ["K00958", "K00394", "K00395", "K11180", "K11181"],
    "M00377_Wood_Ljungdahl": ["K00198", "K01938", "K01491", "K00297", "K15023"],
    "M00122_Cobalamin_B12": ["K02224", "K02232", "K02227", "K02231", "K02233"],
}


def _clean_bin(name):
    s = Path(str(name)).name
    if s.endswith(".gz"):
        s = s[:-3]
    for suf in ("-Pathways.tbl", "_pathways.tsv", ".sbml.xml", ".xml", ".sbml", ".tbl", ".tsv"):
        if s.lower().endswith(suf.lower()):
            s = s[: -len(suf)]
            break
    return bin_stem(s)


def load_module_map(path):
    df = pd.read_csv(path, sep="\t" if str(path).endswith(".tsv") else ",", dtype=str).fillna("")
    if not {"module_id", "ko_list"} <= set(df.columns):
        fail("module_map requires module_id and ko_list columns")
    out = {}
    for _, r in df.iterrows():
        mid = r["module_id"].strip()
        mname = r["module_name"].strip() if "module_name" in df.columns and r["module_name"].strip() else ""
        label = f"{mid}_{mname}" if mname else mid
        kos = [k.strip().removeprefix("ko:") for k in re.split(r"[,;\s]+", r["ko_list"]) if k.strip()]
        if mid and kos:
            out[label] = kos
    if not out:
        fail(f"no valid modules parsed from {path}")
    return out


def score_ko_matrix(ko_matrix_path, modules):
    df = pd.read_csv(ko_matrix_path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    cols_clean = {c.strip().removeprefix("ko:"): c for c in df.columns}
    scores = {}
    for b in df.index:
        present = {k for k, orig in cols_clean.items() if float(df.loc[b, orig]) > 0}
        scores[b] = {
            mod: round(len(present.intersection(kos)) / float(len(kos)), 4)
            for mod, kos in modules.items()
        }
    return pd.DataFrame.from_dict(scores, orient="index")


def parse_gapseq_pathways(paths):
    scores = {}
    counts = {}
    for raw in paths:
        p = Path(raw)
        if not p.is_file():
            fail(f"gapseq pathway file not found: {p}")
        b = _clean_bin(p)
        df = pd.read_csv(p, sep="\t", comment="#", dtype=str).fillna("")
        if "ID" not in df.columns or "Completeness" not in df.columns:
            fail(f"{p} does not look like a gapseq *-Pathways.tbl file (needs ID and Completeness)")
        comp_map = {}
        n_pred = 0
        for _, r in df.iterrows():
            pid = r["ID"].strip().strip("|")
            name = r["Name"].strip() if "Name" in df.columns and r["Name"].strip() else pid
            label = f"{pid}:{name}" if name != pid else pid
            try:
                cval = float(r["Completeness"])
                if cval > 1.0:
                    cval = cval / 100.0
            except ValueError:
                cval = 0.0
            comp_map[label] = round(max(0.0, min(1.0, cval)), 4)
            if str(r.get("Prediction", "")).strip().lower() in ("true", "1", "yes"):
                n_pred += 1
        scores[b] = comp_map
        counts[b] = n_pred
    mat = pd.DataFrame.from_dict(scores, orient="index").fillna(0.0)
    return mat, counts


def parse_sbml_models(paths):
    """Parse SBML XML (.xml or .xml.gz) GEMs from CarveMe / gapseq with stdlib ET."""
    stats = {}
    for raw in paths:
        p = Path(raw)
        if not p.is_file():
            fail(f"SBML model file not found: {p}")
        b = _clean_bin(p)
        opener = gzip.open if str(p).endswith(".gz") else open
        with opener(p, "rb") as fh:
            tree = ET.parse(fh)
        root = tree.getroot()
        n_species = 0
        n_reactions = 0
        n_exchange = 0
        n_genes = 0
        for elem in root.iter():
            tag = elem.tag.split("}")[-1]
            if tag == "species":
                n_species += 1
            elif tag == "reaction":
                n_reactions += 1
                rid = str(elem.attrib.get("id", ""))
                if rid.startswith("R_EX_") or rid.startswith("EX_"):
                    n_exchange += 1
            elif tag == "geneProduct":
                n_genes += 1
        stats[b] = {
            "sbml_metabolites": n_species,
            "sbml_reactions": n_reactions,
            "sbml_exchange_reactions": n_exchange,
            "sbml_genes": n_genes,
        }
    return stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ko_matrix", help="MAG x KO count matrix TSV (e.g. mag-ko-matrix.tsv from mg_mag_function).")
    ap.add_argument("--module_map", help="Optional custom TSV/CSV with module_id, module_name, ko_list columns.")
    ap.add_argument("--gapseq_pathways", nargs="+", help="gapseq *-Pathways.tbl files (one per MAG).")
    ap.add_argument("--sbml", nargs="+", help="CarveMe or gapseq SBML .xml(.gz) files (one per MAG).")
    ap.add_argument("--completeness_cutoff", type=float, default=0.75)
    ap.add_argument("--output_dir", required=True)
    a = ap.parse_args()

    if not (a.ko_matrix or a.gapseq_pathways or a.sbml):
        fail("give at least one of --ko_matrix / --gapseq_pathways / --sbml")

    sources = []
    bins = set()
    mod_frames = []
    gapseq_counts = {}
    sbml_stats = {}

    if a.ko_matrix:
        if not Path(a.ko_matrix).is_file():
            fail(f"ko_matrix not found: {a.ko_matrix}")
        modules = load_module_map(a.module_map) if a.module_map else DEFAULT_MODULES
        ko_comp = score_ko_matrix(a.ko_matrix, modules)
        sources.append("ko_modules")
        bins.update(ko_comp.index.tolist())
        mod_frames.append(ko_comp)

    if a.gapseq_pathways:
        gp_comp, gapseq_counts = parse_gapseq_pathways(a.gapseq_pathways)
        sources.append("gapseq_pathways")
        bins.update(gp_comp.index.tolist())
        mod_frames.append(gp_comp)

    if a.sbml:
        sbml_stats = parse_sbml_models(a.sbml)
        sources.append("sbml_gem")
        bins.update(sbml_stats.keys())

    if not bins:
        fail("no MAGs parsed from inputs")

    ordered_bins = sorted(bins)
    comp_df = None
    if mod_frames:
        comp_df = pd.concat(mod_frames, axis=1).reindex(ordered_bins).fillna(0.0)
        comp_df.index.name = "bin"

    summary_rows = []
    for b in ordered_bins:
        row = {"bin": b, "sources": ",".join(sources)}
        if comp_df is not None:
            vals = comp_df.loc[b].to_numpy(dtype=float)
            row["n_modules_evaluated"] = int(len(vals))
            row["n_complete_modules"] = int((vals >= a.completeness_cutoff).sum())
            row["mean_module_completeness"] = round(float(vals.mean()), 4) if len(vals) else 0.0
        if gapseq_counts:
            row["n_gapseq_predicted"] = int(gapseq_counts.get(b, 0))
        if sbml_stats:
            st = sbml_stats.get(b, {})
            row["sbml_metabolites"] = st.get("sbml_metabolites")
            row["sbml_reactions"] = st.get("sbml_reactions")
            row["sbml_exchange_reactions"] = st.get("sbml_exchange_reactions")
            row["sbml_genes"] = st.get("sbml_genes")
        summary_rows.append(row)

    outdir = Path(a.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = (
        ([a.ko_matrix] if a.ko_matrix else [])
        + ([a.module_map] if a.module_map else [])
        + (a.gapseq_pathways or [])
        + (a.sbml or [])
    )
    params = {"sources": sources, "completeness_cutoff": a.completeness_cutoff}
    artifacts = []
    outputs = {}

    sum_df = pd.DataFrame(summary_rows)
    sum_path = outdir / "mag-metabolism-summary.tsv"
    sum_df.to_csv(sum_path, sep="\t", index=False)
    outputs["metabolism_summary"] = str(sum_path)
    artifacts.append(
        write_artifact(
            sum_path,
            "mg_mag_metabolism",
            "Per-MAG metabolic module completeness and GEM summary",
            skill="mag-workflow",
            inputs=inputs,
            parameters=params,
            shape=sum_df.shape,
        )
    )

    if comp_df is not None:
        comp_path = outdir / "mag-module-completeness.tsv"
        comp_df.to_csv(comp_path, sep="\t")
        outputs["module_completeness"] = str(comp_path)
        artifacts.append(
            write_artifact(
                comp_path,
                "mg_mag_metabolism",
                "MAG x metabolic module completeness matrix",
                skill="mag-workflow",
                inputs=inputs,
                parameters=params,
                shape=comp_df.shape,
            )
        )

        # Render completeness heatmap (cap display at top 30 variable/complete modules)
        show_df = comp_df
        if show_df.shape[1] > 30:
            top_cols = show_df.mean(axis=0).sort_values(ascending=False).index[:30]
            show_df = show_df[top_cols]
        plt = setup_matplotlib()
        fig, ax = plt.subplots(
            figsize=(max(6, 0.45 * show_df.shape[1] + 2), max(4, 0.45 * show_df.shape[0] + 1.5))
        )
        im = ax.imshow(show_df.to_numpy(dtype=float), aspect="auto", vmin=0.0, vmax=1.0, cmap="YlGnBu")
        ax.set_xticks(range(show_df.shape[1]))
        ax.set_xticklabels([c[:28] for c in show_df.columns], rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(show_df.shape[0]))
        ax.set_yticklabels(show_df.index, fontsize=8)
        ax.set_title("MAG metabolic pathway/module completeness")
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("Completeness (0-1)")
        fig.tight_layout()
        fig_path = outdir / "mag-metabolism-heatmap.png"
        fig.savefig(fig_path, dpi=180)
        plt.close(fig)
        outputs["metabolism_figure"] = str(fig_path)
        artifacts.append(
            write_artifact(
                fig_path,
                "mg_mag_metabolism",
                "MAG metabolic module completeness heatmap",
                skill="mag-workflow",
                inputs=inputs,
                parameters=params,
            )
        )

    emit(
        {
            "tool": "mg_mag_metabolism",
            "n_bins": int(len(ordered_bins)),
            "sources": sources,
            "n_modules_evaluated": int(comp_df.shape[1]) if comp_df is not None else 0,
            "completeness_cutoff": float(a.completeness_cutoff),
            "n_sbml_models": int(len(sbml_stats)),
            "outputs": outputs,
            "artifacts": artifacts,
            "caveat": "Default KEGG module completeness scores presence of curated KO sets "
                      "(not full Boolean pathway logic); SBML counts summarize CarveMe/gapseq "
                      "reconstructed reactions without running flux balance analysis (FBA).",
        }
    )


if __name__ == "__main__":
    main()
