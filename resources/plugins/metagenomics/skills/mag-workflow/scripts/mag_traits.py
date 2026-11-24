#!/usr/bin/env python3
"""mg_mag_traits: summarize MAG AMR, virulence, and secondary-metabolite BGC
annotations produced by upstream wrappers (AMRFinderPlus, ABRicate, RGI/CARD,
and antiSMASH).
"""
import argparse
import gzip
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd

from _common import bin_stem, emit, fail, setup_matplotlib, write_artifact

FILE_SUFFIXES = (
    ".amrfinder.tsv", ".amrfinder", ".abricate.tsv", ".abricate",
    ".rgi.txt", ".rgi.tsv", ".antismash.json", ".json", ".tsv", ".txt",
    ".fna", ".fasta", ".fa",
)


def open_text(path):
    p = str(path)
    return gzip.open(p, "rt", encoding="utf-8", errors="replace") if p.endswith(".gz") \
        else open(p, encoding="utf-8", errors="replace")


def clean_bin_id(name):
    s = Path(str(name)).name
    if s.endswith(".gz"):
        s = s[:-3]
    lower = s.lower()
    for suf in FILE_SUFFIXES:
        if lower.endswith(suf):
            s = s[: -len(suf)]
            break
    s = re.sub(r"_(amrfinder|abricate|rgi|antismash)$", "", s, flags=re.IGNORECASE)
    return bin_stem(s)


def load_bin_map(path):
    df = pd.read_csv(path, sep="\t" if str(path).endswith(".tsv") else ",")
    if not {"file", "bin"} <= set(df.columns):
        fail("bin_map needs file,bin columns")
    return {Path(str(f)).name: str(b) for f, b in zip(df["file"], df["bin"])}


def attribute(name_or_path, bin_map):
    name = Path(str(name_or_path)).name
    if name in bin_map:
        return bin_map[name]
    return clean_bin_id(name)


def _to_float(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def _ensure_bin(entries, bin_id, source):
    row = entries.setdefault(
        bin_id,
        {
            "bin": bin_id,
            "sources": set(),
            "amr_hits": 0,
            "virulence_hits": 0,
            "stress_hits": 0,
            "bgc_regions": 0,
            "amr_classes": Counter(),
            "bgc_types": Counter(),
        },
    )
    row["sources"].add(source)
    return row


def parse_amrfinder(paths, bin_map, entries, min_id, min_cov):
    for raw_path in paths:
        p = Path(raw_path)
        if not p.is_file():
            fail(f"amrfinder file not found: {p}")
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        cols = {c.strip(): c for c in df.columns}
        elem_col = cols.get("Element type") or cols.get("Type")
        class_col = cols.get("Class") or cols.get("Subclass")
        id_col = cols.get("% Identity to reference sequence") or cols.get("% Identity to reference")
        cov_col = cols.get("% Coverage of reference sequence") or cols.get("% Coverage of reference")
        name_col = cols.get("Name") or cols.get("bin") or cols.get("Genome")
        if elem_col is None and class_col is None:
            fail(f"{p} does not look like an AMRFinderPlus TSV (missing Element type / Class)")
        default_bin = attribute(p, bin_map)
        _ensure_bin(entries, default_bin, "amrfinderplus")
        for _, r in df.iterrows():
            ident = _to_float(r[id_col]) if id_col else None
            cov = _to_float(r[cov_col]) if cov_col else None
            if ident is not None and ident < min_id:
                continue
            if cov is not None and cov < min_cov:
                continue
            b = attribute(r[name_col], bin_map) if (name_col and r[name_col].strip()) else default_bin
            rec = _ensure_bin(entries, b, "amrfinderplus")
            etype = (r[elem_col].strip().upper() if elem_col else "AMR") or "AMR"
            if etype == "VIRULENCE":
                rec["virulence_hits"] += 1
            elif etype == "STRESS":
                rec["stress_hits"] += 1
            else:
                rec["amr_hits"] += 1
                cls = (r[class_col].strip() if class_col else "") or "UNCLASSIFIED"
                for part in re.split(r"[/;,]+", cls):
                    part = part.strip()
                    if part:
                        rec["amr_classes"][part] += 1


def parse_abricate(paths, bin_map, entries, min_id, min_cov):
    for raw_path in paths:
        p = Path(raw_path)
        if not p.is_file():
            fail(f"abricate file not found: {p}")
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        df.columns = [c.lstrip("#").strip() for c in df.columns]
        if "GENE" not in df.columns:
            fail(f"{p} does not look like an ABRicate TSV (missing GENE column)")
        default_bin = attribute(p, bin_map)
        for _, r in df.iterrows():
            ident = _to_float(r.get("%IDENTITY"))
            cov = _to_float(r.get("%COVERAGE"))
            if ident is not None and ident < min_id:
                continue
            if cov is not None and cov < min_cov:
                continue
            file_val = str(r.get("FILE", "")).strip()
            b = attribute(file_val, bin_map) if file_val else default_bin
            rec = _ensure_bin(entries, b, "abricate")
            db = str(r.get("DATABASE", "")).strip().lower()
            if "vfdb" in db or "virulence" in db:
                rec["virulence_hits"] += 1
            else:
                rec["amr_hits"] += 1
                res = str(r.get("RESISTANCE", "")).strip() or str(r.get("GENE", "")).strip() or "UNCLASSIFIED"
                for part in re.split(r"[;,/]+", res):
                    part = part.strip()
                    if part:
                        rec["amr_classes"][part] += 1


def parse_rgi(paths, bin_map, entries, min_id, min_cov):
    for raw_path in paths:
        p = Path(raw_path)
        if not p.is_file():
            fail(f"rgi file not found: {p}")
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        if "Best_Hit_ARO" not in df.columns and "Drug Class" not in df.columns:
            fail(f"{p} does not look like an RGI TSV (missing Best_Hit_ARO / Drug Class)")
        b = attribute(p, bin_map)
        rec = _ensure_bin(entries, b, "rgi")
        for _, r in df.iterrows():
            ident = _to_float(r.get("Best_Identities"))
            cov = _to_float(r.get("Percentage Length of Reference Sequence"))
            if ident is not None and ident < min_id:
                continue
            if cov is not None and cov < min_cov:
                continue
            rec["amr_hits"] += 1
            dc = str(r.get("Drug Class", "")).strip() or "UNCLASSIFIED"
            for part in re.split(r"[;]+", dc):
                part = part.strip()
                if part:
                    rec["amr_classes"][part] += 1


def _extract_antismash_json(json_path, bin_id, entries):
    data = json.loads(Path(json_path).read_text(encoding="utf-8", errors="replace"))
    rec = _ensure_bin(entries, bin_id, "antismash")
    records = data.get("records", []) if isinstance(data, dict) else []
    for record in records:
        for feat in record.get("features", []):
            if feat.get("type") != "region":
                continue
            rec["bgc_regions"] += 1
            quals = feat.get("qualifiers", {})
            products = quals.get("product", [])
            if isinstance(products, str):
                products = [products]
            if not products:
                rec["bgc_types"]["unknown"] += 1
            else:
                for prod in products:
                    p_str = str(prod).strip()
                    if p_str:
                        rec["bgc_types"][p_str] += 1


def parse_antismash(paths, bin_map, entries):
    for raw_path in paths:
        p = Path(raw_path)
        if p.is_dir():
            # Either per-bin subdirectories containing *.json or direct *.json files
            found = False
            for sub in sorted(d for d in p.iterdir() if d.is_dir()):
                jsons = sorted(sub.glob("*.json"))
                if jsons:
                    found = True
                    b = bin_map.get(sub.name, bin_stem(sub.name))
                    _extract_antismash_json(jsons[0], b, entries)
            if not found:
                for jf in sorted(p.glob("*.json")):
                    found = True
                    _extract_antismash_json(jf, attribute(jf, bin_map), entries)
            if not found:
                fail(f"no antiSMASH *.json files found under {p}")
        elif p.is_file():
            _extract_antismash_json(p, attribute(p, bin_map), entries)
        else:
            fail(f"antismash path not found: {p}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amrfinder", nargs="+")
    ap.add_argument("--abricate", nargs="+")
    ap.add_argument("--rgi", nargs="+")
    ap.add_argument("--antismash", nargs="+")
    ap.add_argument("--bin_map")
    ap.add_argument("--min_identity", type=float, default=0.0)
    ap.add_argument("--min_coverage", type=float, default=0.0)
    ap.add_argument("--output_dir", required=True)
    a = ap.parse_args()

    if not (a.amrfinder or a.abricate or a.rgi or a.antismash):
        fail("give at least one of --amrfinder / --abricate / --rgi / --antismash")

    bin_map = load_bin_map(a.bin_map) if a.bin_map else {}
    entries = {}

    if a.amrfinder:
        parse_amrfinder(a.amrfinder, bin_map, entries, a.min_identity, a.min_coverage)
    if a.abricate:
        parse_abricate(a.abricate, bin_map, entries, a.min_identity, a.min_coverage)
    if a.rgi:
        parse_rgi(a.rgi, bin_map, entries, a.min_identity, a.min_coverage)
    if a.antismash:
        parse_antismash(a.antismash, bin_map, entries)

    if not entries:
        fail("no MAG trait records parsed from the given inputs")

    outdir = Path(a.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = (a.amrfinder or []) + (a.abricate or []) + (a.rgi or []) + (a.antismash or [])
    sources = sorted({s for rec in entries.values() for s in rec["sources"]})
    params = {
        "sources": sources,
        "min_identity": a.min_identity,
        "min_coverage": a.min_coverage,
    }
    artifacts = []
    outputs = {}

    summary_rows = []
    amr_dict = {}
    bgc_dict = {}
    for b in sorted(entries.keys()):
        rec = entries[b]
        summary_rows.append(
            {
                "bin": b,
                "sources": ",".join(sorted(rec["sources"])),
                "n_amr_hits": rec["amr_hits"],
                "n_amr_classes": len(rec["amr_classes"]),
                "amr_classes": ";".join(sorted(rec["amr_classes"].keys())),
                "n_virulence_hits": rec["virulence_hits"],
                "n_stress_hits": rec["stress_hits"],
                "n_bgc_regions": rec["bgc_regions"],
                "bgc_types": ";".join(sorted(rec["bgc_types"].keys())),
            }
        )
        if rec["amr_classes"]:
            amr_dict[b] = dict(rec["amr_classes"])
        if rec["bgc_types"]:
            bgc_dict[b] = dict(rec["bgc_types"])

    summary_df = pd.DataFrame(summary_rows)
    summary_path = outdir / "mag-traits-summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    outputs["traits_summary"] = str(summary_path)
    artifacts.append(
        write_artifact(
            summary_path,
            "mg_mag_traits",
            "Per-MAG AMR, virulence, stress, and BGC trait summary",
            skill="mag-workflow",
            inputs=inputs,
            parameters=params,
            shape=summary_df.shape,
        )
    )

    n_amr_classes = 0
    if amr_dict:
        amr_df = (
            pd.DataFrame.from_dict(amr_dict, orient="index")
            .reindex(sorted(entries.keys()))
            .fillna(0)
            .astype(int)
        )
        amr_df = amr_df.reindex(sorted(amr_df.columns), axis=1)
        amr_df.index.name = "bin"
        n_amr_classes = int(len(amr_df.columns))
        amr_path = outdir / "mag-amr-matrix.tsv"
        amr_df.to_csv(amr_path, sep="\t")
        outputs["amr_matrix"] = str(amr_path)
        artifacts.append(
            write_artifact(
                amr_path,
                "mg_mag_traits",
                "MAG x AMR drug class count matrix",
                skill="mag-workflow",
                inputs=inputs,
                parameters=params,
                shape=amr_df.shape,
            )
        )

    n_bgc_types = 0
    if bgc_dict:
        bgc_df = (
            pd.DataFrame.from_dict(bgc_dict, orient="index")
            .reindex(sorted(entries.keys()))
            .fillna(0)
            .astype(int)
        )
        bgc_df = bgc_df.reindex(sorted(bgc_df.columns), axis=1)
        bgc_df.index.name = "bin"
        n_bgc_types = int(len(bgc_df.columns))
        bgc_path = outdir / "mag-bgc-matrix.tsv"
        bgc_df.to_csv(bgc_path, sep="\t")
        outputs["bgc_matrix"] = str(bgc_path)
        artifacts.append(
            write_artifact(
                bgc_path,
                "mg_mag_traits",
                "MAG x BGC product type count matrix (antiSMASH)",
                skill="mag-workflow",
                inputs=inputs,
                parameters=params,
                shape=bgc_df.shape,
            )
        )

    # Stacked overview barplot
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(max(6, 0.55 * len(summary_df)), 4.5))
    x_pos = range(len(summary_df))
    amr_vals = summary_df["n_amr_hits"].to_numpy()
    vf_vals = summary_df["n_virulence_hits"].to_numpy()
    bgc_vals = summary_df["n_bgc_regions"].to_numpy()
    ax.bar(x_pos, amr_vals, label="AMR hits", color="#E11D48", width=0.72)
    ax.bar(x_pos, vf_vals, bottom=amr_vals, label="Virulence hits", color="#F59E0B", width=0.72)
    ax.bar(x_pos, bgc_vals, bottom=amr_vals + vf_vals, label="BGC regions", color="#0EA5E9", width=0.72)
    ax.set_xticks(list(x_pos))
    ax.set_xticklabels(summary_df["bin"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Annotated trait count")
    ax.set_title("AMR, virulence, and BGC traits per MAG")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig_path = outdir / "mag-traits-overview.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)
    outputs["traits_figure"] = str(fig_path)
    artifacts.append(
        write_artifact(
            fig_path,
            "mg_mag_traits",
            "Per-MAG AMR, virulence, and BGC counts",
            skill="mag-workflow",
            inputs=inputs,
            parameters=params,
        )
    )

    emit(
        {
            "tool": "mg_mag_traits",
            "n_bins": int(len(summary_df)),
            "sources": sources,
            "total_amr_hits": int(summary_df["n_amr_hits"].sum()),
            "total_virulence_hits": int(summary_df["n_virulence_hits"].sum()),
            "total_stress_hits": int(summary_df["n_stress_hits"].sum()),
            "total_bgc_regions": int(summary_df["n_bgc_regions"].sum()),
            "n_amr_classes": n_amr_classes,
            "n_bgc_types": n_bgc_types,
            "outputs": outputs,
            "artifacts": artifacts,
            "caveat": "Gene/BGC presence reflects sequence-homology predictions from upstream "
                      "wrappers (AMRFinderPlus/ABRicate/RGI/antiSMASH), not phenotypic resistance "
                      "or active metabolite production.",
        }
    )


if __name__ == "__main__":
    main()
