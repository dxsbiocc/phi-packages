"""Summarize MAG gene prediction and functional annotation (prokka/prodigal/eggNOG)."""
import argparse
import gzip
import re
from collections import Counter
from pathlib import Path

import pandas as pd

from _common import bin_stem, emit, fail, setup_matplotlib, write_artifact

PROKKA_KEYS = {"contigs": "n_contigs", "bases": "n_bases", "CDS": "n_cds",
               "gene": "n_genes", "rRNA": "n_rrna", "tRNA": "n_trna",
               "tmRNA": "n_tmrna", "misc_RNA": "n_misc_rna"}
RRNA_PRODUCTS = {"5S": "has_5s", "16S": "has_16s", "23S": "has_23s"}
FILE_SUFFIXES = (".emapper.annotations", ".faa", ".gff", ".fna", ".fasta", ".fa")


def open_text(path):
    p = str(path)
    return gzip.open(p, "rt", encoding="utf-8", errors="replace") if p.endswith(".gz") \
        else open(p, encoding="utf-8", errors="replace")


def clean_bin_id(name):
    s = Path(str(name)).name
    if s.endswith(".gz"):
        s = s[:-3]
    for suf in FILE_SUFFIXES:
        if s.endswith(suf):
            s = s[: -len(suf)]
            break
    if s.endswith("_prodigal"):
        s = s[: -len("_prodigal")]
    return re.sub(r"^(megahit|spades|metaspades)[-_]", "", s, flags=re.IGNORECASE)


def load_bin_map(path):
    df = pd.read_csv(path, sep="\t" if str(path).endswith(".tsv") else ",")
    if not {"file", "bin"} <= set(df.columns):
        fail("bin_map needs file,bin columns")
    return {Path(str(f)).name: str(b) for f, b in zip(df["file"], df["bin"])}


def attribute(path, bin_map):
    name = Path(path).name
    if name in bin_map:
        return bin_map[name]
    return clean_bin_id(name)


def parse_prokka_dir(root, bin_map):
    """Per-bin subdirectories (nf-core/prokka layout) or flat *.txt files."""
    root = Path(root)
    if not root.is_dir():
        fail(f"prokka_dir not found: {root}")
    entries = {}
    subdirs = sorted(p for p in root.iterdir() if p.is_dir())
    targets = []
    for d in subdirs:
        txts = sorted(d.glob("*.txt"))
        if txts:
            targets.append((bin_map.get(d.name, bin_stem(d.name)), txts[0], d))
    if not targets:
        for txt in sorted(root.glob("*.txt")):
            targets.append((attribute(txt, bin_map), txt, root))
    if not targets:
        fail(f"no prokka *.txt summaries under {root}")
    for bin_id, txt, folder in targets:
        row = entries.setdefault(bin_id, {"bin": bin_id, "sources": set()})
        row["sources"].add("prokka")
        with open(txt, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if ":" not in line:
                    continue
                key, _, val = line.partition(":")
                key = key.strip()
                if key in PROKKA_KEYS:
                    try:
                        row[PROKKA_KEYS[key]] = int(val.strip())
                    except ValueError:
                        pass
        tsvs = sorted(folder.glob("*.tsv"))
        if tsvs:
            try:
                feat = pd.read_csv(tsvs[0], sep="\t")
            except Exception:
                feat = None
            if feat is not None and {"ftype", "product"} <= set(feat.columns):
                products = feat.loc[feat["ftype"] == "rRNA", "product"].astype(str)
                for label, col in RRNA_PRODUCTS.items():
                    row[col] = bool(products.str.contains(
                        rf"\b{label}\b.*ribosomal RNA", case=False, regex=True).any())
                if "n_trna" not in row:
                    row["n_trna"] = int((feat["ftype"] == "tRNA").sum())
    return entries


def parse_prodigal(paths, bin_map):
    entries = {}
    for path in paths:
        p = Path(path)
        if not p.is_file():
            fail(f"not a file: {p}")
        bin_id = attribute(p, bin_map)
        row = entries.setdefault(bin_id, {"bin": bin_id, "sources": set()})
        row["sources"].add("prodigal")
        name = p.name.lower()
        n = 0
        with open_text(p) as fh:
            if ".gff" in name:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.split("\t")
                    if len(f) >= 3 and f[2] == "CDS":
                        n += 1
            else:
                for line in fh:
                    if line.startswith(">"):
                        n += 1
        if n == 0:
            fail(f"no genes parsed from {p} (expected prodigal .faa/.gff, optionally .gz)")
        row["n_cds"] = max(row.get("n_cds", 0), n)
    return entries


def parse_eggnog(paths, bin_map):
    ko = {}
    cog = {}
    stats = {}
    for path in paths:
        p = Path(path)
        if not p.is_file():
            fail(f"not a file: {p}")
        bin_id = attribute(p, bin_map)
        with open_text(p) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            if not header or header[0].lstrip("#") != "query":
                fail(f"{p} is not an eggNOG-mapper .emapper.annotations table "
                     "(first column must be #query)")
            idx = {c.lstrip("#"): i for i, c in enumerate(header)}
            ko_i = idx.get("KEGG_ko")
            cog_i = idx.get("COG_category")
            name_i = idx.get("Preferred_name")
            if ko_i is None and cog_i is None:
                fail(f"{p} lacks KEGG_ko and COG_category columns")
            s = stats.setdefault(bin_id, {"n_queries": 0, "n_annotated": 0})
            k = ko.setdefault(bin_id, Counter())
            c = cog.setdefault(bin_id, Counter())
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                s["n_queries"] += 1
                annotated = False
                if ko_i is not None and ko_i < len(f):
                    for term in f[ko_i].split(","):
                        term = term.strip()
                        if term and term != "-":
                            k[term.removeprefix("ko:")] += 1
                            annotated = True
                if cog_i is not None and cog_i < len(f):
                    for letter in f[cog_i].strip():
                        if letter.isalpha() and letter.isupper():
                            c[letter] += 1
                            annotated = True
                if name_i is not None and name_i < len(f) and f[name_i].strip() not in ("", "-"):
                    annotated = True
                if annotated:
                    s["n_annotated"] += 1
    return ko, cog, stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prokka_dir")
    ap.add_argument("--prodigal", nargs="+")
    ap.add_argument("--eggnog", nargs="+")
    ap.add_argument("--bin_map")
    ap.add_argument("--output_dir", required=True)
    a = ap.parse_args()
    if not (a.prokka_dir or a.prodigal or a.eggnog):
        fail("give at least one of --prokka_dir / --prodigal / --eggnog")
    bin_map = load_bin_map(a.bin_map) if a.bin_map else {}

    entries = {}
    if a.prokka_dir:
        for b, row in parse_prokka_dir(a.prokka_dir, bin_map).items():
            entries.setdefault(b, {"bin": b, "sources": set()}).update(row)
            entries[b]["sources"] |= row["sources"]
    if a.prodigal:
        for b, row in parse_prodigal(a.prodigal, bin_map).items():
            target = entries.setdefault(b, {"bin": b, "sources": set()})
            target["sources"] |= row["sources"]
            if "n_cds" not in target:
                target["n_cds"] = row["n_cds"]
    ko, cog, egg_stats = ({}, {}, {})
    if a.eggnog:
        ko, cog, egg_stats = parse_eggnog(a.eggnog, bin_map)
        for b, s in egg_stats.items():
            target = entries.setdefault(b, {"bin": b, "sources": set()})
            target["sources"].add("eggnog")
            target["n_queries"] = s["n_queries"]
            target["n_annotated"] = s["n_annotated"]
    if not entries:
        fail("no bins parsed from the given inputs")

    for row in entries.values():
        if {"has_5s", "has_16s", "has_23s"} <= set(row) and "n_trna" in row:
            row["mimag_structural_pass"] = bool(
                row["has_5s"] and row["has_16s"] and row["has_23s"] and row["n_trna"] >= 18)
        row["sources"] = ",".join(sorted(row["sources"]))

    outdir = Path(a.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = ([a.prokka_dir] if a.prokka_dir else []) + (a.prodigal or []) + (a.eggnog or [])
    params = {"sources": sorted({s for r in entries.values() for s in r["sources"].split(",")})}
    artifacts = []
    outputs = {}

    cols = ["bin", "sources", "n_contigs", "n_bases", "n_cds", "n_genes", "n_rrna",
            "n_trna", "n_tmrna", "n_misc_rna", "has_5s", "has_16s", "has_23s",
            "mimag_structural_pass", "n_queries", "n_annotated"]
    counts = pd.DataFrame(sorted(entries.values(), key=lambda r: r["bin"]))
    counts = counts.reindex(columns=cols)
    counts_path = outdir / "mag-gene-counts.tsv"
    counts.to_csv(counts_path, sep="\t", index=False)
    outputs["gene_counts"] = str(counts_path)
    artifacts.append(write_artifact(counts_path, "mg_mag_function",
                                    "MAG gene counts and MIMAG structural evidence",
                                    skill="mag-workflow", inputs=inputs,
                                    parameters=params, shape=counts.shape))

    if ko:
        ko_df = pd.DataFrame.from_dict({b: dict(c) for b, c in ko.items()},
                                       orient="index").fillna(0).astype(int)
        ko_df = ko_df.reindex(sorted(ko_df.columns), axis=1)
        ko_df.index.name = "bin"
        ko_path = outdir / "mag-ko-matrix.tsv"
        ko_df.to_csv(ko_path, sep="\t")
        outputs["ko_matrix"] = str(ko_path)
        artifacts.append(write_artifact(ko_path, "mg_mag_function",
                                        "MAG x KO count matrix (eggNOG-mapper)",
                                        skill="mag-workflow", inputs=inputs,
                                        parameters=params, shape=ko_df.shape))
    if cog:
        cog_df = pd.DataFrame.from_dict({b: dict(c) for b, c in cog.items()},
                                        orient="index").fillna(0).astype(int)
        cog_df = cog_df.reindex(sorted(cog_df.columns), axis=1)
        cog_df.index.name = "bin"
        cog_path = outdir / "mag-cog-matrix.tsv"
        cog_df.to_csv(cog_path, sep="\t")
        outputs["cog_matrix"] = str(cog_path)
        artifacts.append(write_artifact(cog_path, "mg_mag_function",
                                        "MAG x COG category count matrix (eggNOG-mapper)",
                                        skill="mag-workflow", inputs=inputs,
                                        parameters=params, shape=cog_df.shape))
        plt = setup_matplotlib()
        fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(cog_df)), 5))
        bottom = None
        cmap = plt.get_cmap("tab20")
        for i, letter in enumerate(cog_df.columns):
            vals = cog_df[letter].values
            ax.bar(range(len(cog_df)), vals, bottom=bottom,
                   color=cmap(i % 20), label=letter, width=0.8)
            bottom = vals if bottom is None else bottom + vals
        ax.set_xticks(range(len(cog_df)))
        ax.set_xticklabels(cog_df.index, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Annotated genes (COG hits)")
        ax.set_title("COG functional categories per MAG")
        ax.legend(fontsize=6, ncol=4, loc="upper right")
        fig.tight_layout()
        fig_path = outdir / "mag-cog-categories.png"
        fig.savefig(fig_path, dpi=200)
        plt.close(fig)
        outputs["cog_figure"] = str(fig_path)
        artifacts.append(write_artifact(fig_path, "mg_mag_function",
                                        "COG functional categories per MAG",
                                        skill="mag-workflow", inputs=inputs, parameters=params))

    structural = counts["mimag_structural_pass"]
    emit({
        "tool": "mg_mag_function",
        "n_bins": int(len(counts)),
        "sources": params["sources"],
        "n_kos": int(len(ko_df.columns)) if ko else 0,
        "n_cog_categories": int(len(cog_df.columns)) if cog else 0,
        "n_annotated_queries": int(sum(s["n_annotated"] for s in egg_stats.values())),
        "mimag_structural": {
            "pass": int((structural == True).sum()),  # noqa: E712
            "fail": int((structural == False).sum()),  # noqa: E712
            "unknown": int(structural.isna().sum()),
        },
        "outputs": outputs,
        "artifacts": artifacts,
        "caveat": "rRNA/tRNA evidence is prokka's prediction, not barrnap/tRNAscan-SE "
                  "confirmation; KO/COG counts are eggNOG-mapper hit counts, not "
                  "expression or abundance.",
    })


if __name__ == "__main__":
    main()