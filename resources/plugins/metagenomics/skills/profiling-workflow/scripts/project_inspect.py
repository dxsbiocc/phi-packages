#!/usr/bin/env python3
"""mg_inspect: identify metagenomics profiling inputs and report what they support."""
import argparse
import re
from pathlib import Path

from _common import emit, fail, read_table, sample_stem

RANK_CODES = set("URDKPCOFGS")
KRAKEN_RANKS = ["U", "R", "D", "K", "P", "C", "O", "F", "G", "S"]
MPA_RANKS = ["d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__", "t__"]


def classify(path):
    """Return (kind, detail) for one file, using headers and row shape only."""
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            lines = []
            for i, line in enumerate(fh):
                lines.append(line.rstrip("\n"))
                if i >= 40:
                    break
    except OSError as exc:
        return "unreadable", {"error": str(exc)}
    if not lines:
        return "empty", {}
    head = lines[0]
    fields = head.split("\t")
    low = head.lower()
    # MetaPhlAn/MPA files lead with comment lines (#mpa_v..., #SampleID,
    # #clade_name header); look past them for the real shape.
    data_lines = [l for l in lines if l and not l.startswith("#")]
    first_data = data_lines[0] if data_lines else ""
    data_fields = first_data.split("\t")
    comment_low = " ".join(l.lower() for l in lines[:5] if l.startswith("#"))

    if ("clade_name" in low and "relative_abundance" in low) or (
        "clade_name" in comment_low and "relative_abundance" in comment_low
    ):
        ranks = set()
        n = 0
        for line in data_lines:
            clade = line.split("\t")[0]
            n += 1
            for marker in MPA_RANKS:
                if marker in clade:
                    ranks.add(marker.strip("_"))
        return "metaphlan_profile", {"rows": n, "ranks": sorted(ranks)}

    if "taxonomy_lvl" in low and "new_est_reads" in low:
        return "bracken_report", {"header": True}

    if low.startswith("# gene family") or low.startswith("gene family") or (
        "humann" in low and "\t" in head
    ):
        strat = sum(1 for l in lines[1:] if "|" in l.split("\t")[0])
        return "humann_genefamilies", {"stratified_rows": strat}

    if re.match(r"^(PWY-|UNMAPPED|UNINTEGRATED)", fields[0]) and len(fields) >= 2:
        return "humann_pathabundance", {}

    if head.startswith("#Classification") or head.startswith("#SampleID") or (
        len(data_fields) >= 2 and re.match(r"^[a-z]__", data_fields[0]) and "|" in data_fields[0]
    ):
        return "mpa_table", {}

    if "completeness" in low and "contamination" in low and (
        "name" in low or "bin id" in low
    ):
        return "bin_quality_table", {"note": "belongs to the mag-workflow route"}

    if fields[0].strip() in ("Assembly", "N50", "# contigs") or re.match(
        r"^(N50|L50|Total length|# contigs)", fields[0]
    ):
        return "quast_report", {"note": "belongs to the mag-workflow route"}

    if len(fields) == 6 and fields[3].strip()[:1] in RANK_CODES:
        try:
            float(fields[0])
            int(fields[1])
            int(fields[2])
            int(fields[4])
        except ValueError:
            pass
        else:
            ranks = set()
            n = 0
            for line in lines:
                f = line.split("\t")
                if len(f) == 6:
                    ranks.add(f[3].strip())
                    n += 1
            return "kraken2_report", {"rows": n, "ranks": sorted(ranks)}

    if len(fields) == 6 and fields[2].strip() in RANK_CODES:
        try:
            int(fields[1])
            int(fields[3])
            float(fields[4])
        except ValueError:
            pass
        else:
            return "bracken_report", {"header": False}

    if len(data_fields) >= 3 and re.match(r"^[a-z]__", data_fields[0]):
        return "mpa_table", {}

    if len(fields) >= 2:
        try:
            df = read_table(path)
        except Exception:
            return "unknown", {}
        numeric = df.apply(lambda c: c.astype(str).str.match(r"^-?\d+(\.\d+)?$").all())
        if numeric.all() and len(df.columns) >= 1:
            if {"sample", "fastq_1", "fastq_2"} & {c.lower() for c in df.columns}:
                return "samplesheet", {}
            return "abundance_matrix", {
                "samples": df.columns.tolist(),
                "taxa": int(len(df)),
            }
    return "unknown", {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    args = ap.parse_args()
    root = Path(args.input)
    if not root.exists():
        fail(f"input does not exist: {root}")
    paths = (
        sorted(p for p in root.rglob("*") if p.is_file() and not p.name.startswith("."))
        if root.is_dir()
        else [root]
    )
    if not paths:
        fail(f"no files found under {root}")
    files = []
    for p in paths[:200]:
        kind, detail = classify(p)
        entry = {"path": str(p), "kind": kind, "sample": sample_stem(p)}
        entry.update(detail)
        files.append(entry)
    kinds = {f["kind"] for f in files}
    profile_kinds = {
        "kraken2_report", "bracken_report", "metaphlan_profile", "mpa_table",
        "humann_genefamilies", "humann_pathabundance", "abundance_matrix",
    }
    ready = []
    if kinds & {"kraken2_report", "bracken_report", "metaphlan_profile", "mpa_table"}:
        ready += ["mg_merge"]
        if kinds & {"kraken2_report", "bracken_report", "mpa_table"}:
            ready += ["mg_rarefaction"]
        ready += ["mg_diversity", "mg_differential", "mg_plot"]
    if "abundance_matrix" in kinds:
        ready += ["mg_diversity", "mg_differential", "mg_plot"]
    if kinds & {"humann_genefamilies", "humann_pathabundance"}:
        ready += ["mg_function", "mg_plot"]
    if "bin_quality_table" in kinds:
        ready += ["mg_bin_qc (mag-workflow)", "mg_mag_table (mag-workflow)"]
    if "quast_report" in kinds:
        ready += ["mg_assembly_stats (mag-workflow)"]
    emit({
        "input": str(root),
        "input_type": sorted(kinds)[0] if len(kinds) == 1 else "mixed",
        "n_files": len(files),
        "files": files,
        "samples": sorted({f["sample"] for f in files if f["kind"] in profile_kinds}),
        "ready_steps": list(dict.fromkeys(ready)),
        "notes": [] if ready else ["no recognized profiling input; see references/step-notes/inspect.md"],
    })


if __name__ == "__main__":
    main()