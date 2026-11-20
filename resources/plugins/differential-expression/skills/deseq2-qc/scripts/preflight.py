"""Preflight: validate user-described design against sample metadata.

Never infers design, never echoes sample/donor ids or row content.
Compact JSON on stdout; exit 0 for both ok and needs_input.
"""

import argparse
import csv
import io
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deqclib as lib

MISSING_TOKENS = {"", "na", "n/a", "nan", "null"}

DESIGN_OPTIONS = [
    {"label": "Paired samples", "description": "The same donors contributed samples to both groups. Provide the sample-to-donor mapping."},
    {"label": "Independent samples", "description": "Each biological replicate comes from a different subject."},
]
DESCRIPTIVE = {"label": "Descriptive summary only", "description": "Report the available data and stop before differential testing or enrichment."}


def ask(reasons, header, question=None, options=None):
    if header in ("Donor mapping", "Donor check", "Donor conflict", "Unbalanced pairs"):
        question = "Can you supply or correct the sample-to-donor mapping for this comparison?"
        options = [{"label": "Provide donor mapping", "description": "Supply the metadata column or corrected table and recheck the design before testing."}]
    elif header == "Replication":
        question = "How would you like to address the replication required by this design?"
        options = [{"label": "Provide replicates", "description": "Supply the independent biological replicates or complete donor pairs required for this contrast."}]
    else:
        options = [{"label": o["label"].replace("_", " ").capitalize(), "description": o["description"]} for o in options]
    lib.needs_input(reasons, {"header": header, "question": question, "options": options[:3] + [DESCRIPTIVE]})


def is_missing(value):
    return value.strip().lower() in MISSING_TOKENS


def main():
    ap = argparse.ArgumentParser(description="Design-aware metadata preflight")
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--design", required=True, choices=["independent", "paired", "unknown"])
    ap.add_argument("--sample_id_col", required=True)
    ap.add_argument("--condition_col", required=True)
    ap.add_argument("--reference_level", required=True)
    ap.add_argument("--target_level", required=True)
    ap.add_argument("--donor_col", default="")
    args = ap.parse_args()

    try:
        with open(args.metadata, "rb") as fh:
            data = fh.read()
    except OSError:
        ask(["metadata file not found or unreadable"], "Metadata",
            "Where is the sample metadata file?",
            [{"label": "provide_path", "description": "Provide the correct path to the sample metadata file."},
             {"label": "reexport", "description": "Re-export the metadata as UTF-8 TSV or CSV and share the new path."}])
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        ask(["metadata is not valid UTF-8"], "Metadata",
            "Can you re-export the metadata as UTF-8?",
            [{"label": "reexport", "description": "Re-export the metadata as UTF-8 TSV or CSV."}])

    ext = os.path.splitext(args.metadata)[1].lower()
    first_line = text.split("\n", 1)[0]
    if ext == ".csv":
        delim = ","
    elif ext in (".tsv", ".tab"):
        delim = "\t"
    elif "\t" in first_line:
        delim = "\t"
    elif "," in first_line:
        delim = ","
    else:
        ask(["cannot detect TSV/CSV delimiter in metadata"], "Metadata",
            "Is the metadata tab-separated or comma-separated?",
            [{"label": "tsv", "description": "Tab-separated; re-export with a .tsv extension."},
             {"label": "csv", "description": "Comma-separated; re-export with a .csv extension."}])

    try:
        rows = list(csv.reader(io.StringIO(text), delimiter=delim, strict=True))
    except csv.Error:
        ask(["metadata is malformed (CSV quoting error)"], "Metadata",
            "Can you provide a cleanly exported metadata file?",
            [{"label": "reexport", "description": "Re-export the metadata with standard CSV/TSV quoting."}])
    if len(rows) < 2:
        ask(["metadata needs a header and at least one sample record"], "Metadata",
            "Can you provide complete sample metadata?",
            [{"label": "provide_complete", "description": "Provide metadata with a header row and one record per sample."}])

    header = rows[0]
    if any(h == "" for h in header):
        ask(["metadata header contains an unnamed column"], "Metadata",
            "Can you name every metadata column?",
            [{"label": "name_columns", "description": "Give every column a unique, non-empty name."}])
    if len(set(header)) != len(header):
        ask(["duplicate column names in metadata header"], "Metadata",
            "Can you rename columns so each name is unique?",
            [{"label": "rename_columns", "description": "Rename duplicated columns and re-export the metadata."}])
    for col in (args.sample_id_col, args.condition_col):
        if col not in header:
            ask(["requested column %r not found in metadata header" % col], "Columns",
                "Which metadata columns hold the sample ids and group labels?",
                [{"label": "fix_columns", "description": "Re-check the exact column names and pass them as "
                                                         "--sample_id_col / --condition_col."}])
    si = header.index(args.sample_id_col)
    ci = header.index(args.condition_col)

    n_bad = n_missing = n_dup = 0
    seen = set()
    samples = []
    for row in rows[1:]:
        if len(row) != len(header):
            n_bad += 1
            continue
        sid = row[si].strip()
        cond = row[ci]
        if is_missing(sid) or is_missing(cond):
            n_missing += 1
            continue
        if sid in seen:
            n_dup += 1
            continue
        seen.add(sid)
        samples.append((sid, cond, row))
    if n_bad:
        ask(["%d metadata record(s) have a field count that does not match the header" % n_bad], "Metadata",
            "Can you fix the malformed records?",
            [{"label": "fix_rows", "description": "Correct the records so every row has exactly the header's field count."}])
    if n_missing:
        ask(["missing sample id or group label in %d record(s)" % n_missing], "Metadata",
            "Can you complete the missing sample ids or group labels?",
            [{"label": "fill_values", "description": "Fill in every sample id and group label."},
             {"label": "remove_incomplete", "description": "Remove incomplete records only if they are truly out of scope."}])
    if n_dup:
        ask(["duplicate sample ids in %d record(s)" % n_dup], "Metadata",
            "Are the repeated sample ids technical replicates or labeling errors?",
            [{"label": "fix_ids", "description": "Make sample ids unique."},
             {"label": "collapse_replicates", "description": "Collapse technical replicates only with explicit user confirmation."}])

    reference = args.reference_level
    target = args.target_level
    if not reference or not target or reference == target:
        ask(["reference_level and target_level must be two distinct non-empty group labels"], "Contrast",
            "Which two distinct groups define the contrast?",
            [{"label": "fix_levels", "description": "Provide distinct --reference_level and --target_level exactly as written in the metadata."}])

    if args.design == "unknown":
        ask(["study design has not been described by the user; it is never inferred from data"],
            "Study design",
            "Did the same biological donors contribute samples to both groups, or are the groups independent subjects?",
            DESIGN_OPTIONS)

    group_counts = Counter(cond for _, cond, _ in samples)
    for role, level in (("reference", reference), ("target", target)):
        if level not in group_counts:
            ask(["%s level not found in the condition column (exact string match required)" % role], "Contrast",
                "What are the exact group labels as written in the metadata?",
                [{"label": "fix_levels", "description": "Re-check spelling, case and whitespace; numeric-looking labels are matched as exact strings."},
                 {"label": "describe_groups", "description": "Describe the intended groups and how they are labeled in the metadata."}])

    contrast_samples = [(sid, cond, row) for sid, cond, row in samples if cond in (reference, target)]
    warnings = []

    if args.design == "paired":
        if not args.donor_col:
            ask(["paired design requires a donor mapping (donor_col)"], "Donor mapping")
        if args.donor_col not in header:
            ask(["donor column %r not found in metadata header" % args.donor_col], "Donor mapping")
        di = header.index(args.donor_col)
        n_missing_donor = sum(1 for _, _, row in contrast_samples if is_missing(row[di]))
        if n_missing_donor:
            ask(["missing donor ids for %d contrast sample(s)" % n_missing_donor], "Donor mapping")
        per_donor = {}
        for _, cond, row in contrast_samples:
            per_donor.setdefault(row[di].strip(), Counter())[cond] += 1
        n_within_dup = sum(1 for c in per_donor.values() if c[target] > 1 or c[reference] > 1)
        if n_within_dup:
            ask(["%d donor(s) appear more than once within a contrast group" % n_within_dup], "Donor mapping")
        n_pairs = sum(1 for c in per_donor.values() if c[target] == 1 and c[reference] == 1)
        n_unpaired = len(per_donor) - n_pairs
        if n_unpaired:
            ask(["unbalanced donor groups: %d complete pair(s), %d donor(s) present in only one group; "
                 "no donors are excluded silently" % (n_pairs, n_unpaired)], "Unbalanced pairs")
        if n_pairs < 2:
            ask(["insufficient pairs: %d complete pair(s); paired DESeq2 needs at least 2 (recommend >= 3)" % n_pairs],
                "Replication")
        if n_pairs < 3:
            warnings.append("only %d complete pairs; dispersion estimates will be unstable (recommend >= 3)" % n_pairs)
        params = {"study_design": "paired", "sample_id_col": args.sample_id_col,
                  "contrast_variable": args.condition_col, "reference_level": reference,
                  "target_level": target, "donor_col": args.donor_col}
        contrast_counts = {target: n_pairs, reference: n_pairs}
        extra = {"n_pairs": n_pairs}
    else:
        for role, level in (("reference", reference), ("target", target)):
            if group_counts[level] < 2:
                ask(["insufficient distinct independent replicates in the %s group: %d (need >= 2, recommend >= 3)"
                     % (role, group_counts[level])], "Replication")
        for role, level in (("reference", reference), ("target", target)):
            if group_counts[level] < 3:
                warnings.append("%s group has only %d replicates (recommend >= 3)" % (role, group_counts[level]))
        if args.donor_col:
            if args.donor_col not in header:
                ask(["donor column %r not found in metadata header" % args.donor_col], "Donor check")
            di = header.index(args.donor_col)
            donor_values = [row[di].strip() for _, _, row in contrast_samples]
            if any(is_missing(v) for v in donor_values):
                ask(["missing donor ids in contrast samples; donor uniqueness cannot be verified"], "Donor check")
            if len(set(donor_values)) != len(donor_values):
                ask(["donor values repeat across contrast samples; an independent design requires unique donors across both groups"],
                    "Donor conflict")
        params = {"study_design": "independent", "sample_id_col": args.sample_id_col,
                  "contrast_variable": args.condition_col, "reference_level": reference,
                  "target_level": target}
        contrast_counts = {target: group_counts[target], reference: group_counts[reference]}
        extra = {}

    result = {"status": "ok", "study_design": args.design, "wrapper_params": params,
              "contrast_counts": contrast_counts, "warnings": warnings}
    result.update(extra)
    lib.emit(result, 0)


if __name__ == "__main__":
    main()