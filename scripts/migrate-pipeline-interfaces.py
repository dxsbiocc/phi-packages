#!/usr/bin/env python3
"""Adapt vendored pipelines to the shared HEAD component interfaces.

After dedupe_subworkflows() retargets pipeline includes to the shared tree, the
pinned pipeline code still references the old component interfaces. This script
applies the mechanical, behaviour-preserving fixes:

  1. Drop ``.mix(X.out.versions)`` segments where X is a shared module/subworkflow:
     HEAD components publish versions via ``topic: versions`` (collected by the
     pipeline's softwareVersionsToYAML/channel.topic collation), not a ``versions``
     emit. Legacy local-module versions mixes are kept.
  2. Rename ``.out.bai`` / ``.out.tbi`` to ``.out.index`` where the shared target
     renamed its index emit.

Usage:
  python3 scripts/migrate-pipeline-interfaces.py --dest resources/wrappers [--dry-run]
"""
import argparse
import re
from pathlib import Path

INC = re.compile(r"include\s*\{([^}]*)\}\s*from\s*['\"]([^'\"]+)['\"]")
COMMENT = re.compile(r"//[^\n]*")


def resolve(base, p):
    b = base / p
    for c in (b, b.with_name(b.name + ".nf"), b / "main.nf"):
        if c.is_file():
            return c.resolve()
    return None


def output_names(f):
    return set(re.findall(r"emit:\s*(\w+)", COMMENT.sub("", f.read_text())))


def emit_names(f):
    text = COMMENT.sub("", f.read_text())
    names = set()
    for m in re.finditer(r"\bworkflow\s+(\w+)\s*\{", text):
        i, d = m.end(), 1
        while i < len(text) and d:
            if text[i] == "{":
                d += 1
            elif text[i] == "}":
                d -= 1
            i += 1
        em = re.search(r"\bemit\s*:(.*)", text[m.end():i - 1], re.S)
        if em:
            for line in em.group(1).splitlines():
                nm = re.match(r"\s*(\w+)\s*(?:=|$)", line)
                if nm and nm.group(1) != "}":
                    names.add(nm.group(1))
    return names


def alias_targets(f, tree):
    """alias -> (kind, names) for shared modules/subworkflows only."""
    out = {}
    for m in INC.finditer(f.read_text()):
        tgt = resolve(f.parent, m.group(2))
        if tgt is None:
            continue
        rel = str(tgt.relative_to(tree))
        if rel.startswith("modules/nf-core"):
            kind, names = "module", output_names(tgt)
        elif rel.startswith("subworkflows/nf-core"):
            kind, names = "sw", emit_names(tgt)
        else:
            continue
        for a in m.group(1).split(","):
            alias = a.strip().split(" as ")[-1].strip()
            out[alias] = (kind, names)
    return out


def migrate_file(f, tree, report):
    targets = alias_targets(f, tree)
    text = f.read_text()
    orig = text
    no_versions = [a for a, (k, n) in targets.items() if "versions" not in n]
    rename_index = [a for a, (k, n) in targets.items()
                    if "index" in n and "bai" not in n and "tbi" not in n]

    for alias in no_versions:
        # .mix(ALIAS.out.versions[.first()|.unique{...}])
        pat = re.compile(r"\.mix\(\s*" + re.escape(alias) +
                         r"\.out\.versions(?:\.first\(\)|\.unique\{[^}]*\})?\s*\)")
        text, n = pat.subn("", text)
        report["versions_mixes_removed"] += n
        # versions = ALIAS.out.versions  (direct emit assignment)
        pat2 = re.compile(r"(versions\s*=\s*)" + re.escape(alias) + r"\.out\.versions\b")
        text, n2 = pat2.subn(r"\1channel.empty()", text)
        report["versions_assigns_fixed"] += n2

    for alias in rename_index:
        for old in ("bai", "tbi"):
            pat = re.compile(r"\b" + re.escape(alias) + r"\.out\." + old + r"\b")
            text, n = pat.subn(alias + ".out.index", text)
            report["index_renames"] += n

    # drop no-op lines left behind in workflow bodies: `ch_versions = ch_versions`
    # or bare `ch_versions`. Never touch take:/emit: declarations (a bare
    # `ch_versions` there is a channel declaration, not a no-op).
    lines = []
    section = None
    for line in text.splitlines(keepends=True):
        sec = re.match(r"\s*(take|emit|main|script|shell|exec|input|output|when|stub)\s*:", line)
        if sec:
            section = sec.group(1)
        stripped = COMMENT.sub("", line).strip()
        if section == "main" and stripped in ("ch_versions = ch_versions", "ch_versions"):
            report["noop_lines_removed"] += 1
            continue
        lines.append(line)
    text = "".join(lines)

    if text != orig:
        report["files"].append(str(f.relative_to(tree)))
        if not report["dry_run"]:
            f.write_text(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", type=Path, required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    tree = args.dest.resolve()
    report = {"files": [], "versions_mixes_removed": 0, "versions_assigns_fixed": 0,
              "index_renames": 0, "noop_lines_removed": 0, "dry_run": args.dry_run}
    for pipeline in sorted((tree / "workflows" / "nf-core").iterdir()):
        for nf in sorted(pipeline.rglob("*.nf")):
            if "/tests/" in str(nf):
                continue
            migrate_file(nf, tree, report)
    print(f"files changed: {len(report['files'])}")
    print(f"versions mixes removed: {report['versions_mixes_removed']}, "
          f"assigns fixed: {report['versions_assigns_fixed']}, "
          f"index renames: {report['index_renames']}, "
          f"no-op lines removed: {report['noop_lines_removed']}")
    for x in report["files"]:
        print("  ", x)


if __name__ == "__main__":
    main()