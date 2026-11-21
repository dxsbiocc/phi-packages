#!/usr/bin/env python3
"""Vendor complete nf-core pipelines into the phi wrapper tree.

Follows the recipe proven by workflows/rna-seq:
  - keep the pipeline skeleton (main.nf, conf/, assets/, bin/, workflows/,
    subworkflows/, modules/local/, schema, configs) intact;
  - by default SELECTIVELY dedupe: a pipeline-pinned module/subworkflow whose
    files are identical to the shared tree's copy is retargeted to the shared
    tree and its private copy dropped; on ANY version conflict the pipeline's
    pinned copy wins and stays private (subworkflows only dedupe when their
    whole include closure also dedupes, so no HEAD drift leaks in);
  - with --dedupe, dedupe ALL modules/nf-core AND subworkflows/nf-core into the
    shared tree regardless of version (requires interface migration, see git
    history);
  - patch $projectDir references in .nf files one level up, because the
    wrapper entrypoint runs from wrapper/ (projectDir = wrapper/);
  - wrapper/main.nf is the pipeline main.nf with './' includes lifted to
    '../'; wrapper/nextflow.config includeConfig's the real pipeline config;
  - wrapper/params.json explicitly overrides schema defaults that embed
    ${projectDir} (nf-schema re-derives those from the schema literal).

Usage:
  python3 scripts/vendor-pipelines.py --source <dir-with-pipeline-folders> --dest resources/wrappers [--only atacseq] [--dry-run]
"""

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import yaml

KEEP_FILES = ["main.nf", "nextflow.config", "nextflow_schema.json", "modules.json",
              "README.md", "LICENSE", "CITATIONS.md", "CHANGELOG.md"]
KEEP_DIRS = ["conf", "assets", "bin", "workflows", "subworkflows"]
INCLUDE_RE = re.compile(r"(include\s*\{[^}]*\}\s*from\s*['\"])([^'\"]+)(['\"])")
PROJECTDIR_RE = re.compile(r"(\$\{projectDir\}|\$projectDir)/")
GROOVY_PROJECTDIR_RE = re.compile(r"(\$\{projectDir\}|\$projectDir|\$\{workflow\.projectDir\})/")


def rel_include(from_file, target_abs, tree_root):
    import os
    return os.path.relpath(target_abs, from_file.parent)


def rewrite_includes(dest_root, src_root, tree_root, report, retarget=True):
    """Patch $projectDir references one level up; with retarget=True also rewrite
    includes that pointed at src_root/modules/nf-core to the shared tree."""
    shared_modules = tree_root / "modules" / "nf-core"
    for nf in sorted(dest_root.rglob("*.nf")):
        rel = nf.relative_to(dest_root)
        orig = src_root / rel
        text = nf.read_text()

        def sub(m):
            if not retarget:
                return m.group(0)
            pre, path, post = m.groups()
            resolved = (orig.parent / path).resolve()
            try:
                inside = resolved.relative_to((src_root / "modules" / "nf-core").resolve())
            except ValueError:
                return m.group(0)
            target = shared_modules / inside
            if not target.exists() and not target.with_name(target.name + ".nf").exists():
                report["missing_modules"].append(f"{rel}: {inside}")
            new = rel_include(nf, target, tree_root)
            return f"{pre}{new}{post}"

        new_text = INCLUDE_RE.sub(sub, text)
        new_text = PROJECTDIR_RE.sub(r"\1/../", new_text)
        if new_text != text:
            nf.write_text(new_text)
            report["rewritten_files"].append(str(rel))


def fill_missing_modules(src_root, tree_root, report):
    """Copy pipeline-pinned module versions the shared tree lacks (additive only)."""
    shared = tree_root / "modules" / "nf-core"
    for entry in sorted(set(report["missing_modules"])):
        rel_mod = entry.split(": ", 1)[1]
        src = src_root / "modules" / "nf-core" / rel_mod
        dst = shared / rel_mod
        if src.exists() and not dst.exists():
            shutil.copytree(src, dst)
            report["filled_modules"].append(str(rel_mod))


def dedupe_subworkflows(dest_root, src_root, tree_root, report):
    """Retarget includes that pointed at the pipeline's private subworkflows/nf-core
    copy to the shared tree, then drop the private copy. Idempotent: only touches
    include paths, never projectDir patches. Fills shared gaps additively."""
    shared = tree_root / "subworkflows" / "nf-core"
    private_src = (src_root / "subworkflows" / "nf-core").resolve()
    missing = []
    for nf in sorted(dest_root.rglob("*.nf")):
        rel = nf.relative_to(dest_root)
        if rel.parts[:2] == ("subworkflows", "nf-core"):
            continue
        orig = src_root / rel
        text = nf.read_text()

        def sub(m):
            pre, path, post = m.groups()
            resolved = (orig.parent / path).resolve()
            try:
                inside = resolved.relative_to(private_src)
            except ValueError:
                return m.group(0)
            target = shared / inside
            if not target.exists() and not target.with_name(target.name + ".nf").exists():
                missing.append(f"{rel}: {inside}")
            new = rel_include(nf, target, tree_root)
            return f"{pre}{new}{post}"

        new_text = INCLUDE_RE.sub(sub, text)
        if new_text != text:
            nf.write_text(new_text)
            report["rewritten_files"].append(str(rel))
    for entry in sorted(set(missing)):
        rel_sw = entry.split(": ", 1)[1]
        src = src_root / "subworkflows" / "nf-core" / rel_sw
        dst = shared / rel_sw
        if src.exists() and not dst.exists():
            shutil.copytree(src, dst)
            report["filled_subworkflows"].append(str(rel_sw))
        else:
            report["missing_subworkflows"].append(entry)
    private_dst = dest_root / "subworkflows" / "nf-core"
    if private_dst.exists():
        report["dropped_subworkflows"] = sorted(p.name for p in private_dst.iterdir() if p.is_dir())
        shutil.rmtree(private_dst)


def component_fingerprint(comp_dir):
    """Files that define component behaviour: everything except tests/, meta.yml
    (documentation), wrapper/ (phi artifacts, absent from pipeline copies) and
    .DS_Store. Includes main.nf, environment.yml, templates/, bin/."""
    out = {}
    for f in sorted(comp_dir.rglob("*")):
        if not f.is_file() or f.name == ".DS_Store" or f.name == "meta.yml":
            continue
        rel = f.relative_to(comp_dir)
        if rel.parts[0] in ("tests", "wrapper"):
            continue
        out[str(rel)] = f.read_bytes()
    return out


def resolve_nf(base, path):
    p = (base / path).resolve()
    if p.is_dir():
        p = p / "main.nf"
    if not p.exists() and Path(str(p) + ".nf").exists():
        p = Path(str(p) + ".nf")
    return p


def selective_dedupe(dest_root, tree_root, report):
    """Dedupe only version-identical components; conflicts keep the pipeline copy.
    Modules are leaves. A subworkflow dedupes only if its files match the shared
    copy AND every include resolves to an already-deduped component (fixpoint),
    so a deduped subworkflow can never pull in a HEAD module the pipeline conflicts with.
    Note: comparison runs AFTER the $projectDir one-level-up patch, so components
    containing $projectDir references stay private — the shared copy is unpatched
    and would resolve those paths wrong under wrapper execution."""
    shared_m = (tree_root / "modules" / "nf-core").resolve()
    shared_s = (tree_root / "subworkflows" / "nf-core").resolve()
    priv_m = dest_root / "modules" / "nf-core"
    priv_s = dest_root / "subworkflows" / "nf-core"

    dedup_m = set()
    if priv_m.exists():
        for main in sorted(priv_m.rglob("main.nf")):
            rel = main.parent.relative_to(priv_m)
            shared = shared_m / rel
            if shared.exists() and component_fingerprint(main.parent) == component_fingerprint(shared):
                dedup_m.add(str(rel))
            else:
                report["kept_modules"].append(str(rel))

    dedup_s = set()
    if priv_s.exists():
        candidates = {}
        for main in sorted(priv_s.rglob("main.nf")):
            rel = main.parent.relative_to(priv_s)
            shared = shared_s / rel
            if shared.exists() and component_fingerprint(main.parent) == component_fingerprint(shared):
                candidates[str(rel)] = main
            else:
                report["kept_subworkflows"].append(str(rel))
        while True:
            progressed = False
            for rel, main in sorted(candidates.items()):
                if rel in dedup_s:
                    continue
                text = re.sub(r"//[^\n]*", "", main.read_text())
                ok = True
                for _pre, path, _post in INCLUDE_RE.findall(text):
                    if path.startswith("plugin/"):
                        continue
                    tgt = resolve_nf(main.parent, path)
                    for root, marks in ((priv_m.resolve(), dedup_m), (priv_s.resolve(), dedup_s)):
                        try:
                            inside = tgt.relative_to(root)
                        except ValueError:
                            continue
                        comp = str(inside.parent) if inside.name == "main.nf" else None
                        if comp is None or comp not in marks:
                            ok = False
                        break
                    else:
                        ok = False  # local module / outside the dedupable tree
                    if not ok:
                        break
                if ok:
                    dedup_s.add(rel)
                    progressed = True
            if not progressed:
                break
        report["kept_subworkflows"].extend(sorted(set(candidates) - dedup_s))

    if not dedup_m and not dedup_s:
        return

    for nf in sorted(dest_root.rglob("*.nf")):
        text = nf.read_text()

        def sub(m):
            pre, path, post = m.groups()
            if path.startswith("plugin/"):
                return m.group(0)
            tgt = resolve_nf(nf.parent, path)
            for priv, shared, marks in ((priv_m, shared_m, dedup_m), (priv_s, shared_s, dedup_s)):
                try:
                    inside = tgt.relative_to(priv.resolve())
                except ValueError:
                    continue
                if inside.name != "main.nf" or str(inside.parent) not in marks:
                    return m.group(0)
                new = rel_include(nf, shared / inside.parent / "main", tree_root)
                return f"{pre}{new}{post}"
            return m.group(0)

        new_text = INCLUDE_RE.sub(sub, text)
        if new_text != text:
            nf.write_text(new_text)
            report["rewritten_files"].append(str(nf.relative_to(dest_root)))

    for rel in sorted(dedup_m):
        shutil.rmtree(priv_m / rel)
    for rel in sorted(dedup_s):
        shutil.rmtree(priv_s / rel)
    for root in (priv_m, priv_s):
        if root.exists():
            for d in sorted((p for p in root.rglob("*") if p.is_dir()), key=lambda p: len(p.parts), reverse=True):
                if not any(d.iterdir()):
                    d.rmdir()
            if not any(root.iterdir()):
                root.rmdir()
    report["deduped_modules"] = sorted(dedup_m)
    report["deduped_subworkflows"] = sorted(dedup_s)


def schema_info(dest_root):
    schema_path = dest_root / "nextflow_schema.json"
    if not schema_path.exists():
        return {}, []
    schema = json.loads(schema_path.read_text())
    props, required = {}, []
    for section in schema.get("$defs", schema.get("definitions", {})).values():
        props.update(section.get("properties", {}))
        required.extend(section.get("required", []))
    props.update(schema.get("properties", {}))
    required.extend(schema.get("required", []))
    return props, sorted(set(required))


def vendor_lib(src_root, dest_root, report):
    """Copy lib/ into wrapper/lib (projectDir = wrapper/ at runtime, and Nextflow
    auto-loads projectDir/lib). Groovy asset paths are patched one level up,
    mirroring the .nf projectDir patch, so they resolve against the pipeline root."""
    src_lib = src_root / "lib"
    if not src_lib.exists():
        return
    dst_lib = dest_root / "wrapper" / "lib"
    shutil.copytree(src_lib, dst_lib, dirs_exist_ok=True, ignore=shutil.ignore_patterns(".DS_Store"))
    for groovy in sorted(dst_lib.glob("*.groovy")):
        text = groovy.read_text()
        new_text = GROOVY_PROJECTDIR_RE.sub(r"\1/../", text)
        if new_text != text:
            groovy.write_text(new_text)
    report["lib_files"] = [p.name for p in sorted(dst_lib.iterdir())]


def make_wrapper(dest_root, name, comp_id, report):
    props, required = schema_info(dest_root)
    wrapper_dir = dest_root / "wrapper"
    wrapper_dir.mkdir(exist_ok=True)

    main_text = (dest_root / "main.nf").read_text()

    def lift(m):
        pre, path, post = m.groups()
        if path.startswith("plugin/") or path.startswith("/"):
            return m.group(0)
        if path.startswith("./"):
            path = path[2:]
        return f"{pre}../{path}{post}"

    main_text = INCLUDE_RE.sub(lift, main_text)
    header = (
        "#!/usr/bin/env nextflow\n"
        f"// Wrapper entrypoint: ../main.nf verbatim with './' include paths lifted\n"
        "// to '../' (this copy lives one directory deeper). projectDir references\n"
        "// in the vendored .nf files were patched one level up for the same reason.\n"
        "// Generated by scripts/vendor-pipelines.py; regenerate, do not hand-edit.\n"
    )
    if not main_text.startswith("#!"):
        main_text = header + main_text
    (wrapper_dir / "main.nf").write_text(main_text)

    (wrapper_dir / "nextflow.config").write_text(
        "// Launcher config: pull in the real pipeline config (profiles, resource\n"
        "// labels, containers). includeConfig resolves its own relative includes\n"
        "// against ../, not this file. Generated by scripts/vendor-pipelines.py.\n"
        "includeConfig '../nextflow.config'\n\n"
        "// Conservative caps so a smoke run fits a laptop; remove for production.\n"
        "process {\n"
        "    resourceLimits = [\n"
        "        cpus: 4,\n"
        "        memory: '15.GB',\n"
        "        time: '4.h'\n"
        "    ]\n"
        "}\n"
    )

    defaults = {"outdir": "results"}
    for key, spec in props.items():
        d = spec.get("default")
        if isinstance(d, str) and "${projectDir}" in d:
            defaults[key] = d.replace("${projectDir}/", "")
            report["projectdir_overrides"].append(f"{name}: {key} = {defaults[key]}")
    test_cfg = dest_root / "conf" / "test.config"
    if test_cfg.exists():
        text = test_cfg.read_text()
        pm = re.search(r"\bparams\s*\{", text)
        if pm:
            i = text.find("{", pm.start())
            depth, j = 0, i
            while j < len(text):
                if text[j] == "{":
                    depth += 1
                elif text[j] == "}":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            for k, v in re.findall(r"^\s*(\w+)\s*=\s*['\"]([^'\"]+)['\"]", text[i + 1:j], re.M):
                if "${projectDir}" in v:
                    v = v.replace("${projectDir}/", "")
                if defaults.get(k) is None or k not in defaults:
                    defaults[k] = v

    params = {}
    summary = (props.get("input", {}) or {}).get("description") or f"nf-core/{name} pipeline."
    if "input" in props:
        params["input"] = {"kind": "input", "type": "csv", "required": True,
                           "description": (props["input"].get("description") or "Samplesheet CSV.")[:300]}
        defaults.setdefault("input", None)
    for key in ("fasta", "gtf", "gff", "genome"):
        if key in props:
            ptype = {"fasta": "fasta", "gtf": "gtf", "gff": "gff"}.get(key, "string")
            params[key] = {"kind": "input", "type": ptype, "required": key in required,
                           "description": (props[key].get("description") or key)[:300]}
            defaults.setdefault(key, props[key].get("default"))
    for key in required:
        if key in params or key not in props:
            continue
        spec = props[key]
        jtype = {"integer": "integer", "number": "number", "boolean": "boolean"}.get(spec.get("type"), "string")
        params[key] = {"kind": "option", "type": jtype, "required": True,
                       "description": (spec.get("description") or key)[:300]}
        defaults.setdefault(key, spec.get("default"))
    params["outdir"] = {"kind": "output", "type": "path", "required": True,
                        "description": "Output directory for all pipeline results."}
    for key, value in defaults.items():
        if key not in params and re.match(r"^[a-z][a-z0-9_]*$", key):
            spec = props.get(key, {})
            jtype = {True: "boolean", False: "boolean"}.get(value, None)
            if jtype is None:
                jtype = {"integer": "integer", "number": "number", "boolean": "boolean"}.get(spec.get("type"), "string")
            params[key] = {"kind": "option", "type": jtype, "required": False,
                           "description": (spec.get("description") or f"Pipeline parameter {key}.")[:300]}
    schema_desc = json.loads((dest_root / "nextflow_schema.json").read_text()).get("description", "") if (dest_root / "nextflow_schema.json").exists() else ""
    wrapper_yaml = {
        "id": comp_id,
        "name": f"nf-core/{name}"[:80],
        "summary": (schema_desc or f"Full nf-core/{name} pipeline.")[:300],
        "params": params,
        "outputs": {"results": {"type": "directory", "path": "${outdir}", "primary": True}},
    }
    comment = (
        f"# Agent-facing wrapper adapter for the vendored nf-core/{name} pipeline at ../.\n"
        "# Follows docs/contracts/wrapper.md. Params are a starter subset of\n"
        "# nextflow_schema.json; REVIEW and curate before publishing.\n"
    )
    (wrapper_dir / "wrapper.yaml").write_text(comment + yaml.safe_dump(wrapper_yaml, sort_keys=False, allow_unicode=True, width=100))
    (wrapper_dir / "params.json").write_text(json.dumps(defaults, indent=2) + "\n")

    bin_dir = dest_root / "bin"
    if bin_dir.exists():
        shutil.copytree(bin_dir, wrapper_dir / "bin", dirs_exist_ok=True)


def vendor_pipeline(src, tree_root, name, dest_rel, dry_run=False, dedupe=False):
    report = {"rewritten_files": [], "missing_modules": [], "filled_modules": [], "projectdir_overrides": [], "lib_files": [],
              "filled_subworkflows": [], "missing_subworkflows": [], "dropped_subworkflows": [],
              "deduped_modules": [], "deduped_subworkflows": [], "kept_modules": [], "kept_subworkflows": []}
    dest = tree_root / dest_rel
    if dry_run:
        report["note"] = "dry-run: nothing written"
        return report
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    for f in KEEP_FILES:
        if (src / f).exists():
            shutil.copy2(src / f, dest / f)
    for d in KEEP_DIRS:
        if (src / d).exists():
            shutil.copytree(src / d, dest / d, ignore=shutil.ignore_patterns(".DS_Store"))
    if (src / "modules" / "local").exists():
        shutil.copytree(src / "modules" / "local", dest / "modules" / "local", ignore=shutil.ignore_patterns(".DS_Store"))
    if not dedupe and (src / "modules" / "nf-core").exists():
        shutil.copytree(src / "modules" / "nf-core", dest / "modules" / "nf-core", ignore=shutil.ignore_patterns(".DS_Store"))
    rewrite_includes(dest, src, tree_root, report, retarget=dedupe)
    if dedupe:
        fill_missing_modules(src, tree_root, report)
        dedupe_subworkflows(dest, src, tree_root, report)
    else:
        selective_dedupe(dest, tree_root, report)
    slug = "rna-seq" if name == "rnaseq" else name
    comp_id = f"nf-core/workflows/{slug}"
    make_wrapper(dest, name, comp_id, report)
    vendor_lib(src, dest, report)
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--dest", type=Path, required=True, help="phi wrapper tree (resources/wrappers)")
    ap.add_argument("--only", help="pipeline name substring, e.g. atacseq")
    ap.add_argument("--dedupe", action="store_true",
                    help="retarget includes to the shared module/subworkflow tree and drop private copies (default: keep pipeline-pinned versions)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    pipelines = sorted(p for p in args.source.iterdir() if p.is_dir() and (p / "main.nf").exists())
    for src in pipelines:
        name = re.sub(r"-\d+(\.\d+)*$", "", src.name)
        if args.only and args.only not in name:
            continue
        dest_rel = f"workflows/nf-core/{'rna-seq' if name == 'rnaseq' else name}"
        report = vendor_pipeline(src, args.dest, name, dest_rel, dry_run=args.dry_run, dedupe=args.dedupe)
        print(f"== {name} -> {dest_rel}")
        print(f"   rewritten files: {len(report['rewritten_files'])}, missing modules: {len(report['missing_modules'])}, filled: {len(report['filled_modules'])}, projectDir overrides: {len(report['projectdir_overrides'])}, lib files: {len(report['lib_files'])}")
        print(f"   subworkflows dropped: {len(report['dropped_subworkflows'])}, filled: {len(report['filled_subworkflows'])}, missing: {len(report['missing_subworkflows'])}")
        print(f"   selective dedupe: modules shared={len(report['deduped_modules'])} kept={len(report['kept_modules'])}, subworkflows shared={len(report['deduped_subworkflows'])} kept={len(report['kept_subworkflows'])}")
        for m in report["missing_modules"][:10]:
            print(f"   [missing] {m}")
        for m in report["filled_modules"][:10]:
            print(f"   [filled] {m}")


if __name__ == "__main__":
    main()