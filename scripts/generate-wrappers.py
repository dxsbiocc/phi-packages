#!/usr/bin/env python3
"""Generate wrapper/ adapters for vendored nf-core modules.

Reads nf-core module sources (main.nf + meta.yml) and emits the adapter layer
defined by docs/contracts/wrapper.md: wrapper.yaml, main.nf, params.json,
nextflow.config. Existing wrapper/ directories are never touched.

Usage:
  python3 scripts/generate-wrappers.py --source <nf-core-modules-root> --dest <wrappers-tree> [--vendor] [--only tool/subcmd] [--dry-run]
  python3 scripts/generate-wrappers.py --dest <wrappers-tree> --validate

Modules whose input shape cannot be wired safely are reported for manual
review instead of being guessed at.
"""

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import yaml

PARAM_NAME_RE = re.compile(r"^[a-z][a-z0-9_]*$")
SECTION_END_RE = re.compile(r"^\s*(output|when|script|shell|exec):\s*$", re.M)
INPUT_BLOCK_RE = re.compile(r"^\s*input:\s*\n(.*?)(?=^\s*(?:output|when|script|shell|exec):\s*$)", re.M | re.S)
OUTPUT_BLOCK_RE = re.compile(r"^\s*output:\s*\n(.*?)(?=^\s*(?:when|script|shell|exec):\s*$)", re.M | re.S)
ELEM_RE = re.compile(r"(val|path|file|env|stdin)\s*\(?\s*([A-Za-z_]\w*)|(?:path|file)\s*\(?\s*[\"']([^\"']+)[\"']")
EMIT_RE = re.compile(r"emit:\s*(\w+)")
GLOB_RE = re.compile(r'path\(\s*"([^"]+)"')
LABEL_RE = re.compile(r"label\s+'(process_\w+)'")

RESOURCES = {
    "process_single": (1, "2.GB", "1.h"),
    "process_low": (2, "4.GB", "2.h"),
    "process_medium": (4, "8.GB", "4.h"),
    "process_high": (8, "16.GB", "8.h"),
    "process_long": (2, "4.GB", "20.h"),
    "process_high_memory": (8, "64.GB", "8.h"),
}

FASTQ_NAME_RE = re.compile(r"(reads|fastq|fq|fastqgz)", re.I)
REF_NAME_RE = re.compile(r"(fasta|fa$|index|indices|db$|database|reference|genome|gtf|gff|bed$|sizes|fai$|dict$|bwt|sa$|ann$|amb$|pac$|hisat2|star|salmon|kallisto)", re.I)
BOOL_NAME_RE = re.compile(r"^(save_|skip_|enable_|use_|discard_|no_|keep_|allow_|is_|has_|with_|write_|emit_|do_)")


class ReviewNeeded(Exception):
    pass


def snake(name):
    s = re.sub(r"(?<!^)(?=[A-Z])", "_", name).lower()
    s = re.sub(r"[^a-z0-9_]", "_", s)
    return s.lstrip("_0123456789") or "input"


def parse_input_lines(text):
    m = INPUT_BLOCK_RE.search(text)
    if not m:
        return []
    lines = []
    for raw in m.group(1).splitlines():
        line = raw.strip()
        if not line or line.startswith("//"):
            continue
        lines.append(line)
    return lines


def parse_channel(line):
    """Return (kind, elements) where kind is 'tuple' or 'single'."""
    kind = "tuple" if line.startswith("tuple") else "single"
    elems = []
    matches = list(ELEM_RE.finditer(line))
    for i, em in enumerate(matches):
        seg_end = matches[i + 1].start() if i + 1 < len(matches) else len(line)
        sa_m = re.search(r"stageAs:\s*['\"]([^'\"]+)['\"]", line[em.end():seg_end])
        sa = sa_m.group(1) if sa_m else None
        if em.group(1):
            etype, ename = em.group(1), em.group(2)
            if etype in ("env", "stdin"):
                raise ReviewNeeded(f"unsupported input qualifier '{etype}'")
            if etype == "file":
                etype = "path"
            elems.append({"type": etype, "name": ename, "stageAs": sa if etype == "path" else None})
        else:
            glob = em.group(3)
            seg = glob.split("/")[0]
            if not seg or "*" in seg:
                seg = glob.rsplit("/", 1)[-1]
            name = snake(re.sub(r"^[\*.]+", "", seg).replace("*", "")) or "staged_files"
            elems.append({"type": "path", "name": name, "glob": glob, "stageAs": sa})
    if not elems:
        raise ReviewNeeded(f"cannot parse input line: {line}")
    return kind, elems


def parse_outputs(text):
    m = OUTPUT_BLOCK_RE.search(text)
    if not m:
        return []
    outs = []
    for raw in m.group(1).splitlines():
        line = raw.strip()
        if not line or line.startswith("//"):
            continue
        emit = EMIT_RE.search(line)
        if not emit:
            continue
        name = emit.group(1)
        if name.startswith("versions") or "topic: versions" in line:
            continue
        globs = GLOB_RE.findall(line)
        outs.append({"name": name, "globs": globs, "optional": "optional:true" in line.replace(" ", "")})
    return outs


def meta_io_map(meta, section):
    """Flatten meta.yml input/output entries to {variable_name: info}."""
    result = {}

    def walk(node):
        if isinstance(node, list):
            for item in node:
                walk(item)
        elif isinstance(node, dict):
            for key, val in node.items():
                if isinstance(val, dict) and ("type" in val or "description" in val):
                    result[key] = val
                else:
                    walk(val)

    walk(meta.get(section, []))
    return result


def first_sentence(desc, limit=300):
    desc = " ".join((desc or "").split())
    if not desc:
        return ""
    cut = re.split(r"(?<=\.)\s", desc, maxsplit=1)[0]
    return cut[:limit]


def test_dataset_defaults(module_dir):
    """Pull https test-dataset URLs used by the module's nf-test file."""
    urls = {}
    test_file = module_dir / "tests" / "main.nf.test"
    if not test_file.exists():
        return urls
    text = test_file.read_text()
    for m in re.finditer(r'(?:path|input|reads|fasta|bam|vcf|gtf)\s*[:=]?\s*\(?\s*["\'](https://[^"\']+)["\']', text):
        urls.setdefault(m.group(1), None)
    return list(urls)


def is_fastq_input(elem, meta_info, process_text):
    name = elem["name"]
    pattern = (meta_info.get(name, {}) or {}).get("pattern", "")
    if FASTQ_NAME_RE.search(name) and ("single_end" in process_text or "fastq" in pattern or "fq" in pattern):
        return True
    return "fastq" in pattern or "fq.gz" in pattern


def yaml_dump(data):
    return yaml.safe_dump(data, sort_keys=False, allow_unicode=True, width=100)


def generate_module(src_dir, dest_dir, rel, dry_run=False):
    """rel is like 'macs2/callpeak'. Returns (status, detail)."""
    main_nf = src_dir / "main.nf"
    meta_yml = src_dir / "meta.yml"
    text = main_nf.read_text()
    if re.search(r"assert\s+false", text):
        return "skipped-deprecated", "module asserts false (deprecated upstream)"
    if not meta_yml.exists():
        return "review", "no meta.yml"
    meta = yaml.safe_load(meta_yml.read_text()) or {}
    pm = re.search(r"process\s+(\w+)\s*\{", text)
    if not pm:
        return "review", "no process block found"
    process = pm.group(1)
    proc_dir = process.lower()

    try:
        lines = parse_input_lines(text)
        channels = [parse_channel(l) for l in lines]
    except ReviewNeeded as e:
        return "review", str(e)

    meta_in = meta_io_map(meta, "input")
    outputs = parse_outputs(text)
    if not outputs:
        return "review", "no declared outputs (emit)"

    flags = []
    params = {}       # wrapper.yaml params

    def claim(base):
        name, n = base, 2
        while name in params:
            name = f"{base}_{n}"
            n += 1
        return name

    defaults = {}     # params.json
    wiring = []       # main.nf channel construction statements
    call_args = []    # process invocation arguments
    test_urls = test_dataset_defaults(src_dir)

    tuple_channels = [(i, c) for i, (k, c) in enumerate(channels) if k == "tuple"]
    if not tuple_channels:
        main_idx = None
    else:
        main_idx = tuple_channels[0][0]
        if len(tuple_channels) > 1:
            flags.append("multiple tuple channels: secondary channels wired as single references (.first()); verify pairing semantics")

    for idx, (kind, elems) in enumerate(channels):
        for e in elems:
            sa = e.get("stageAs")
            if not sa:
                continue
            if e["type"] == "path" and re.match(r"^[\w-]+/\*$", sa):
                e["glob"] = sa
                e["staged_dir"] = True
            elif e["type"] == "path" and re.match(r"^[\w./-]+/$", sa):
                e["glob"] = sa + "*"
                e["staged_dir"] = True
            elif sa in ("?/*", "*") or re.match(r"^\*\.\w+$", sa):
                e["stageAs"] = None  # per-file staging; plain fromPath is equivalent
            elif FASTQ_NAME_RE.search(e["name"]) and re.match(r"^input\*?/\*$", sa):
                e["stageAs"] = None  # fromFilePairs staging is equivalent
            else:
                return "review", f"stageAs '{sa}' needs manual wiring"
        if kind == "tuple":
            has_meta = elems[0]["type"] == "val" and elems[0]["name"].startswith("meta")
            if not has_meta:
                flags.append("tuple has no leading val(meta); sample meta synthesized from file name")
            file_elems = [e for e in elems if e["type"] == "path"]
            val_elems = [e for e in elems[1:] if e["type"] == "val"] if has_meta else [e for e in elems if e["type"] == "val"]
            if not file_elems:
                if not val_elems:
                    return "review", "tuple channel without path or val elements"
                vnames = []
                for ve in val_elems:
                    vname = claim(snake(ve["name"]))
                    add_val_param(vname, ve, meta_in, params, defaults, flags)
                    vnames.append(vname)
                vals_expr = ", ".join(f"params.{v}" for v in vnames)
                call_args.append(f"Channel.value([[id: 'input'], {vals_expr}])")
                flags.append("val-only tuple wired as a single sample; verify multi-sample semantics")
                continue
            is_main = idx == main_idx
            primary_elem = file_elems[0]
            pname = claim(snake(primary_elem["name"]))
            fastq = is_fastq_input(primary_elem, meta_in, text) and is_main
            desc = first_sentence((meta_in.get(primary_elem["name"], {}) or {}).get("description", "")) or f"Input {primary_elem['name']}."
            if primary_elem.get("staged_dir"):
                ptype = "path"
                desc = desc or f"Directory whose contents are staged as {primary_elem['glob']}."
            elif primary_elem.get("glob"):
                ptype = "path_glob"
                desc = desc or f"Files matching {primary_elem['glob']}."
            elif fastq:
                ptype = "fastq_glob"
                desc = desc or "FASTQ files (glob matching paired- or single-end reads)."
            else:
                ptype = "path"
            params[pname] = {"kind": "input", "type": ptype, "required": True, "description": desc}
            defaults[pname] = test_urls[0] if test_urls and is_main else None

            ch_name = f"{pname}_ch"
            if primary_elem.get("glob"):
                src = f'"${{params.{pname}}}/*"' if primary_elem.get("staged_dir") else f"params.{pname}"
                wiring.append(
                    f"    {ch_name} = Channel\n"
                    f"        .fromPath({src}, checkIfExists: true)\n"
                    f"        .collect()\n"
                    f"        .map {{ files -> [[id: 'input'], files] }}"
                )
                flags.append(f"glob input '{pname}' collected into one staged list; verify staging layout")
            elif fastq:
                wiring.append(
                    f"    {ch_name} = Channel\n"
                    f"        .fromFilePairs(params.{pname}, size: -1)\n"
                    f"        .map {{ sample, reads -> [[id: sample, single_end: reads instanceof List ? reads.size() == 1 : true], reads] }}"
                )
            else:
                wiring.append(
                    f"    {ch_name} = Channel\n"
                    f"        .fromPath(params.{pname}, checkIfExists: true)\n"
                    f"        .map {{ f -> [[id: f.simpleName], f] }}"
                )
            tuple_expr = ch_name
            # companion path elements inside the same tuple: optional, combined
            for ce in file_elems[1:]:
                if ce.get("glob"):
                    return "review", f"companion input '{ce['name']}' with glob/stageAs staging needs manual wiring"
                cpname = claim(snake(ce["name"]))
                cdesc = first_sentence((meta_in.get(ce["name"], {}) or {}).get("description", "")) or f"Optional companion input {ce['name']}."
                params[cpname] = {"kind": "input", "type": "path", "required": False, "description": cdesc}
                defaults[cpname] = None
                opt_ch = f"{cpname}_ch"
                wiring.append(
                    f"    {opt_ch} = params.{cpname} ? Channel.fromPath(params.{cpname}, checkIfExists: true).map {{ f -> [f] }} : Channel.value([[]])"
                )
                tuple_expr = f"{tuple_expr}.combine({opt_ch})"
                flags.append(f"companion input '{cpname}' combined cartesian-style; verify if per-sample pairing is required")
            for ve in val_elems:
                vname = claim(snake(ve["name"]))
                add_val_param(vname, ve, meta_in, params, defaults, flags)
                tuple_expr = f"{tuple_expr}.map {{ t -> t + [params.{vname}] }}"
            if not is_main:
                tuple_expr = f"{tuple_expr}.first()"
            call_args.append(tuple_expr)
        else:
            elem = elems[0]
            ename = claim(snake(elem["name"]))
            if elem["type"] == "val":
                add_val_param(ename, elem, meta_in, params, defaults, flags)
                call_args.append(f"params.{ename}")
            else:
                desc = first_sentence((meta_in.get(elem["name"], {}) or {}).get("description", "")) or f"Input {elem['name']}."
                if elem.get("glob"):
                    if elem.get("staged_dir"):
                        params[ename] = {"kind": "input", "type": "path", "required": True, "description": desc or f"Directory whose contents are staged as {elem['glob']}."}
                        src = f'"${{params.{ename}}}/*"'
                    else:
                        params[ename] = {"kind": "input", "type": "path_glob", "required": True, "description": desc or f"Files matching {elem['glob']}."}
                        src = f"params.{ename}"
                    defaults[ename] = None
                    wiring.append(
                        f"    {ename}_ch = Channel.fromPath({src}, checkIfExists: true).collect()"
                    )
                    flags.append(f"glob input '{ename}' collected into one staged list; verify staging layout")
                else:
                    required = not REF_NAME_RE.search(ename)
                    params[ename] = {"kind": "input", "type": "path", "required": required, "description": desc}
                    defaults[ename] = None
                    if REF_NAME_RE.search(ename):
                        wiring.append(
                            f"    {ename}_ch = params.{ename} ? Channel.fromPath(params.{ename}, checkIfExists: true).first() : Channel.value([])"
                        )
                        flags.append(f"standalone path '{ename}' treated as single reference (.first()); verify if a list/collection is expected")
                    else:
                        wiring.append(
                            f"    {ename}_ch = Channel.fromPath(params.{ename}, checkIfExists: true)"
                        )
                call_args.append(f"{ename}_ch")

    # generic args option when the module reads task.ext.args
    uses_ext_args = "task.ext.args" in text
    if uses_ext_args:
        params["args"] = {
            "kind": "option", "type": "string", "required": False,
            "description": f"Extra command-line arguments passed to {process.split('_')[0].lower()}.",
        }
        defaults["args"] = ""

    params["outdir"] = {
        "kind": "output", "type": "path", "required": True,
        "description": "Output directory for the module's results.",
    }
    defaults["outdir"] = "results"

    tool, _, subcmd = rel.partition("/")
    comp_id = f"nf-core/modules/{snake(tool).replace('_', '-')}" + (f"-{snake(subcmd).replace('_', '-')}" if subcmd else "")
    name = tool + (f" {subcmd}" if subcmd else "")
    summary = first_sentence(meta.get("description", "")) or f"Run {name}."
    primary_out = outputs[0]["name"]
    out_globs = sorted({g for o in outputs for g in o["globs"] if "*" in g})
    pattern = ",".join(out_globs) if out_globs else "*"
    if len(out_globs) > 1:
        pattern = "{" + ",".join(out_globs) + "}" if all("." in g for g in out_globs) else "*"

    wrapper_yaml = {
        "id": comp_id,
        "name": name[:80],
        "summary": summary[:300],
        "params": params,
        "outputs": {
            primary_out: {"type": "directory", "path": "${outdir}/" + proc_dir, "primary": True},
        },
    }

    header = (
        "#!/usr/bin/env nextflow\n"
        "// Thin agent-facing adapter over the vendored module at ../main.nf.\n"
        "// See docs/contracts/wrapper.md. Generated by scripts/generate-wrappers.py;\n"
        "// review flagged semantics before publishing.\n"
        "nextflow.enable.dsl = 2\n\n"
        f"include {{ {process} }} from '../main.nf'\n\n"
    )
    param_decls = "\n".join(f"params.{k} = {json.dumps(v)}" if v is not None else f"params.{k} = null" for k, v in defaults.items())
    body = "workflow {\n" + "\n".join(wiring) + f"\n\n    {process}({', '.join(call_args)})\n}}\n"
    main_nf_out = header + param_decls + "\n\n" + body

    label = (LABEL_RE.search(text) or [None, "process_medium"])[1] if LABEL_RE.search(text) else "process_medium"
    cpus, memory, time = RESOURCES.get(label, RESOURCES["process_medium"])
    ext_args_line = "\n        ext.args = { params.args ?: '' }" if uses_ext_args else ""
    config = f"""// Launcher config for the {name} wrapper adapter. Not part of the fixed
// wrapper/{{main.nf,params.json,wrapper.yaml}} contract — Nextflow picks this
// up automatically because it sits beside wrapper/main.nf.

process {{
    cpus   = {cpus}
    memory = '{memory}'
    time   = '{time}'

    withName: '{process}' {{
        publishDir = [
            path: {{ "${{params.outdir}}/{proc_dir}" }},
            mode: 'copy',
            pattern: '{pattern}'
        ]{ext_args_line}
    }}
}}

profiles {{
    docker {{
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
        conda.enabled = false
        singularity.enabled = false
    }}
    singularity {{
        singularity.enabled = true
        singularity.autoMounts = true
        docker.enabled = false
        conda.enabled = false
    }}
    conda {{
        conda.enabled = true
        docker.enabled = false
        singularity.enabled = false
    }}
}}
"""

    if dry_run:
        return ("review" if flags else "ok"), "; ".join(flags) or "dry-run"

    wrapper_dir = dest_dir / "wrapper"
    wrapper_dir.mkdir(parents=True, exist_ok=True)
    comment = "# Agent-facing wrapper adapter for the vendored module at ../main.nf.\n# Follows docs/contracts/wrapper.md. Generated; review flagged semantics.\n"
    if flags:
        comment += "# REVIEW: " + "; ".join(flags) + "\n"
    (wrapper_dir / "wrapper.yaml").write_text(comment + yaml_dump(wrapper_yaml))
    (wrapper_dir / "main.nf").write_text(main_nf_out)
    (wrapper_dir / "params.json").write_text(json.dumps(defaults, indent=2) + "\n")
    (wrapper_dir / "nextflow.config").write_text(config)
    return ("review" if flags else "ok"), "; ".join(flags)



def split_top_level(s):
    toks, depth, cur = [], 0, ""
    for ch in s:
        if ch in "([":
            depth += 1
        if ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            toks.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        toks.append(cur.strip())
    return toks


def struct_to_elems(struct_str, take_name):
    """Parse a channel structure comment into elements, or None if complex."""
    s = struct_str.strip()
    s = re.sub(r"\[\s*([A-Za-z_][\w.]*)\s*\]", r"path(\1)", s)
    if s.startswith("[") and s.endswith("]"):
        s = s[1:-1]
    if re.search(r"\[[^\]]*,", s):
        return None  # multi-file nested list
    elems = []
    for t in split_top_level(s):
        m = re.match(r"(val|path)\(\s*([\w.]+)\s*\)", t)
        if m:
            elems.append({"type": m.group(1), "name": m.group(2)})
        elif re.match(r"^meta\w*$", t):
            elems.append({"type": "val", "name": t})
        elif re.match(r"^[A-Za-z_]\w*$", t):
            elems.append({"type": "path", "name": t})
        elif re.match(r"^[A-Za-z_][\w ]*$", t):
            elems.append({"type": "val", "name": snake(t), "vtype": "string"})
        else:
            return None
    return elems or None


def take_structures(text, meta, sw_dir):
    """Return [(take_name, elems|None)] using comments, meta.yml, then module-call inference."""
    m = re.search(r"take:\s*\n(.*?)(?=\n\s*main:)", text, re.S)
    if not m:
        return []
    takes = []
    for line in m.group(1).splitlines():
        line = line.strip()
        if not line or line.startswith("//"):
            continue
        name, _, comment = line.partition("//")
        takes.append((name.strip(), comment))
    mstruct = {}
    for item in (meta.get("input", []) or []):
        if isinstance(item, dict):
            for k, v in item.items():
                if isinstance(v, dict):
                    sm = re.search(r"Structure:\s*(\[.*\])", v.get("description") or "", re.S)
                    if sm:
                        mstruct[k] = sm.group(1)
    # module-call inference map: take name -> module input channel signature
    includes = dict(re.findall(r"include\s*\{\s*(\w+)(?:\s+as\s+\w+)?\s*\}\s*from\s*'([^']+)'", text))
    main_m = re.search(r"\n\s*main:\s*\n(.*?)(?=\n\s*emit:|\Z)", text, re.S)
    inferred = {}
    if main_m:
        for call in re.finditer(r"(\w+)\(([^)]*)\)", main_m.group(1)):
            alias, args = call.group(1), call.group(2)
            if alias not in includes:
                continue
            mod_path = (sw_dir / (includes[alias] + ".nf")).resolve()
            if not mod_path.exists():
                continue
            try:
                mod_channels = [parse_channel(l) for l in parse_input_lines(mod_path.read_text())]
            except ReviewNeeded:
                continue
            for pos, arg in enumerate(split_top_level(args)):
                arg = arg.strip()
                if pos < len(mod_channels) and re.match(r"^[A-Za-z_]\w*$", arg):
                    inferred.setdefault(arg, mod_channels[pos][1])
    result = []
    for name, comment in takes:
        elems = None
        cm = re.search(r"channel:\s*(\[.*\]|(?:path|val)\([^)]*\))", comment) or re.search(r"^\s*(\[.*\])\s*$", comment)
        if cm:
            elems = struct_to_elems(cm.group(1), name)
        vm = re.search(r"val:\s*(\w+)", comment)
        if elems is None and vm:
            elems = [{"type": "val", "name": name, "vtype": vm.group(1)}]
        if elems is None and name in mstruct:
            elems = struct_to_elems(mstruct[name], name)
        if elems is None and name in inferred:
            elems = inferred[name]
        result.append((name, elems))
    return result


def generate_subworkflow(src_dir, dest_dir, rel, dry_run=False):
    main_nf = src_dir / "main.nf"
    meta_yml = src_dir / "meta.yml"
    text = main_nf.read_text()
    meta = yaml.safe_load(meta_yml.read_text()) if meta_yml.exists() else {}
    wm = re.search(r"workflow\s+(\w+)\s*\{", text)
    if not wm:
        return "review", "no workflow block found"
    workflow = wm.group(1)
    proc_dir = workflow.lower()
    takes = take_structures(text, meta or {}, src_dir)
    if not takes:
        return "review", "no take channels found"
    flags = ["subworkflow wiring is best-effort; verify companion/reference semantics and internal module prerequisites"]
    params, defaults, wiring, call_args = {}, {}, [], []

    def claim(base):
        name, n = base, 2
        while name in params:
            name = f"{base}_{n}"
            n += 1
        return name

    for take_name, elems in takes:
        if not elems:
            return "review", f"take '{take_name}' structure unknown; needs manual wiring"
        has_meta = elems[0]["type"] == "val" and elems[0]["name"].startswith("meta")
        paths = [e for e in elems if e["type"] == "path"]
        vals = [e for e in elems[1:] if e["type"] == "val"] if has_meta else [e for e in elems if e["type"] == "val"]
        if not paths:
            if not vals:
                return "review", f"take '{take_name}' has no usable elements"
            vnames = []
            for ve in vals:
                vname = claim(snake(ve["name"]))
                add_val_param(vname, ve, {}, params, defaults, flags)
                vnames.append(vname)
            if len(vnames) == 1:
                call_args.append(f"params.{vnames[0]}")
            else:
                call_args.append("Channel.value([" + ", ".join(f"params.{v}" for v in vnames) + "])")
            continue
        primary = paths[0]
        pname = claim(snake(primary["name"]))
        params[pname] = {"kind": "input", "type": "path", "required": True,
                         "description": f"Input {primary['name']} for {take_name}."}
        defaults[pname] = None
        ch = f"{pname}_ch"
        wiring.append(
            f"    {ch} = Channel\n"
            f"        .fromPath(params.{pname}, checkIfExists: true)\n"
            f"        .map {{ f -> [[id: f.simpleName], f] }}"
        )
        expr = ch
        for ce in paths[1:]:
            cpname = claim(snake(ce["name"]))
            params[cpname] = {"kind": "input", "type": "path", "required": False,
                              "description": f"Optional companion {ce['name']} for {take_name}."}
            defaults[cpname] = None
            wiring.append(
                f"    {cpname}_ch = params.{cpname} ? Channel.fromPath(params.{cpname}, checkIfExists: true).map {{ f -> [f] }} : Channel.value([[]])"
            )
            expr = f"{expr}.combine({cpname}_ch)"
        for ve in vals:
            vname = claim(snake(ve["name"]))
            add_val_param(vname, ve, {}, params, defaults, flags)
            expr = f"{expr}.map {{ t -> t + [params.{vname}] }}"
        if has_meta and elems[0]["name"] != "meta":
            expr = f"{expr}.first()"
        call_args.append(expr)

    params["outdir"] = {"kind": "output", "type": "path", "required": True,
                        "description": "Output directory for the subworkflow's results."}
    defaults["outdir"] = "results"
    comp_id = f"nf-core/subworkflows/{rel.replace('_', '-')}"
    summary = first_sentence((meta or {}).get("description", "")) or f"Run the {rel} subworkflow."
    wrapper_yaml = {
        "id": comp_id,
        "name": rel.replace("_", " ")[:80],
        "summary": summary[:300],
        "params": params,
        "outputs": {"results": {"type": "directory", "path": "${outdir}/" + proc_dir, "primary": True}},
    }
    header = (
        "#!/usr/bin/env nextflow\n"
        "// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.\n"
        "// See docs/contracts/wrapper.md. Generated by scripts/generate-wrappers.py;\n"
        "// review flagged semantics before publishing.\n"
        "nextflow.enable.dsl = 2\n\n"
        f"include {{ {workflow} }} from '../main.nf'\n\n"
    )
    param_decls = "\n".join(f"params.{k} = {json.dumps(v)}" if v is not None else f"params.{k} = null" for k, v in defaults.items())
    body = "workflow {\n" + "\n".join(wiring) + f"\n\n    {workflow}({', '.join(call_args)})\n}}\n"
    config = f"""// Launcher config for the {rel} subworkflow wrapper adapter. Not part of the
// fixed wrapper/{{main.nf,params.json,wrapper.yaml}} contract — Nextflow picks
// this up automatically because it sits beside wrapper/main.nf.

process {{
    cpus   = 2
    memory = '4.GB'
    time   = '2.h'

    publishDir = [
        path: {{ "${{params.outdir}}/{proc_dir}" }},
        mode: 'copy'
    ]
}}

profiles {{
    docker {{
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
        conda.enabled = false
        singularity.enabled = false
    }}
    singularity {{
        singularity.enabled = true
        singularity.autoMounts = true
        docker.enabled = false
        conda.enabled = false
    }}
    conda {{
        conda.enabled = true
        docker.enabled = false
        singularity.enabled = false
    }}
}}
"""
    if dry_run:
        return "review", "; ".join(flags)
    wrapper_dir = dest_dir / "wrapper"
    wrapper_dir.mkdir(parents=True, exist_ok=True)
    comment = "# Agent-facing wrapper adapter for the vendored subworkflow at ../main.nf.\n# Follows docs/contracts/wrapper.md. Generated; review flagged semantics.\n# REVIEW: " + "; ".join(flags) + "\n"
    (wrapper_dir / "wrapper.yaml").write_text(comment + yaml_dump(wrapper_yaml))
    (wrapper_dir / "main.nf").write_text(header + param_decls + "\n\n" + body)
    (wrapper_dir / "params.json").write_text(json.dumps(defaults, indent=2) + "\n")
    (wrapper_dir / "nextflow.config").write_text(config)
    return "review", "; ".join(flags)


def add_val_param(vname, elem, meta_in, params, defaults, flags):
    info = meta_in.get(elem["name"], {}) or {}
    mtype = info.get("type", "string")
    desc = first_sentence(info.get("description", "")) or f"Option {elem['name']}."
    if mtype == "boolean" or BOOL_NAME_RE.search(vname):
        params[vname] = {"kind": "option", "type": "boolean", "required": False, "description": desc}
        defaults[vname] = False
    else:
        jtype = {"integer": "integer", "number": "number", "float": "number", "double": "number"}.get(mtype, "string")
        params[vname] = {"kind": "option", "type": jtype, "required": True, "description": desc}
        defaults[vname] = None
        flags.append(f"val input '{vname}' exposed as required option; confirm whether a safe default exists")


def vendor_and_generate(source, dest, only=None, dry_run=False, vendor=False, kind="modules"):
    modules_src = source / kind / "nf-core"
    modules_dest = dest / kind / "nf-core"
    report = {"ok": [], "review": [], "skipped-deprecated": [], "existing": [], "error": []}
    candidates = sorted(p for p in modules_src.rglob("main.nf"))
    for mn in candidates:
        rel = str(mn.parent.relative_to(modules_src))
        if only and rel != only and not rel.startswith(only + "/"):
            continue
        dest_mod = modules_dest / rel
        if (dest_mod / "wrapper").exists():
            report["existing"].append(rel)
            continue
        if vendor and not dry_run and not dest_mod.exists():
            shutil.copytree(mn.parent, dest_mod)
        if not dest_mod.exists() and not dry_run:
            report["error"].append(f"{rel}: not vendored at destination (use --vendor)")
            continue
        try:
            gen = generate_subworkflow if kind == "subworkflows" else generate_module
            status, detail = gen(mn.parent, dest_mod, rel, dry_run=dry_run)
            report[status].append(f"{rel}: {detail}" if detail else rel)
        except Exception as e:  # noqa: BLE001 - report and continue batch
            report["error"].append(f"{rel}: {type(e).__name__}: {e}")
    return report


def validate_tree(dest):
    problems = []
    warnings = []
    ids = {}
    for wy in sorted(dest.rglob("wrapper/wrapper.yaml")):
        rel = str(wy.parent.parent.relative_to(dest))
        try:
            data = yaml.safe_load(wy.read_text())
        except Exception as e:  # noqa: BLE001
            problems.append(f"{rel}: YAML parse error: {e}")
            continue
        for key in ("id", "name", "summary", "params", "outputs"):
            if key not in data:
                problems.append(f"{rel}: missing required key '{key}'")
        if not isinstance(data, dict):
            continue
        allowed = {"id", "name", "summary", "params", "outputs"}
        for key in data:
            if key not in allowed:
                problems.append(f"{rel}: unknown key '{key}'")
        cid = data.get("id", "")
        if cid in ids:
            problems.append(f"{rel}: duplicate id '{cid}' (also {ids[cid]})")
        ids[cid] = rel
        if not re.match(r"^(nf-core|local)/(modules|subworkflows|workflows)/[a-z0-9-]+$", cid):
            problems.append(f"{rel}: id '{cid}' does not match <provider>/<kind>/<name>")
        if not (1 <= len(data.get("name", "")) <= 80):
            problems.append(f"{rel}: name length out of range")
        if not (1 <= len(data.get("summary", "")) <= 300):
            problems.append(f"{rel}: summary length out of range")
        params = data.get("params", {}) or {}
        for pname, p in params.items():
            if not PARAM_NAME_RE.match(pname):
                problems.append(f"{rel}: bad param name '{pname}'")
            if p.get("kind") not in ("input", "output", "option"):
                problems.append(f"{rel}: param '{pname}' bad kind")
            if not p.get("type"):
                problems.append(f"{rel}: param '{pname}' missing type")
            if "default" in p:
                problems.append(f"{rel}: param '{pname}' carries a default (defaults belong in params.json)")
        outputs = data.get("outputs", {}) or {}
        if not any(o.get("primary") for o in outputs.values()):
            problems.append(f"{rel}: no primary output")
        for oname, o in outputs.items():
            for ref in re.findall(r"\$\{(\w+)\}", o.get("path", "")):
                if ref not in params:
                    problems.append(f"{rel}: output '{oname}' references undeclared param '{ref}'")
        pj = wy.parent / "params.json"
        if pj.exists():
            try:
                defaults = json.loads(pj.read_text())
                for k in defaults:
                    if k not in params and PARAM_NAME_RE.match(k):
                        problems.append(f"{rel}: params.json key '{k}' not declared in wrapper.yaml")
                    elif k not in params:
                        warnings.append(f"{rel}: params.json passthrough default '{k}' uses a pipeline-native name (not agent-facing)")
            except Exception as e:  # noqa: BLE001
                problems.append(f"{rel}: params.json parse error: {e}")
        else:
            problems.append(f"{rel}: missing params.json")
        mnf = wy.parent / "main.nf"
        if not mnf.exists():
            problems.append(f"{rel}: missing wrapper/main.nf")
        else:
            for inc in re.findall(r"include\s*\{[^}]*\}\s*from\s*'([^']+)'", mnf.read_text()):
                if inc.startswith("plugin/"):
                    continue
                target = (wy.parent / inc).resolve()
                if not target.exists() and not target.with_name(target.name + ".nf").exists():
                    problems.append(f"{rel}: include target missing: {inc}")
        comp_mnf = wy.parent.parent / "main.nf"
        if comp_mnf.exists():
            for inc in re.findall(r"include\s*\{[^}]*\}\s*from\s*['\"]([^'\"]+)['\"]", comp_mnf.read_text()):
                if inc.startswith("plugin/"):
                    continue
                target = (comp_mnf.parent / inc).resolve()
                if not target.exists() and not target.with_name(target.name + ".nf").exists():
                    problems.append(f"{rel}: component include target missing: {inc}")
    for w in warnings:
        print(f"[warn] {w}")
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, help="nf-core modules repo root (contains modules/nf-core)")
    ap.add_argument("--dest", type=Path, required=True, help="phi wrapper tree (resources/wrappers)")
    ap.add_argument("--vendor", action="store_true", help="copy missing module sources into the tree first")
    ap.add_argument("--only", help="restrict to one module path, e.g. macs2/callpeak")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--validate", action="store_true", help="validate the destination tree instead of generating")
    ap.add_argument("--report", type=Path, help="write JSON report here")
    ap.add_argument("--kind", choices=["modules", "subworkflows"], default="modules")
    args = ap.parse_args()

    if args.validate:
        problems = validate_tree(args.dest)
        for p in problems:
            print(p)
        print(f"\n{len(problems)} problem(s)")
        sys.exit(1 if problems else 0)

    if not args.source:
        ap.error("--source is required unless --validate")
    report = vendor_and_generate(args.source, args.dest, only=args.only, dry_run=args.dry_run, vendor=args.vendor, kind=args.kind)
    for status in ("ok", "review", "skipped-deprecated", "existing", "error"):
        print(f"{status}: {len(report[status])}")
    for status in ("review", "error"):
        for item in report[status]:
            print(f"  [{status}] {item}")
    if args.report:
        args.report.write_text(json.dumps(report, indent=2))
    sys.exit(1 if report["error"] else 0)


if __name__ == "__main__":
    main()