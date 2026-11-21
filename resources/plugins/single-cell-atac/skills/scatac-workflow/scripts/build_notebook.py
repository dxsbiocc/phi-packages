#!/usr/bin/env python3
"""Render a runnable analysis notebook into the user's working directory.

The template is a fragment library, not a stored .ipynb: analysis code in
assets/fragments/, presets and parameter defaults in assets/steps.json,
per-step explanations in references/step-notes/.

* A package may not execute anything under `assets/` (skill contract § 1), so
  fragments ship as data and are rendered into the user's project, where Phi's
  notebook tools and the managed kernel can run them.
* The renderer fills every placeholder from the user's intent (which steps)
  and the input data (paths, column names). The deliverable is one finished
  notebook, not a template directory to copy and edit.
* All tunable values live in a single config cell at the top; fragment code
  references those variables, so changing a threshold means editing one cell
  and re-running the cells below it.
* `--inject` adds an external R or Python script as a code cell; `--merge`
  splices in another notebook, whole or a cell range (`other.ipynb[0:6]`).
  A notebook has one kernel, so every source must match the target language.

The renderer is language-neutral machinery. Everything specific to a toolkit
(kernel names, footer, attribution, presets, parameter defaults) lives in
assets/steps.json and the teaching text in references/step-notes/, so the same
file ships in every workflow skill of the plugin and behaves according to the
library next to it.

Raw-read processing (Cell Ranger, alignment) is deliberately absent: that is
the wrapper tree's job. The notebook starts from a count matrix or a saved
object.

With --self-check every code cell is executed in order in the output directory
before this script returns, so a broken cell is reported here rather than
discovered by a learner.
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import uuid
from datetime import datetime, timezone

RENDERER_VERSION = 3
SKILL_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
ASSET_ROOT = os.path.join(SKILL_ROOT, "assets")
NOTES_ROOT = os.path.join(SKILL_ROOT, "references", "step-notes")

LANGUAGES = {
    "r": {"info": "R", "exts": {".r"}, "runner": ["Rscript"], "script": "_selfcheck.R"},
    "python": {"info": "python", "exts": {".py"}, "runner": [sys.executable], "script": "_selfcheck.py"},
}

HOW_TO = (
    "## How to use this notebook\n\n"
    "Run the cells in order with `Shift+Enter`. Each cell states what it does and why "
    "before it runs.\n\n"
    "**Edit the config cell.** Every threshold, path, and column name lives in the first "
    "code cell; the analysis cells reference those variables. Change a value there and "
    "re-run the cells below it — nothing else needs editing. Each value's comment says "
    "where it came from (default, config file, or an explicit flag).\n\n"
    "**Read the output, not just the absence of errors.** A clustering run that "
    "produces clusters but no biology is still a failed analysis."
)


def die(message):
    print(json.dumps({"error": message}, indent=2, ensure_ascii=False))
    sys.exit(1)


def info(message):
    print(f"[render] {message}", file=sys.stderr)


# ---------------------------------------------------------------------------
# fragment library
# ---------------------------------------------------------------------------

def load_library():
    path = os.path.join(ASSET_ROOT, "steps.json")
    if not os.path.exists(path):
        die(f"fragment library missing: {path}")
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def resolve_steps(registry, requested):
    """Close the requested steps over their dependencies, in registry order."""
    by_id = {step["id"]: step for step in registry["steps"]}
    unknown = [name for name in requested if name not in by_id]
    if unknown:
        die(f"unknown steps: {', '.join(unknown)}. Available: {', '.join(by_id)}")
    wanted = set()
    stack = list(requested)
    while stack:
        name = stack.pop()
        if name in wanted:
            continue
        wanted.add(name)
        stack.extend(by_id[name].get("requires", []))
    return [step["id"] for step in registry["steps"] if step["id"] in wanted]


def read_fragment(language, step_id):
    base = os.path.join(ASSET_ROOT, "fragments", step_id)
    code_path = base + "." + language_ext(language)
    md_path = os.path.join(NOTES_ROOT, step_id + ".md")
    if not os.path.exists(code_path):
        die(f"fragment code missing: {code_path}")
    with open(code_path, encoding="utf-8") as handle:
        code = handle.read().rstrip("\n")
    description = ""
    if os.path.exists(md_path):
        with open(md_path, encoding="utf-8") as handle:
            description = handle.read().strip()
    digest = hashlib.sha256(code.encode("utf-8")).hexdigest()[:12]
    return code, description, digest


def language_ext(language):
    return "R" if language == "r" else "py"


# ---------------------------------------------------------------------------
# parameters and the config cell
# ---------------------------------------------------------------------------

def resolve_parameters(registry, args):
    """Flag > config file > library default; each value keeps its source."""
    defaults = registry["parameters"] if registry else {}
    config = {}
    if args.config:
        with open(args.config, encoding="utf-8") as handle:
            config = json.load(handle)
        unknown = [key for key in config if key not in defaults]
        if unknown:
            die(f"unknown keys in --config: {', '.join(unknown)}. Known: {', '.join(defaults)}")
    resolved = {}
    for name, spec in defaults.items():
        flag_value = getattr(args, name, None)
        if flag_value is not None:
            if spec.get("type") == "list":
                flag_value = flag_list(flag_value)
            resolved[name] = {"value": flag_value, "source": "flag"}
        elif name in config:
            resolved[name] = {"value": config[name], "source": "config"}
        else:
            resolved[name] = {"value": spec["default"], "source": "default"}
    return resolved


def flag_bool(value):
    return str(value).lower() in ("true", "1", "yes", "on")


def flag_list(values):
    out = []
    for value in values:
        try:
            out.append(float(value))
        except ValueError:
            out.append(value)
    return out


def literal(language, value):
    if value is None:
        return "NULL" if language == "r" else "None"
    if isinstance(value, bool):
        return ("TRUE" if value else "FALSE") if language == "r" else repr(value)
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=False)
    if isinstance(value, list):
        if language == "r":
            return "c(" + ", ".join(literal(language, item) for item in value) + ")"
        return repr(value)
    return repr(value)


def config_cell(language, structural, parameters, registry):
    assign = "<-" if language == "r" else "="
    lines = [
        "# " + "=" * 62,
        "# Analysis parameters — the single place to tune this notebook.",
        "# Edit a value here, then re-run the cells below. The tag in each",
        "# comment says where the value came from at render time.",
        "# " + "=" * 62,
        "",
        "# Input and output (paths are relative to this notebook's folder)",
    ]
    width = max(len(name) for name in list(structural) + list(parameters))
    for name, value in structural.items():
        lines.append(f"{name.ljust(width)} {assign} {literal(language, value)}")
    lines.append("")
    lines.append("# QC and analysis parameters")
    comments = {name: spec.get("comment", "") for name, spec in (registry["parameters"].items() if registry else [])}
    for name, entry in parameters.items():
        comment = comments.get(name, "")
        lines.append(f"{name.ljust(width)} {assign} {literal(language, entry['value'])}  # [{entry['source']}] {comment}".rstrip())
    return "\n".join(lines)


def substitution_map(structural, parameters, footer):
    tokens = {
        "{{DATA}}": structural["data_path"],
        "{{PROJECT}}": structural["project_name"],
        "{{OUTPUT}}": structural["output_object"],
        "{{SAMPLE_COL}}": structural["sample_column"],
        "{{FOOTER}}": footer,
    }
    for name, entry in parameters.items():
        tokens["{{" + name.upper() + "}}"] = str(entry["value"])
    return tokens


def expand(text, tokens):
    for key, value in tokens.items():
        text = text.replace(key, value)
    return text


# ---------------------------------------------------------------------------
# external sources: injected scripts and merged notebooks
# ---------------------------------------------------------------------------

def detect_script_language(path):
    ext = os.path.splitext(path)[1].lower()
    for language, spec in LANGUAGES.items():
        if ext in spec["exts"]:
            return language
    die(f"cannot tell the language of '{path}': use a .R or .py file")


def parse_merge_spec(spec):
    match = re.match(r"^(?P<path>.+?)(?:\[(?P<start>\d+)(?::(?P<end>\d+))?\])?$", spec)
    if not match or not match.group("path"):
        die(f"cannot parse --merge '{spec}'; use path.ipynb or path.ipynb[start:end]")
    start = match.group("start")
    end = match.group("end")
    selection = None
    if start is not None:
        selection = (int(start), int(end) if end is not None else int(start) + 1)
    return match.group("path"), selection


def notebook_language(nb, path):
    meta = nb.get("metadata", {})
    name = meta.get("language_info", {}).get("name") or meta.get("kernelspec", {}).get("language")
    if not name:
        die(f"cannot determine the language of '{path}': no language_info or kernelspec")
    name = name.lower()
    if name in LANGUAGES:
        return name
    die(f"'{path}' declares language '{name}', which this renderer does not support")


def load_merge(spec, language):
    path, selection = parse_merge_spec(spec)
    if not os.path.exists(path):
        die(f"merge source not found: {path}")
    with open(path, encoding="utf-8") as handle:
        nb = json.load(handle)
    source_language = notebook_language(nb, path)
    if source_language != language:
        die(f"'{path}' is a {source_language} notebook; this render targets {language}. "
            "A notebook has one kernel — render a separate notebook for the other language.")
    cells = nb.get("cells", [])
    if selection:
        cells = cells[selection[0]:selection[1]]
    kept = []
    for cell in cells:
        kind = cell.get("cell_type")
        if kind not in ("code", "markdown"):
            continue
        source = "".join(cell.get("source", []))
        kept.append((kind, source, {"phi_source": f"merge:{path}"}))
    label = f"{path}[{selection[0]}:{selection[1]}]" if selection else path
    return kept, label


def load_inject(path, language):
    if not os.path.exists(path):
        die(f"inject source not found: {path}")
    script_language = detect_script_language(path)
    if script_language != language:
        die(f"'{path}' is a {script_language} script; this render targets {language}. "
            "A notebook has one kernel — inject it into a notebook of its own language.")
    with open(path, encoding="utf-8") as handle:
        code = handle.read().rstrip("\n")
    return code


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------

def build_cells(args, language, registry, library, structural, parameters):
    tokens = substitution_map(structural, parameters, library["footer"])
    cells = []
    sources = []

    analysis_title = library.get("analysis_title", "Analysis")
    if args.module is not None:
        preset = registry["presets"][str(args.module)]
        title, intro = preset["title"], preset["intro"]
        step_ids = preset["steps"]
    elif args.steps:
        step_ids = resolve_steps(registry, args.steps)
        title = f"{analysis_title} — rendered from the {library['library_title']}"
        intro = (
            f"Rendered for this project from the {library['skill']} fragment library. Steps: "
            + ", ".join(step_ids)
            + ".\n\n"
            + library["attribution"]
        )
    else:
        step_ids = []
        title = f"{analysis_title} — assembled from your sources"
        intro = (
            "Rendered for this project from the injected scripts and merged notebooks "
            "listed below. The config cell holds the shared paths and parameters."
        )
    cells.append(("markdown", f"# {title}\n\n{intro}", {}))
    cells.append(("markdown", HOW_TO, {}))
    cells.append(("markdown", "## Configuration", {}))
    cells.append(("code", config_cell(language, structural, parameters, registry or library), {"phi_source": "config"}))

    for index, step_id in enumerate(step_ids, start=1):
        step = next(s for s in registry["steps"] if s["id"] == step_id)
        code, description, digest = read_fragment(language, step_id)
        heading = f"## {index}. {step['title']}"
        if description:
            heading += "\n\n" + description
        cells.append(("markdown", expand(heading, tokens), {}))
        cells.append(("code", expand(code, tokens), {"phi_source": f"fragment:{language}/{step_id}"}))
        sources.append({"kind": "fragment", "id": step_id, "sha256": digest})

    for path in args.inject or []:
        code = load_inject(path, language)
        cells.append(("markdown", f"## Injected script: `{os.path.basename(path)}`\n\nAdded verbatim from `{path}`.", {}))
        cells.append(("code", code, {"phi_source": f"inject:{path}"}))
        sources.append({"kind": "inject", "path": path})

    for spec in args.merge or []:
        merged, label = load_merge(spec, language)
        cells.append(("markdown", f"## Merged from `{label}`", {}))
        cells.extend(merged)
        sources.append({"kind": "merge", "spec": label, "cells": len(merged)})

    cells.append(("markdown", library["footer"], {}))
    return cells, step_ids, sources


def make_notebook(cells, language, kernel, display, provenance):
    spec = LANGUAGES[language]
    return {
        "cells": [
            {
                "cell_type": kind,
                "id": uuid.uuid4().hex[:8],
                "metadata": meta,
                "source": text.splitlines(keepends=True),
                **({"outputs": [], "execution_count": None} if kind == "code" else {}),
            }
            for kind, text, meta in cells
        ],
        "metadata": {
            "kernelspec": {"display_name": display, "language": spec["info"], "name": kernel},
            "language_info": {"name": spec["info"]},
            "phi_render": provenance,
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def data_fingerprint(path):
    if os.path.isdir(path):
        entries = sorted(os.listdir(path))
        return {"kind": "directory", "entries": entries[:20], "entry_count": len(entries)}
    stat = os.stat(path)
    return {"kind": "file", "bytes": stat.st_size, "mtime": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat()}


# ---------------------------------------------------------------------------
# self-check
# ---------------------------------------------------------------------------

def run_cells(cells, workdir, language):
    spec = LANGUAGES[language]
    parts = []
    if language == "python":
        parts.append('import matplotlib; matplotlib.use("Agg"); import scanpy as _sc; _sc.settings.autoshow = False')
    for index, (kind, text, _meta) in enumerate(cells):
        if kind != "code":
            continue
        parts.append(f"cat('\\n=== cell {index} ===\\n')" if language == "r" else f"print('\\n=== cell {index} ===')")
        parts.append(text)
    path = os.path.join(workdir, spec["script"])
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n\n".join(parts))
    try:
        # Self-check runs outside Jupyter: without a headless backend macOS
        # matplotlib opens a native window per plot (and plt.show blocks).
        result = subprocess.run(
            spec["runner"] + [path],
            capture_output=True,
            text=True,
            cwd=workdir,
            timeout=900,
            env={**os.environ, "MPLBACKEND": "Agg"},
        )
    except subprocess.TimeoutExpired:
        return False, "self-check timed out after 900 seconds"
    finally:
        if os.path.exists(path):
            os.remove(path)
    output = (result.stdout or "") + (result.stderr or "")
    return result.returncode == 0, output


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------

def main():
    library = load_library()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", type=int, help="library preset (see steps.json presets)")
    parser.add_argument("--steps", action="append", help="step id to render (repeatable); declared dependencies are added automatically")
    parser.add_argument("--output", required=True, help="notebook path inside the project")
    parser.add_argument("--data", required=True, help="input matrix directory or saved object")
    parser.add_argument("--project", default=None, help="project name recorded on the object")
    parser.add_argument("--output-object", default=None, help="where the notebook saves its object")
    parser.add_argument("--sample-column", default="sample")
    parser.add_argument("--config", help="JSON file of parameter overrides (keys mirror flag names)")
    parser.add_argument("--inject", action="append", help="R or Python script to add as a code cell (repeatable)")
    parser.add_argument("--merge", action="append", help="notebook to splice in, whole or path.ipynb[start:end] (repeatable)")
    parser.add_argument("--language", choices=sorted(LANGUAGES), default=None)
    parser.add_argument("--kernel", default=None, help="kernel name; defaults to the library's kernel for the language")
    parser.add_argument("--self-check", action="store_true", help="execute every cell in the output directory now")
    # Every library parameter is also a flag, so the engine can pass any of them.
    for name, spec in library["parameters"].items():
        if any(a.dest == name for a in parser._actions):
            continue
        kind = spec.get("type", "string")
        if kind == "int":
            parser.add_argument(f"--{name}", type=int, default=None)
        elif kind == "number":
            parser.add_argument(f"--{name}", type=float, default=None)
        elif kind == "bool":
            parser.add_argument(f"--{name}", type=flag_bool, default=None)
        elif kind == "list":
            parser.add_argument(f"--{name}", action="append", default=None)
        else:
            parser.add_argument(f"--{name}", default=None)
    args = parser.parse_args()

    language = args.language or library["language"]
    if language not in library["kernels"]:
        die(f"the {library['skill']} library supports {', '.join(sorted(library['kernels']))}, not {language}; "
            "use the sibling workflow skill for the other language")
    registry = library if library["language"] == language else None

    if language == library["language"]:
        if (args.module is None) == (args.steps is None):
            die("pass exactly one of --module (a preset) or --steps (a comma-separated step list)")
    elif args.module is not None or args.steps:
        die(f"the {library['skill']} fragment library renders {library['language']} steps; "
            f"a {language} notebook is assembled from --inject/--merge sources only")
    if args.module is not None and str(args.module) not in (registry or {}).get("presets", {}):
        die(f"unknown preset {args.module}. Available: {', '.join(sorted((registry or {}).get('presets', {})))}")
    args.steps = [name.strip() for item in (args.steps or []) for name in item.split(",") if name.strip()] or None

    parameters = resolve_parameters(registry or library, args)

    outdir = os.path.dirname(os.path.abspath(args.output)) or "."
    data_abs = os.path.abspath(args.data)
    if not os.path.exists(data_abs):
        die(f"data path does not exist: {args.data}")
    output_object = args.output_object or library.get("output_object_default", "analyzed.h5ad")
    if not os.path.isabs(output_object):
        output_object = os.path.join(outdir, output_object)
    structural = {
        "data_path": os.path.relpath(data_abs, outdir),
        "project_name": args.project or os.path.basename(data_abs.rstrip("/")) or "scatac",
        "output_object": os.path.relpath(output_object, outdir),
        "sample_column": args.sample_column,
    }

    cells, step_ids, sources = build_cells(args, language, registry, library, structural, parameters)
    kernel_spec = library["kernels"][language]
    kernel = args.kernel or kernel_spec["name"]
    provenance = {
        "renderer": library.get("renderer_id", "single-cell-atac/build_notebook.py"),
        "skill": library["skill"],
        "version": RENDERER_VERSION,
        "generated": datetime.now(timezone.utc).isoformat(),
        "language": language,
        "preset": str(args.module) if args.module is not None else None,
        "steps": step_ids,
        "parameters": parameters,
        "structural": structural,
        "data": data_fingerprint(data_abs),
        "sources": sources,
    }
    notebook = make_notebook(cells, language, kernel, kernel_spec["display"], provenance)

    report = {
        "notebook": os.path.abspath(args.output),
        "language": language,
        "steps": step_ids,
        "cells": len(cells),
        "code_cells": sum(1 for kind, _, _ in cells if kind == "code"),
        "sources": sources,
    }

    os.makedirs(outdir, exist_ok=True)
    if args.self_check:
        info(f"self-check: executing {report['code_cells']} code cells in {outdir}")
        ok, output = run_cells(cells, outdir, language)
        report["self_check"] = "passed" if ok else "failed"
        if not ok:
            report["error"] = "A cell failed when executed in the output directory"
            report["output_tail"] = output[-4000:]
            print(json.dumps(report, indent=2, ensure_ascii=False))
            return 1
        report["check_workdir"] = outdir

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(notebook, handle, indent=1, ensure_ascii=False)
        handle.write("\n")
    info(f"wrote {args.output}")

    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())