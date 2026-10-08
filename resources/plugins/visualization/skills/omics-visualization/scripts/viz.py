#!/usr/bin/env python3
"""Omics visualization tools: examples, route, prepare, and render.

Success writes one JSON object to stdout and exits 0. Failure writes
{"error": "..."} to stdout and exits 1. Diagnostics go to stderr.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any, NoReturn

# The skill directory is read-only when bundled and must stay free of run leftovers:
# never write __pycache__ for the helper modules imported below.
sys.dont_write_bytecode = True
_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)

from viz_common import (
    COMMON_R,
    QA_SCRIPT,
    ROUTER,
    SKILL_ROOT,
    as_string,
    clamp,
    fail,
    json_number,
    note,
    r_environment,
    r_failure,
    node_resolve,
    project_relative,
    resolve_inside_project,
    run_process,
    string_map,
    strings,
    succeed,
    summarize_qa,
)
from viz_templates import (
    MAX_EXAMPLES,
    MAX_WHY,
    closest_names,
    copy_assets,
    installed_preview,
    is_installed_preview,
    line_section,
    list_template_sources,
    parse_template_source,
    patch_bootstrap,
    refresh_project_common_r,
    score_template,
    script_input_names,
    terms_for,
    write_descriptor,
)

DEFAULT_ROUTE_TOP = 4
MAX_ROUTE_TOP = 6
MAX_SUGGESTIONS = 4
ROUTER_TIMEOUT_SECONDS = 60
QA_TIMEOUT_SECONDS = 60
DEFAULT_RENDER_TIMEOUT_SECONDS = 180
MAX_RENDER_TIMEOUT_SECONDS = 900
ROUTER_ERROR_CHARS = 800
MESSAGE_CHARS = 400
OUTPUT_EXTENSIONS = {".png", ".pdf", ".svg"}

def compact_router_output(raw: object) -> dict[str, Any]:
    root = raw if isinstance(raw, dict) else {}
    profile = root.get("input_profile") if isinstance(root.get("input_profile"), dict) else {}
    numeric = set(strings(profile.get("numeric_columns")))
    sidecars_value = profile.get("sidecars")
    sidecars = list(sidecars_value.keys()) if isinstance(sidecars_value, dict) else []
    alignment_value = profile.get("sidecar_alignment")
    alignment = alignment_value.get("status") if isinstance(alignment_value, dict) else None
    recommendations = root.get("recommendations")
    if not isinstance(recommendations, list):
        recommendations = []
    candidates: list[dict[str, Any]] = []
    for rec in recommendations:
        if not isinstance(rec, dict):
            continue
        template_id = as_string(rec.get("id"))
        preview_value = rec.get("preview")
        preview = os.path.join(SKILL_ROOT, preview_value) if isinstance(preview_value, str) else ""
        candidates.append(
            {
                "template_id": template_id,
                "confidence": as_string(rec.get("confidence"), "unknown"),
                "score": json_number(rec.get("score")),
                "title": as_string(rec.get("title")),
                "why": strings(rec.get("rationale"))[:MAX_WHY],
                "risks": strings(rec.get("risks")),
                "use_when": as_string(rec.get("use_when")),
                "avoid_when": as_string(rec.get("avoid_when")),
                "required_roles": strings(rec.get("required_roles")),
                "role_mapping": string_map(rec.get("role_mapping")),
                "preview": preview,
                "preview_markdown": f"![{template_id}]({preview})" if preview else "",
            }
        )
    payload: dict[str, Any] = {
        "rows": json_number(profile.get("row_count")),
        "columns": [
            {"name": name, "type": "number" if name in numeric else "text"}
            for name in strings(profile.get("columns"))
        ],
        "roles": string_map(profile.get("role_mapping")),
    }
    if sidecars:
        payload["sidecars"] = sidecars
    if isinstance(alignment, str) and alignment != "not_checked":
        payload["sidecar_alignment"] = alignment
    return {"input": payload, "candidates": candidates}

class JsonArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("allow_abbrev", False)
        kwargs.setdefault("add_help", False)
        super().__init__(*args, **kwargs)

    def error(self, message: str) -> NoReturn:
        fail(message)

def command_examples(purpose: str, top: int | None) -> None:
    purpose = purpose.strip()
    if not purpose:
        fail("purpose is required.")
    catalog_path = os.path.join(SKILL_ROOT, "references", "template_contracts.json")
    with open(catalog_path, encoding="utf-8") as handle:
        catalog = json.load(handle)
    templates = catalog.get("templates") if isinstance(catalog, dict) else None
    if not isinstance(templates, list):
        templates = []
    terms = terms_for(purpose)
    matches: list[tuple[int, int, dict[str, Any]]] = []
    for order, template in enumerate(templates):
        if not isinstance(template, dict):
            continue
        score = score_template(template, terms)
        if score > 0:
            matches.append((score, order, template))
    matches.sort(key=lambda item: (-item[0], item[1]))
    limit = clamp(top if top is not None else MAX_EXAMPLES, 1, MAX_EXAMPLES)
    candidates: list[dict[str, Any]] = []
    for _score, _order, template in matches:
        preview = installed_preview(template)
        if preview is None:
            continue
        template_id = str(template.get("id", ""))
        roles = template.get("roles") if isinstance(template.get("roles"), dict) else {}
        required = roles.get("required") if isinstance(roles, dict) else None
        if isinstance(required, list):
            required_roles = [item for item in required if isinstance(item, str)]
        else:
            required_roles = []
        candidates.append(
            {
                "template_id": template_id,
                "title": str(template.get("title", "")),
                "family": str(template.get("family", "")),
                "preview": preview,
                "preview_markdown": f"![{template_id}]({preview})",
                "use_when": str(template.get("use_when") or ""),
                "avoid_when": str(template.get("avoid_when") or ""),
                "required_roles": required_roles,
            }
        )
        if len(candidates) >= limit:
            break
    succeed(
        {
            "purpose": purpose,
            "source": "bundled-template-preview",
            "dataFitted": False,
            "totalMatches": len(matches),
            "candidates": candidates,
        }
    )

def command_route(
    data_path: str, purpose: str, mode: str, top: int | None, sidecar_dir: str | None
) -> None:
    data_path = data_path.strip()
    purpose = purpose.strip()
    if not data_path or not purpose:
        fail("data_path and purpose are required.")
    if not os.path.isfile(data_path):
        fail(f"The data table was not found: {data_path}")
    limit = clamp(top if top is not None else DEFAULT_ROUTE_TOP, 1, MAX_ROUTE_TOP)
    argv = [
        sys.executable,
        ROUTER,
        "--input",
        data_path,
        "--query",
        purpose,
        "--mode",
        mode,
        "--top",
        str(limit),
    ]
    if sidecar_dir and sidecar_dir.strip():
        argv.extend(["--sidecar-dir", sidecar_dir.strip()])
    argv.append("--json")
    outcome = run_process(argv, timeout=ROUTER_TIMEOUT_SECONDS)
    if outcome.spawn_error:
        fail(
            f"python3 was not found ({outcome.spawn_error}). "
            "The template router needs Python 3 on PATH; report this as a missing dependency."
        )
    if outcome.timed_out or outcome.code != 0:
        note(outcome.stderr)
        detail = "timed out" if outcome.timed_out else outcome.stderr.strip()[-ROUTER_ERROR_CHARS:]
        fail(f"The template router failed: {detail}")
    try:
        parsed = json.loads(outcome.stdout)
    except json.JSONDecodeError:
        note(outcome.stderr)
        fail("The template router returned output that is not JSON.")
    if outcome.stderr:
        note(outcome.stderr)
    succeed(compact_router_output(parsed))

def command_prepare(template_id: str, workdir: str, reset: bool) -> None:
    template_id = template_id.strip()
    workdir = workdir.strip()
    if not template_id or not workdir:
        fail("template_id and workdir are required.")
    templates = list_template_sources()
    template = next((candidate for candidate in templates if candidate["id"] == template_id), None)
    if template is None:
        close = closest_names(template_id, [candidate["id"] for candidate in templates], MAX_SUGGESTIONS)
        suggestion = f" Closest: {', '.join(close)}." if close else ""
        fail(
            f'There is no template "{template_id}".{suggestion}'
            " Use viz_route to find one for the data; do not guess template ids."
        )
    directory = resolve_inside_project(os.getcwd(), workdir, "The working directory")
    script = os.path.join(directory, "plot.R")
    resolve_inside_project(os.getcwd(), os.path.realpath(script), "The script")
    existing = os.path.exists(script) and not reset
    assets: list[str] = []
    if not existing:
        os.makedirs(directory, exist_ok=True)
        original = open(template["path"], encoding="utf-8").read()
        try:
            patched = patch_bootstrap(original, COMMON_R)
        except ValueError as exc:
            fail(str(exc))
        with open(script, "w", encoding="utf-8") as handle:
            handle.write(patched)
        assets = copy_assets(template["path"], directory)
    current = parse_template_source(refresh_project_common_r(script, COMMON_R))
    names = template["input_names"]
    result: dict[str, Any] = {"template_id": template["id"], "script": script, "common_r": COMMON_R}
    if existing:
        result["existing"] = True
    if assets:
        result["assets"] = assets
    result["purpose"] = template["purpose"]
    result["inputs"] = names
    joined = " ".join(f"<{name}>" for name in names)
    result["run"] = f"Rscript {script} {joined} <output>"
    result["dependencies"] = template["dependencies"]
    if template.get("adaptation"):
        result["adaptation"] = template["adaptation"]
    if template.get("assumptions"):
        result["assumptions"] = template["assumptions"]
    config = line_section(current.get("config"))
    data_preparation = line_section(current.get("data_preparation"))
    if config:
        result["config"] = config
    if data_preparation:
        result["data_preparation"] = data_preparation
    plot = current.get("plot")
    if isinstance(plot, dict):
        result["plot"] = {"start_line": plot["start_line"], "end_line": plot["end_line"]}
    succeed(result)

def command_render(script: str, inputs: list[str], output: str, timeout_seconds: int | None) -> None:
    script = script.strip()
    output = output.strip()
    if not script or not output:
        fail("script and output are required.")
    cleaned = [item.strip() for item in inputs if item.strip()]
    if not cleaned:
        fail("At least one input table is required.")
    cwd = os.getcwd()
    script_path = resolve_inside_project(cwd, script, "The script")
    output_path = resolve_inside_project(cwd, output, "The output file")
    extension = os.path.splitext(output_path)[1].lower()
    if extension not in OUTPUT_EXTENSIONS:
        fail("The output file must be a .png, .pdf or .svg file.")
    if not os.path.exists(script_path):
        fail(f"The script was not found: {script_path}. Create it with viz_prepare.")
    resolved_inputs = [node_resolve(cwd, item) for item in cleaned]
    missing = [item for item in resolved_inputs if not os.path.exists(item)]
    if missing:
        fail(f"Input file(s) not found: {', '.join(missing)}")
    expected = script_input_names(refresh_project_common_r(script_path, COMMON_R))
    if expected is not None and len(expected) != len(resolved_inputs):
        fail(
            f"This template expects {len(expected)} input(s): {', '.join(expected)}; "
            f"got {len(resolved_inputs)}."
        )
    parent = os.path.dirname(output_path) or cwd
    os.makedirs(parent, exist_ok=True)
    seconds = clamp(
        timeout_seconds if timeout_seconds is not None else DEFAULT_RENDER_TIMEOUT_SECONDS,
        1,
        MAX_RENDER_TIMEOUT_SECONDS,
    )
    ran = run_process(
        ["Rscript", script_path, *resolved_inputs, output_path],
        cwd=os.path.dirname(script_path),
        env=r_environment(),
        timeout=seconds,
    )
    if ran.spawn_error:
        fail(
            f"Rscript was not found ({ran.spawn_error}). "
            "Rendering needs R on PATH; report this as a missing dependency and do not install it."
        )
    if ran.timed_out:
        note(ran.stderr)
        fail(f"R timed out after {seconds} s.")
    if ran.code != 0:
        note(ran.stderr)
        fail(r_failure(ran.stderr))
    if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
        note(ran.stderr)
        fail(
            f"R finished but did not write {output_path}. "
            f"Check the output path the script uses ({os.path.basename(script_path)} "
            "takes it as its last argument)."
        )
    qa = run_process(
        [sys.executable, QA_SCRIPT, output_path, "--json"],
        timeout=QA_TIMEOUT_SECONDS,
    )
    if qa.spawn_error or qa.timed_out:
        summary: dict[str, Any] = {"qa": {"ok": False, "failed": ["qa_unavailable"]}}
    else:
        summary = summarize_qa(qa.stdout)
    if ran.stderr:
        note(ran.stderr)
    if qa.stderr:
        note(qa.stderr)
    try:
        write_descriptor(
            cwd=cwd,
            script_path=script_path,
            inputs=resolved_inputs,
            output_path=output_path,
            timeout_seconds=seconds if timeout_seconds is not None else None,
        )
    except Exception as exc:
        fail(f"Could not write the artifact descriptor: {exc}")
    messages = ran.stderr.strip()[-MESSAGE_CHARS:]
    result: dict[str, Any] = {
        "ok": True,
        "output": output_path,
        "bytes": os.path.getsize(output_path),
        "qa": summary["qa"],
    }
    if "format" in summary:
        result["format"] = summary["format"]
    if "width" in summary:
        result["width"] = summary["width"]
    if "height" in summary:
        result["height"] = summary["height"]
    if messages:
        result["messages"] = messages
    result["artifacts"] = [project_relative(cwd, output_path)]
    succeed(result)

def build_parser() -> JsonArgumentParser:
    parser = JsonArgumentParser(prog="viz.py")
    commands = parser.add_subparsers(dest="command", required=True, parser_class=JsonArgumentParser)
    examples = commands.add_parser("examples")
    examples.add_argument("--purpose", required=True)
    examples.add_argument("--top", type=int)
    route = commands.add_parser("route")
    route.add_argument("--data_path", required=True)
    route.add_argument("--purpose", required=True)
    route.add_argument("--mode", choices=("preview", "publication"), default="preview")
    route.add_argument("--top", type=int)
    route.add_argument("--sidecar_dir")
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--template_id", required=True)
    prepare.add_argument("--workdir", required=True)
    prepare.add_argument("--reset", action="store_true")
    render = commands.add_parser("render")
    render.add_argument("--script", required=True)
    render.add_argument("--inputs", action="append", required=True)
    render.add_argument("--output", required=True)
    render.add_argument("--timeout_seconds", type=int)
    return parser

def reject_unknown_flags(argv: list[str]) -> None:
    """Report an unknown option before a missing required one.

    Subparsers validate required flags before leftover options, so `--help`
    would otherwise be swallowed. The message matches argparse.
    """
    if not argv or argv[0].startswith("-"):
        return
    known = {
        "examples": {"--purpose", "--top"},
        "route": {"--data_path", "--purpose", "--mode", "--top", "--sidecar_dir"},
        "prepare": {"--template_id", "--workdir", "--reset"},
        "render": {"--script", "--inputs", "--output", "--timeout_seconds"},
    }.get(argv[0])
    if known is None:
        return
    for token in argv[1:]:
        if token == "--":
            break
        if token.startswith("-") and token.split("=", 1)[0] not in known:
            fail(f"unrecognized arguments: {token}")

def main(argv: list[str] | None = None) -> None:
    args_list = list(sys.argv[1:] if argv is None else argv)
    reject_unknown_flags(args_list)
    args = build_parser().parse_args(args_list)
    if args.command == "examples":
        command_examples(args.purpose, args.top)
        return
    if args.command == "route":
        command_route(args.data_path, args.purpose, args.mode, args.top, args.sidecar_dir)
        return
    if args.command == "prepare":
        command_prepare(args.template_id, args.workdir, args.reset)
        return
    if args.command == "render":
        command_render(args.script, args.inputs if args.inputs else [], args.output, args.timeout_seconds)
        return
    fail(f"unknown command {args.command}")

if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:
        fail(str(exc))
