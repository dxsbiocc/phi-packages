#!/usr/bin/env python3
"""Omics visualization tools: catalog, profile, prepare, and render.

Success writes one JSON object to stdout and exits 0. Failure writes
{"error": "..."} to stdout and exits 1. Diagnostics go to stderr.

The tools look things up, count, copy, check, and run. Choosing a template
for a purpose is the agent's job; nothing here ranks templates.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any, NoReturn

# The skill directory is read-only when bundled and must stay free of run leftovers:
# never write __pycache__ for the helper modules imported below.
sys.dont_write_bytecode = True
_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)

from viz_catalog import (
    entry,
    load_catalog,
    load_contracts,
    load_families,
    lookup_alias,
    summary,
)
from viz_common import (
    COMMON_R,
    QA_SCRIPT,
    clamp,
    fail,
    node_resolve,
    note,
    project_relative,
    r_environment,
    r_failure,
    resolve_inside_project,
    run_process,
    succeed,
    summarize_qa,
)
from viz_profile import check_columns, profile
from viz_templates import (
    closest_names,
    copy_assets,
    line_section,
    list_template_sources,
    parse_template_source,
    patch_bootstrap,
    refresh_project_common_r,
    script_input_names,
    write_descriptor,
)

MAX_SUGGESTIONS = 4
QA_TIMEOUT_SECONDS = 60
R_PACKAGE_TIMEOUT_SECONDS = 60
DEFAULT_RENDER_TIMEOUT_SECONDS = 180
MAX_RENDER_TIMEOUT_SECONDS = 900
MESSAGE_CHARS = 400
OUTPUT_EXTENSIONS = {".png", ".pdf", ".svg"}
R_PACKAGE_CHECK = (
    "pk <- commandArgs(TRUE); "
    "cat(pk[!nzchar(vapply(pk, function(p) system.file(package = p), character(1)))], sep = '\\n')"
)


class JsonArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("allow_abbrev", False)
        kwargs.setdefault("add_help", False)
        super().__init__(*args, **kwargs)

    def error(self, message: str) -> NoReturn:
        fail(message)


def unknown_ids(ids: list[str], known: list[str]) -> NoReturn:
    parts = []
    for template_id in ids:
        close = closest_names(template_id, known, MAX_SUGGESTIONS)
        hint = f" (closest: {', '.join(close)})" if close else ""
        parts.append(f'"{template_id}"{hint}')
    fail(f"There is no template {'; '.join(parts)}. Call viz_catalog without arguments to list every id.")


def command_catalog(family: str | None, ids: list[str], alias: str | None) -> None:
    family = (family or "").strip()
    ids = [item.strip() for item in ids if item.strip()]
    alias = (alias or "").strip()
    if sum(1 for given in (family, ids, alias) if given) > 1:
        fail("Give at most one of family, ids, or alias.")
    templates = load_catalog()
    if family:
        families = [item["family"] for item in load_families()]
        if family not in families:
            fail(f'There is no family "{family}". Families: {", ".join(families)}.')
        chosen = [summary(template) for template in templates if template["family"] == family]
        succeed({"mode": "family", "family": family, "templates": chosen})
        return
    if ids:
        by_id = {template["template_id"]: template for template in templates}
        missing = [template_id for template_id in ids if template_id not in by_id]
        if missing:
            unknown_ids(missing, list(by_id))
        succeed({"mode": "ids", "templates": [entry(by_id[template_id]) for template_id in dict.fromkeys(ids)]})
        return
    if alias:
        succeed({"mode": "alias", "alias": alias, "templates": [entry(template) for template in lookup_alias(alias)]})
        return
    index = []
    for item in load_families():
        index.append(
            {
                "family": item["family"],
                "title": item["title"],
                "use_when": item["use_when"],
                "templates": [
                    {"template_id": template["template_id"], "title": template["title"]}
                    for template in templates
                    if template["family"] == item["family"]
                ],
            }
        )
    succeed({"mode": "index", "template_count": len(templates), "families": index})


def command_profile(data_path: str, tables: list[str], joins: list[str]) -> None:
    data_path = data_path.strip()
    if not data_path:
        fail("data_path is required.")
    paths = [data_path, *[item.strip() for item in tables if item.strip()]]
    paths = list(dict.fromkeys(paths))
    missing = [path for path in paths if not os.path.isfile(path)]
    if missing:
        fail(f"File(s) not found: {', '.join(missing)}")
    try:
        succeed(profile(paths, [item.strip() for item in joins if item.strip()]))
    except ValueError as exc:
        fail(str(exc))


def r_packages(dependencies: list[str]) -> dict[str, Any]:
    if not dependencies:
        return {"checked": True, "missing": []}
    ran = run_process(
        ["Rscript", "-e", R_PACKAGE_CHECK, *dependencies],
        env=r_environment(),
        timeout=R_PACKAGE_TIMEOUT_SECONDS,
    )
    if ran.spawn_error or ran.timed_out or ran.code != 0:
        note(ran.stderr)
        return {"checked": False, "missing": []}
    return {"checked": True, "missing": [line.strip() for line in ran.stdout.splitlines() if line.strip()]}


def command_prepare(
    template_id: str, workdir: str, reset: bool, inputs: list[str], columns: list[str]
) -> None:
    template_id = template_id.strip()
    workdir = workdir.strip()
    inputs = [item.strip() for item in inputs if item.strip()]
    columns = [item.strip() for item in columns if item.strip()]
    if not template_id or not workdir:
        fail("template_id and workdir are required.")
    if bool(inputs) != bool(columns):
        fail("Give inputs and columns together: the column check needs both.")
    templates = list_template_sources()
    template = next((candidate for candidate in templates if candidate["id"] == template_id), None)
    if template is None:
        unknown_ids([template_id], [candidate["id"] for candidate in templates])
    missing_inputs = [path for path in inputs if not os.path.isfile(path)]
    if missing_inputs:
        fail(f"Input file(s) not found: {', '.join(missing_inputs)}")
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
    result["r_packages"] = r_packages(template["dependencies"])
    if columns:
        catalog = next((item for item in load_catalog() if item["template_id"] == template_id), None)
        roles = catalog["roles"] if catalog else {"required": [], "optional": []}
        role_types = load_contracts().get("role_types", {})
        try:
            result["column_check"] = check_columns(roles, role_types, inputs, columns)
        except ValueError as exc:
            fail(str(exc))
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


KNOWN_FLAGS = {
    "catalog": {"--family", "--ids", "--alias"},
    "profile": {"--data_path", "--tables", "--joins"},
    "prepare": {"--template_id", "--workdir", "--reset", "--inputs", "--columns"},
    "render": {"--script", "--inputs", "--output", "--timeout_seconds"},
}


def build_parser() -> JsonArgumentParser:
    parser = JsonArgumentParser(prog="viz.py")
    commands = parser.add_subparsers(dest="command", required=True, parser_class=JsonArgumentParser)
    catalog = commands.add_parser("catalog")
    catalog.add_argument("--family")
    catalog.add_argument("--ids", action="append", default=[])
    catalog.add_argument("--alias")
    profile_parser = commands.add_parser("profile")
    profile_parser.add_argument("--data_path", required=True)
    profile_parser.add_argument("--tables", action="append", default=[])
    profile_parser.add_argument("--joins", action="append", default=[])
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--template_id", required=True)
    prepare.add_argument("--workdir", required=True)
    prepare.add_argument("--reset", action="store_true")
    prepare.add_argument("--inputs", action="append", default=[])
    prepare.add_argument("--columns", action="append", default=[])
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
    known = KNOWN_FLAGS.get(argv[0])
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
    if args.command == "catalog":
        command_catalog(args.family, args.ids, args.alias)
        return
    if args.command == "profile":
        command_profile(args.data_path, args.tables, args.joins)
        return
    if args.command == "prepare":
        command_prepare(args.template_id, args.workdir, args.reset, args.inputs, args.columns)
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
