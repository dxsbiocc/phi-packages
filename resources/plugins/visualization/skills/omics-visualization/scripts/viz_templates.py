"""Template parsing, listing, copying, and scoring for the omics visualization CLI."""

from __future__ import annotations

import ast
import json
import ntpath
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from typing import Any

sys.dont_write_bytecode = True

from viz_common import (
    MEDIA_TYPES,
    SKILL_ROOT,
    note,
    one_line,
    png_size,
    project_relative,
    sha256_file,
)

MAX_EXAMPLES = 4
MAX_WHY = 3
MAX_SUGGESTIONS = 4
MAX_SCAN_DEPTH = 4
TITLE_MAX = 200

GENERIC_TERMS = frozenset(
    {
        "example",
        "examples",
        "preview",
        "previews",
        "show",
        "template",
        "templates",
        "plot",
        "figure",
        "示例",
        "模板",
        "预览",
        "图片",
        "图",
    }
)
EXAMPLE_DATA_EXTENSIONS = {".tsv", ".csv"}
NOT_ASSETS = {"plot.R", "preview.png"}

HEADER_LABEL = re.compile(r"^# ([A-Z][A-Za-z-]*(?: [a-z]+)*):[ \t]*(.*)$")
SECTION_MARKER = re.compile(r"^# (CONFIG|DATA PREPARATION|PLOT|SAVE)\b")
RULER = re.compile(r"^# -{10,}")
BOOTSTRAP = re.compile(r"^local\(\{\n[\s\S]*?^\}\)\n", re.MULTILINE)
INPUT_NAMES = re.compile(r"input_names\s*=\s*c\(([^)]*)\)")
QUOTED = re.compile(r"""["']([^"']+)["']""")
PARSE_IO = re.compile(r"parse_io_args\s*\(")
DEPENDENCY_TOKEN = re.compile(r"^[A-Za-z][A-Za-z0-9.]*$")

_TEMPLATES: list[dict[str, Any]] | None = None


def header_fields(lines: list[str]) -> dict[str, str]:
    fields: dict[str, str] = {}
    label: str | None = None
    parts: list[str] = []

    def flush() -> None:
        nonlocal label
        if label is not None:
            fields[label] = "\n".join(parts).strip()

    for line in lines:
        if line.startswith("#!"):
            continue
        if not line.startswith("#"):
            if line.strip() == "" and len(fields) + len(parts) > 0:
                break
            if line.strip() != "" and not line.startswith("#!"):
                break
            continue
        heading = HEADER_LABEL.match(line)
        if heading:
            flush()
            label = heading.group(1)
            parts = [heading.group(2)] if heading.group(2) else []
        elif label is not None:
            parts.append(re.sub(r"^#\s{0,2}", "", line, count=1))
    flush()
    return fields


def input_names_of(text: str) -> list[str]:
    match = INPUT_NAMES.search(text)
    if not match:
        return ["input"]
    quoted = QUOTED.findall(match.group(1))
    return quoted if quoted else ["input"]

def script_input_names(text: str) -> list[str] | None:
    if PARSE_IO.search(text):
        return input_names_of(text)
    return None


def section_map(lines: list[str]) -> dict[str, dict[str, Any]]:
    markers: list[tuple[str, int]] = []
    for index, line in enumerate(lines):
        match = SECTION_MARKER.match(line)
        if match is None:
            continue
        following = lines[index + 1] if index + 1 < len(lines) else ""
        if not RULER.match(following):
            continue
        markers.append((match.group(1), index))
    found: dict[str, dict[str, Any]] = {}
    for position, (name, index) in enumerate(markers):
        first = index + 2
        next_index = markers[position + 1][1] if position + 1 < len(markers) else None
        last = next_index - 2 if next_index is not None else len(lines) - 1
        while last > first and lines[last].strip() == "":
            last -= 1
        if last < first:
            continue
        found[name] = {
            "start_line": first + 1,
            "end_line": last + 1,
            "text": "\n".join(lines[first : last + 1]),
        }
    return found


def dependencies_of(header: dict[str, str]) -> list[str]:
    return [
        token
        for token in re.split(r"[\s,]+", header.get("Dependencies", ""))
        if DEPENDENCY_TOKEN.fullmatch(token) and token != "and"
    ]


def parse_template_source(text: str) -> dict[str, Any]:
    lines = text.split("\n")
    header = header_fields(lines)
    found = section_map(lines)
    source: dict[str, Any] = {
        "id": header.get("Template-ID", ""),
        "purpose": one_line(header.get("Purpose", "")),
        "input_names": input_names_of(text),
        "dependencies": dependencies_of(header),
        "title": one_line(header["Title"]) if header.get("Title") else "",
    }
    adaptation = header.get("Agent adaptation")
    if adaptation:
        source["adaptation"] = one_line(adaptation)
    assumptions = header.get("Scientific assumptions")
    if assumptions:
        source["assumptions"] = one_line(assumptions)
    if "CONFIG" in found:
        source["config"] = found["CONFIG"]
    if "DATA PREPARATION" in found:
        source["data_preparation"] = found["DATA PREPARATION"]
    if "PLOT" in found:
        source["plot"] = found["PLOT"]
    return source


def r_string(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def patch_bootstrap(text: str, common_r: str) -> str:
    match = BOOTSTRAP.search(text)
    if match is None or "common.R" not in match.group(0):
        raise ValueError("The template has no common.R bootstrap block to replace.")
    replacement = f"source({r_string(common_r)})\n"
    return text[: match.start()] + replacement + text[match.end() :]


def refresh_project_common_r(script_path: str, common_r: str) -> str:
    """Rebind a missing installed helper without resetting an adapted plot."""
    with open(script_path, encoding="utf-8", newline="") as handle:
        text = handle.read()
    if not os.path.isfile(common_r):
        return text
    source = re.compile(
        r"""^[ \t]*source\([ \t]*(?:file[ \t]*=[ \t]*)?(?P<path>"(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*')(?=[ \t]*[,\)])""",
        re.MULTILINE,
    )

    def rebind(match: re.Match[str]) -> str:
        try:
            path = ast.literal_eval(match.group("path"))
        except (SyntaxError, ValueError):
            return match.group(0)
        absolute = os.path.isabs(path) or ntpath.isabs(path)
        installed = absolute and path.replace("\\", "/").endswith(
            "/skills/omics-visualization/scripts/lib/common.R"
        )
        if not installed or os.path.exists(path):
            return match.group(0)
        prefix = match.group(0)[: match.start("path") - match.start()]
        return prefix + r_string(common_r)

    updated = source.sub(rebind, text)
    if updated != text:
        with open(script_path, "w", encoding="utf-8", newline="") as handle:
            handle.write(updated)
        note("Updated a missing installed common.R reference in the project plotting script.")
    return updated


def find_plot_scripts(directory: str, depth: int) -> list[str]:
    if depth > MAX_SCAN_DEPTH:
        return []
    found: list[str] = []
    try:
        names = sorted(os.listdir(directory))
    except OSError:
        return []
    for name in names:
        full = os.path.join(directory, name)
        try:
            if os.path.isdir(full):
                found.extend(find_plot_scripts(full, depth + 1))
            elif name == "plot.R" and os.path.isfile(full):
                found.append(full)
        except OSError:
            continue
    return found


def list_template_sources() -> list[dict[str, Any]]:
    global _TEMPLATES
    if _TEMPLATES is not None:
        return _TEMPLATES
    scripts = os.path.join(SKILL_ROOT, "scripts")
    sources: list[dict[str, Any]] = []
    if os.path.isdir(scripts):
        for path in find_plot_scripts(scripts, 0):
            try:
                text = open(path, encoding="utf-8").read()
            except OSError:
                continue
            source = parse_template_source(text)
            if not source["id"]:
                continue
            source["path"] = path
            sources.append(source)
    _TEMPLATES = sources
    return sources


def words(text: str) -> list[str]:
    return [part for part in re.split(r"[^a-z0-9]+", text.lower()) if part]


def edit_distance(left: str, right: str) -> int:
    previous = list(range(len(right) + 1))
    for row in range(1, len(left) + 1):
        current = [row]
        for column in range(1, len(right) + 1):
            current.append(
                min(
                    previous[column] + 1,
                    current[column - 1] + 1,
                    previous[column - 1] + (0 if left[row - 1] == right[column - 1] else 1),
                )
            )
        previous = current
    return previous[len(right)]


def closest_names(typed: str, candidates: list[str], limit: int) -> list[str]:
    threshold = max(2, int(len(typed) * 0.4))
    input_words = set(words(typed))
    scored: list[tuple[int, bool, int, str]] = []
    for name in candidates:
        lower_input = typed.lower()
        lower_name = name.lower()
        shared = sum(1 for word in words(name) if word in input_words)
        contains = lower_input in lower_name or lower_name in lower_input
        distance = edit_distance(lower_input, lower_name)
        if shared > 0 or contains or distance <= threshold:
            scored.append((shared, contains, distance, name))
    scored.sort(key=lambda item: (-item[0], 0 if item[1] else 1, item[2], item[3]))
    return [item[3] for item in scored[:limit]]


def copy_entry(src: str, dst: str) -> None:
    if os.path.lexists(dst):
        if os.path.isdir(dst) and not os.path.islink(dst):
            shutil.rmtree(dst)
        else:
            os.remove(dst)
    if os.path.islink(src):
        os.symlink(os.readlink(src), dst)
        return
    if os.path.isdir(src):
        shutil.copytree(src, dst, symlinks=True)
        return
    shutil.copy2(src, dst, follow_symlinks=False)


def copy_assets(template_script: str, target_dir: str) -> list[str]:
    source_dir = os.path.dirname(template_script)
    names = sorted(
        name
        for name in os.listdir(source_dir)
        if name not in NOT_ASSETS and os.path.splitext(name)[1].lower() not in EXAMPLE_DATA_EXTENSIONS
    )
    for name in names:
        copy_entry(os.path.join(source_dir, name), os.path.join(target_dir, name))
    return names


def line_section(section: dict[str, Any] | None) -> dict[str, Any] | None:
    if section is None:
        return None
    return {"start_line": section["start_line"], "end_line": section["end_line"], "text": section["text"]}


def terms_for(purpose: str) -> list[str]:
    return [term for term in re.split(r"[\s,;，；。:：/]+", purpose.lower()) if term and term not in GENERIC_TERMS]


def score_template(template: dict[str, Any], terms: list[str]) -> int:
    if not terms:
        return 1
    labels = [
        str(template.get("id", "")),
        str(template.get("title", "")),
        str(template.get("family", "")),
    ]
    keywords = template.get("intent_keywords") or []
    labels.extend(str(keyword) for keyword in keywords)
    lowered = [label.lower() for label in labels]
    use_when = str(template.get("use_when") or "").lower()
    score = 0
    for term in terms:
        if any(term in label or label in term for label in lowered):
            score += 4
        elif term in use_when:
            score += 1
    return score


def is_installed_preview(path: str) -> bool:
    if not os.path.isabs(path) or os.path.basename(path) != "preview.png":
        return False
    try:
        root = os.path.realpath(SKILL_ROOT)
        real = os.path.realpath(path)
        parts = os.path.relpath(real, root).split(os.sep)
        return (
            len(parts) == 4
            and parts[0] == "scripts"
            and parts[3] == "preview.png"
            and os.path.isfile(real)
        )
    except OSError:
        return False


def installed_preview(template: dict[str, Any]) -> str | None:
    preview = template.get("preview")
    if not isinstance(preview, str) or not preview or os.path.isabs(preview):
        return None
    path = os.path.normpath(os.path.join(SKILL_ROOT, preview))
    if not is_installed_preview(path):
        return None
    return os.path.realpath(path)


def figure_title(script_path: str, output_path: str) -> str:
    try:
        source = parse_template_source(open(script_path, encoding="utf-8").read())
        title = str(source.get("title") or "").strip()
    except OSError:
        title = ""
    if not title:
        title = os.path.basename(output_path)
    return title[:TITLE_MAX]


def write_descriptor(
    *,
    cwd: str,
    script_path: str,
    inputs: list[str],
    output_path: str,
    timeout_seconds: int | None,
) -> None:
    extension = os.path.splitext(output_path)[1].lower()
    figure: dict[str, Any] = {"format": extension[1:]}
    if extension == ".png":
        size = png_size(output_path)
        if size is not None:
            figure["widthPx"] = size[0]
            figure["heightPx"] = size[1]
    provenance: dict[str, Any] = {
        "createdAt": datetime.now(timezone.utc).isoformat(),
        "tool": "viz_render",
        "skill": "omics-visualization",
        "inputs": [
            {"path": project_relative(cwd, path), "sha256": sha256_file(path)} for path in inputs
        ],
        "script": project_relative(cwd, script_path),
    }
    if timeout_seconds is not None:
        provenance["parameters"] = {"timeout_seconds": timeout_seconds}
    descriptor = {
        "contractVersion": "1.0.0",
        "kind": "figure",
        "file": os.path.basename(output_path),
        "mediaType": MEDIA_TYPES[extension],
        "title": figure_title(script_path, output_path),
        "figure": figure,
        "provenance": provenance,
    }
    target = output_path + ".phi-artifact.json"
    temporary = target + ".tmp"
    try:
        with open(temporary, "w", encoding="utf-8") as handle:
            json.dump(descriptor, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        os.replace(temporary, target)
    except Exception:
        if os.path.exists(temporary):
            try:
                os.remove(temporary)
            except OSError:
                pass
        raise
