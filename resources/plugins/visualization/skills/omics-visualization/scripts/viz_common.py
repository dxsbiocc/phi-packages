"""JSON, path, and process helpers for the omics visualization CLI."""

from __future__ import annotations

import errno
import hashlib
import json
import os
import re
import struct
import subprocess
import sys
from typing import Any, NoReturn

# The skill directory is read-only when bundled and must stay free of run leftovers.
sys.dont_write_bytecode = True

SKILL_ROOT = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
COMMON_R = os.path.join(SKILL_ROOT, "scripts", "lib", "common.R")
ROUTER = os.path.join(SKILL_ROOT, "scripts", "route_template.py")
QA_SCRIPT = os.path.join(SKILL_ROOT, "scripts", "qa_single_plot.py")

MAX_CAPTURE_CHARS = 200_000
ERROR_TAIL_CHARS = 1200
MEDIA_TYPES = {".png": "image/png", ".pdf": "application/pdf", ".svg": "image/svg+xml"}
MISSING_PACKAGE = re.compile(r"there is no package called [‘'\"]([^’'\"]+)[’'\"]")
UTF8_LOCALE = re.compile(r"utf-?8", re.IGNORECASE)
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"

def dump_json(payload: object) -> str:
    """Compact JSON, matching JSON.stringify so stdout round-trips in the tests."""
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"))


def fail(message: str) -> NoReturn:
    sys.stdout.write(dump_json({"error": message}) + "\n")
    sys.stdout.flush()
    raise SystemExit(1)


def succeed(payload: dict[str, Any]) -> None:
    sys.stdout.write(dump_json(payload) + "\n")
    sys.stdout.flush()


def note(text: str) -> None:
    if text:
        sys.stderr.write(text if text.endswith("\n") else text + "\n")


def clamp(value: int, low: int, high: int) -> int:
    return min(max(value, low), high)


def is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def json_number(value: object, default: int | float = 0) -> int | float:
    """JSON numbers that are whole stay ints, so stdout matches JSON.stringify."""
    if not is_number(value):
        return default
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def as_string(value: object, default: str = "") -> str:
    if value is None:
        return default
    return str(value)


def strings(value: object) -> list[str]:
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, str)]


def string_map(value: object) -> dict[str, str]:
    if not isinstance(value, dict):
        return {}
    return {key: item for key, item in value.items() if isinstance(item, str)}


def one_line(value: str) -> str:
    return re.sub(r"\s*\n\s*", " ", value).strip()


def node_resolve(cwd: str, requested: str) -> str:
    if os.path.isabs(requested):
        return os.path.normpath(requested)
    return os.path.normpath(os.path.join(cwd, requested))


def nearest_existing(path: str) -> str:
    current = path
    while not os.path.exists(current):
        parent = os.path.dirname(current)
        if parent == current:
            return current
        current = parent
    return current


def resolve_inside_project(cwd: str, requested: str, what: str) -> str:
    target = node_resolve(cwd, requested)
    real_cwd = os.path.realpath(cwd)
    anchor = nearest_existing(target)
    relative_tail = os.path.relpath(target, anchor)
    real_target = os.path.normpath(os.path.join(os.path.realpath(anchor), relative_tail))
    inside = os.path.relpath(real_target, real_cwd)
    if inside == ".." or inside.startswith(".." + os.sep) or os.path.isabs(inside):
        fail(
            f"{what} must be inside the project directory ({cwd}); got {requested}. "
            "Data files elsewhere can be read, but everything written goes in the project."
        )
    return target


def project_relative(cwd: str, path: str) -> str:
    relative = os.path.relpath(os.path.realpath(path), os.path.realpath(cwd))
    return relative.replace(os.sep, "/")


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decode_output(data: bytes | str | None) -> str:
    if data is None:
        return ""
    text = data if isinstance(data, str) else data.decode("utf-8", "replace")
    if len(text) > MAX_CAPTURE_CHARS:
        return text[-MAX_CAPTURE_CHARS:]
    return text


class ProcessResult:
    def __init__(
        self,
        code: int | None,
        stdout: str,
        stderr: str,
        timed_out: bool,
        spawn_error: str | None,
    ) -> None:
        self.code = code
        self.stdout = stdout
        self.stderr = stderr
        self.timed_out = timed_out
        self.spawn_error = spawn_error


def run_process(
    argv: list[str],
    *,
    cwd: str | None = None,
    env: dict[str, str] | None = None,
    timeout: float,
) -> ProcessResult:
    try:
        completed = subprocess.run(
            argv,
            cwd=cwd,
            env=env,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        return ProcessResult(None, decode_output(exc.stdout), decode_output(exc.stderr), True, None)
    except OSError as exc:
        spawn = "ENOENT" if exc.errno == errno.ENOENT else (exc.strerror or str(exc))
        return ProcessResult(None, "", "", False, spawn)
    return ProcessResult(
        completed.returncode,
        decode_output(completed.stdout),
        decode_output(completed.stderr),
        False,
        None,
    )


def utf8_locale_env() -> dict[str, str]:
    effective = os.environ.get("LC_ALL") or os.environ.get("LC_CTYPE") or os.environ.get("LANG") or ""
    if UTF8_LOCALE.search(effective):
        return {}
    locale = "en_US.UTF-8" if sys.platform == "darwin" else "C.UTF-8"
    return {"LC_ALL": locale}


def r_environment() -> dict[str, str]:
    env = dict(os.environ)
    env.update(utf8_locale_env())
    env["OMICS_VISUALIZATION_SKILL_ROOT"] = SKILL_ROOT
    return env

def png_size(path: str) -> tuple[int, int] | None:
    try:
        with open(path, "rb") as handle:
            data = handle.read(24)
    except OSError:
        return None
    if len(data) < 24 or not data.startswith(PNG_SIGNATURE) or data[12:16] != b"IHDR":
        return None
    width, height = struct.unpack(">II", data[16:24])
    if width < 1 or height < 1:
        return None
    return width, height


def summarize_qa(stdout: str) -> dict[str, Any]:
    try:
        parsed = json.loads(stdout)
    except json.JSONDecodeError:
        return {"qa": {"ok": False, "failed": ["qa_unavailable"]}}
    checks = parsed.get("checks") if isinstance(parsed, dict) else None
    if isinstance(checks, list):
        records = [check for check in checks if isinstance(check, dict)]
    else:
        records = []
    failed = [str(check.get("name")) for check in records if check.get("ok") is False]
    dimensions = next((check for check in records if isinstance(check.get("format"), str)), None)
    summary: dict[str, Any] = {
        "qa": {"ok": isinstance(parsed, dict) and parsed.get("ok") is True, "failed": failed}
    }
    if dimensions is not None:
        summary["format"] = str(dimensions["format"])
        width = dimensions.get("width")
        height = dimensions.get("height")
        if is_number(width):
            summary["width"] = json_number(width)
        if is_number(height):
            summary["height"] = json_number(height)
    return summary


def r_failure(stderr: str) -> str:
    missing = list(dict.fromkeys(MISSING_PACKAGE.findall(stderr)))
    tail = stderr.strip()[-ERROR_TAIL_CHARS:]
    if missing:
        head = (
            f"R package(s) not installed: {', '.join(missing)}. "
            "Do not install packages; report the missing dependency.\n"
        )
    else:
        head = "R failed.\n"
    return head + tail


