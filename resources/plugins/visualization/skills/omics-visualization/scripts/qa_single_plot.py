#!/usr/bin/env python3
"""Lightweight single-plot artifact QA for omics visualization outputs."""

from __future__ import annotations

import argparse
import json
import re
import struct
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def parse_svg_dimension(value: str | None) -> float | None:
    if not value:
        return None
    match = re.match(r"^\s*([0-9]*\.?[0-9]+)", value)
    if not match:
        return None
    return float(match.group(1))


def inspect_png(path: Path) -> dict[str, Any]:
    data = path.read_bytes()
    result: dict[str, Any] = {"format": "png"}
    if len(data) < 24 or not data.startswith(PNG_SIGNATURE):
        result.update({"ok": False, "error": "invalid PNG signature or missing IHDR"})
        return result
    chunk_type = data[12:16]
    if chunk_type != b"IHDR":
        result.update({"ok": False, "error": "first PNG chunk is not IHDR"})
        return result
    width, height = struct.unpack(">II", data[16:24])
    result.update({"ok": width > 0 and height > 0, "width": width, "height": height})
    if width <= 0 or height <= 0:
        result["error"] = "PNG width and height must be positive"
    return result


def inspect_svg(path: Path) -> dict[str, Any]:
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        return {"format": "svg", "ok": False, "error": f"invalid SVG XML: {exc}"}
    if not root.tag.lower().endswith("svg"):
        return {"format": "svg", "ok": False, "error": "root element is not svg"}
    width = parse_svg_dimension(root.attrib.get("width"))
    height = parse_svg_dimension(root.attrib.get("height"))
    view_box = root.attrib.get("viewBox") or root.attrib.get("viewbox")
    if (width is None or height is None) and view_box:
        parts = [parse_svg_dimension(part) for part in re.split(r"[\s,]+", view_box.strip())]
        if len(parts) == 4 and parts[2] is not None and parts[3] is not None:
            width = parts[2]
            height = parts[3]
    ok = bool(width and height and width > 0 and height > 0)
    result: dict[str, Any] = {
        "format": "svg",
        "ok": ok,
        "width": width,
        "height": height,
        "viewBox": view_box,
    }
    if not ok:
        result["error"] = "SVG needs positive width/height or a valid viewBox"
    return result


def inspect_file(path: Path, min_width: float, min_height: float) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    if not path.exists():
        return {"path": str(path), "ok": False, "checks": [{"ok": False, "error": "file does not exist"}]}
    size = path.stat().st_size
    checks.append({"name": "exists", "ok": True})
    checks.append({"name": "non_empty", "ok": size > 0, "size_bytes": size})
    if size <= 0:
        return {"path": str(path), "ok": False, "checks": checks}

    suffix = path.suffix.lower()
    if suffix == ".png":
        format_check = inspect_png(path)
    elif suffix == ".svg":
        format_check = inspect_svg(path)
    else:
        format_check = {"format": suffix.lstrip(".") or "unknown", "ok": True, "warning": "dimension check skipped"}
    format_check["name"] = "format"
    checks.append(format_check)

    width = format_check.get("width")
    height = format_check.get("height")
    if width is not None and height is not None:
        checks.append(
            {
                "name": "min_dimensions",
                "ok": width >= min_width and height >= min_height,
                "min_width": min_width,
                "min_height": min_height,
                "width": width,
                "height": height,
            }
        )

    return {"path": str(path), "ok": all(check.get("ok", False) for check in checks), "checks": checks}


def render_text(payload: dict[str, Any]) -> str:
    status = "PASS" if payload["ok"] else "FAIL"
    lines = [f"{status}: {payload['path']}"]
    for check in payload["checks"]:
        check_status = "ok" if check.get("ok") else "failed"
        detail = ""
        if "width" in check and "height" in check:
            detail = f" ({check['width']} x {check['height']})"
        elif "size_bytes" in check:
            detail = f" ({check['size_bytes']} bytes)"
        elif "error" in check:
            detail = f" ({check['error']})"
        lines.append(f"- {check.get('name', 'check')}: {check_status}{detail}")
    return "\n".join(lines)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run lightweight QA on a rendered single-plot artifact.")
    parser.add_argument("path", help="Rendered PNG, SVG, PDF, or other artifact path.")
    parser.add_argument("--min-width", type=float, default=32, help="Minimum expected width for dimension-aware formats.")
    parser.add_argument("--min-height", type=float, default=32, help="Minimum expected height for dimension-aware formats.")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    payload = inspect_file(Path(args.path).expanduser().resolve(), args.min_width, args.min_height)
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print(render_text(payload))
    return 0 if payload["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
