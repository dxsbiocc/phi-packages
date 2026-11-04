#!/usr/bin/env python3
"""Optional MCP *client* for phi-literature search_publications / related_publications.

This plugin does NOT ship an MCP server. Default search remains providers_http.
When --backend mcp or LIT_SEARCH_BACKEND=mcp, call an external MCP bridge
(command configured via LIT_MCP_CMD) or a mock fixture for offline tests.

Soft-fail: if MCP is unavailable and force=False, return ok=False so caller
can fall back to HTTP. If force=True, return actionable error.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, Optional


DEFAULT_TOOLS = {
    "search": "search_publications",
    "related": "related_publications",
}


def resolve_backend(cli_backend: Optional[str] = None) -> str:
    raw = (cli_backend or os.environ.get("LIT_SEARCH_BACKEND") or "http").strip().lower()
    if raw in ("mcp", "http", "auto"):
        return raw
    return "http"


def load_mcp_mock_fixture(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _hits_from_fixture(fixture: dict[str, Any], *, tool: str) -> list[dict[str, Any]]:
    if tool == "related":
        return list(fixture.get("related_hits") or fixture.get("hits") or [])
    return list(fixture.get("search_hits") or fixture.get("hits") or [])


def call_mcp_tool(
    tool: str,
    arguments: dict[str, Any],
    *,
    mcp_cmd: Optional[str] = None,
    mock_fixture: Optional[Path] = None,
    timeout: float = 60.0,
) -> dict[str, Any]:
    """Invoke one MCP tool via mock fixture or LIT_MCP_CMD JSON-lines bridge.

    Bridge protocol (stdio): caller writes one JSON line
      {"tool": "...", "arguments": {...}}
    and reads one JSON line
      {"ok": true, "hits": [...], ...}  or {"ok": false, "error": "..."}
    """
    if mock_fixture:
        fixture = load_mcp_mock_fixture(Path(mock_fixture))
        hits = _hits_from_fixture(fixture, tool="related" if "related" in tool else "search")
        # optional filter by query key presence
        return {
            "ok": True,
            "provider": "mcp_mock",
            "tool": tool,
            "hits": hits,
            "n": len(hits),
            "fixture": str(mock_fixture),
            "arguments": arguments,
        }

    cmd = (mcp_cmd or os.environ.get("LIT_MCP_CMD") or "").strip()
    if not cmd:
        return {
            "ok": False,
            "provider": "mcp",
            "tool": tool,
            "error": "mcp_unavailable",
            "message": (
                "LIT_MCP_CMD not set and no --mcp_mock_fixture. "
                "Configure a stdio JSON-lines bridge to phi-literature, "
                "or use --backend http (default)."
            ),
            "hits": [],
        }

    payload = json.dumps({"tool": tool, "arguments": arguments}, ensure_ascii=False) + "\n"
    try:
        proc = subprocess.run(
            cmd,
            input=payload.encode("utf-8"),
            shell=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as e:
        return {
            "ok": False,
            "provider": "mcp",
            "tool": tool,
            "error": "mcp_exec_failed",
            "message": str(e),
            "hits": [],
        }

    raw = (proc.stdout or b"").decode("utf-8", errors="replace").strip()
    if not raw:
        err = (proc.stderr or b"").decode("utf-8", errors="replace")[:500]
        return {
            "ok": False,
            "provider": "mcp",
            "tool": tool,
            "error": "mcp_empty_response",
            "message": err or f"exit {proc.returncode}",
            "hits": [],
        }
    # last non-empty line
    line = [ln for ln in raw.splitlines() if ln.strip()][-1]
    try:
        data = json.loads(line)
    except json.JSONDecodeError as e:
        return {
            "ok": False,
            "provider": "mcp",
            "tool": tool,
            "error": "mcp_bad_json",
            "message": str(e),
            "hits": [],
            "raw_preview": line[:300],
        }
    hits = list(data.get("hits") or data.get("result") or [])
    ok = bool(data.get("ok", True)) and not data.get("error")
    return {
        "ok": ok,
        "provider": "mcp",
        "tool": tool,
        "hits": hits,
        "n": len(hits),
        "error": data.get("error"),
        "message": data.get("message"),
        "raw": data,
    }


def search_publications_mcp(
    query: str,
    *,
    sources: Optional[list[str]] = None,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
    mcp_cmd: Optional[str] = None,
    mock_fixture: Optional[Path] = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    tool = DEFAULT_TOOLS["search"]
    args = {
        "query": query,
        "sources": sources or [],
        "limit": limit,
        "date_from": date_from,
        "date_to": date_to,
    }
    result = call_mcp_tool(
        tool, args, mcp_cmd=mcp_cmd, mock_fixture=mock_fixture
    )
    hits = list(result.get("hits") or [])
    # tag source
    for h in hits:
        h.setdefault("source", h.get("source") or "mcp")
    meta = {k: v for k, v in result.items() if k != "hits"}
    meta["source"] = "mcp"
    return hits, meta


def related_publications_mcp(
    *,
    relation: str = "similar",
    limit: int = 10,
    pmid: Optional[str] = None,
    doi: Optional[str] = None,
    openalex: Optional[str] = None,
    pmcid: Optional[str] = None,
    arxiv: Optional[str] = None,
    mcp_cmd: Optional[str] = None,
    mock_fixture: Optional[Path] = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    tool = DEFAULT_TOOLS["related"]
    args = {
        "relation": relation,
        "limit": limit,
        "pmid": pmid,
        "doi": doi,
        "openalex": openalex,
        "pmcid": pmcid,
        "arxiv": arxiv,
    }
    result = call_mcp_tool(
        tool, args, mcp_cmd=mcp_cmd, mock_fixture=mock_fixture
    )
    hits = list(result.get("hits") or [])
    for h in hits:
        h.setdefault("source", h.get("source") or "mcp")
    meta = {k: v for k, v in result.items() if k != "hits"}
    meta["relation"] = relation
    meta["source"] = "mcp"
    return hits, meta


def mcp_available(*, mcp_cmd: Optional[str] = None, mock_fixture: Optional[Path] = None) -> bool:
    if mock_fixture and Path(mock_fixture).is_file():
        return True
    cmd = (mcp_cmd or os.environ.get("LIT_MCP_CMD") or "").strip()
    if not cmd:
        return False
    # crude: if first token is a path/binary, check which
    first = cmd.split()[0]
    if first in ("npx", "node", "python", "python3"):
        return shutil.which(first) is not None
    return True
