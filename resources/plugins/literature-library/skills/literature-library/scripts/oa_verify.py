#!/usr/bin/env python3
"""Soft OA / PDF URL verification (stdlib only).

Candidate OA URLs from Unpaywall-like fields or provider ``url`` are NOT treated
as verified storable PDFs until checked. Never downloads a paywalled body into
the library as owned full text.
"""
from __future__ import annotations

import json
import re
import urllib.error
import urllib.request
from typing import Any, Callable, Optional

OA_STATUSES = (
    "verified_oa",
    "unknown",
    "paywall_or_blocked",
    "not_pdf",
)

_PAYWALL_HINTS = (
    "login",
    "signin",
    "sign-in",
    "subscribe",
    "purchase",
    "paywall",
    "checkout",
    "entitlements",
    "getaccess",
    "get-access",
)

_HTML_TYPES = ("text/html", "application/xhtml", "text/xml")
_PDF_TYPES = ("application/pdf", "application/x-pdf")

USER_AGENT = "literature-library/0.5.1-draft (OA verify; soft HEAD/GET; mailto:research@example.com)"


def candidate_oa_urls(record: dict[str, Any]) -> list[str]:
    """Collect candidate OA/PDF URLs from a hit/draft without trusting them."""
    urls: list[str] = []
    for key in (
        "url",
        "pdf_url",
        "oa_url",
        "best_oa_url",
        "open_access_url",
        "landing_page_url",
    ):
        val = record.get(key)
        if isinstance(val, str) and val.startswith("http"):
            urls.append(val.strip())
    oa = record.get("open_access") or record.get("oa") or {}
    if isinstance(oa, dict):
        for key in ("url", "pdf_url", "oa_url", "best_oa_location"):
            val = oa.get(key)
            if isinstance(val, str) and val.startswith("http"):
                urls.append(val.strip())
            elif isinstance(val, dict):
                u = val.get("url") or val.get("pdf_url")
                if isinstance(u, str) and u.startswith("http"):
                    urls.append(u.strip())
    best = record.get("best_oa_location")
    if isinstance(best, dict):
        u = best.get("url") or best.get("pdf_url")
        if isinstance(u, str) and u.startswith("http"):
            urls.append(u.strip())
    # dedupe preserve order
    seen: set[str] = set()
    out: list[str] = []
    for u in urls:
        if u not in seen:
            seen.add(u)
            out.append(u)
    return out


def _classify_from_response(
    *,
    url: str,
    final_url: str,
    status: int,
    content_type: str,
    content_disp: str = "",
) -> str:
    ct = (content_type or "").lower()
    cd = (content_disp or "").lower()
    low_url = (final_url or url or "").lower()

    if status in (401, 402, 403, 407):
        return "paywall_or_blocked"
    if status == 404:
        return "unknown"
    if any(h in low_url for h in _PAYWALL_HINTS):
        return "paywall_or_blocked"

    looks_pdf = (
        any(p in ct for p in _PDF_TYPES)
        or ".pdf" in low_url.split("?")[0]
        or "filename=" in cd and ".pdf" in cd
    )
    looks_html = any(h in ct for h in _HTML_TYPES) or ct.startswith("text/")

    if 200 <= status < 300:
        if looks_pdf and not looks_html:
            return "verified_oa"
        if looks_html:
            return "not_pdf"
        if looks_pdf:
            return "verified_oa"
        return "unknown"
    if 300 <= status < 400:
        return "unknown"
    if status >= 500:
        return "unknown"
    return "unknown"


def verify_oa_url(
    url: str,
    *,
    timeout: float = 8.0,
    mock: Optional[dict[str, Any]] = None,
    opener: Optional[Callable[..., Any]] = None,
) -> dict[str, Any]:
    """Soft-verify a candidate OA URL via HEAD (then tiny GET).

    Returns dict with ``oa_status`` in OA_STATUSES. Never returns PDF bytes.
    ``mock`` short-circuits with a fixture dict (must include oa_status or status/code).
    """
    if not url or not str(url).startswith("http"):
        return {
            "ok": False,
            "url": url,
            "oa_status": "unknown",
            "message": "missing or non-http url",
        }

    if mock is not None:
        out = dict(mock)
        out.setdefault("url", url)
        status = out.get("oa_status") or out.get("status")
        if status in OA_STATUSES:
            out["oa_status"] = status
        elif isinstance(out.get("http_status"), int):
            out["oa_status"] = _classify_from_response(
                url=url,
                final_url=str(out.get("final_url") or url),
                status=int(out["http_status"]),
                content_type=str(out.get("content_type") or ""),
                content_disp=str(out.get("content_disposition") or ""),
            )
        else:
            out["oa_status"] = "unknown"
        out["ok"] = True
        out["mock"] = True
        out["stored_pdf_bytes"] = False
        return out

    headers = {"User-Agent": USER_AGENT, "Accept": "application/pdf,*/*;q=0.8"}
    http_status = 0
    content_type = ""
    content_disp = ""
    final_url = url
    method_used = "HEAD"
    err_msg = None

    def _do(req: urllib.request.Request) -> tuple[int, str, str, str]:
        open_fn = opener or urllib.request.urlopen
        with open_fn(req, timeout=timeout) as resp:  # noqa: S310 — intentional soft probe
            st = getattr(resp, "status", None) or resp.getcode()
            ct = resp.headers.get("Content-Type", "") if resp.headers else ""
            cd = resp.headers.get("Content-Disposition", "") if resp.headers else ""
            fu = getattr(resp, "geturl", lambda: url)()
            # Do NOT read body for HEAD; for GET read at most 0 bytes intent —
            # we only need headers. If body starts streaming, close immediately.
            return int(st), str(ct), str(cd), str(fu)

    try:
        req = urllib.request.Request(url, method="HEAD", headers=headers)
        http_status, content_type, content_disp, final_url = _do(req)
    except urllib.error.HTTPError as e:
        http_status = int(e.code)
        content_type = e.headers.get("Content-Type", "") if e.headers else ""
        content_disp = e.headers.get("Content-Disposition", "") if e.headers else ""
        final_url = getattr(e, "url", None) or url
        err_msg = str(e.reason)
    except Exception as e:  # noqa: BLE001
        # Fall through to GET — some servers reject HEAD
        err_msg = str(e)
        http_status = 0

    need_get = http_status == 0 or (
        200 <= http_status < 300
        and not content_type
    ) or http_status in (405, 501)

    if need_get:
        method_used = "GET"
        try:
            # Range request: ask for first byte only so we never pull a full PDF
            get_headers = dict(headers)
            get_headers["Range"] = "bytes=0-0"
            req = urllib.request.Request(url, method="GET", headers=get_headers)
            http_status, content_type, content_disp, final_url = _do(req)
            err_msg = None
        except urllib.error.HTTPError as e:
            http_status = int(e.code)
            content_type = e.headers.get("Content-Type", "") if e.headers else ""
            content_disp = e.headers.get("Content-Disposition", "") if e.headers else ""
            final_url = getattr(e, "url", None) or url
            err_msg = str(e.reason)
        except Exception as e:  # noqa: BLE001
            return {
                "ok": False,
                "url": url,
                "oa_status": "unknown",
                "message": f"request failed: {e}",
                "method": method_used,
            }

    oa_status = _classify_from_response(
        url=url,
        final_url=final_url,
        status=http_status,
        content_type=content_type,
        content_disp=content_disp,
    )
    return {
        "ok": True,
        "url": url,
        "final_url": final_url,
        "http_status": http_status,
        "content_type": content_type,
        "content_disposition": content_disp,
        "method": method_used,
        "oa_status": oa_status,
        "message": err_msg,
        # Explicit: verification never yields storable PDF bytes
        "stored_pdf_bytes": False,
    }


def verify_record_oa(
    record: dict[str, Any],
    *,
    mock_map: Optional[dict[str, dict[str, Any]]] = None,
    timeout: float = 8.0,
) -> dict[str, Any]:
    """Verify first candidate URL on a record; return aggregated result."""
    urls = candidate_oa_urls(record)
    if not urls:
        return {
            "ok": True,
            "oa_status": "unknown",
            "url": None,
            "message": "no candidate oa url",
            "candidates": [],
        }
    url = urls[0]
    mock = None
    if mock_map:
        mock = mock_map.get(url) or mock_map.get("*")
    result = verify_oa_url(url, timeout=timeout, mock=mock)
    result["candidates"] = urls
    return result


def dumps_oa(result: dict[str, Any]) -> str:
    return json.dumps(result, ensure_ascii=False)
