#!/usr/bin/env python3
"""Campus / SSO authorized full-text pathway (scaffold).

Requires explicit --institutional_access. Own-account only:
  - user_provided_file: mount an already-authorized local PDF
  - vpn_proxy / session_reuse: optional fetch via LIT_INSTITUTIONAL_PROXY
    (HTTP(S) proxy or gateway base URL). Credentials stay in the user vault —
    never stored in SQLite / plugin config.
  - saml_browser: record method only; Playwright/SAML is a host/agent concern.

No paywall bypass. No borrowed campus accounts. No commercial gateways.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Optional

_DB = Path(__file__).resolve().parent
_INGEST = Path(__file__).resolve().parent
sys.path.insert(0, str(_DB))
sys.path.insert(0, str(_INGEST))

from db import link_to_collection, open_db, record_pub_file, resolve_or_create, utc_now  # noqa: E402
from ingest import (  # noqa: E402
    extract_ids_from_filename,
    extract_ids_from_text,
    file_content_hash,
    load_records,
)

ACCESS_METHODS = ("vpn_proxy", "session_reuse", "saml_browser", "user_provided_file")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch_via_institutional_proxy(
    url: str,
    *,
    proxy_base: Optional[str] = None,
    mock_body: Optional[bytes] = None,
    timeout: float = 30.0,
) -> dict[str, Any]:
    """Fetch URL through user-configured proxy. Credentials never passed here.

    proxy_base: from LIT_INSTITUTIONAL_PROXY — either a full proxy URL
    (http://host:port) used as HTTPS/HTTP proxy, or empty for direct
    (session_reuse assumes ambient VPN / cookies outside this process).
    """
    if mock_body is not None:
        return {
            "ok": True,
            "body": mock_body,
            "content_type": "application/pdf",
            "mode": "mock",
            "url": url,
        }
    proxy = (proxy_base or os.environ.get("LIT_INSTITUTIONAL_PROXY") or "").strip()
    handlers: list[Any] = []
    if proxy:
        handlers.append(
            urllib.request.ProxyHandler({"http": proxy, "https": proxy})
        )
    opener = urllib.request.build_opener(*handlers) if handlers else urllib.request.build_opener()
    req = urllib.request.Request(url, headers={"User-Agent": "phi-literature-library/0.6"})
    try:
        with opener.open(req, timeout=timeout) as resp:
            body = resp.read()
            ctype = resp.headers.get("Content-Type") or ""
            return {
                "ok": True,
                "body": body,
                "content_type": ctype,
                "mode": "proxy" if proxy else "direct",
                "url": url,
                "proxy_configured": bool(proxy),
            }
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
        return {
            "ok": False,
            "error": str(e),
            "url": url,
            "proxy_configured": bool(proxy),
            "mode": "proxy" if proxy else "direct",
        }


def mount_authorized_file(
    *,
    db_path: Path,
    collection_id: str,
    file_path: Path,
    access_method: str,
    access_origin: str = "institutional",
    doi: Optional[str] = None,
    pmid: Optional[str] = None,
    title: Optional[str] = None,
    try_ocr: bool = True,
    ocr_cmd: Optional[str] = None,
) -> dict[str, Any]:
    """Ingest/mount a user-authorized local PDF with institutional metadata."""
    if access_method not in ACCESS_METHODS:
        return {"ok": False, "error": f"unsupported access_method {access_method}"}
    if not file_path.is_file():
        return {"ok": False, "error": f"missing file {file_path}"}

    records = load_records(file_path, fmt="pdf", try_ocr=try_ocr, ocr_cmd=ocr_cmd)
    if not records:
        return {"ok": False, "error": "could not parse local file"}
    rec = dict(records[0])
    if doi:
        rec["doi"] = doi
    if pmid:
        rec["pmid"] = pmid
    if title:
        rec["title"] = title
    rec["access_origin"] = access_origin
    rec["access_method"] = access_method
    rec["ingest_origin"] = "institutional"
    rec["_local_file"] = str(file_path.resolve())

    conn = open_db(db_path)
    try:
        col = conn.execute(
            "SELECT id FROM collections WHERE id=?", (collection_id,)
        ).fetchone()
        if not col:
            return {"ok": False, "error": f"unknown collection_id {collection_id}"}
        result = resolve_or_create(conn, rec)
        pub_id = result["publication_id"]
        link_to_collection(conn, collection_id, pub_id)
        ch = None
        try:
            ch = file_content_hash(file_path)
        except OSError:
            ch = None
        meta = record_pub_file(
            conn,
            publication_id=pub_id,
            file_path=str(file_path.resolve()),
            content_hash=ch,
            ingest_origin="institutional",
            access_origin=access_origin,
            access_method=access_method,
        )
        conn.commit()
        return {
            "ok": True,
            "action": "institutional_mount",
            "publication_id": pub_id,
            "resolve": result,
            "file": meta,
            "access_origin": access_origin,
            "access_method": access_method,
            "extract_method": rec.get("extract_method"),
            "extract_warnings": rec.get("extract_warnings") or [],
        }
    finally:
        conn.close()


def institutional_fetch_and_mount(
    *,
    db_path: Path,
    collection_id: str,
    url: str,
    dest_dir: Path,
    access_method: str,
    access_origin: str = "institutional",
    doi: Optional[str] = None,
    pmid: Optional[str] = None,
    title: Optional[str] = None,
    mock_fixture: Optional[Path] = None,
    proxy_base: Optional[str] = None,
) -> dict[str, Any]:
    """Fetch via institutional proxy (or mock) → save locally → mount.

    Never stores credentials. For saml_browser without a local file, returns
    actionable guidance instead of launching a browser.
    """
    if access_method == "saml_browser" and not mock_fixture:
        return {
            "ok": False,
            "action": "institutional_fetch",
            "error": "saml_browser_host_concern",
            "message": (
                "SAML/CAS web SSO is a host/agent concern (Playwright). "
                "Complete login in the user's own session, save the authorized PDF "
                "locally, then re-run with --method user_provided_file --file_path …"
            ),
            "access_method": access_method,
        }
    if access_method not in ("vpn_proxy", "session_reuse") and not mock_fixture:
        return {
            "ok": False,
            "error": f"fetch not applicable for method {access_method}; use --file_path",
        }

    mock_body = None
    if mock_fixture:
        mock_body = Path(mock_fixture).read_bytes()

    fetched = fetch_via_institutional_proxy(
        url, proxy_base=proxy_base, mock_body=mock_body
    )
    if not fetched.get("ok"):
        return {
            "ok": False,
            "action": "institutional_fetch",
            "error": fetched.get("error") or "fetch_failed",
            "fetch": {k: v for k, v in fetched.items() if k != "body"},
            "hint": (
                "Ensure campus VPN/proxy is up, set LIT_INSTITUTIONAL_PROXY, "
                "or download the PDF yourself and use --method user_provided_file."
            ),
        }

    body: bytes = fetched["body"]
    dest_dir.mkdir(parents=True, exist_ok=True)
    # Prefer DOI-ish name when known
    name_bits = []
    if doi:
        name_bits.append(doi.replace("/", "_"))
    elif pmid:
        name_bits.append(f"PMID{pmid}")
    else:
        name_bits.append(_sha256_bytes(body)[:12])
    dest = dest_dir / ("".join(name_bits) + ".pdf")
    dest.write_bytes(body)

    mount = mount_authorized_file(
        db_path=db_path,
        collection_id=collection_id,
        file_path=dest,
        access_method=access_method,
        access_origin=access_origin,
        doi=doi,
        pmid=pmid,
        title=title,
        try_ocr=False,  # binary may be mock fixture without real PDF structure
    )
    mount["fetch"] = {k: v for k, v in fetched.items() if k != "body"}
    mount["saved_path"] = str(dest)
    mount["action"] = "institutional_fetch_and_mount"
    return mount


def main() -> int:
    p = argparse.ArgumentParser(
        description="Mount campus/SSO-authorized full text (own account only; explicit flag required)."
    )
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", required=True)
    p.add_argument(
        "--institutional_access",
        action="store_true",
        required=True,
        help="Required gate: engage institutional pathway (own account only)",
    )
    p.add_argument(
        "--method",
        required=True,
        choices=list(ACCESS_METHODS),
        help="Authorization method recorded on pub_files.access_method",
    )
    p.add_argument("--file_path", default="", help="Already-authorized local PDF path")
    p.add_argument("--url", default="", help="Publisher/fulltext URL (proxy/session fetch)")
    p.add_argument("--dest_dir", default="", help="Where to save fetched PDF (default: project papers/)")
    p.add_argument("--doi", default="")
    p.add_argument("--pmid", default="")
    p.add_argument("--title", default="")
    p.add_argument(
        "--mock_fixture",
        default="",
        help="Offline bytes fixture for fetch (CI); skips live network",
    )
    p.add_argument(
        "--proxy",
        default="",
        help="Override LIT_INSTITUTIONAL_PROXY for this call (URL only; no credentials)",
    )
    p.add_argument("--no_ocr", action="store_true")
    p.add_argument("--ocr_cmd", default="")
    args = p.parse_args()

    if not args.institutional_access:
        print(json.dumps({"ok": False, "error": "refusing without --institutional_access"}))
        return 1

    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    access_origin = "institutional"

    if args.file_path:
        out = mount_authorized_file(
            db_path=db_path,
            collection_id=args.collection_id,
            file_path=Path(args.file_path),
            access_method=args.method,
            access_origin=access_origin,
            doi=args.doi or None,
            pmid=args.pmid or None,
            title=args.title or None,
            try_ocr=not args.no_ocr,
            ocr_cmd=args.ocr_cmd or None,
        )
    elif args.url or args.mock_fixture:
        dest = Path(args.dest_dir) if args.dest_dir else (root / "papers" / "institutional")
        out = institutional_fetch_and_mount(
            db_path=db_path,
            collection_id=args.collection_id,
            url=args.url or "https://example.org/institutional/mock.pdf",
            dest_dir=dest,
            access_method=args.method,
            access_origin=access_origin,
            doi=args.doi or None,
            pmid=args.pmid or None,
            title=args.title or None,
            mock_fixture=Path(args.mock_fixture) if args.mock_fixture else None,
            proxy_base=args.proxy or None,
        )
    else:
        out = {
            "ok": False,
            "error": "provide --file_path (preferred) or --url / --mock_fixture",
            "access_methods": list(ACCESS_METHODS),
            "docs": (
                "Own campus VPN/proxy or already-authorized file only. "
                "Credentials remain in user vault — never passed to this script."
            ),
        }

    out["checked_at"] = utc_now()
    out["institutional_access"] = True
    print(json.dumps(out, ensure_ascii=False))
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
