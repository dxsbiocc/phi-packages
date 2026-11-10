"""Read text from the local files saved for a paper (PDF or text sidecar)."""
from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

from pdf_extract import extract_pdf_text

TEXT_SUFFIXES = {".txt", ".md"}


def saved_files(conn: sqlite3.Connection, publication_id: str) -> list[dict[str, Any]]:
    rows = conn.execute(
        "SELECT file_path, content_hash, ingest_origin, access_origin, access_method "
        "FROM pub_files WHERE publication_id = ? ORDER BY id",
        (publication_id,),
    ).fetchall()
    out = []
    for r in rows:
        d = dict(r)
        d["exists"] = Path(d["file_path"]).is_file()
        out.append(d)
    return out


def local_text(conn: sqlite3.Connection, publication_id: str, *, max_chars: int = 6000) -> dict[str, Any]:
    """Text of the first readable saved file. Returns {"text", "path", "method", "truncated"}."""
    for f in saved_files(conn, publication_id):
        if not f["exists"]:
            continue
        path = Path(f["file_path"])
        suffix = path.suffix.lower()
        if suffix == ".pdf":
            res = extract_pdf_text(path)
            text, method = (res.get("text") or ""), res.get("method") or "pdf"
        elif suffix in TEXT_SUFFIXES:
            text, method = path.read_text(encoding="utf-8", errors="replace"), "text"
        else:
            continue
        text = text.strip()
        if not text:
            continue
        return {
            "text": text[:max_chars],
            "path": str(path),
            "method": method,
            "truncated": len(text) > max_chars,
            "n_chars": len(text),
        }
    return {"text": "", "path": None, "method": None, "truncated": False, "n_chars": 0}
