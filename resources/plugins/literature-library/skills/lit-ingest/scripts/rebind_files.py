#!/usr/bin/env python3
"""Rebind pub_files.file_path by matching content_hash after files move.

Dry-run by default; pass --confirm to write. Stdlib only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Optional

_DB = Path(__file__).resolve().parents[2] / "lit-db" / "scripts"
sys.path.insert(0, str(_DB))
from db import open_db, utc_now  # noqa: E402


def file_content_hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _scan_hashes(scan_dir: Path) -> dict[str, list[Path]]:
    """Map content_hash → list of paths under scan_dir (files only)."""
    out: dict[str, list[Path]] = {}
    if not scan_dir.is_dir():
        return out
    for path in sorted(scan_dir.rglob("*")):
        if not path.is_file():
            continue
        # skip huge / obvious non-paper noise lightly: still hash; caller scopes dir
        try:
            digest = file_content_hash(path)
        except OSError:
            continue
        out.setdefault(digest, []).append(path.resolve())
    return out


def rebind_pub_files(
    conn,
    *,
    scan_dir: Path,
    collection_id: Optional[str] = None,
    confirm: bool = False,
    sample_limit: int = 10,
) -> dict[str, Any]:
    """Match scanned hashes to pub_files; update file_path when missing/different."""
    hash_map = _scan_hashes(scan_dir)
    n_scanned = sum(len(v) for v in hash_map.values())

    sql = "SELECT id, publication_id, file_path, content_hash, ingest_origin FROM pub_files"
    params: list[Any] = []
    if collection_id:
        sql += (
            " WHERE publication_id IN "
            "(SELECT publication_id FROM collection_pubs WHERE collection_id=?)"
        )
        params.append(collection_id)
    rows = [dict(r) for r in conn.execute(sql, params).fetchall()]

    matched: list[dict[str, Any]] = []
    updated: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []

    for row in rows:
        ch = (row.get("content_hash") or "").strip()
        if not ch:
            continue
        candidates = hash_map.get(ch) or []
        if not candidates:
            continue
        old_path = str(row.get("file_path") or "")
        old_exists = bool(old_path) and Path(old_path).is_file()
        # Prefer a path different from missing/stale current; if current still exists and
        # hash matches, no update needed.
        if old_exists:
            try:
                if file_content_hash(Path(old_path)) == ch:
                    skipped.append(
                        {
                            "publication_id": row["publication_id"],
                            "file_path": old_path,
                            "reason": "current_path_ok",
                        }
                    )
                    continue
            except OSError:
                pass
        # Choose first candidate not equal to old_path (or first if old missing)
        new_path = None
        for cand in candidates:
            if str(cand) != old_path:
                new_path = cand
                break
        if new_path is None:
            new_path = candidates[0]
        if str(new_path) == old_path and old_exists:
            continue

        entry = {
            "pub_files_id": row["id"],
            "publication_id": row["publication_id"],
            "old_path": old_path,
            "new_path": str(new_path),
            "content_hash": ch,
            "old_exists": old_exists,
        }
        matched.append(entry)

        if not confirm:
            continue

        # UNIQUE(publication_id, file_path): if new path already registered for same pub, drop old
        clash = conn.execute(
            "SELECT id FROM pub_files WHERE publication_id=? AND file_path=? AND id!=?",
            (row["publication_id"], str(new_path), row["id"]),
        ).fetchone()
        if clash:
            conn.execute("DELETE FROM pub_files WHERE id=?", (row["id"],))
            entry["action"] = "dropped_duplicate_kept_existing"
        else:
            conn.execute(
                "UPDATE pub_files SET file_path=? WHERE id=?",
                (str(new_path), row["id"]),
            )
            entry["action"] = "updated"
        updated.append(entry)

    if confirm and updated:
        conn.commit()

    return {
        "ok": True,
        "action": "rebind_files",
        "dry_run": not confirm,
        "scan_dir": str(scan_dir.resolve()),
        "collection_id": collection_id,
        "n_scanned": n_scanned,
        "n_pub_file_rows": len(rows),
        "n_matched": len(matched),
        "n_updated": len(updated) if confirm else 0,
        "n_would_update": len(matched) if not confirm else len(updated),
        "samples": (updated if confirm else matched)[:sample_limit],
        "n_skipped_ok": len(skipped),
        "checked_at": utc_now(),
    }


def main() -> int:
    p = argparse.ArgumentParser(
        description="Rebind pub_files paths by content_hash after files move (dry-run default)."
    )
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    p.add_argument("--scan_dir", required=True, help="Directory to scan recursively for files")
    p.add_argument("--confirm", action="store_true", help="Write path updates (default: dry-run)")
    p.add_argument("--sample_limit", type=int, default=10)
    args = p.parse_args()

    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    if not db_path.exists():
        print(json.dumps({"ok": False, "error": f"missing {db_path}"}))
        return 1
    scan_dir = Path(args.scan_dir).resolve()
    if not scan_dir.is_dir():
        print(json.dumps({"ok": False, "error": f"scan_dir not a directory: {scan_dir}"}))
        return 1

    conn = open_db(db_path)
    try:
        out = rebind_pub_files(
            conn,
            scan_dir=scan_dir,
            collection_id=args.collection_id or None,
            confirm=bool(args.confirm),
            sample_limit=args.sample_limit,
        )
        out["db_path"] = str(db_path)
        out["project_root"] = str(root)
    finally:
        conn.close()

    print(json.dumps(out, ensure_ascii=False))
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
