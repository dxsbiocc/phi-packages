"""Small helpers shared by the tool entry scripts."""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any, Optional

from db import open_db


def emit(obj: dict[str, Any]) -> None:
    """Print the single JSON object a tool must write to stdout."""
    print(json.dumps(obj, ensure_ascii=False))


def fail(message: str, **extra: Any) -> int:
    emit({"ok": False, "error": message, **extra})
    return 1


def db_path_for(project_root: str) -> Path:
    return Path(project_root).resolve() / "library.sqlite"


def open_project_db(project_root: str) -> tuple[Optional[sqlite3.Connection], Optional[str]]:
    """Open the project's library.sqlite. Returns (conn, error)."""
    path = db_path_for(project_root)
    if not path.exists():
        return None, f"missing {path}; create a collection first (lit_collection)"
    return open_db(path), None


def resolve_collection(
    conn: sqlite3.Connection, collection_id: Optional[str]
) -> tuple[Optional[sqlite3.Row], Optional[str]]:
    """Return the collection row. With no id, use the only collection if there is exactly one."""
    if collection_id:
        row = conn.execute(
            "SELECT id, name, path FROM collections WHERE id = ?", (collection_id,)
        ).fetchone()
        return (row, None) if row else (None, f"unknown collection_id {collection_id}")
    rows = conn.execute("SELECT id, name, path FROM collections ORDER BY created_at").fetchall()
    if len(rows) == 1:
        return rows[0], None
    if not rows:
        return None, "no collection yet; create one with lit_collection"
    names = ", ".join(f"{r['id']} ({r['name']})" for r in rows)
    return None, f"several collections; pass collection_id. Available: {names}"
