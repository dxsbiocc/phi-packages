#!/usr/bin/env python3
"""Create / list literature collections; init library.sqlite under collection path."""
from __future__ import annotations

import argparse
import json
import re
import sys
import uuid
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "lib"))
sys.path.insert(0, str(_HERE))
from db import init_db, open_db, utc_now  # noqa: E402
from wiki import init_vault  # noqa: E402


def _slug(name: str) -> str:
    s = re.sub(r"[^\w\-]+", "-", name.strip(), flags=re.UNICODE).strip("-").lower()
    return s[:64] or uuid.uuid4().hex[:8]


def create_collection(
    project_root: Path,
    name: str,
    description: str = "",
    *,
    init_wiki: bool = True,
) -> dict:
    project_root = project_root.resolve()
    project_root.mkdir(parents=True, exist_ok=True)
    db_path = project_root / "library.sqlite"
    init = init_db(db_path)
    cid = f"col-{_slug(name)}-{uuid.uuid4().hex[:6]}"
    col_dir = project_root / "collections" / cid
    col_dir.mkdir(parents=True, exist_ok=True)
    (col_dir / "runs").mkdir(exist_ok=True)
    (col_dir / "files").mkdir(exist_ok=True)
    now = utc_now()
    conn = open_db(db_path)
    try:
        conn.execute(
            "INSERT INTO collections(id, name, description, created_at, path) VALUES (?, ?, ?, ?, ?)",
            (cid, name, description or None, now, str(col_dir)),
        )
        conn.commit()
    finally:
        conn.close()

    wiki_out = None
    if init_wiki:
        wiki_out = init_vault(col_dir, collection_id=cid, collection_name=name)

    return {
        "ok": True,
        "mode": "live",
        "action": "create",
        "collection_id": cid,
        "name": name,
        "path": str(col_dir),
        "db_path": str(db_path),
        "fts": init.get("fts"),
        "wiki_init": wiki_out,
        "wiki_vault": wiki_out.get("vault_root") if wiki_out else None,
    }


def list_collections(project_root: Path) -> dict:
    db_path = project_root.resolve() / "library.sqlite"
    if not db_path.exists():
        return {"ok": True, "mode": "live", "action": "list", "collections": [], "db_path": str(db_path)}
    conn = open_db(db_path)
    try:
        rows = conn.execute(
            "SELECT c.id, c.name, c.description, c.created_at, c.path, "
            "(SELECT COUNT(*) FROM collection_pubs cp WHERE cp.collection_id = c.id) AS n_papers "
            "FROM collections c ORDER BY c.created_at"
        ).fetchall()
        cols = [dict(r) for r in rows]
    finally:
        conn.close()
    return {"ok": True, "mode": "live", "action": "list", "collections": cols, "db_path": str(db_path)}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--action", choices=["create", "list"], required=True)
    p.add_argument("--project_root", required=True)
    p.add_argument("--name", default="")
    p.add_argument("--description", default="")
    p.add_argument(
        "--no_wiki",
        action="store_true",
        help="Skip automatic wiki-vault init on create (default: init wiki)",
    )
    args = p.parse_args()
    root = Path(args.project_root)
    if args.action == "create":
        if not args.name.strip():
            print(json.dumps({"ok": False, "error": "name required for create"}, ensure_ascii=False))
            return 1
        out = create_collection(
            root,
            args.name.strip(),
            args.description,
            init_wiki=not args.no_wiki,
        )
    else:
        out = list_collections(root)
    print(json.dumps(out, ensure_ascii=False))
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
