#!/usr/bin/env python3
"""Shared SQLite library helpers for literature-library (v0.6.0).

Publication primary key (natural key on publications.id, not UUID):
  doi:… > pmid:… > pmcid:… > arxiv:… > openalex:… > soft:… > tmp:…

Schema v6: v5 + pub_files.access_origin / access_method (institutional / OA / local).
Schema v5: v4 + pub_files (local authorized files) + publications.oa_status/oa_json.
Schema v4: runs.status / coverage_path / gate_json. Schema v3: publication_versions + review_reasons_json.

Identity policies:
  - Soft fingerprint never auto-merges → needs_review.
  - True conflicting strong IDs → needs_review (+ optional #conflict- suffix).
  - Complementary multi-match clusters auto-merge to DOI-first.
  - Single strong-ID match that later gains a higher-priority ID (e.g. pmid-only
    row then DOI+PMID) reassigns publications.id via _reassign_publication_id.
  - bioRxiv/medRxiv use DOI (10.1101/…); S2 stays out of STRONG_ID_KEYS.
  - Preprint↔journal version links: propose → confirm; high-confidence same_work
    may merge into journal DOI canonical (preprint IDs kept as aliases).
  - Missing opposite DOI may be preview-fetched as metadata-only shell + pending link.
"""
from __future__ import annotations

import hashlib
import json
import re
import sqlite3
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

SCHEMA_VERSION = 6

STRONG_ID_KEYS = ("doi", "pmid", "pmcid", "arxiv", "openalex")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


# sqlite3.Connection may lack __dict__ / weakref support — key by id(conn)
_ID_REMAPS: dict[int, list[dict[str, str]]] = {}


def _record_id_remap(conn: sqlite3.Connection, old_id: str, new_id: str) -> None:
    """Track publication id renames / absorbs for wiki_repoint callers."""
    if not old_id or not new_id or old_id == new_id:
        return
    key = id(conn)
    remaps = _ID_REMAPS.get(key)
    if remaps is None:
        remaps = []
        _ID_REMAPS[key] = remaps
    remaps.append({"old_id": old_id, "new_id": new_id})


def drain_id_remaps(conn: sqlite3.Connection) -> list[dict[str, str]]:
    """Return and clear remaps recorded on this connection."""
    key = id(conn)
    remaps = list(_ID_REMAPS.pop(key, []) or [])
    return remaps


def normalize_doi(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    s = str(value).strip()
    s = re.sub(r"^https?://(dx\.)?doi\.org/", "", s, flags=re.I)
    s = re.sub(r"^doi:\s*", "", s, flags=re.I)
    s = s.strip().lower().rstrip(".")
    if not s or "/" not in s:
        return None
    return s


def normalize_pmid(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    s = str(value).strip()
    s = re.sub(r"^(pmid:|pubmed:)\s*", "", s, flags=re.I)
    digits = re.sub(r"\D", "", s)
    return digits or None


def normalize_pmcid(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    s = str(value).strip().upper()
    s = re.sub(r"^PMCID:\s*", "", s, flags=re.I)
    if s.startswith("PMC"):
        body = s[3:]
    else:
        body = s
    body = re.sub(r"\D", "", body)
    return f"PMC{body}" if body else None


def normalize_arxiv(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    s = str(value).strip()
    s = re.sub(r"^https?://arxiv\.org/(abs|pdf)/", "", s, flags=re.I)
    s = re.sub(r"^arxiv:\s*", "", s, flags=re.I)
    s = s.replace(".pdf", "").strip()
    return s.lower() or None


def normalize_openalex(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    s = str(value).strip()
    s = re.sub(r"^https?://openalex\.org/", "", s, flags=re.I)
    if s.lower().startswith("openalex:"):
        s = s.split(":", 1)[1]
    s = s.strip()
    if not s:
        return None
    if not s.upper().startswith("W"):
        s = f"W{s}" if s.isdigit() else s
    return s.upper() if s.upper().startswith("W") else s


def normalize_title(title: Optional[str]) -> str:
    if not title:
        return ""
    t = title.lower()
    t = re.sub(r"[^\w\s]", " ", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


def first_author_key(authors: Any) -> str:
    if not authors:
        return ""
    if isinstance(authors, str):
        try:
            authors = json.loads(authors)
        except json.JSONDecodeError:
            return authors.split(",")[0].strip().lower()
    if isinstance(authors, list) and authors:
        a0 = authors[0]
        if isinstance(a0, dict):
            return str(a0.get("family") or a0.get("name") or a0.get("last") or "").strip().lower()
        return str(a0).strip().lower()
    return ""


def soft_fingerprint(title: Optional[str], authors: Any, year: Any) -> Optional[str]:
    nt = normalize_title(title)
    fa = first_author_key(authors)
    y = str(year or "").strip()
    if not nt or not y:
        return None
    raw = f"{nt}|{fa}|{y}"
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()[:16]


def canonical_id_from_ids(
    *,
    doi: Optional[str] = None,
    pmid: Optional[str] = None,
    pmcid: Optional[str] = None,
    arxiv: Optional[str] = None,
    openalex: Optional[str] = None,
    soft: Optional[str] = None,
) -> str:
    if doi:
        return f"doi:{doi}"
    if pmid:
        return f"pmid:{pmid}"
    if pmcid:
        return f"pmcid:{pmcid}"
    if arxiv:
        return f"arxiv:{arxiv}"
    if openalex:
        return f"openalex:{openalex}"
    if soft:
        return f"soft:{soft}"
    return f"tmp:{uuid.uuid4().hex[:12]}"


def _connect(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(path))
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def init_db(path: str | Path) -> dict[str, Any]:
    """Create schema (idempotent). Returns {ok, path, fts}."""
    db_path = Path(path)
    conn = _connect(db_path)
    try:
        conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS meta (
              key TEXT PRIMARY KEY,
              value TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS collections (
              id TEXT PRIMARY KEY,
              name TEXT NOT NULL,
              description TEXT,
              created_at TEXT NOT NULL,
              path TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS publications (
              id TEXT PRIMARY KEY,
              doi TEXT,
              pmid TEXT,
              pmcid TEXT,
              arxiv TEXT,
              openalex TEXT,
              title TEXT,
              abstract TEXT,
              authors_json TEXT,
              year INTEGER,
              venue TEXT,
              type TEXT,
              tags_json TEXT,
              needs_review INTEGER NOT NULL DEFAULT 0,
              review_reasons_json TEXT,
              created_at TEXT NOT NULL,
              updated_at TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS publication_versions (
              id INTEGER PRIMARY KEY AUTOINCREMENT,
              from_id TEXT NOT NULL REFERENCES publications(id),
              to_id TEXT NOT NULL REFERENCES publications(id),
              relation TEXT NOT NULL,
              confidence REAL NOT NULL DEFAULT 0.5,
              needs_confirm INTEGER NOT NULL DEFAULT 1,
              evidence_json TEXT,
              created_at TEXT NOT NULL,
              updated_at TEXT NOT NULL,
              UNIQUE(from_id, to_id, relation)
            );

            CREATE TABLE IF NOT EXISTS aliases (
              alias TEXT PRIMARY KEY,
              publication_id TEXT NOT NULL REFERENCES publications(id)
            );

            CREATE TABLE IF NOT EXISTS collection_pubs (
              collection_id TEXT NOT NULL REFERENCES collections(id),
              publication_id TEXT NOT NULL REFERENCES publications(id),
              added_at TEXT NOT NULL,
              source_run_id TEXT,
              PRIMARY KEY (collection_id, publication_id)
            );

            CREATE TABLE IF NOT EXISTS runs (
              id TEXT PRIMARY KEY,
              collection_id TEXT NOT NULL REFERENCES collections(id),
              query TEXT,
              sources_json TEXT,
              config_json TEXT,
              started_at TEXT,
              finished_at TEXT,
              hits_path TEXT,
              status TEXT,
              coverage_path TEXT,
              gate_json TEXT
            );

            CREATE TABLE IF NOT EXISTS run_hits (
              id INTEGER PRIMARY KEY AUTOINCREMENT,
              run_id TEXT NOT NULL REFERENCES runs(id),
              publication_id TEXT,
              raw_json TEXT,
              source TEXT,
              fetched_at TEXT,
              dedupe_action TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS theme_nodes (
              id TEXT PRIMARY KEY,
              collection_id TEXT NOT NULL REFERENCES collections(id),
              parent_id TEXT,
              name TEXT NOT NULL,
              sort INTEGER NOT NULL DEFAULT 0
            );

            CREATE TABLE IF NOT EXISTS theme_pubs (
              node_id TEXT NOT NULL REFERENCES theme_nodes(id),
              publication_id TEXT NOT NULL REFERENCES publications(id),
              PRIMARY KEY (node_id, publication_id)
            );

            CREATE TABLE IF NOT EXISTS pub_files (
              id INTEGER PRIMARY KEY AUTOINCREMENT,
              publication_id TEXT NOT NULL REFERENCES publications(id),
              file_path TEXT NOT NULL,
              content_hash TEXT,
              ingest_origin TEXT NOT NULL DEFAULT 'local_file',
              access_origin TEXT,
              access_method TEXT,
              created_at TEXT NOT NULL,
              UNIQUE(publication_id, file_path)
            );

            CREATE INDEX IF NOT EXISTS idx_pubs_year ON publications(year);
            CREATE INDEX IF NOT EXISTS idx_collection_pubs_pub ON collection_pubs(publication_id);
            CREATE INDEX IF NOT EXISTS idx_run_hits_run ON run_hits(run_id);
            CREATE INDEX IF NOT EXISTS idx_pub_versions_confirm ON publication_versions(needs_confirm);
            """
        )
        _ensure_unique_id_indexes(conn)
        _ensure_review_reasons_column(conn)
        _ensure_publication_versions(conn)
        _ensure_runs_gate_columns(conn)
        _ensure_pub_files(conn)
        _ensure_oa_columns(conn)
        conn.execute(
            "INSERT OR REPLACE INTO meta(key, value) VALUES ('schema_version', ?)",
            (str(SCHEMA_VERSION),),
        )
        fts_ok = False
        try:
            conn.execute(
                """
                CREATE VIRTUAL TABLE IF NOT EXISTS pubs_fts USING fts5(
                  publication_id UNINDEXED,
                  title,
                  abstract,
                  tags
                )
                """
            )
            fts_ok = True
        except sqlite3.OperationalError:
            fts_ok = False
        conn.commit()
        return {"ok": True, "path": str(db_path.resolve()), "fts": fts_ok, "schema_version": SCHEMA_VERSION}
    finally:
        conn.close()


def open_db(path: str | Path) -> sqlite3.Connection:
    p = Path(path)
    if not p.exists():
        init_db(p)
        return _connect(p)
    conn = _connect(p)
    _migrate_schema(conn)
    return conn



def _ensure_unique_id_indexes(conn: sqlite3.Connection) -> None:
    """Idempotent UNIQUE indexes on strong ID columns.

    SQLite UNIQUE already allows multiple NULLs. Partial predicate excludes
    `#conflict-` review rows so true ID clashes can still store both values
    on a needs_review shell without breaking uniqueness of canonical rows.
    """
    for col in STRONG_ID_KEYS:
        conn.execute(
            f"CREATE UNIQUE INDEX IF NOT EXISTS uq_pubs_{col} "
            f"ON publications({col}) "
            f"WHERE {col} IS NOT NULL AND instr(id, '#conflict-') = 0"
        )


def _ensure_review_reasons_column(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(publications)").fetchall()}
    if "review_reasons_json" not in cols:
        conn.execute("ALTER TABLE publications ADD COLUMN review_reasons_json TEXT")


def _ensure_publication_versions(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS publication_versions (
          id INTEGER PRIMARY KEY AUTOINCREMENT,
          from_id TEXT NOT NULL REFERENCES publications(id),
          to_id TEXT NOT NULL REFERENCES publications(id),
          relation TEXT NOT NULL,
          confidence REAL NOT NULL DEFAULT 0.5,
          needs_confirm INTEGER NOT NULL DEFAULT 1,
          evidence_json TEXT,
          created_at TEXT NOT NULL,
          updated_at TEXT NOT NULL,
          UNIQUE(from_id, to_id, relation)
        )
        """
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_pub_versions_confirm "
        "ON publication_versions(needs_confirm)"
    )


def _ensure_runs_gate_columns(conn: sqlite3.Connection) -> None:
    """Idempotent runs.status / coverage_path / gate_json for protocol close."""
    cols = {r[1] for r in conn.execute("PRAGMA table_info(runs)").fetchall()}
    if not cols:
        return
    if "status" not in cols:
        conn.execute("ALTER TABLE runs ADD COLUMN status TEXT")
    if "coverage_path" not in cols:
        conn.execute("ALTER TABLE runs ADD COLUMN coverage_path TEXT")
    if "gate_json" not in cols:
        conn.execute("ALTER TABLE runs ADD COLUMN gate_json TEXT")




def _ensure_pub_files(conn: sqlite3.Connection) -> None:
    """Idempotent pub_files table for authorized local file sidecars / PDFs."""
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS pub_files (
          id INTEGER PRIMARY KEY AUTOINCREMENT,
          publication_id TEXT NOT NULL REFERENCES publications(id),
          file_path TEXT NOT NULL,
          content_hash TEXT,
          ingest_origin TEXT NOT NULL DEFAULT 'local_file',
          access_origin TEXT,
          access_method TEXT,
          created_at TEXT NOT NULL,
          UNIQUE(publication_id, file_path)
        )
        """
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_pub_files_pub ON pub_files(publication_id)"
    )
    _ensure_pub_files_access_columns(conn)


def _ensure_pub_files_access_columns(conn: sqlite3.Connection) -> None:
    """Idempotent pub_files.access_origin / access_method (schema v6)."""
    cols = {r[1] for r in conn.execute("PRAGMA table_info(pub_files)").fetchall()}
    if not cols:
        return
    if "access_origin" not in cols:
        conn.execute("ALTER TABLE pub_files ADD COLUMN access_origin TEXT")
    if "access_method" not in cols:
        conn.execute("ALTER TABLE pub_files ADD COLUMN access_method TEXT")


def _ensure_oa_columns(conn: sqlite3.Connection) -> None:
    """Idempotent publications.oa_status / oa_json for candidate URL verification."""
    cols = {r[1] for r in conn.execute("PRAGMA table_info(publications)").fetchall()}
    if "oa_status" not in cols:
        conn.execute("ALTER TABLE publications ADD COLUMN oa_status TEXT")
    if "oa_json" not in cols:
        conn.execute("ALTER TABLE publications ADD COLUMN oa_json TEXT")


def record_pub_file(
    conn: sqlite3.Connection,
    *,
    publication_id: str,
    file_path: str,
    content_hash: Optional[str] = None,
    ingest_origin: str = "local_file",
    access_origin: Optional[str] = None,
    access_method: Optional[str] = None,
) -> dict[str, Any]:
    """Register an authorized local file against a publication (no PDF bytes stored)."""
    now = utc_now()
    path = str(file_path)
    origin = ingest_origin or "local_file"
    conn.execute(
        """
        INSERT INTO pub_files(
          publication_id, file_path, content_hash, ingest_origin,
          access_origin, access_method, created_at
        )
        VALUES (?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(publication_id, file_path) DO UPDATE SET
          content_hash=COALESCE(excluded.content_hash, pub_files.content_hash),
          ingest_origin=excluded.ingest_origin,
          access_origin=COALESCE(excluded.access_origin, pub_files.access_origin),
          access_method=COALESCE(excluded.access_method, pub_files.access_method)
        """,
        (publication_id, path, content_hash, origin, access_origin, access_method, now),
    )
    return {
        "publication_id": publication_id,
        "file_path": path,
        "content_hash": content_hash,
        "ingest_origin": origin,
        "access_origin": access_origin,
        "access_method": access_method,
    }


def set_publication_oa(
    conn: sqlite3.Connection,
    publication_id: str,
    oa_result: dict[str, Any],
) -> None:
    """Store OA URL verification result (status + JSON); never stores PDF bytes."""
    status = oa_result.get("oa_status") or oa_result.get("status")
    payload = dict(oa_result)
    payload.setdefault("checked_at", utc_now())
    conn.execute(
        "UPDATE publications SET oa_status=?, oa_json=?, updated_at=? WHERE id=?",
        (status, json.dumps(payload, ensure_ascii=False), utc_now(), publication_id),
    )


def _repoint_pub_files(conn: sqlite3.Connection, old_id: str, new_id: str) -> None:
    """Move pub_files rows from old_id → new_id; drop duplicates on conflict."""
    if old_id == new_id:
        return
    rows = conn.execute(
        "SELECT id, file_path, content_hash, ingest_origin, created_at FROM pub_files WHERE publication_id=?",
        (old_id,),
    ).fetchall()
    for r in rows:
        exists = conn.execute(
            "SELECT 1 FROM pub_files WHERE publication_id=? AND file_path=?",
            (new_id, r["file_path"]),
        ).fetchone()
        if exists:
            conn.execute("DELETE FROM pub_files WHERE id=?", (r["id"],))
        else:
            conn.execute(
                "UPDATE pub_files SET publication_id=? WHERE id=?",
                (new_id, r["id"]),
            )


def _migrate_schema(conn: sqlite3.Connection) -> None:
    """Bring an existing DB up to SCHEMA_VERSION (idempotent)."""
    row = conn.execute(
        "SELECT value FROM meta WHERE key = 'schema_version'"
    ).fetchone()
    current = int(row["value"]) if row and str(row["value"]).isdigit() else 0
    if current < 2:
        try:
            _ensure_unique_id_indexes(conn)
        except sqlite3.OperationalError:
            # Duplicate strong IDs already in DB — indexes skipped; still bump version
            # so resolve path remains the gate. Repair via merge_ids_cli.
            pass
    if current < 3:
        _ensure_review_reasons_column(conn)
        _ensure_publication_versions(conn)
    if current < 4:
        _ensure_runs_gate_columns(conn)
    if current < 5:
        _ensure_pub_files(conn)
        _ensure_oa_columns(conn)
    if current < 6:
        _ensure_pub_files(conn)  # creates table if missing
        _ensure_pub_files_access_columns(conn)
    if current < SCHEMA_VERSION:
        conn.execute(
            "INSERT OR REPLACE INTO meta(key, value) VALUES ('schema_version', ?)",
            (str(SCHEMA_VERSION),),
        )
        conn.commit()


def _id_prefix_rank(pub_id: str) -> int:
    """Lower = higher priority. Unknown prefixes sort last."""
    order = ("doi:", "pmid:", "pmcid:", "arxiv:", "openalex:", "soft:", "tmp:")
    for i, prefix in enumerate(order):
        if pub_id.startswith(prefix):
            return i
    return len(order)


def _alias_keys(draft: dict[str, Any]) -> dict[str, Optional[str]]:
    return {
        "doi": normalize_doi(draft.get("doi")),
        "pmid": normalize_pmid(draft.get("pmid")),
        "pmcid": normalize_pmcid(draft.get("pmcid")),
        "arxiv": normalize_arxiv(draft.get("arxiv")),
        "openalex": normalize_openalex(draft.get("openalex")),
    }


def _strong_aliases(ids: dict[str, Optional[str]]) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for key in STRONG_ID_KEYS:
        val = ids.get(key)
        if val:
            out.append((f"{key}:{val}", key))
    return out


def _merge_field(existing: Any, incoming: Any) -> Any:
    if existing is None or existing == "" or existing == "null":
        return incoming
    return existing


def _fts_upsert(conn: sqlite3.Connection, pub_id: str, title: str, abstract: str, tags: str) -> None:
    try:
        conn.execute("DELETE FROM pubs_fts WHERE publication_id = ?", (pub_id,))
        conn.execute(
            "INSERT INTO pubs_fts(publication_id, title, abstract, tags) VALUES (?, ?, ?, ?)",
            (pub_id, title or "", abstract or "", tags or ""),
        )
    except sqlite3.OperationalError:
        pass


def resolve_or_create(conn: sqlite3.Connection, publication_draft: dict[str, Any]) -> dict[str, Any]:
    """Insert or merge a publication. Soft-only matches → needs_review, no auto-merge.

    Returns {publication_id, action: inserted|merged|needs_review, matched_via?, soft_candidates?}
    """
    ids = _alias_keys(publication_draft)
    soft = soft_fingerprint(
        publication_draft.get("title"),
        publication_draft.get("authors") or publication_draft.get("authors_json"),
        publication_draft.get("year"),
    )
    strong = _strong_aliases(ids)
    now = utc_now()

    matched_pub_ids: dict[str, str] = {}
    for alias, _kind in strong:
        row = conn.execute(
            "SELECT publication_id FROM aliases WHERE alias = ?", (alias,)
        ).fetchone()
        if row:
            matched_pub_ids[alias] = row["publication_id"]

    distinct = set(matched_pub_ids.values())
    if len(distinct) > 1:
        # Multiple strong aliases hit different pubs. Complementary IDs (no
        # contradictory values per ID kind) → high-confidence merge into DOI-first
        # canonical. True clashes → needs_review conflict row.
        cluster_ids = sorted(distinct)
        rows = []
        for pid in cluster_ids:
            row = conn.execute("SELECT * FROM publications WHERE id = ?", (pid,)).fetchone()
            if row:
                rows.append(dict(row))
        # Synthetic incoming view for complementarity check
        incoming_view = {
            "id": "__incoming__",
            "doi": ids.get("doi"),
            "pmid": ids.get("pmid"),
            "pmcid": ids.get("pmcid"),
            "arxiv": ids.get("arxiv"),
            "openalex": ids.get("openalex"),
            "title": publication_draft.get("title"),
            "abstract": publication_draft.get("abstract"),
            "authors_json": _authors_json(publication_draft),
            "year": _year(publication_draft),
            "venue": publication_draft.get("venue") or publication_draft.get("journal"),
            "type": publication_draft.get("type"),
            "tags_json": _tags_json(publication_draft),
            "needs_review": 0,
        }
        if cluster_high_confidence(rows + [incoming_view]):
            merged = merge_publication_cluster(
                conn,
                cluster_ids,
                preferred_id=choose_canonical_id(rows, ids),
                incoming=incoming_view,
                clear_review=True,
            )
            pub_id = merged["publication_id"]
            # Fold incoming fields
            _merge_pub_richest(conn, pub_id, ids, publication_draft, soft, now)
            for alias, _ in strong:
                conn.execute(
                    "INSERT OR REPLACE INTO aliases(alias, publication_id) VALUES (?, ?)",
                    (alias, pub_id),
                )
            if soft:
                conn.execute(
                    "INSERT OR REPLACE INTO aliases(alias, publication_id) VALUES (?, ?)",
                    (f"soft:{soft}", pub_id),
                )
            conn.execute(
                "UPDATE publications SET needs_review=0, updated_at=? WHERE id=?",
                (now, pub_id),
            )
            return {
                "publication_id": pub_id,
                "action": "merged",
                "matched_via": sorted(matched_pub_ids.keys()),
                "reason": "complementary_strong_ids",
                "merged_ids": merged.get("absorbed_ids") or cluster_ids,
            }

        canon = canonical_id_from_ids(
            doi=ids["doi"], pmid=ids["pmid"], pmcid=ids["pmcid"],
            arxiv=ids["arxiv"], openalex=ids["openalex"], soft=soft,
        )
        # Ensure unique if conflict
        if conn.execute("SELECT 1 FROM publications WHERE id=?", (canon,)).fetchone():
            canon = f"{canon}#conflict-{uuid.uuid4().hex[:8]}"
        _insert_pub(
            conn, canon, ids, publication_draft, soft, now,
            needs_review=1,
            review_reasons=["conflicting_strong_ids"],
        )
        for alias, _ in strong:
            conn.execute(
                "INSERT OR IGNORE INTO aliases(alias, publication_id) VALUES (?, ?)",
                (alias, canon),
            )
        if soft:
            conn.execute(
                "INSERT OR IGNORE INTO aliases(alias, publication_id) VALUES (?, ?)",
                (f"soft:{soft}", canon),
            )
        return {
            "publication_id": canon,
            "action": "needs_review",
            "reason": "conflicting_strong_ids",
            "matched_ids": sorted(distinct),
        }

    if len(distinct) == 1:
        pub_id = next(iter(distinct))
        previous_id: Optional[str] = None
        # Prefer non-shell title / richer abstract when folding bib onto ids-only shells
        _merge_pub_richest(conn, pub_id, ids, publication_draft, soft, now)
        # Canonical upgrade: pmid-only row + incoming DOI → reassign id to doi:…
        merged_row = conn.execute(
            "SELECT * FROM publications WHERE id = ?", (pub_id,)
        ).fetchone()
        if merged_row:
            md = dict(merged_row)
            desired = canonical_id_from_ids(
                doi=md.get("doi"),
                pmid=md.get("pmid"),
                pmcid=md.get("pmcid"),
                arxiv=md.get("arxiv"),
                openalex=md.get("openalex"),
                soft=soft if not any(md.get(k) for k in STRONG_ID_KEYS) else None,
            )
            if (
                desired
                and desired != pub_id
                and not desired.startswith("tmp:")
                and "#conflict-" not in desired
            ):
                # Prefer upgrade when desired is higher priority, or any preferred form
                if _id_prefix_rank(desired) <= _id_prefix_rank(pub_id):
                    target = conn.execute(
                        "SELECT 1 FROM publications WHERE id = ?", (desired,)
                    ).fetchone()
                    if target:
                        clustered = merge_publication_cluster(
                            conn,
                            [pub_id, desired],
                            preferred_id=desired,
                            clear_review=True,
                        )
                        previous_id = pub_id
                        pub_id = clustered["publication_id"]
                    else:
                        _reassign_publication_id(conn, pub_id, desired)
                        previous_id = pub_id
                        pub_id = desired
        for alias, _ in strong:
            conn.execute(
                "INSERT OR REPLACE INTO aliases(alias, publication_id) VALUES (?, ?)",
                (alias, pub_id),
            )
        if soft:
            conn.execute(
                "INSERT OR REPLACE INTO aliases(alias, publication_id) VALUES (?, ?)",
                (f"soft:{soft}", pub_id),
            )
        out: dict[str, Any] = {
            "publication_id": pub_id,
            "action": "merged",
            "matched_via": sorted(matched_pub_ids.keys()),
        }
        if previous_id:
            out["previous_id"] = previous_id
            out["reason"] = "canonical_id_upgrade"
        return out

    # No strong match — check soft fingerprint
    soft_candidates: list[str] = []
    if soft:
        row = conn.execute(
            "SELECT publication_id FROM aliases WHERE alias = ?", (f"soft:{soft}",)
        ).fetchone()
        if row:
            soft_candidates.append(row["publication_id"])

    needs_review = 1 if soft_candidates else 0
    soft_reasons = ["soft_fingerprint_conflict"] if soft_candidates else None
    canon = canonical_id_from_ids(
        doi=ids["doi"], pmid=ids["pmid"], pmcid=ids["pmcid"],
        arxiv=ids["arxiv"], openalex=ids["openalex"], soft=soft,
    )
    if conn.execute("SELECT 1 FROM publications WHERE id=?", (canon,)).fetchone():
        # Extremely rare collision on soft/tmp — merge into existing id
        _merge_pub_richest(conn, canon, ids, publication_draft, soft, now)
        action = "merged"
    else:
        _insert_pub(
            conn, canon, ids, publication_draft, soft, now,
            needs_review=needs_review,
            review_reasons=soft_reasons,
        )
        action = "needs_review" if needs_review else "inserted"

    for alias, _ in strong:
        conn.execute(
            "INSERT OR IGNORE INTO aliases(alias, publication_id) VALUES (?, ?)",
            (alias, canon),
        )
    if soft and action != "needs_review":
        # For soft-candidate conflicts we still record soft alias on the NEW record
        # but leave the old one alone (no auto-merge).
        conn.execute(
            "INSERT OR IGNORE INTO aliases(alias, publication_id) VALUES (?, ?)",
            (f"soft:{soft}", canon),
        )
    elif soft and not soft_candidates:
        conn.execute(
            "INSERT OR IGNORE INTO aliases(alias, publication_id) VALUES (?, ?)",
            (f"soft:{soft}", canon),
        )
    elif soft and soft_candidates:
        # Soft conflict: do not claim the soft alias (already owned); flag review
        pass

    result: dict[str, Any] = {"publication_id": canon, "action": action}
    if soft_candidates:
        result["soft_candidates"] = soft_candidates
        result["reason"] = "soft_fingerprint_conflict"
        for cand in soft_candidates:
            _append_review_reason(conn, cand, "soft_fingerprint_conflict", related=[canon])
        _append_review_reason(conn, canon, "soft_fingerprint_conflict", related=soft_candidates)
    return result


def _authors_json(draft: dict[str, Any]) -> Optional[str]:
    authors = draft.get("authors_json")
    if authors is not None:
        if isinstance(authors, str):
            return authors
        return json.dumps(authors, ensure_ascii=False)
    authors = draft.get("authors")
    if authors is None:
        return None
    if isinstance(authors, str):
        return json.dumps([a.strip() for a in authors.split(" and ")], ensure_ascii=False)
    return json.dumps(authors, ensure_ascii=False)


def _tags_json(draft: dict[str, Any]) -> Optional[str]:
    tags = draft.get("tags_json", draft.get("tags"))
    if tags is None:
        return None
    if isinstance(tags, str):
        try:
            json.loads(tags)
            return tags
        except json.JSONDecodeError:
            return json.dumps([t.strip() for t in tags.split(",") if t.strip()], ensure_ascii=False)
    return json.dumps(tags, ensure_ascii=False)


def _year(draft: dict[str, Any]) -> Optional[int]:
    y = draft.get("year")
    if y is None or y == "":
        return None
    try:
        return int(str(y)[:4])
    except ValueError:
        return None


def _insert_pub(
    conn: sqlite3.Connection,
    pub_id: str,
    ids: dict[str, Optional[str]],
    draft: dict[str, Any],
    soft: Optional[str],
    now: str,
    *,
    needs_review: int,
    review_reasons: Optional[list[str]] = None,
) -> None:
    title = draft.get("title")
    abstract = draft.get("abstract")
    tags = _tags_json(draft)
    reasons = review_reasons or draft.get("review_reasons")
    reasons_json = json.dumps(reasons, ensure_ascii=False) if reasons else None
    conn.execute(
        """
        INSERT INTO publications(
          id, doi, pmid, pmcid, arxiv, openalex, title, abstract,
          authors_json, year, venue, type, tags_json, needs_review,
          review_reasons_json, created_at, updated_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            pub_id,
            ids.get("doi"),
            ids.get("pmid"),
            ids.get("pmcid"),
            ids.get("arxiv"),
            ids.get("openalex"),
            title,
            abstract,
            _authors_json(draft),
            _year(draft),
            draft.get("venue") or draft.get("journal"),
            draft.get("type"),
            tags,
            int(needs_review),
            reasons_json,
            now,
            now,
        ),
    )
    _fts_upsert(conn, pub_id, title or "", abstract or "", tags or "")


def _merge_pub(
    conn: sqlite3.Connection,
    pub_id: str,
    ids: dict[str, Optional[str]],
    draft: dict[str, Any],
    soft: Optional[str],
    now: str,
) -> None:
    row = conn.execute("SELECT * FROM publications WHERE id = ?", (pub_id,)).fetchone()
    if not row:
        _insert_pub(conn, pub_id, ids, draft, soft, now, needs_review=0)
        return
    d = dict(row)
    new_doi = _merge_field(d.get("doi"), ids.get("doi"))
    new_pmid = _merge_field(d.get("pmid"), ids.get("pmid"))
    new_pmcid = _merge_field(d.get("pmcid"), ids.get("pmcid"))
    new_arxiv = _merge_field(d.get("arxiv"), ids.get("arxiv"))
    new_openalex = _merge_field(d.get("openalex"), ids.get("openalex"))
    new_title = _merge_field(d.get("title"), draft.get("title"))
    new_abstract = _merge_field(d.get("abstract"), draft.get("abstract"))
    new_authors = _merge_field(d.get("authors_json"), _authors_json(draft))
    new_year = _merge_field(d.get("year"), _year(draft))
    new_venue = _merge_field(d.get("venue"), draft.get("venue") or draft.get("journal"))
    new_type = _merge_field(d.get("type"), draft.get("type"))
    new_tags = _merge_field(d.get("tags_json"), _tags_json(draft))
    conn.execute(
        """
        UPDATE publications SET
          doi=?, pmid=?, pmcid=?, arxiv=?, openalex=?,
          title=?, abstract=?, authors_json=?, year=?, venue=?, type=?,
          tags_json=?, updated_at=?
        WHERE id=?
        """,
        (
            new_doi, new_pmid, new_pmcid, new_arxiv, new_openalex,
            new_title, new_abstract, new_authors, new_year, new_venue, new_type,
            new_tags, now, pub_id,
        ),
    )
    _fts_upsert(conn, pub_id, new_title or "", new_abstract or "", new_tags or "")


def link_to_collection(
    conn: sqlite3.Connection,
    collection_id: str,
    publication_id: str,
    source_run_id: Optional[str] = None,
) -> None:
    conn.execute(
        """
        INSERT OR IGNORE INTO collection_pubs(collection_id, publication_id, added_at, source_run_id)
        VALUES (?, ?, ?, ?)
        """,
        (collection_id, publication_id, utc_now(), source_run_id),
    )


def record_run_hit(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    publication_id: Optional[str],
    raw: Any,
    source: str,
    dedupe_action: str,
) -> None:
    conn.execute(
        """
        INSERT INTO run_hits(run_id, publication_id, raw_json, source, fetched_at, dedupe_action)
        VALUES (?, ?, ?, ?, ?, ?)
        """,
        (
            run_id,
            publication_id,
            json.dumps(raw, ensure_ascii=False) if not isinstance(raw, str) else raw,
            source,
            utc_now(),
            dedupe_action,
        ),
    )


def fts_available(conn: sqlite3.Connection) -> bool:
    try:
        conn.execute("SELECT 1 FROM pubs_fts LIMIT 1")
        return True
    except sqlite3.OperationalError:
        return False



def _is_shell_title(row: dict[str, Any]) -> bool:
    """True when title is just an id placeholder (doi/pmid string or canonical id)."""
    title = (row.get("title") or "").strip()
    if not title:
        return True
    low = title.lower()
    pid = str(row.get("id") or "")
    if title == pid or low == pid.lower():
        return True
    for key in STRONG_ID_KEYS:
        val = row.get(key)
        if not val:
            continue
        if low == str(val).lower() or low == f"{key}:{val}".lower():
            return True
    if low.startswith("doi:") or low.startswith("pmid:") or low.startswith("pmcid:"):
        return True
    return False


def _field_richness(row: dict[str, Any]) -> tuple:
    abs_len = len(row.get("abstract") or "")
    title_len = 0 if _is_shell_title(row) else len(row.get("title") or "")
    n_ids = sum(1 for k in STRONG_ID_KEYS if row.get(k))
    authors = row.get("authors_json") or ""
    return (title_len > 0, abs_len, n_ids, len(authors), title_len)


def ids_are_complementary(rows: list[dict[str, Any]]) -> bool:
    """True when no ID kind has two different non-null values across rows."""
    if len(rows) < 2:
        return True
    for key in STRONG_ID_KEYS:
        vals = set()
        for r in rows:
            v = r.get(key)
            if v is None or v == "":
                continue
            # normalize lightly
            if key == "doi":
                v = normalize_doi(str(v)) or v
            elif key == "pmid":
                v = normalize_pmid(str(v)) or v
            elif key == "pmcid":
                v = normalize_pmcid(str(v)) or v
            elif key == "arxiv":
                v = normalize_arxiv(str(v)) or v
            elif key == "openalex":
                v = normalize_openalex(str(v)) or v
            vals.add(str(v))
        if len(vals) > 1:
            return False
    return True



def _titles_compatible(rows: list[dict[str, Any]]) -> bool:
    """Non-shell titles agree (equal normalized, or one empty/shell)."""
    real = []
    for r in rows:
        if _is_shell_title(r):
            continue
        nt = normalize_title(r.get("title"))
        if nt:
            real.append(nt)
    if len(real) <= 1:
        return True
    # Allow near-equal / substring (reprint meta suffixes)
    base = min(real, key=len)
    for t in real:
        if t == base or base in t or t in base:
            continue
        # token overlap
        ba, bb = set(base.split()), set(t.split())
        if not ba or not bb:
            return False
        overlap = len(ba & bb) / max(len(ba), len(bb))
        if overlap < 0.6:
            return False
    return True


def cluster_high_confidence(rows: list[dict[str, Any]]) -> bool:
    """Complementary strong IDs + (shells / matching titles / conflict marker)."""
    if len(rows) < 2:
        return False
    if not ids_are_complementary(rows):
        return False
    if any("#conflict-" in str(r.get("id") or "") for r in rows):
        return True
    n_shell = sum(1 for r in rows if _is_shell_title(r))
    if n_shell >= 1 and _titles_compatible(rows):
        return True
    if _titles_compatible(rows):
        # All real titles compatible and IDs complementary
        return True
    return False


def choose_canonical_id(rows: list[dict[str, Any]], incoming_ids: Optional[dict[str, Optional[str]]] = None) -> str:
    """DOI-first canonical id for a complementary cluster."""
    ids = {k: None for k in STRONG_ID_KEYS}
    for r in rows:
        for k in STRONG_ID_KEYS:
            if ids[k] is None and r.get(k):
                ids[k] = r.get(k)
    if incoming_ids:
        for k in STRONG_ID_KEYS:
            if ids[k] is None and incoming_ids.get(k):
                ids[k] = incoming_ids.get(k)
    desired = canonical_id_from_ids(
        doi=ids.get("doi"),
        pmid=ids.get("pmid"),
        pmcid=ids.get("pmcid"),
        arxiv=ids.get("arxiv"),
        openalex=ids.get("openalex"),
    )
    # Prefer an existing non-conflict row that already has this id
    for r in rows:
        rid = str(r.get("id") or "")
        if rid == desired:
            return desired
    for r in sorted(rows, key=_field_richness, reverse=True):
        rid = str(r.get("id") or "")
        if "#conflict-" not in rid and not rid.startswith("tmp:"):
            # If desired is doi:X and this row is doi:X already handled; else keep desired
            break
    return desired


def _reassign_publication_id(conn: sqlite3.Connection, old_id: str, new_id: str) -> None:
    """Rename publication primary key and repoint dependents. new_id must not exist."""
    if old_id == new_id:
        return
    if conn.execute("SELECT 1 FROM publications WHERE id=?", (new_id,)).fetchone():
        raise ValueError(f"cannot rename {old_id} → {new_id}: target exists")
    row = conn.execute("SELECT * FROM publications WHERE id=?", (old_id,)).fetchone()
    if not row:
        raise ValueError(f"missing publication {old_id}")
    _record_id_remap(conn, old_id, new_id)
    cols = [c for c in row.keys()]
    data = dict(row)
    data["id"] = new_id
    # Clear UNIQUE strong-ID columns on the old row so the INSERT does not collide
    conn.execute(
        "UPDATE publications SET doi=NULL, pmid=NULL, pmcid=NULL, arxiv=NULL, openalex=NULL WHERE id=?",
        (old_id,),
    )
    placeholders = ",".join("?" * len(cols))
    col_names = ",".join(cols)
    conn.execute(
        f"INSERT INTO publications({col_names}) VALUES ({placeholders})",
        tuple(data[c] for c in cols),
    )
    conn.execute("UPDATE aliases SET publication_id=? WHERE publication_id=?", (new_id, old_id))
    # collection_pubs PK includes publication_id
    for cp in conn.execute(
        "SELECT collection_id, added_at, source_run_id FROM collection_pubs WHERE publication_id=?",
        (old_id,),
    ).fetchall():
        conn.execute(
            "INSERT OR IGNORE INTO collection_pubs(collection_id, publication_id, added_at, source_run_id) VALUES (?,?,?,?)",
            (cp["collection_id"], new_id, cp["added_at"], cp["source_run_id"]),
        )
    conn.execute("DELETE FROM collection_pubs WHERE publication_id=?", (old_id,))
    conn.execute("UPDATE run_hits SET publication_id=? WHERE publication_id=?", (new_id, old_id))
    for tp in conn.execute(
        "SELECT node_id FROM theme_pubs WHERE publication_id=?", (old_id,)
    ).fetchall():
        conn.execute(
            "INSERT OR IGNORE INTO theme_pubs(node_id, publication_id) VALUES (?, ?)",
            (tp["node_id"], new_id),
        )
    conn.execute("DELETE FROM theme_pubs WHERE publication_id=?", (old_id,))
    _repoint_publication_versions(conn, old_id, new_id)
    _repoint_pub_files(conn, old_id, new_id)
    try:
        conn.execute("DELETE FROM pubs_fts WHERE publication_id=?", (old_id,))
    except sqlite3.OperationalError:
        pass
    conn.execute("DELETE FROM publications WHERE id=?", (old_id,))
    _fts_upsert(
        conn,
        new_id,
        data.get("title") or "",
        data.get("abstract") or "",
        data.get("tags_json") or "",
    )


def _pick_richest_value(rows: list[dict[str, Any]], key: str) -> Any:
    best = None
    best_score = None
    for r in sorted(rows, key=_field_richness, reverse=True):
        v = r.get(key)
        if v is None or v == "":
            continue
        if key == "title" and _is_shell_title({**r, "title": v}):
            score = (0, 0)
        elif key == "abstract":
            score = (1, len(str(v)))
        else:
            score = (1, len(str(v)))
        if best_score is None or score > best_score:
            best, best_score = v, score
    return best


def merge_publication_cluster(
    conn: sqlite3.Connection,
    pub_ids: list[str],
    *,
    preferred_id: Optional[str] = None,
    incoming: Optional[dict[str, Any]] = None,
    clear_review: bool = True,
) -> dict[str, Any]:
    """Merge complementary publications into one canonical row.

    Repoints aliases / collection_pubs / run_hits / theme_pubs / FTS.
    Deletes absorbed publication rows.
    """
    ids = sorted({p for p in pub_ids if p})
    if len(ids) < 1:
        return {"ok": False, "error": "empty cluster"}
    rows: list[dict[str, Any]] = []
    for pid in ids:
        row = conn.execute("SELECT * FROM publications WHERE id=?", (pid,)).fetchone()
        if row:
            rows.append(dict(row))
    if not rows:
        return {"ok": False, "error": "no rows in cluster"}
    if incoming:
        # incoming may be synthetic; only use for field picking / canonical choice
        pass

    id_map = {k: None for k in STRONG_ID_KEYS}
    for r in rows + ([incoming] if incoming else []):
        if not r:
            continue
        for k in STRONG_ID_KEYS:
            if id_map[k] is None and r.get(k):
                id_map[k] = r.get(k)

    desired = preferred_id or choose_canonical_id(rows, id_map)
    # Strip conflict suffix if someone passed a conflict id as preferred incorrectly
    if "#conflict-" in desired:
        desired = desired.split("#conflict-", 1)[0]

    # Ensure a keeper row exists under desired id
    if not conn.execute("SELECT 1 FROM publications WHERE id=?", (desired,)).fetchone():
        # Rename the richest non-conflict row → desired
        candidates = [r for r in rows if "#conflict-" not in str(r["id"])]
        if not candidates:
            candidates = list(rows)
        candidates.sort(key=_field_richness, reverse=True)
        _reassign_publication_id(conn, candidates[0]["id"], desired)
        # refresh rows list ids
        rows = []
        for pid in set(ids) | {desired}:
            if pid == candidates[0]["id"]:
                continue
            row = conn.execute("SELECT * FROM publications WHERE id=?", (pid,)).fetchone()
            if row:
                rows.append(dict(row))
        keeper_row = conn.execute("SELECT * FROM publications WHERE id=?", (desired,)).fetchone()
        if keeper_row:
            rows.append(dict(keeper_row))
        ids = [desired] + [r["id"] for r in rows if r["id"] != desired]

    keeper = desired
    absorb = [r["id"] for r in rows if r["id"] != keeper]
    # Also absorb any leftover cluster ids still present
    for pid in ids:
        if pid != keeper and pid not in absorb:
            if conn.execute("SELECT 1 FROM publications WHERE id=?", (pid,)).fetchone():
                absorb.append(pid)

    all_rows: list[dict[str, Any]] = []
    for pid in [keeper] + absorb:
        row = conn.execute("SELECT * FROM publications WHERE id=?", (pid,)).fetchone()
        if row:
            all_rows.append(dict(row))
    if incoming:
        all_rows.append(incoming)

    now = utc_now()
    new_doi = _pick_richest_value(all_rows, "doi") or id_map.get("doi")
    new_pmid = _pick_richest_value(all_rows, "pmid") or id_map.get("pmid")
    new_pmcid = _pick_richest_value(all_rows, "pmcid") or id_map.get("pmcid")
    new_arxiv = _pick_richest_value(all_rows, "arxiv") or id_map.get("arxiv")
    new_openalex = _pick_richest_value(all_rows, "openalex") or id_map.get("openalex")
    new_title = _pick_richest_value(all_rows, "title")
    new_abstract = _pick_richest_value(all_rows, "abstract")
    new_authors = _pick_richest_value(all_rows, "authors_json")
    new_year = _pick_richest_value(all_rows, "year")
    new_venue = _pick_richest_value(all_rows, "venue")
    new_type = _pick_richest_value(all_rows, "type")
    new_tags = _pick_richest_value(all_rows, "tags_json")
    needs = 0 if clear_review else 1

    # Absorb rows may still hold the same strong IDs — clear before keeper UPDATE
    for aid in absorb:
        conn.execute(
            "UPDATE publications SET doi=NULL, pmid=NULL, pmcid=NULL, arxiv=NULL, openalex=NULL WHERE id=?",
            (aid,),
        )

    conn.execute(
        """
        UPDATE publications SET
          doi=?, pmid=?, pmcid=?, arxiv=?, openalex=?,
          title=?, abstract=?, authors_json=?, year=?, venue=?, type=?,
          tags_json=?, needs_review=?, updated_at=?
        WHERE id=?
        """,
        (
            normalize_doi(new_doi) if new_doi else new_doi,
            normalize_pmid(str(new_pmid)) if new_pmid else new_pmid,
            normalize_pmcid(str(new_pmcid)) if new_pmcid else new_pmcid,
            normalize_arxiv(str(new_arxiv)) if new_arxiv else new_arxiv,
            normalize_openalex(str(new_openalex)) if new_openalex else new_openalex,
            new_title, new_abstract, new_authors, new_year, new_venue, new_type,
            new_tags, needs, now, keeper,
        ),
    )
    _fts_upsert(conn, keeper, new_title or "", new_abstract or "", new_tags or "")

    for aid in absorb:
        _record_id_remap(conn, aid, keeper)
        # repoint dependents
        conn.execute("UPDATE aliases SET publication_id=? WHERE publication_id=?", (keeper, aid))
        for cp in conn.execute(
            "SELECT collection_id, added_at, source_run_id FROM collection_pubs WHERE publication_id=?",
            (aid,),
        ).fetchall():
            conn.execute(
                "INSERT OR IGNORE INTO collection_pubs(collection_id, publication_id, added_at, source_run_id) VALUES (?,?,?,?)",
                (cp["collection_id"], keeper, cp["added_at"], cp["source_run_id"]),
            )
        conn.execute("DELETE FROM collection_pubs WHERE publication_id=?", (aid,))
        conn.execute("UPDATE run_hits SET publication_id=? WHERE publication_id=?", (keeper, aid))
        for tp in conn.execute(
            "SELECT node_id FROM theme_pubs WHERE publication_id=?", (aid,)
        ).fetchall():
            conn.execute(
                "INSERT OR IGNORE INTO theme_pubs(node_id, publication_id) VALUES (?, ?)",
                (tp["node_id"], keeper),
            )
        conn.execute("DELETE FROM theme_pubs WHERE publication_id=?", (aid,))
        _repoint_publication_versions(conn, aid, keeper)
        _repoint_pub_files(conn, aid, keeper)
        try:
            conn.execute("DELETE FROM pubs_fts WHERE publication_id=?", (aid,))
        except sqlite3.OperationalError:
            pass
        conn.execute("DELETE FROM publications WHERE id=?", (aid,))

    # Ensure strong aliases point at keeper
    for key in STRONG_ID_KEYS:
        val = {"doi": new_doi, "pmid": new_pmid, "pmcid": new_pmcid, "arxiv": new_arxiv, "openalex": new_openalex}.get(key)
        if key == "doi":
            val = normalize_doi(val) if val else None
        elif key == "pmid":
            val = normalize_pmid(str(val)) if val else None
        elif key == "pmcid":
            val = normalize_pmcid(str(val)) if val else None
        elif key == "arxiv":
            val = normalize_arxiv(str(val)) if val else None
        elif key == "openalex":
            val = normalize_openalex(str(val)) if val else None
        if val:
            conn.execute(
                "INSERT OR REPLACE INTO aliases(alias, publication_id) VALUES (?, ?)",
                (f"{key}:{val}", keeper),
            )

    return {
        "ok": True,
        "publication_id": keeper,
        "absorbed_ids": absorb,
        "doi": new_doi,
        "pmid": new_pmid,
        "needs_review": needs,
        "title": new_title,
    }


def _merge_pub_richest(
    conn: sqlite3.Connection,
    pub_id: str,
    ids: dict[str, Optional[str]],
    draft: dict[str, Any],
    soft: Optional[str],
    now: str,
) -> None:
    """Like _merge_pub but prefers longer abstract / non-shell title from incoming."""
    row = conn.execute("SELECT * FROM publications WHERE id = ?", (pub_id,)).fetchone()
    if not row:
        _insert_pub(conn, pub_id, ids, draft, soft, now, needs_review=0)
        return
    d = dict(row)
    incoming = {
        "id": pub_id,
        "doi": ids.get("doi"),
        "pmid": ids.get("pmid"),
        "pmcid": ids.get("pmcid"),
        "arxiv": ids.get("arxiv"),
        "openalex": ids.get("openalex"),
        "title": draft.get("title"),
        "abstract": draft.get("abstract"),
        "authors_json": _authors_json(draft),
        "year": _year(draft),
        "venue": draft.get("venue") or draft.get("journal"),
        "type": draft.get("type"),
        "tags_json": _tags_json(draft),
    }
    rows = [d, incoming]
    new_doi = _pick_richest_value(rows, "doi") or _merge_field(d.get("doi"), ids.get("doi"))
    new_pmid = _pick_richest_value(rows, "pmid") or _merge_field(d.get("pmid"), ids.get("pmid"))
    new_pmcid = _pick_richest_value(rows, "pmcid") or _merge_field(d.get("pmcid"), ids.get("pmcid"))
    new_arxiv = _pick_richest_value(rows, "arxiv") or _merge_field(d.get("arxiv"), ids.get("arxiv"))
    new_openalex = _pick_richest_value(rows, "openalex") or _merge_field(d.get("openalex"), ids.get("openalex"))
    new_title = _pick_richest_value(rows, "title") or d.get("title")
    new_abstract = _pick_richest_value(rows, "abstract") or d.get("abstract")
    new_authors = _pick_richest_value(rows, "authors_json") or d.get("authors_json")
    new_year = _pick_richest_value(rows, "year") or d.get("year")
    new_venue = _pick_richest_value(rows, "venue") or d.get("venue")
    new_type = _pick_richest_value(rows, "type") or d.get("type")
    new_tags = _pick_richest_value(rows, "tags_json") or d.get("tags_json")
    conn.execute(
        """
        UPDATE publications SET
          doi=?, pmid=?, pmcid=?, arxiv=?, openalex=?,
          title=?, abstract=?, authors_json=?, year=?, venue=?, type=?,
          tags_json=?, updated_at=?
        WHERE id=?
        """,
        (
            normalize_doi(new_doi) if new_doi else new_doi,
            normalize_pmid(str(new_pmid)) if new_pmid else new_pmid,
            normalize_pmcid(str(new_pmcid)) if new_pmcid else new_pmcid,
            normalize_arxiv(str(new_arxiv)) if new_arxiv else new_arxiv,
            normalize_openalex(str(new_openalex)) if new_openalex else new_openalex,
            new_title, new_abstract, new_authors, new_year, new_venue, new_type,
            new_tags, now, pub_id,
        ),
    )
    _fts_upsert(conn, pub_id, new_title or "", new_abstract or "", new_tags or "")


def find_merge_clusters(conn: sqlite3.Connection) -> list[dict[str, Any]]:
    """Find high-confidence same-paper clusters via shared strong ID values.

    Publications that share a DOI or PMID (etc.) value form a graph component.
    Complementary components (no contradictory ID values) are merge candidates.
    Soft-only conflicts are not included.
    """
    rows = [dict(r) for r in conn.execute("SELECT * FROM publications").fetchall()]
    by_id = {r["id"]: r for r in rows}
    parent: dict[str, str] = {r["id"]: r["id"] for r in rows}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    # Index by each strong id value
    for key in STRONG_ID_KEYS:
        buckets: dict[str, list[str]] = {}
        for r in rows:
            v = r.get(key)
            if not v:
                continue
            buckets.setdefault(str(v), []).append(r["id"])
        for ids in buckets.values():
            for other in ids[1:]:
                union(ids[0], other)

    # Also: conflict rows connect via their embedded doi/pmid even if aliases diverge
    for r in rows:
        rid = r["id"]
        if "#conflict-" in rid:
            base = rid.split("#conflict-", 1)[0]
            if base in by_id:
                union(rid, base)

    components: dict[str, list[str]] = {}
    for r in rows:
        root = find(r["id"])
        components.setdefault(root, []).append(r["id"])

    proposals: list[dict[str, Any]] = []
    for members in components.values():
        if len(members) < 2:
            continue
        cluster_rows = [by_id[m] for m in members if m in by_id]
        complementary = ids_are_complementary(cluster_rows)
        high = cluster_high_confidence(cluster_rows)
        preferred = choose_canonical_id(cluster_rows)
        proposals.append(
            {
                "publication_ids": sorted(members),
                "complementary": complementary,
                "preferred_id": preferred,
                "high_confidence": high,
                "n": len(members),
                "titles": [
                    {
                        "id": r["id"],
                        "title": (r.get("title") or "")[:80],
                        "doi": r.get("doi"),
                        "pmid": r.get("pmid"),
                        "needs_review": r.get("needs_review"),
                        "abs_len": len(r.get("abstract") or ""),
                    }
                    for r in cluster_rows
                ],
            }
        )
    proposals.sort(key=lambda p: (-p["n"], p["preferred_id"]))
    return proposals


def plan_id_merges(
    conn: sqlite3.Connection,
    *,
    only_high_confidence: bool = True,
) -> dict[str, Any]:
    """Dry-run plan for repairing complementary ID clusters."""
    clusters = find_merge_clusters(conn)
    actionable = [c for c in clusters if c.get("high_confidence")] if only_high_confidence else clusters
    skipped = [c for c in clusters if not c.get("high_confidence")]
    return {
        "ok": True,
        "n_clusters": len(clusters),
        "n_actionable": len(actionable),
        "n_skipped_conflicts": len(skipped),
        "merges": [
            {
                "preferred_id": c["preferred_id"],
                "absorb": [i for i in c["publication_ids"] if i != c["preferred_id"]],
                "publication_ids": c["publication_ids"],
                "titles": c["titles"],
            }
            for c in actionable
        ],
        "skipped": skipped,
    }


def apply_id_merges(
    conn: sqlite3.Connection,
    *,
    only_high_confidence: bool = True,
) -> dict[str, Any]:
    """Apply high-confidence complementary merges. Caller commits."""
    plan = plan_id_merges(conn, only_high_confidence=only_high_confidence)
    applied: list[dict[str, Any]] = []
    for m in plan["merges"]:
        # Re-check rows still exist (prior merge may have absorbed)
        existing = [
            pid
            for pid in m["publication_ids"]
            if conn.execute("SELECT 1 FROM publications WHERE id=?", (pid,)).fetchone()
        ]
        if len(existing) < 2:
            continue
        result = merge_publication_cluster(
            conn,
            existing,
            preferred_id=m["preferred_id"],
            clear_review=True,
        )
        applied.append(result)
    return {
        "ok": True,
        "n_planned": len(plan["merges"]),
        "n_applied": len(applied),
        "applied": applied,
        "skipped_conflicts": plan.get("n_skipped_conflicts", 0),
    }


def query_publications(
    conn: sqlite3.Connection,
    *,
    collection_id: Optional[str] = None,
    q: Optional[str] = None,
    tag: Optional[str] = None,
    needs_review: Optional[int] = None,
    limit: int = 50,
) -> list[dict[str, Any]]:
    limit = max(1, min(int(limit), 500))
    params: list[Any] = []
    where: list[str] = []
    join = ""
    if collection_id:
        join = "JOIN collection_pubs cp ON cp.publication_id = p.id"
        where.append("cp.collection_id = ?")
        params.append(collection_id)
    if needs_review is not None:
        where.append("p.needs_review = ?")
        params.append(int(needs_review))
    if tag:
        where.append("CAST(p.tags_json AS TEXT) LIKE ?")
        params.append(f"%{tag}%")

    use_fts = bool(q) and fts_available(conn)
    if q and use_fts:
        # Match FTS then filter
        fts_rows = conn.execute(
            "SELECT publication_id FROM pubs_fts WHERE pubs_fts MATCH ? LIMIT ?",
            (q, limit * 5),
        ).fetchall()
        ids = [r["publication_id"] for r in fts_rows]
        if not ids:
            return []
        placeholders = ",".join("?" * len(ids))
        where.append(f"p.id IN ({placeholders})")
        params.extend(ids)
    elif q:
        where.append("(p.title LIKE ? OR IFNULL(p.abstract,'') LIKE ? OR IFNULL(p.tags_json,'') LIKE ?)")
        like = f"%{q}%"
        params.extend([like, like, like])

    sql = f"SELECT p.* FROM publications p {join}"
    if where:
        sql += " WHERE " + " AND ".join(where)
    sql += " ORDER BY p.year DESC NULLS LAST, p.updated_at DESC LIMIT ?"
    # SQLite may not support NULLS LAST on older builds — fallback
    try:
        params_final = list(params) + [limit]
        rows = conn.execute(sql, params_final).fetchall()
    except sqlite3.OperationalError:
        sql2 = f"SELECT p.* FROM publications p {join}"
        if where:
            sql2 += " WHERE " + " AND ".join(where)
        sql2 += " ORDER BY p.year DESC, p.updated_at DESC LIMIT ?"
        rows = conn.execute(sql2, list(params) + [limit]).fetchall()
    return [dict(r) for r in rows]



# ---------------------------------------------------------------------------
# 0.4.7: review reasons, human_queue, coverage, publication_versions
# ---------------------------------------------------------------------------

REVIEW_SOFT = "soft_fingerprint_conflict"
REVIEW_CONFLICT = "conflicting_strong_ids"
REVIEW_SHELL = "empty_unenriched_shell"
REVIEW_VERSION = "pending_version_link"
REVIEW_PREVIEW = "preview_fetched_opposite"

VERSION_RELATIONS = ("preprint_of", "published_as", "same_work")

PROTOCOL_VERSION = "0.5.2"
GOLD_SET_MIN_GUIDANCE = 5  # guidance only; <5 → evidence_weak, not hard-fail

PROVIDER_PROTOCOL_NOTES = {
    "pubmed": "NCBI E-utilities esearch+esummary/efetch",
    "europepmc": "Europe PMC search API",
    "openalex": "OpenAlex /works",
    "crossref": "Crossref /works",
    "semantic_scholar": "S2 Graph API (soft-fail without S2_API_KEY)",
    "arxiv": "arXiv Atom API",
    "biorxiv": "bioRxiv details/date browse + Europe PMC fallback",
    "medrxiv": "medRxiv details/date browse + Europe PMC fallback",
    "pmc": "NCBI E-utilities db=pmc",
    "doaj": "DOAJ /api/search/articles",
    "datacite": "DataCite REST /dois",
    "zenodo": "Zenodo /api/records",
}

RUN_STATUS_CONFIGURED = "configured"
RUN_STATUS_SEARCHED = "searched"
RUN_STATUS_PASSED = "passed"
RUN_STATUS_NEEDS_ATTENTION = "needs_attention"

SUGGESTED_ACTIONS = {
    REVIEW_CONFLICT: "Inspect aliases / #conflict- row; merge complementary IDs via merge_ids_cli --apply or keep distinct",
    REVIEW_SOFT: "Compare related_ids soft-fingerprint peers; adjudicate same work vs distinct; never auto-merge soft",
    REVIEW_VERSION: "Confirm with lit-versions --action confirm --link_id <id> (merges to journal DOI by default)",
    REVIEW_PREVIEW: "Review preview-fetched metadata shell; confirm version link to merge into journal DOI preferred",
    REVIEW_SHELL: "Run lit-enrich on publication_id (metadata only; --mock offline)",
    "needs_review": "Open item in lit-review-queue --dig and adjudicate identity",
}


def _parse_reasons(raw: Any) -> list[str]:
    if not raw:
        return []
    if isinstance(raw, list):
        return [str(x) for x in raw if x]
    if isinstance(raw, str):
        try:
            val = json.loads(raw)
            if isinstance(val, list):
                return [str(x) for x in val if x]
        except json.JSONDecodeError:
            return [raw] if raw.strip() else []
    return []


def _append_review_reason(
    conn: sqlite3.Connection,
    pub_id: str,
    reason: str,
    *,
    related: Optional[list[str]] = None,
    set_flag: bool = True,
) -> None:
    row = conn.execute(
        "SELECT needs_review, review_reasons_json FROM publications WHERE id=?",
        (pub_id,),
    ).fetchone()
    if not row:
        return
    reasons = _parse_reasons(row["review_reasons_json"] if "review_reasons_json" in row.keys() else None)
    if reason not in reasons:
        reasons.append(reason)
    payload = reasons
    # Optionally stash related ids inside a tagged entry is overkill; keep flat reasons.
    conn.execute(
        "UPDATE publications SET needs_review=?, review_reasons_json=?, updated_at=? WHERE id=?",
        (
            1 if set_flag else int(row["needs_review"] or 0),
            json.dumps(payload, ensure_ascii=False),
            utc_now(),
            pub_id,
        ),
    )


def _repoint_publication_versions(conn: sqlite3.Connection, old_id: str, new_id: str) -> None:
    """Repoint or collapse version links when a publication id is renamed/absorbed."""
    if old_id == new_id:
        return
    try:
        rows = conn.execute(
            "SELECT id, from_id, to_id, relation FROM publication_versions "
            "WHERE from_id=? OR to_id=?",
            (old_id, old_id),
        ).fetchall()
    except sqlite3.OperationalError:
        return
    for r in rows:
        new_from = new_id if r["from_id"] == old_id else r["from_id"]
        new_to = new_id if r["to_id"] == old_id else r["to_id"]
        if new_from == new_to:
            conn.execute("DELETE FROM publication_versions WHERE id=?", (r["id"],))
            continue
        # UNIQUE(from_id,to_id,relation) — delete self if collision, else update
        clash = conn.execute(
            "SELECT id FROM publication_versions WHERE from_id=? AND to_id=? AND relation=? AND id!=?",
            (new_from, new_to, r["relation"], r["id"]),
        ).fetchone()
        if clash:
            conn.execute("DELETE FROM publication_versions WHERE id=?", (r["id"],))
        else:
            conn.execute(
                "UPDATE publication_versions SET from_id=?, to_id=?, updated_at=? WHERE id=?",
                (new_from, new_to, utc_now(), r["id"]),
            )


def is_preprint_doi(doi: Optional[str]) -> bool:
    d = normalize_doi(doi) or ""
    return d.startswith("10.1101/")


def is_preprint_pub(row: dict[str, Any]) -> bool:
    if row.get("arxiv"):
        return True
    if is_preprint_doi(row.get("doi")):
        return True
    t = str(row.get("type") or "").lower()
    if t in ("preprint", "posted-content", "posted_content"):
        return True
    return False


def is_journal_pub(row: dict[str, Any]) -> bool:
    """Heuristic: non-preprint DOI and/or PMID without being a preprint row."""
    if is_preprint_pub(row):
        # Still journal-ish if it also carries a non-preprint DOI (after merge)
        doi = normalize_doi(row.get("doi"))
        if doi and not is_preprint_doi(doi):
            return True
        return False
    doi = normalize_doi(row.get("doi"))
    if doi and not is_preprint_doi(doi):
        return True
    if row.get("pmid"):
        return True
    return False


def _collection_pub_ids(conn: sqlite3.Connection, collection_id: Optional[str]) -> Optional[set[str]]:
    if not collection_id:
        return None
    rows = conn.execute(
        "SELECT publication_id FROM collection_pubs WHERE collection_id=?",
        (collection_id,),
    ).fetchall()
    return {r["publication_id"] for r in rows}


def _is_unenriched_shell(row: dict[str, Any]) -> bool:
    """ID shell waiting for enrich: shell-like title and empty abstract."""
    abstract = (row.get("abstract") or "").strip()
    if abstract:
        return False
    return _is_shell_title(row)


def derive_review_reasons(conn: sqlite3.Connection, row: dict[str, Any]) -> list[str]:
    """Combine stored reasons with structural signals."""
    reasons = _parse_reasons(row.get("review_reasons_json"))
    pid = str(row.get("id") or "")
    if "#conflict-" in pid and REVIEW_CONFLICT not in reasons:
        reasons.append(REVIEW_CONFLICT)
    if int(row.get("needs_review") or 0) == 1 and not reasons:
        # Soft conflict often leaves reason on row after 0.4.7; older rows: sniff soft
        soft = soft_fingerprint(row.get("title"), row.get("authors_json"), row.get("year"))
        if soft:
            owner = conn.execute(
                "SELECT publication_id FROM aliases WHERE alias=?",
                (f"soft:{soft}",),
            ).fetchone()
            if owner and owner["publication_id"] != pid:
                reasons.append(REVIEW_SOFT)
            else:
                # Same soft owned by this row — look for other pubs with same fingerprint
                for other in conn.execute("SELECT id, title, authors_json, year FROM publications").fetchall():
                    if other["id"] == pid:
                        continue
                    if soft_fingerprint(other["title"], other["authors_json"], other["year"]) == soft:
                        reasons.append(REVIEW_SOFT)
                        break
        if not reasons:
            reasons.append("needs_review")
    if _is_unenriched_shell(row) and REVIEW_SHELL not in reasons:
        # Only flag shells that have some identity handle
        if any(row.get(k) for k in STRONG_ID_KEYS) or pid.startswith(("soft:", "tmp:", "doi:", "pmid:", "pmcid:", "arxiv:", "openalex:")):
            reasons.append(REVIEW_SHELL)
    return reasons


def list_human_queue(
    conn: sqlite3.Connection,
    *,
    collection_id: Optional[str] = None,
    limit: int = 200,
    include_shells: bool = True,
    include_versions: bool = True,
) -> dict[str, Any]:
    """Publications needing human *identity* adjudication (not full-text review)."""
    limit = max(1, min(int(limit), 2000))
    allowed = _collection_pub_ids(conn, collection_id)
    items: dict[str, dict[str, Any]] = {}

    def _touch(pid: str, reasons: list[str], **extra: Any) -> None:
        if allowed is not None and pid not in allowed:
            return
        row = conn.execute("SELECT * FROM publications WHERE id=?", (pid,)).fetchone()
        if not row:
            return
        d = dict(row)
        cur = items.get(pid)
        if not cur:
            cur = {
                "publication_id": pid,
                "reasons": [],
                "related_ids": [],
                "title": d.get("title"),
                "doi": d.get("doi"),
                "pmid": d.get("pmid"),
                "pmcid": d.get("pmcid"),
                "arxiv": d.get("arxiv"),
                "year": d.get("year"),
                "needs_review": d.get("needs_review"),
                "is_shell": _is_unenriched_shell(d),
                "is_preprint": is_preprint_pub(d),
            }
            items[pid] = cur
        for r in reasons:
            if r not in cur["reasons"]:
                cur["reasons"].append(r)
        for rid in extra.get("related_ids") or []:
            if rid and rid not in cur["related_ids"] and rid != pid:
                cur["related_ids"].append(rid)
        if extra.get("version_link_id") is not None:
            links = cur.setdefault("version_link_ids", [])
            if extra["version_link_id"] not in links:
                links.append(extra["version_link_id"])

    # 1) needs_review rows + stored/derived reasons
    q = "SELECT * FROM publications WHERE needs_review = 1"
    for row in conn.execute(q).fetchall():
        d = dict(row)
        if allowed is not None and d["id"] not in allowed:
            continue
        reasons = derive_review_reasons(conn, d)
        related: list[str] = []
        if REVIEW_SOFT in reasons:
            soft = soft_fingerprint(d.get("title"), d.get("authors_json"), d.get("year"))
            if soft:
                for other in conn.execute(
                    "SELECT id, title, authors_json, year FROM publications"
                ).fetchall():
                    if other["id"] == d["id"]:
                        continue
                    if soft_fingerprint(other["title"], other["authors_json"], other["year"]) == soft:
                        related.append(other["id"])
        if REVIEW_CONFLICT in reasons or "#conflict-" in d["id"]:
            # matched peers via shared strong aliases pointing elsewhere is hard; expose matched aliases
            for key in STRONG_ID_KEYS:
                val = d.get(key)
                if not val:
                    continue
                peer = conn.execute(
                    "SELECT publication_id FROM aliases WHERE alias=?",
                    (f"{key}:{val}",),
                ).fetchone()
                if peer and peer["publication_id"] != d["id"]:
                    related.append(peer["publication_id"])
        _touch(d["id"], reasons, related_ids=related)

    # 2) empty unenriched shells (even if needs_review=0)
    if include_shells:
        for row in conn.execute("SELECT * FROM publications").fetchall():
            d = dict(row)
            if allowed is not None and d["id"] not in allowed:
                continue
            if _is_unenriched_shell(d):
                _touch(d["id"], [REVIEW_SHELL])

    # 3) pending version links
    version_pending = 0
    if include_versions:
        try:
            vrows = conn.execute(
                "SELECT * FROM publication_versions WHERE needs_confirm = 1"
            ).fetchall()
        except sqlite3.OperationalError:
            vrows = []
        for vr in vrows:
            version_pending += 1
            _touch(
                vr["from_id"],
                [REVIEW_VERSION],
                related_ids=[vr["to_id"]],
                version_link_id=vr["id"],
            )
            _touch(
                vr["to_id"],
                [REVIEW_VERSION],
                related_ids=[vr["from_id"]],
                version_link_id=vr["id"],
            )

    ordered = sorted(
        items.values(),
        key=lambda x: (
            0 if REVIEW_CONFLICT in x["reasons"] else 1,
            0 if REVIEW_SOFT in x["reasons"] else 1,
            0 if REVIEW_VERSION in x["reasons"] else 1,
            str(x.get("publication_id") or ""),
        ),
    )
    truncated = ordered[:limit]
    by_reason: dict[str, int] = {}
    for it in truncated:
        for r in it["reasons"]:
            by_reason[r] = by_reason.get(r, 0) + 1

    # Attach suggested next action (priority: conflict > soft > version > preview > shell)
    priority = (
        REVIEW_CONFLICT,
        REVIEW_SOFT,
        REVIEW_VERSION,
        REVIEW_PREVIEW,
        REVIEW_SHELL,
        "needs_review",
    )
    for it in truncated:
        reasons = it.get("reasons") or []
        action = None
        for code in priority:
            if code in reasons:
                action = SUGGESTED_ACTIONS.get(code)
                break
        if action is None and reasons:
            action = SUGGESTED_ACTIONS.get(reasons[0]) or SUGGESTED_ACTIONS["needs_review"]
        it["suggested_next_action"] = action or SUGGESTED_ACTIONS["needs_review"]

    return {
        "ok": True,
        "n": len(truncated),
        "n_total": len(ordered),
        "by_reason": by_reason,
        "version_links_pending": version_pending,
        "collection_id": collection_id,
        "items": truncated,
        "backlog_empty": len(ordered) == 0,
    }



def coverage_actionable_hints(report: dict[str, Any]) -> list[str]:
    """Honest, conditional next-step hints when recall/coverage looks incomplete."""
    hints: list[str] = []
    revid = report.get("recall_evidence") or {}
    gold = revid.get("gold_set_recall") or {}
    oow = revid.get("out_of_window") or {}
    unc = revid.get("declared_but_uncontacted_sources") or {}
    coll = report.get("collection") or {}
    hits = report.get("hits") or {}

    recall = gold.get("recall")
    gold_n = int(gold.get("gold_set_n") or 0)
    n_missed = int(gold.get("n_missed") or 0)
    if gold_n > 0 and recall is not None and float(recall) < 1.0:
        hints.append(
            f"gold_set recall={recall} ({gold.get('n_found', 0)}/{gold_n}) — "
            "increase lit-search --limit and/or broaden query terms"
        )
        if n_missed:
            hints.append(
                "ingest missed gold IDs via lit-ingest --format ids (seed file), then re-run coverage"
            )
            hints.append(
                "run lit-expand snowball from found seeds (--confirm_write) to pull related work"
            )
        if gold.get("evidence_weak") or report.get("evidence_weak"):
            hints.append(
                "gold_set is empty or below guidance size — evidence_weak; enlarge gold_set in run-config"
            )
    if oow.get("n_controls") and oow.get("filter_ok") is False:
        hints.append(
            "out-of-window controls leaked into hit corpus — tighten year/window filters or query"
        )
    if int(unc.get("n_uncontacted") or 0) > 0:
        hints.append(
            "declared sources were not contacted — check provider errors / keys / --sources list"
        )
    if int(coll.get("n_unenriched_shells") or 0) > 0:
        hints.append(
            f"{coll.get('n_unenriched_shells')} unenriched shells — run lit-enrich before wiki compile"
        )
    if int(coll.get("n_needs_review") or 0) > 0:
        hints.append(
            f"{coll.get('n_needs_review')} needs_review — dig lit-review-queue before treating corpus as closed"
        )
    if int(hits.get("n") or 0) == 0 and gold_n > 0:
        hints.append("run has zero hits while gold_set is set — re-run lit-search or seed ingest")
    return hints


def coverage_report(
    conn: sqlite3.Connection,
    *,
    run_id: Optional[str] = None,
    collection_id: Optional[str] = None,
) -> dict[str, Any]:
    """Search hits vs unique ingested pubs + dedupe / review / shell stats."""
    if not run_id and not collection_id:
        return {"ok": False, "error": "provide run_id and/or collection_id"}

    run_meta = None
    if run_id:
        row = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
        if not row:
            return {"ok": False, "error": f"unknown run_id {run_id}"}
        run_meta = dict(row)
        if not collection_id:
            collection_id = run_meta.get("collection_id")

    # Hits
    hit_where = []
    hit_params: list[Any] = []
    if run_id:
        hit_where.append("run_id = ?")
        hit_params.append(run_id)
    elif collection_id:
        hit_where.append("run_id IN (SELECT id FROM runs WHERE collection_id = ?)")
        hit_params.append(collection_id)
    hit_sql = "SELECT dedupe_action, publication_id, source FROM run_hits"
    if hit_where:
        hit_sql += " WHERE " + " AND ".join(hit_where)
    hits = conn.execute(hit_sql, hit_params).fetchall()
    n_hits = len(hits)
    dedupe_counts: dict[str, int] = {}
    source_counts: dict[str, int] = {}
    hit_pub_ids: set[str] = set()
    for h in hits:
        act = h["dedupe_action"] or "unknown"
        dedupe_counts[act] = dedupe_counts.get(act, 0) + 1
        src = h["source"] or "unknown"
        source_counts[src] = source_counts.get(src, 0) + 1
        if h["publication_id"]:
            hit_pub_ids.add(h["publication_id"])

    # Collection pubs
    if collection_id:
        cp = conn.execute(
            "SELECT publication_id FROM collection_pubs WHERE collection_id=?",
            (collection_id,),
        ).fetchall()
        coll_ids = {r["publication_id"] for r in cp}
    else:
        coll_ids = set(hit_pub_ids)

    pubs = []
    for pid in sorted(coll_ids):
        row = conn.execute("SELECT * FROM publications WHERE id=?", (pid,)).fetchone()
        if row:
            pubs.append(dict(row))

    n_pubs = len(pubs)
    n_needs_review = sum(1 for p in pubs if int(p.get("needs_review") or 0) == 1)
    n_shells = sum(1 for p in pubs if _is_unenriched_shell(p))
    n_missing_title = sum(1 for p in pubs if not (p.get("title") or "").strip() or _is_shell_title(p))
    n_missing_abstract = sum(1 for p in pubs if not (p.get("abstract") or "").strip())
    missing_strong = []
    for p in pubs:
        if not any(p.get(k) for k in STRONG_ID_KEYS):
            missing_strong.append(p["id"])
    n_preprint = sum(1 for p in pubs if is_preprint_pub(p))

    v_pending = 0
    if coll_ids:
        try:
            placeholders = ",".join("?" * len(coll_ids))
            v_pending = conn.execute(
                f"SELECT COUNT(*) c FROM publication_versions WHERE needs_confirm=1 "
                f"AND (from_id IN ({placeholders}) OR to_id IN ({placeholders}))",
                tuple(coll_ids) + tuple(coll_ids),
            ).fetchone()["c"]
        except sqlite3.OperationalError:
            v_pending = 0

    unique_from_hits = len(hit_pub_ids)

    # --- 0.5.2 recall evidence (gold / OOW / uncontacted) ---
    cfg: dict[str, Any] = {}
    if run_meta and run_meta.get("config_json"):
        try:
            cfg = json.loads(run_meta["config_json"]) or {}
        except json.JSONDecodeError:
            cfg = {}
    # Prefer corpus = hit pubs when run-scoped; else collection pubs
    recall_corpus = set(hit_pub_ids) if run_id else set(coll_ids)
    # Also accept collection membership for gold found when run has hits empty
    if run_id:
        recall_corpus |= set(coll_ids) & set(hit_pub_ids)
        # Gold recall: present in ingested collection OR in this run's hits
        gold_corpus = set(coll_ids) | set(hit_pub_ids)
        oow_corpus = set(hit_pub_ids)  # windowed = this run's hits
    else:
        gold_corpus = set(coll_ids)
        oow_corpus = set(coll_ids)

    gold_set = cfg.get("gold_set") or []
    oow = cfg.get("out_of_window_controls") or []
    configured_sources: list[str] = []
    if cfg.get("sources"):
        configured_sources = list(cfg.get("sources") or [])
    elif run_meta and run_meta.get("sources_json"):
        try:
            configured_sources = json.loads(run_meta["sources_json"]) or []
        except json.JSONDecodeError:
            configured_sources = []

    sources_executed = cfg.get("sources_executed")
    providers = cfg.get("providers")
    # Disk fallback: hits.meta.json next to hits_path
    if run_meta and run_meta.get("hits_path"):
        meta_path = Path(run_meta["hits_path"]).parent / "hits.meta.json"
        if meta_path.is_file():
            try:
                disk_meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                disk_meta = {}
            if sources_executed is None:
                sources_executed = disk_meta.get("sources")
            if not providers:
                providers = disk_meta.get("providers")

    gold_recall = evaluate_gold_set_recall(conn, list(gold_set), corpus_pub_ids=gold_corpus)
    oow_eval = evaluate_out_of_window_controls(conn, list(oow), corpus_pub_ids=oow_corpus)
    uncontacted = evaluate_uncontacted_sources(
        configured_sources=configured_sources,
        sources_executed=sources_executed,
        providers=providers,
        hit_source_counts=source_counts,
    )
    recall_evidence = {
        "gold_set_recall": gold_recall,
        "out_of_window": oow_eval,
        "declared_but_uncontacted_sources": uncontacted,
        "evidence_weak": bool(gold_recall.get("evidence_weak")),
    }

    out = {
        "ok": True,
        "run_id": run_id,
        "collection_id": collection_id,
        "run": {
            "query": (run_meta or {}).get("query"),
            "sources_json": (run_meta or {}).get("sources_json"),
            "started_at": (run_meta or {}).get("started_at"),
            "finished_at": (run_meta or {}).get("finished_at"),
            "status": (run_meta or {}).get("status"),
        } if run_meta else None,
        "hits": {
            "n": n_hits,
            "unique_publication_ids": unique_from_hits,
            "dedupe_action_counts": dedupe_counts,
            "source_counts": source_counts,
        },
        "collection": {
            "n_publications": n_pubs,
            "n_needs_review": n_needs_review,
            "n_unenriched_shells": n_shells,
            "n_missing_title_or_shell": n_missing_title,
            "n_missing_abstract": n_missing_abstract,
            "n_missing_strong_ids": len(missing_strong),
            "missing_strong_id_sample": missing_strong[:20],
            "n_preprints": n_preprint,
            "n_version_links_pending": v_pending,
        },
        "ratios": {
            "hit_to_unique_pub_ratio": (unique_from_hits / n_hits) if n_hits else None,
            "collection_vs_hit_pubs": {
                "in_collection_not_in_hits": len(coll_ids - hit_pub_ids) if run_id else None,
                "in_hits_not_in_collection": len(hit_pub_ids - coll_ids) if run_id else None,
            },
        },
        "recall_evidence": recall_evidence,
        "evidence_weak": bool(gold_recall.get("evidence_weak")),
    }
    out["hints"] = coverage_actionable_hints(out)
    return out


def propose_version_link(
    conn: sqlite3.Connection,
    from_id: str,
    to_id: str,
    *,
    relation: str = "preprint_of",
    confidence: float = 0.5,
    needs_confirm: int = 1,
    evidence: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Insert or refresh a publication_versions row."""
    if relation not in VERSION_RELATIONS:
        return {"ok": False, "error": f"relation must be one of {VERSION_RELATIONS}"}
    if from_id == to_id:
        return {"ok": False, "error": "from_id and to_id must differ"}
    for pid in (from_id, to_id):
        if not conn.execute("SELECT 1 FROM publications WHERE id=?", (pid,)).fetchone():
            return {"ok": False, "error": f"missing publication {pid}"}
    now = utc_now()
    evidence_json = json.dumps(evidence or {}, ensure_ascii=False)
    existing = conn.execute(
        "SELECT id FROM publication_versions WHERE from_id=? AND to_id=? AND relation=?",
        (from_id, to_id, relation),
    ).fetchone()
    if existing:
        conn.execute(
            "UPDATE publication_versions SET confidence=?, needs_confirm=?, evidence_json=?, updated_at=? WHERE id=?",
            (float(confidence), int(needs_confirm), evidence_json, now, existing["id"]),
        )
        link_id = existing["id"]
        action = "updated"
    else:
        cur = conn.execute(
            """
            INSERT INTO publication_versions(
              from_id, to_id, relation, confidence, needs_confirm, evidence_json, created_at, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (from_id, to_id, relation, float(confidence), int(needs_confirm), evidence_json, now, now),
        )
        link_id = cur.lastrowid
        action = "inserted"
    if needs_confirm:
        _append_review_reason(conn, from_id, REVIEW_VERSION, related=[to_id])
        _append_review_reason(conn, to_id, REVIEW_VERSION, related=[from_id])
    return {
        "ok": True,
        "action": action,
        "link_id": link_id,
        "from_id": from_id,
        "to_id": to_id,
        "relation": relation,
        "confidence": float(confidence),
        "needs_confirm": int(needs_confirm),
    }


def detect_version_candidates(
    conn: sqlite3.Connection,
    *,
    collection_id: Optional[str] = None,
    min_title_overlap: float = 0.6,
) -> list[dict[str, Any]]:
    """Offline heuristic: preprint vs journal rows with matching soft fingerprint / titles.

    Does not call the network. High confidence when soft fingerprints match exactly
    and one side is preprint (arxiv / 10.1101) while the other has a journal DOI.
    """
    allowed = _collection_pub_ids(conn, collection_id)
    rows = [dict(r) for r in conn.execute("SELECT * FROM publications").fetchall()]
    if allowed is not None:
        rows = [r for r in rows if r["id"] in allowed]

    preprints = [r for r in rows if is_preprint_pub(r)]
    journals = [r for r in rows if is_journal_pub(r) and not is_preprint_pub(r)]
    # Also allow journal rows that are simply non-preprint with DOI/PMID
    if not journals:
        journals = [
            r for r in rows
            if (normalize_doi(r.get("doi")) and not is_preprint_doi(r.get("doi"))) or r.get("pmid")
        ]
        journals = [r for r in journals if r["id"] not in {p["id"] for p in preprints}]

    proposals: list[dict[str, Any]] = []
    seen_pairs: set[tuple[str, str]] = set()

    def _title_overlap(a: dict[str, Any], b: dict[str, Any]) -> float:
        if _is_shell_title(a) or _is_shell_title(b):
            return 0.0
        ta, tb = normalize_title(a.get("title")), normalize_title(b.get("title"))
        if not ta or not tb:
            return 0.0
        if ta == tb or ta in tb or tb in ta:
            return 1.0
        aa, bb = set(ta.split()), set(tb.split())
        if not aa or not bb:
            return 0.0
        return len(aa & bb) / max(len(aa), len(bb))

    for pre in preprints:
        soft_p = soft_fingerprint(pre.get("title"), pre.get("authors_json"), pre.get("year"))
        for j in journals:
            if pre["id"] == j["id"]:
                continue
            # Skip if already same row via alias identity
            pair = (pre["id"], j["id"])
            if pair in seen_pairs:
                continue
            # Different DOI values expected
            pd, jd = normalize_doi(pre.get("doi")), normalize_doi(j.get("doi"))
            if pd and jd and pd == jd:
                continue

            soft_j = soft_fingerprint(j.get("title"), j.get("authors_json"), j.get("year"))
            overlap = _title_overlap(pre, j)
            same_soft = bool(soft_p and soft_j and soft_p == soft_j)
            same_year = (
                pre.get("year") is not None
                and j.get("year") is not None
                and int(pre["year"]) == int(j["year"])
            )
            shared_pmid = (
                pre.get("pmid")
                and j.get("pmid")
                and str(pre["pmid"]) == str(j["pmid"])
            )

            if not (same_soft or (overlap >= min_title_overlap and same_year) or shared_pmid):
                continue

            confidence = 0.5
            evidence: dict[str, Any] = {
                "title_overlap": round(overlap, 3),
                "same_soft_fingerprint": same_soft,
                "same_year": same_year,
                "shared_pmid": bool(shared_pmid),
                "preprint_doi": pd,
                "journal_doi": jd,
                "preprint_arxiv": pre.get("arxiv"),
            }
            if shared_pmid:
                confidence = 0.95
            elif same_soft and (pd or pre.get("arxiv")) and jd and not is_preprint_doi(jd):
                confidence = 0.9
            elif same_soft:
                confidence = 0.8
            elif overlap >= 0.85 and same_year:
                confidence = 0.7

            # Prefer orientation: preprint --preprint_of--> journal
            seen_pairs.add(pair)
            proposals.append(
                {
                    "from_id": pre["id"],
                    "to_id": j["id"],
                    "relation": "preprint_of",
                    "confidence": confidence,
                    "high_confidence": confidence >= 0.85,
                    "evidence": evidence,
                    "preprint_title": (pre.get("title") or "")[:100],
                    "journal_title": (j.get("title") or "")[:100],
                }
            )

    proposals.sort(key=lambda p: (-p["confidence"], p["from_id"], p["to_id"]))
    return proposals




def merge_version_proposals(
    *lists: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Dedupe version proposals by (from_id, to_id, relation); keep highest confidence."""
    best: dict[tuple[str, str, str], dict[str, Any]] = {}
    for lst in lists:
        for p in lst or []:
            key = (p.get("from_id") or "", p.get("to_id") or "", p.get("relation") or "preprint_of")
            if not key[0] or not key[1]:
                continue
            prev = best.get(key)
            if prev is None or float(p.get("confidence") or 0) > float(prev.get("confidence") or 0):
                best[key] = p
    out = list(best.values())
    out.sort(key=lambda p: (-float(p.get("confidence") or 0), p.get("from_id") or "", p.get("to_id") or ""))
    return out


def proposals_from_doi_ties(
    conn: sqlite3.Connection,
    ties: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Map external DOI ties onto in-DB publication ids (via columns or aliases)."""

    def _pub_for_doi(doi: str) -> Optional[str]:
        d = normalize_doi(doi)
        if not d:
            return None
        row = conn.execute("SELECT id FROM publications WHERE doi=?", (d,)).fetchone()
        if row:
            return row["id"]
        alias = conn.execute(
            "SELECT publication_id FROM aliases WHERE alias=?", (f"doi:{d}",)
        ).fetchone()
        return alias["publication_id"] if alias else None

    proposals: list[dict[str, Any]] = []
    for t in ties:
        from_doi = normalize_doi(t.get("from_doi"))
        to_doi = normalize_doi(t.get("to_doi"))
        if not from_doi or not to_doi or from_doi == to_doi:
            continue
        from_id = _pub_for_doi(from_doi)
        to_id = _pub_for_doi(to_doi)
        if not from_id or not to_id or from_id == to_id:
            continue
        evidence = dict(t.get("evidence") or {})
        evidence["from_doi"] = from_doi
        evidence["to_doi"] = to_doi
        evidence["http_provider"] = t.get("provider")
        conf = float(t.get("confidence") or 0.5)
        proposals.append(
            {
                "from_id": from_id,
                "to_id": to_id,
                "relation": t.get("relation") or "preprint_of",
                "confidence": conf,
                "high_confidence": bool(t.get("high_confidence", conf >= 0.85)),
                "evidence": evidence,
                "preprint_title": None,
                "journal_title": None,
                "via": "http",
            }
        )
    return proposals

def merge_preprint_journal(
    conn: sqlite3.Connection,
    preprint_id: str,
    journal_id: str,
    *,
    evidence: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Merge preprint into journal DOI canonical; keep preprint IDs as aliases.

    Unlike complementary strong-ID merge, different DOIs (10.1101 vs journal) are expected.
    Journal DOI wins publications.doi + publications.id; preprint DOI/arxiv stay in aliases.
    """
    pre = conn.execute("SELECT * FROM publications WHERE id=?", (preprint_id,)).fetchone()
    jou = conn.execute("SELECT * FROM publications WHERE id=?", (journal_id,)).fetchone()
    if not pre or not jou:
        return {"ok": False, "error": "missing preprint or journal row"}
    pre_d, jou_d = dict(pre), dict(jou)

    # Determine journal canonical
    jdoi = normalize_doi(jou_d.get("doi"))
    if not jdoi or is_preprint_doi(jdoi):
        # Fallback: maybe roles reversed
        if normalize_doi(pre_d.get("doi")) and not is_preprint_doi(pre_d.get("doi")):
            preprint_id, journal_id = journal_id, preprint_id
            pre_d, jou_d = jou_d, pre_d
            jdoi = normalize_doi(jou_d.get("doi"))
    if not jdoi:
        # Prefer pmid canonical if no journal DOI
        desired = canonical_id_from_ids(
            doi=None,
            pmid=jou_d.get("pmid") or pre_d.get("pmid"),
            pmcid=jou_d.get("pmcid") or pre_d.get("pmcid"),
            arxiv=None,
            openalex=jou_d.get("openalex") or pre_d.get("openalex"),
        )
    else:
        desired = f"doi:{jdoi}"

    preprint_doi = normalize_doi(pre_d.get("doi")) if is_preprint_doi(pre_d.get("doi")) else None
    preprint_arxiv = normalize_arxiv(pre_d.get("arxiv")) if pre_d.get("arxiv") else None

    # Clear preprint DOI from column before cluster merge so UNIQUE/complement path is clean
    # We do a specialized merge: reassign journal to desired if needed, absorb preprint.
    if jou_d["id"] != desired:
        if conn.execute("SELECT 1 FROM publications WHERE id=?", (desired,)).fetchone():
            if desired != preprint_id:
                # desired exists as third row — merge journal into it first
                merge_publication_cluster(conn, [jou_d["id"], desired], preferred_id=desired, clear_review=True)
                jou_d = dict(conn.execute("SELECT * FROM publications WHERE id=?", (desired,)).fetchone())
            else:
                # desired somehow is preprint id — force rename path after clearing
                pass
        else:
            _reassign_publication_id(conn, jou_d["id"], desired)
            jou_d = dict(conn.execute("SELECT * FROM publications WHERE id=?", (desired,)).fetchone())
            if preprint_id == desired:
                # extremely odd; bail
                return {"ok": False, "error": "preprint_id collided with desired journal id"}
            # preprint_id unchanged

    keeper = desired
    # Refresh preprint row (may still exist)
    if not conn.execute("SELECT 1 FROM publications WHERE id=?", (preprint_id,)).fetchone():
        return {"ok": False, "error": f"preprint {preprint_id} disappeared"}

    # Manually fold fields: journal DOI stays; preprint DOI only as alias
    absorb = preprint_id
    all_rows = [
        dict(conn.execute("SELECT * FROM publications WHERE id=?", (keeper,)).fetchone()),
        dict(conn.execute("SELECT * FROM publications WHERE id=?", (absorb,)).fetchone()),
    ]
    now = utc_now()
    # Prefer journal DOI in column
    new_doi = jdoi or _pick_richest_value(all_rows, "doi")
    if new_doi and is_preprint_doi(new_doi) and jdoi:
        new_doi = jdoi
    new_pmid = _pick_richest_value(all_rows, "pmid")
    new_pmcid = _pick_richest_value(all_rows, "pmcid")
    new_arxiv = _pick_richest_value(all_rows, "arxiv") or preprint_arxiv
    new_openalex = _pick_richest_value(all_rows, "openalex")
    new_title = _pick_richest_value(all_rows, "title")
    new_abstract = _pick_richest_value(all_rows, "abstract")
    new_authors = _pick_richest_value(all_rows, "authors_json")
    new_year = _pick_richest_value(all_rows, "year")
    new_venue = _pick_richest_value(all_rows, "venue")
    new_type = _pick_richest_value([all_rows[0]], "type") or "journal-article"
    new_tags = _pick_richest_value(all_rows, "tags_json")

    conn.execute(
        "UPDATE publications SET doi=NULL, pmid=NULL, pmcid=NULL, arxiv=NULL, openalex=NULL WHERE id=?",
        (absorb,),
    )
    conn.execute(
        """
        UPDATE publications SET
          doi=?, pmid=?, pmcid=?, arxiv=?, openalex=?,
          title=?, abstract=?, authors_json=?, year=?, venue=?, type=?,
          tags_json=?, needs_review=0, review_reasons_json=NULL, updated_at=?
        WHERE id=?
        """,
        (
            normalize_doi(new_doi) if new_doi else None,
            normalize_pmid(str(new_pmid)) if new_pmid else None,
            normalize_pmcid(str(new_pmcid)) if new_pmcid else None,
            normalize_arxiv(str(new_arxiv)) if new_arxiv else None,
            normalize_openalex(str(new_openalex)) if new_openalex else None,
            new_title, new_abstract, new_authors, new_year, new_venue, new_type,
            new_tags, now, keeper,
        ),
    )
    _fts_upsert(conn, keeper, new_title or "", new_abstract or "", new_tags or "")

    # Repoint dependents from absorb → keeper
    _record_id_remap(conn, absorb, keeper)
    conn.execute("UPDATE aliases SET publication_id=? WHERE publication_id=?", (keeper, absorb))
    for cp in conn.execute(
        "SELECT collection_id, added_at, source_run_id FROM collection_pubs WHERE publication_id=?",
        (absorb,),
    ).fetchall():
        conn.execute(
            "INSERT OR IGNORE INTO collection_pubs(collection_id, publication_id, added_at, source_run_id) VALUES (?,?,?,?)",
            (cp["collection_id"], keeper, cp["added_at"], cp["source_run_id"]),
        )
    conn.execute("DELETE FROM collection_pubs WHERE publication_id=?", (absorb,))
    conn.execute("UPDATE run_hits SET publication_id=? WHERE publication_id=?", (keeper, absorb))
    for tp in conn.execute(
        "SELECT node_id FROM theme_pubs WHERE publication_id=?", (absorb,)
    ).fetchall():
        conn.execute(
            "INSERT OR IGNORE INTO theme_pubs(node_id, publication_id) VALUES (?, ?)",
            (tp["node_id"], keeper),
        )
    conn.execute("DELETE FROM theme_pubs WHERE publication_id=?", (absorb,))
    _repoint_publication_versions(conn, absorb, keeper)
    try:
        conn.execute("DELETE FROM pubs_fts WHERE publication_id=?", (absorb,))
    except sqlite3.OperationalError:
        pass
    conn.execute("DELETE FROM publications WHERE id=?", (absorb,))

    # Aliases: journal + preprint ids
    if new_doi:
        conn.execute(
            "INSERT OR REPLACE INTO aliases(alias, publication_id) VALUES (?, ?)",
            (f"doi:{normalize_doi(new_doi)}", keeper),
        )
    if preprint_doi:
        conn.execute(
            "INSERT OR REPLACE INTO aliases(alias, publication_id) VALUES (?, ?)",
            (f"doi:{preprint_doi}", keeper),
        )
    if new_arxiv or preprint_arxiv:
        ax = normalize_arxiv(str(new_arxiv or preprint_arxiv))
        if ax:
            conn.execute(
                "INSERT OR REPLACE INTO aliases(alias, publication_id) VALUES (?, ?)",
                (f"arxiv:{ax}", keeper),
            )
    for key, val in (
        ("pmid", new_pmid),
        ("pmcid", new_pmcid),
        ("openalex", new_openalex),
    ):
        if not val:
            continue
        norm = {
            "pmid": normalize_pmid,
            "pmcid": normalize_pmcid,
            "openalex": normalize_openalex,
        }[key](str(val))
        if norm:
            conn.execute(
                "INSERT OR REPLACE INTO aliases(alias, publication_id) VALUES (?, ?)",
                (f"{key}:{norm}", keeper),
            )

    # Record confirmed version link (self-loop avoided — use alias evidence only)
    # from_id was absorbed; store a confirmed same_work note with keeper as to_id
    # using a synthetic from alias marker in evidence instead of FK to deleted id.
    # Mark any pending links involving keeper as confirmed (repoint already collapsed absorb)
    ev = json.dumps(
        evidence or {"merged": True, "preprint_doi": preprint_doi, "absorbed": absorb},
        ensure_ascii=False,
    )
    conn.execute(
        "UPDATE publication_versions SET needs_confirm=0, "
        "confidence=CASE WHEN confidence < 0.95 THEN 0.95 ELSE confidence END, "
        "updated_at=?, evidence_json=COALESCE(?, evidence_json) "
        "WHERE (from_id=? OR to_id=?) AND needs_confirm=1",
        (now, ev, keeper, keeper),
    )

    return {
        "ok": True,
        "publication_id": keeper,
        "absorbed_ids": [absorb],
        "journal_doi": new_doi,
        "preprint_doi_alias": preprint_doi,
        "preprint_arxiv_alias": preprint_arxiv or new_arxiv,
        "title": new_title,
    }


def confirm_version_link(
    conn: sqlite3.Connection,
    link_id: int,
    *,
    merge: bool = True,
) -> dict[str, Any]:
    """Confirm a pending version link; optionally merge into journal DOI canonical."""
    row = conn.execute("SELECT * FROM publication_versions WHERE id=?", (link_id,)).fetchone()
    if not row:
        return {"ok": False, "error": f"unknown link_id {link_id}"}
    d = dict(row)
    now = utc_now()
    conn.execute(
        "UPDATE publication_versions SET needs_confirm=0, updated_at=? WHERE id=?",
        (now, link_id),
    )
    result: dict[str, Any] = {
        "ok": True,
        "link_id": link_id,
        "from_id": d["from_id"],
        "to_id": d["to_id"],
        "relation": d["relation"],
        "merged": False,
    }
    if merge and d["relation"] in ("preprint_of", "published_as", "same_work"):
        # Orient: preprint → journal
        a = dict(conn.execute("SELECT * FROM publications WHERE id=?", (d["from_id"],)).fetchone() or {})
        b = dict(conn.execute("SELECT * FROM publications WHERE id=?", (d["to_id"],)).fetchone() or {})
        if not a or not b:
            result["note"] = "link confirmed; one side missing so merge skipped"
            return result
        if d["relation"] == "published_as":
            # from journal? published_as usually journal→? We treat from as preprint-ish source
            pre_id, jou_id = d["from_id"], d["to_id"]
        elif is_preprint_pub(a) and is_journal_pub(b):
            pre_id, jou_id = d["from_id"], d["to_id"]
        elif is_preprint_pub(b) and is_journal_pub(a):
            pre_id, jou_id = d["to_id"], d["from_id"]
        else:
            # Prefer non-preprint DOI as journal
            if is_preprint_doi(a.get("doi")) and not is_preprint_doi(b.get("doi")):
                pre_id, jou_id = d["from_id"], d["to_id"]
            elif is_preprint_doi(b.get("doi")) and not is_preprint_doi(a.get("doi")):
                pre_id, jou_id = d["to_id"], d["from_id"]
            else:
                result["note"] = "link confirmed; could not orient preprint/journal for merge"
                return result
        evidence = {}
        if d.get("evidence_json"):
            try:
                evidence = json.loads(d["evidence_json"])
            except json.JSONDecodeError:
                evidence = {"raw": d["evidence_json"]}
        merged = merge_preprint_journal(conn, pre_id, jou_id, evidence=evidence)
        result["merged"] = bool(merged.get("ok"))
        result["merge"] = merged
    return result


def apply_version_detection(
    conn: sqlite3.Connection,
    *,
    collection_id: Optional[str] = None,
    auto_merge_high_confidence: bool = False,
    confidence_threshold: float = 0.85,
    extra_proposals: Optional[list[dict[str, Any]]] = None,
) -> dict[str, Any]:
    """Detect preprint↔journal candidates; propose links; optionally auto-merge high-confidence."""
    proposals = merge_version_proposals(
        detect_version_candidates(conn, collection_id=collection_id),
        extra_proposals or [],
    )
    proposed = []
    merged = []
    for p in proposals:
        link = propose_version_link(
            conn,
            p["from_id"],
            p["to_id"],
            relation=p["relation"],
            confidence=p["confidence"],
            needs_confirm=0 if (auto_merge_high_confidence and p["high_confidence"]) else 1,
            evidence=p.get("evidence"),
        )
        proposed.append({**p, "link": link})
        if auto_merge_high_confidence and p["high_confidence"] and p["confidence"] >= confidence_threshold:
            m = merge_preprint_journal(
                conn,
                p["from_id"],
                p["to_id"],
                evidence=p.get("evidence"),
            )
            if m.get("ok"):
                merged.append(m)
    return {
        "ok": True,
        "n_proposals": len(proposals),
        "n_links_written": sum(1 for x in proposed if x.get("link", {}).get("ok")),
        "n_merged": len(merged),
        "auto_merge_high_confidence": auto_merge_high_confidence,
        "proposals": proposed,
        "merged": merged,
    }




# ---------------------------------------------------------------------------
# 0.5.2: gold-set recall + OOW controls + uncontacted sources + saturation
# ---------------------------------------------------------------------------


def parse_strong_id_token(raw: Any) -> dict[str, Any]:
    """Soft-parse a strong-ID token (doi/pmid/pmcid/arxiv/openalex). Format only; no network."""
    s = str(raw or "").strip()
    out: dict[str, Any] = {
        "raw": s,
        "kind": None,
        "value": None,
        "alias": None,
        "format_ok": False,
    }
    if not s:
        out["error"] = "empty"
        return out
    low = s.lower()
    if low.startswith("doi:") or "doi.org/" in low or (s.count("/") >= 1 and not low.startswith(("pmid:", "pmcid:", "arxiv:", "openalex:", "pmc"))):
        doi = normalize_doi(s)
        if doi and re.match(r"^10\.\d{4,9}/\S+$", doi):
            out.update(kind="doi", value=doi, alias=f"doi:{doi}", format_ok=True)
            return out
        out.update(kind="doi", value=doi, alias=f"doi:{doi}" if doi else None, format_ok=False, error="doi_format")
        return out
    if low.startswith("pmid:") or low.startswith("pubmed:") or re.fullmatch(r"\d{1,9}", s):
        pmid = normalize_pmid(s)
        if pmid:
            out.update(kind="pmid", value=pmid, alias=f"pmid:{pmid}", format_ok=True)
            return out
        out.update(kind="pmid", format_ok=False, error="pmid_format")
        return out
    if low.startswith("pmcid:") or re.match(r"^pmc\d+$", low):
        pmcid = normalize_pmcid(s)
        if pmcid:
            out.update(kind="pmcid", value=pmcid, alias=f"pmcid:{pmcid}", format_ok=True)
            return out
        out.update(kind="pmcid", format_ok=False, error="pmcid_format")
        return out
    if low.startswith("arxiv:") or re.match(r"^\d{4}\.\d{4,5}(v\d+)?$", s) or re.match(r"^[a-z\-]+/\d+$", s, re.I):
        arxiv = normalize_arxiv(s)
        if arxiv:
            out.update(kind="arxiv", value=arxiv, alias=f"arxiv:{arxiv}", format_ok=True)
            return out
        out.update(kind="arxiv", format_ok=False, error="arxiv_format")
        return out
    if low.startswith("openalex:") or re.match(r"^w\d+$", low) or "openalex.org/w" in low:
        oa = normalize_openalex(s)
        if oa:
            out.update(kind="openalex", value=oa, alias=f"openalex:{oa}", format_ok=True)
            return out
        out.update(kind="openalex", format_ok=False, error="openalex_format")
        return out
    # Last-chance DOI
    doi = normalize_doi(s)
    if doi and re.match(r"^10\.\d{4,9}/\S+$", doi):
        out.update(kind="doi", value=doi, alias=f"doi:{doi}", format_ok=True)
        return out
    out["error"] = "unrecognized_id"
    return out


def normalize_id_token_list(items: list[Any]) -> tuple[list[str], list[dict[str, Any]]]:
    """Normalize a list of ID tokens; prefer alias form; keep order; drop empties."""
    seen: set[str] = set()
    norms: list[str] = []
    meta: list[dict[str, Any]] = []
    for it in items or []:
        parsed = parse_strong_id_token(it)
        meta.append(parsed)
        key = parsed.get("alias") or (str(it).strip() if it is not None else "")
        if not key or key in seen:
            continue
        seen.add(key)
        norms.append(key)
    return norms, meta


def resolve_id_token_publication(
    conn: sqlite3.Connection, token: str
) -> Optional[str]:
    """Map a strong-ID token / alias to publications.id if present."""
    parsed = parse_strong_id_token(token)
    candidates: list[str] = []
    if parsed.get("alias"):
        candidates.append(parsed["alias"])
    raw = str(token or "").strip()
    if raw and raw not in candidates:
        candidates.append(raw)
    for a in candidates:
        row = conn.execute(
            "SELECT publication_id FROM aliases WHERE alias=?", (a,)
        ).fetchone()
        if row:
            return row["publication_id"]
    # Column lookup by kind
    kind = parsed.get("kind")
    val = parsed.get("value")
    if kind and val and kind in STRONG_ID_KEYS:
        row = conn.execute(
            f"SELECT id FROM publications WHERE {kind}=?", (val,)
        ).fetchone()
        if row:
            return row["id"]
    return None


def _corpus_alias_index(
    conn: sqlite3.Connection, pub_ids: set[str]
) -> set[str]:
    """All aliases + strong-ID aliases for pubs in corpus."""
    out: set[str] = set()
    if not pub_ids:
        return out
    for pid in pub_ids:
        out.add(pid)
        for row in conn.execute(
            "SELECT alias FROM aliases WHERE publication_id=?", (pid,)
        ).fetchall():
            out.add(row["alias"])
        prow = conn.execute("SELECT * FROM publications WHERE id=?", (pid,)).fetchone()
        if not prow:
            continue
        d = dict(prow)
        if d.get("doi"):
            out.add(f"doi:{normalize_doi(d['doi'])}")
        if d.get("pmid"):
            out.add(f"pmid:{normalize_pmid(d['pmid'])}")
        if d.get("pmcid"):
            out.add(f"pmcid:{normalize_pmcid(d['pmcid'])}")
        if d.get("arxiv"):
            out.add(f"arxiv:{normalize_arxiv(d['arxiv'])}")
        if d.get("openalex"):
            out.add(f"openalex:{normalize_openalex(d['openalex'])}")
    return {x for x in out if x}


def evaluate_gold_set_recall(
    conn: sqlite3.Connection,
    gold_set: list[str],
    *,
    corpus_pub_ids: set[str],
) -> dict[str, Any]:
    """How many gold_set IDs appear in the corpus (hits and/or collection pubs)."""
    norms, meta = normalize_id_token_list(gold_set or [])
    n = len(norms)
    evidence_weak = n < GOLD_SET_MIN_GUIDANCE
    empty = n == 0
    index = _corpus_alias_index(conn, corpus_pub_ids)
    found: list[dict[str, Any]] = []
    missed: list[dict[str, Any]] = []
    for token in norms:
        parsed = parse_strong_id_token(token)
        alias = parsed.get("alias") or token
        pub_id = None
        present = alias in index or token in index
        if not present:
            # also try resolve against full DB then check membership
            pub_id = resolve_id_token_publication(conn, token)
            present = bool(pub_id and pub_id in corpus_pub_ids)
        else:
            pub_id = resolve_id_token_publication(conn, token)
        entry = {
            "id": alias,
            "raw": parsed.get("raw") or token,
            "kind": parsed.get("kind"),
            "format_ok": bool(parsed.get("format_ok")),
            "publication_id": pub_id,
        }
        if present:
            found.append(entry)
        else:
            missed.append(entry)
    recall = (len(found) / n) if n else None
    return {
        "gold_set_n": n,
        "gold_set_min_guidance": GOLD_SET_MIN_GUIDANCE,
        "evidence_weak": evidence_weak,
        "gold_set_empty": empty,
        "n_found": len(found),
        "n_missed": len(missed),
        "recall": recall,
        "found": found,
        "missed": missed,
        "note": (
            "gold_set is empty — no recall evidence"
            if empty
            else (
                f"gold_set size {n} < guidance {GOLD_SET_MIN_GUIDANCE} — evidence_weak"
                if evidence_weak
                else "gold-set recall against ingested/hit corpus"
            )
        ),
    }


def evaluate_out_of_window_controls(
    conn: sqlite3.Connection,
    controls: list[str],
    *,
    corpus_pub_ids: set[str],
) -> dict[str, Any]:
    """Controls outside the date window must NOT appear in the windowed corpus."""
    norms, _meta = normalize_id_token_list(controls or [])
    index = _corpus_alias_index(conn, corpus_pub_ids)
    leaked: list[dict[str, Any]] = []
    absent: list[dict[str, Any]] = []
    for token in norms:
        parsed = parse_strong_id_token(token)
        alias = parsed.get("alias") or token
        pub_id = resolve_id_token_publication(conn, token)
        present = (alias in index or token in index) or (
            bool(pub_id and pub_id in corpus_pub_ids)
        )
        entry = {
            "id": alias,
            "raw": parsed.get("raw") or token,
            "kind": parsed.get("kind"),
            "format_ok": bool(parsed.get("format_ok")),
            "publication_id": pub_id,
        }
        if present:
            leaked.append(entry)
        else:
            absent.append(entry)
    return {
        "n_controls": len(norms),
        "n_leaked": len(leaked),
        "n_absent_ok": len(absent),
        "filter_ok": len(leaked) == 0,
        "leaked": leaked,
        "absent_ok": absent,
        "note": (
            "no out_of_window_controls configured"
            if not norms
            else (
                "date filter leakage detected"
                if leaked
                else "out-of-window controls absent from windowed corpus (ok)"
            )
        ),
    }


def evaluate_uncontacted_sources(
    *,
    configured_sources: list[str],
    sources_executed: Optional[list[str]] = None,
    providers: Optional[list[dict[str, Any]]] = None,
    hit_source_counts: Optional[dict[str, int]] = None,
) -> dict[str, Any]:
    """Sources listed in run config but never contacted / skipped / errored / unimplemented."""
    configured = [str(s).strip().lower() for s in (configured_sources or []) if str(s).strip()]
    executed = {
        str(s).strip().lower()
        for s in (sources_executed or [])
        if str(s).strip()
    }
    provider_by_src: dict[str, dict[str, Any]] = {}
    for p in providers or []:
        src = str(p.get("source") or p.get("provider") or "").strip().lower()
        if src:
            provider_by_src[src] = p
    hit_sources = {
        str(k).strip().lower() for k in (hit_source_counts or {}).keys() if k
    }

    items: list[dict[str, Any]] = []
    contacted: list[str] = []
    for src in configured:
        meta = provider_by_src.get(src) or {}
        in_executed = src in executed if sources_executed is not None else None
        has_hits = src in hit_sources
        provider_status = str(meta.get("provider") or meta.get("status") or "").lower()
        err = meta.get("error") or meta.get("message")
        reason = None
        contacted_flag = False
        if sources_executed is not None and not executed and not providers:
            # Search never ran
            reason = "search_not_executed"
        elif sources_executed is not None and src not in executed:
            reason = "not_in_sources_executed"
        elif provider_status in {"error", "soft-fail", "soft_fail", "unimplemented", "skipped"}:
            reason = f"provider_{provider_status}"
            contacted_flag = True  # attempted
        elif err and not has_hits and (in_executed or src in provider_by_src):
            reason = "provider_error_or_empty"
            contacted_flag = True
        elif sources_executed is None and not providers and not hit_source_counts:
            reason = "no_execution_metadata"
        elif src in provider_by_src or has_hits or (in_executed is True):
            contacted_flag = True
            reason = None
        else:
            reason = "declared_but_uncontacted"

        if contacted_flag and reason is None:
            contacted.append(src)
        else:
            items.append(
                {
                    "source": src,
                    "reason": reason or "declared_but_uncontacted",
                    "provider_meta": meta or None,
                    "had_hits": has_hits,
                }
            )

    # Always emit a section — blank must not imply "no gaps"
    return {
        "configured_sources": configured,
        "contacted_sources": contacted,
        "uncontacted": items,
        "n_uncontacted": len(items),
        "note": (
            "all configured sources appear contacted"
            if configured and not items
            else (
                "no sources configured"
                if not configured
                else "one or more configured sources were not contacted / failed / skipped"
            )
        ),
    }


def assess_snowball_saturation(
    rounds: list[dict[str, Any]],
    *,
    max_rounds: int,
    low_yield_ratio: float = 0.05,
    consecutive_low: int = 2,
) -> dict[str, Any]:
    """Honest saturation accounting for expand/snowball rounds.

    saturated when:
      (a) consecutive `consecutive_low` rounds with new_included < low_yield_ratio * n_seeds, OR
      (b) no new eligible seeds (n_new_included == 0 / empty next seed queue)
    If max_rounds hit without (a)/(b): saturated=False, stop_reason=max_rounds_not_saturated.
    Never claim saturated falsely.
    """
    max_rounds = max(1, int(max_rounds))
    low_yield_ratio = float(low_yield_ratio)
    consecutive_low = max(1, int(consecutive_low))
    trailing_low = 0
    stop_reason = "in_progress"
    saturated = False

    if not rounds:
        return {
            "saturated": False,
            "stop_reason": "no_rounds",
            "max_rounds": max_rounds,
            "n_rounds": 0,
            "consecutive_low_yield": 0,
            "low_yield_ratio": low_yield_ratio,
            "rounds": [],
            "closure_claimed": False,
            "note": "no expand rounds run — closure NOT claimed",
        }

    for i, rnd in enumerate(rounds):
        n_seeds = int(rnd.get("n_seeds") or 0)
        n_new = int(rnd.get("n_new_included") or 0)
        if n_seeds <= 0 or n_new == 0:
            trailing_low = consecutive_low  # treat as terminal low / no seeds
            saturated = True
            stop_reason = "no_new_eligible_seeds"
            break
        threshold = low_yield_ratio * n_seeds
        if n_new < threshold:
            trailing_low += 1
        else:
            trailing_low = 0
        if trailing_low >= consecutive_low:
            saturated = True
            stop_reason = "consecutive_low_yield"
            break
    else:
        # completed loop without break
        if len(rounds) >= max_rounds:
            saturated = False
            stop_reason = "max_rounds_not_saturated"
        else:
            stop_reason = "in_progress"

    # If we stopped early due to no_new but still under max — saturated ok
    if stop_reason == "no_new_eligible_seeds":
        saturated = True
    if stop_reason == "consecutive_low_yield":
        saturated = True
    if stop_reason == "max_rounds_not_saturated":
        saturated = False

    return {
        "saturated": bool(saturated),
        "stop_reason": stop_reason,
        "max_rounds": max_rounds,
        "n_rounds": len(rounds),
        "consecutive_low_yield": trailing_low,
        "low_yield_ratio": low_yield_ratio,
        "rounds": rounds,
        "closure_claimed": bool(saturated),
        "note": (
            f"saturated via {stop_reason}"
            if saturated
            else (
                "hit max_rounds without saturation — closure NOT claimed"
                if stop_reason == "max_rounds_not_saturated"
                else f"saturation incomplete ({stop_reason})"
            )
        ),
    }



# ---------------------------------------------------------------------------
# 0.5.0: protocol freeze, coverage gate, preview-fetch opposite DOI
# ---------------------------------------------------------------------------


def build_frozen_protocol(
    *,
    run_id: str,
    collection_id: str,
    query: str,
    sources: list[str],
    date_from: Any = None,
    date_to: Any = None,
    pub_type: str = "any",
    limit: int = 20,
    provider_notes: Optional[dict[str, str]] = None,
    gold_set: Optional[list[str]] = None,
    out_of_window_controls: Optional[list[str]] = None,
    extra: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Canonical frozen search protocol stored in run.json / runs.config_json."""
    notes = dict(PROVIDER_PROTOCOL_NOTES)
    if provider_notes:
        notes.update(provider_notes)
    src_notes = {s: notes.get(s, "open HTTP adapter") for s in sources}
    gold_norm, gold_meta = normalize_id_token_list(gold_set or [])
    oow_norm, oow_meta = normalize_id_token_list(out_of_window_controls or [])
    cfg: dict[str, Any] = {
        "protocol_version": PROTOCOL_VERSION,
        "run_id": run_id,
        "collection_id": collection_id,
        "query": query,
        "sources": list(sources),
        "date_from": date_from or None,
        "date_to": date_to or None,
        "type": pub_type or "any",
        "limit": int(limit),
        "limit_per_source": int(limit),
        "provider_notes": src_notes,
        "gold_set": gold_norm,
        "gold_set_meta": gold_meta,
        "out_of_window_controls": oow_norm,
        "out_of_window_meta": oow_meta,
        "gold_set_min_guidance": GOLD_SET_MIN_GUIDANCE,
        "frozen_at": utc_now(),
        "created_at": utc_now(),
        "no_commercial_sources": True,
        "note_gold_verification": (
            "Caller should verify gold DOIs (e.g. doi.org); plugin soft-checks format only offline."
        ),
    }
    if extra:
        cfg.update(extra)
    return cfg


def evaluate_coverage_gate(
    report: dict[str, Any],
    *,
    min_unique_pubs: Optional[int] = None,
    max_needs_review: Optional[int] = None,
    max_version_pending: Optional[int] = None,
    max_unenriched_shells: Optional[int] = None,
    fail_under_coverage: bool = False,
) -> dict[str, Any]:
    """Decide passed vs needs_attention from a coverage_report payload."""
    hits = report.get("hits") or {}
    coll = report.get("collection") or {}
    unique = int(hits.get("unique_publication_ids") or 0)
    n_review = int(coll.get("n_needs_review") or 0)
    n_ver = int(coll.get("n_version_links_pending") or 0)
    n_shells = int(coll.get("n_unenriched_shells") or 0)
    breaches: list[dict[str, Any]] = []

    def _breach(code: str, actual: int, threshold: int, cmp: str) -> None:
        breaches.append(
            {"code": code, "actual": actual, "threshold": threshold, "cmp": cmp}
        )

    if min_unique_pubs is not None and unique < int(min_unique_pubs):
        _breach("min_unique_pubs", unique, int(min_unique_pubs), "lt")
    if max_needs_review is not None and n_review > int(max_needs_review):
        _breach("max_needs_review", n_review, int(max_needs_review), "gt")
    if max_version_pending is not None and n_ver > int(max_version_pending):
        _breach("max_version_pending", n_ver, int(max_version_pending), "gt")
    if max_unenriched_shells is not None and n_shells > int(max_unenriched_shells):
        _breach("max_unenriched_shells", n_shells, int(max_unenriched_shells), "gt")

    status = RUN_STATUS_PASSED if not breaches else RUN_STATUS_NEEDS_ATTENTION
    return {
        "ok": True,
        "status": status,
        "passed": status == RUN_STATUS_PASSED,
        "breaches": breaches,
        "metrics": {
            "unique_publication_ids": unique,
            "n_needs_review": n_review,
            "n_version_links_pending": n_ver,
            "n_unenriched_shells": n_shells,
            "n_hits": int(hits.get("n") or 0),
        },
        "thresholds": {
            "min_unique_pubs": min_unique_pubs,
            "max_needs_review": max_needs_review,
            "max_version_pending": max_version_pending,
            "max_unenriched_shells": max_unenriched_shells,
            "fail_under_coverage": bool(fail_under_coverage),
        },
        "fail_under_coverage": bool(fail_under_coverage),
        "should_fail_cli": bool(fail_under_coverage and breaches),
    }


def close_run(
    conn: sqlite3.Connection,
    run_id: str,
    *,
    min_unique_pubs: Optional[int] = None,
    max_needs_review: Optional[int] = None,
    max_version_pending: Optional[int] = None,
    max_unenriched_shells: Optional[int] = None,
    fail_under_coverage: bool = False,
    write_coverage_json: Optional[Path] = None,
    write_coverage_md: Optional[Path] = None,
) -> dict[str, Any]:
    """Write coverage artifacts and mark run status passed/needs_attention."""
    report = coverage_report(conn, run_id=run_id)
    if not report.get("ok"):
        return report
    gate = evaluate_coverage_gate(
        report,
        min_unique_pubs=min_unique_pubs,
        max_needs_review=max_needs_review,
        max_version_pending=max_version_pending,
        max_unenriched_shells=max_unenriched_shells,
        fail_under_coverage=fail_under_coverage,
    )
    run_row = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
    if not run_row:
        return {"ok": False, "error": f"unknown run_id {run_id}"}
    run_d = dict(run_row)
    col_id = run_d.get("collection_id")
    coverage_path = None
    md_path = None
    payload = {
        **report,
        "gate": gate,
        "closed_at": utc_now(),
        "protocol_version": PROTOCOL_VERSION,
    }
    # Prefer writing next to hits if path known
    hits_path = run_d.get("hits_path")
    run_dir = Path(hits_path).parent if hits_path else None
    if write_coverage_json:
        coverage_path = Path(write_coverage_json)
    elif run_dir:
        coverage_path = run_dir / "coverage.json"
    if coverage_path:
        coverage_path.parent.mkdir(parents=True, exist_ok=True)
        coverage_path.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
        )
    if write_coverage_md:
        md_path = Path(write_coverage_md)
    elif run_dir:
        md_path = run_dir / "coverage.md"
    if md_path:
        hits = report.get("hits") or {}
        coll = report.get("collection") or {}
        lines = [
            "# Literature coverage report",
            "",
            f"- run_id: `{run_id}`",
            f"- collection_id: `{col_id or '—'}`",
            f"- gate_status: **{gate['status']}**",
            "",
            "## Hits",
            f"- n_hits: **{hits.get('n', 0)}**",
            f"- unique_publication_ids: **{hits.get('unique_publication_ids', 0)}**",
            "",
            "## Collection",
            f"- n_publications: **{coll.get('n_publications', 0)}**",
            f"- n_needs_review: {coll.get('n_needs_review', 0)}",
            f"- n_unenriched_shells: {coll.get('n_unenriched_shells', 0)}",
            f"- n_version_links_pending: {coll.get('n_version_links_pending', 0)}",
            "",
            "## Gate breaches",
        ]
        if gate["breaches"]:
            for b in gate["breaches"]:
                lines.append(
                    f"- `{b['code']}`: actual={b['actual']} threshold={b['threshold']} ({b['cmp']})"
                )
        else:
            lines.append("- none")
        lines.append("")
        revid = report.get("recall_evidence") or {}
        gold = revid.get("gold_set_recall") or {}
        oow = revid.get("out_of_window") or {}
        unc = revid.get("declared_but_uncontacted_sources") or {}
        lines.extend(
            [
                "## Gold-set recall",
                f"- gold_set_n: **{gold.get('gold_set_n', 0)}**",
                f"- found: **{gold.get('n_found', 0)}** / missed: **{gold.get('n_missed', 0)}**",
                f"- recall: `{gold.get('recall')}`",
                f"- evidence_weak: **{gold.get('evidence_weak')}**",
                f"- note: {gold.get('note') or '—'}",
                "",
                "## Out-of-window controls",
                f"- n_controls: {oow.get('n_controls', 0)}",
                f"- n_leaked: **{oow.get('n_leaked', 0)}** (filter_ok={oow.get('filter_ok')})",
                f"- note: {oow.get('note') or '—'}",
                "",
                "## Declared-but-uncontacted sources",
                f"- n_uncontacted: **{unc.get('n_uncontacted', 0)}**",
                f"- note: {unc.get('note') or '—'}",
            ]
        )
        for item in unc.get("uncontacted") or []:
            lines.append(f"- `{item.get('source')}`: {item.get('reason')}")
        if not (unc.get("uncontacted") or []):
            lines.append("- (section present; none listed)")
        lines.append("")
        md_path.parent.mkdir(parents=True, exist_ok=True)
        md_path.write_text("\n".join(lines), encoding="utf-8")

    # Merge gate into config_json
    cfg = {}
    if run_d.get("config_json"):
        try:
            cfg = json.loads(run_d["config_json"])
        except json.JSONDecodeError:
            cfg = {}
    cfg["gate"] = gate
    cfg["closed_at"] = utc_now()
    cfg["status"] = gate["status"]
    now = utc_now()
    conn.execute(
        """
        UPDATE runs SET status=?, coverage_path=?, gate_json=?, config_json=?,
          finished_at=COALESCE(finished_at, ?)
        WHERE id=?
        """,
        (
            gate["status"],
            str(coverage_path) if coverage_path else run_d.get("coverage_path"),
            json.dumps(gate, ensure_ascii=False),
            json.dumps(cfg, ensure_ascii=False),
            now,
            run_id,
        ),
    )
    # Sync run.json on disk when possible
    if run_dir and (run_dir / "run.json").is_file():
        try:
            disk = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            disk = {}
        disk.update(
            {
                "status": gate["status"],
                "gate": gate,
                "closed_at": cfg["closed_at"],
                "coverage_json": str(coverage_path) if coverage_path else None,
                "coverage_md": str(md_path) if md_path else None,
            }
        )
        (run_dir / "run.json").write_text(
            json.dumps(disk, ensure_ascii=False, indent=2), encoding="utf-8"
        )

    return {
        "ok": True,
        "run_id": run_id,
        "collection_id": col_id,
        "status": gate["status"],
        "gate": gate,
        "coverage": report,
        "coverage_json": str(coverage_path) if coverage_path else None,
        "coverage_md": str(md_path) if md_path else None,
        "should_fail_cli": gate["should_fail_cli"],
    }


def set_run_status(
    conn: sqlite3.Connection,
    run_id: str,
    status: str,
    *,
    extra_config: Optional[dict[str, Any]] = None,
    hits_path: Optional[str] = None,
) -> None:
    row = conn.execute("SELECT config_json FROM runs WHERE id=?", (run_id,)).fetchone()
    cfg: dict[str, Any] = {}
    if row and row["config_json"]:
        try:
            cfg = json.loads(row["config_json"])
        except json.JSONDecodeError:
            cfg = {}
    cfg["status"] = status
    if extra_config:
        cfg.update(extra_config)
    if hits_path:
        conn.execute(
            "UPDATE runs SET status=?, config_json=?, hits_path=?, finished_at=? WHERE id=?",
            (status, json.dumps(cfg, ensure_ascii=False), hits_path, utc_now(), run_id),
        )
    else:
        conn.execute(
            "UPDATE runs SET status=?, config_json=? WHERE id=?",
            (status, json.dumps(cfg, ensure_ascii=False), run_id),
        )


def unresolved_doi_ties(
    conn: sqlite3.Connection,
    ties: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """DOI ties where exactly one side is already in the library."""

    def _pub_for_doi(doi: str) -> Optional[str]:
        d = normalize_doi(doi)
        if not d:
            return None
        row = conn.execute("SELECT id FROM publications WHERE doi=?", (d,)).fetchone()
        if row:
            return row["id"]
        alias = conn.execute(
            "SELECT publication_id FROM aliases WHERE alias=?", (f"doi:{d}",)
        ).fetchone()
        return alias["publication_id"] if alias else None

    out: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for t in ties:
        from_doi = normalize_doi(t.get("from_doi"))
        to_doi = normalize_doi(t.get("to_doi"))
        if not from_doi or not to_doi or from_doi == to_doi:
            continue
        from_id = _pub_for_doi(from_doi)
        to_id = _pub_for_doi(to_doi)
        if bool(from_id) == bool(to_id):
            # both present (handled elsewhere) or both missing
            continue
        key = (from_doi, to_doi)
        if key in seen:
            continue
        seen.add(key)
        if from_id and not to_id:
            out.append(
                {
                    "present_id": from_id,
                    "present_doi": from_doi,
                    "missing_doi": to_doi,
                    "missing_role": "to",
                    "relation": t.get("relation") or "preprint_of",
                    "confidence": float(t.get("confidence") or 0.5),
                    "high_confidence": bool(t.get("high_confidence", False)),
                    "evidence": dict(t.get("evidence") or {}),
                    "provider": t.get("provider"),
                }
            )
        elif to_id and not from_id:
            out.append(
                {
                    "present_id": to_id,
                    "present_doi": to_doi,
                    "missing_doi": from_doi,
                    "missing_role": "from",
                    "relation": t.get("relation") or "preprint_of",
                    "confidence": float(t.get("confidence") or 0.5),
                    "high_confidence": bool(t.get("high_confidence", False)),
                    "evidence": dict(t.get("evidence") or {}),
                    "provider": t.get("provider"),
                }
            )
    return out


def preview_fetch_opposite(
    conn: sqlite3.Connection,
    *,
    missing_doi: str,
    present_id: str,
    relation: str = "preprint_of",
    confidence: float = 0.5,
    evidence: Optional[dict[str, Any]] = None,
    collection_id: Optional[str] = None,
    hit: Optional[dict[str, Any]] = None,
    missing_role: str = "to",
) -> dict[str, Any]:
    """Create metadata-only shell for missing opposite DOI + pending version link.

    Never downloads PDFs. Shell is needs_review with preview_fetched_opposite.
    """
    doi = normalize_doi(missing_doi)
    if not doi:
        return {"ok": False, "error": "invalid missing_doi"}
    if not conn.execute("SELECT 1 FROM publications WHERE id=?", (present_id,)).fetchone():
        return {"ok": False, "error": f"missing present publication {present_id}"}

    draft: dict[str, Any] = {"doi": doi}
    if hit:
        for k in (
            "doi",
            "pmid",
            "pmcid",
            "arxiv",
            "openalex",
            "title",
            "abstract",
            "authors",
            "year",
            "venue",
            "type",
            "url",
        ):
            if hit.get(k) is not None:
                draft[k] = hit.get(k)
        draft["doi"] = doi  # force normalized target DOI
    if not draft.get("title"):
        draft["title"] = f"[preview shell] doi:{doi}"
    if not draft.get("type"):
        draft["type"] = "preprint" if is_preprint_doi(doi) else "journal-article"

    created = resolve_or_create(conn, draft)
    new_id = created.get("publication_id") or created.get("id")
    if not new_id:
        return {"ok": False, "error": "resolve_or_create failed", "created": created}

    _append_review_reason(conn, new_id, REVIEW_PREVIEW, related=[present_id])
    _append_review_reason(conn, new_id, REVIEW_VERSION, related=[present_id])
    if _is_unenriched_shell(dict(conn.execute("SELECT * FROM publications WHERE id=?", (new_id,)).fetchone())):
        _append_review_reason(conn, new_id, REVIEW_SHELL)

    if collection_id:
        link_to_collection(conn, collection_id, new_id, source_run_id=None)

    ev = dict(evidence or {})
    ev["preview_fetch"] = True
    ev["missing_doi"] = doi
    ev["present_id"] = present_id
    ev["missing_role"] = missing_role
    ev["metadata_only"] = True
    ev["no_pdf"] = True

    # Orient: preprint --preprint_of--> journal when possible
    present = dict(conn.execute("SELECT * FROM publications WHERE id=?", (present_id,)).fetchone())
    new_row = dict(conn.execute("SELECT * FROM publications WHERE id=?", (new_id,)).fetchone())
    if relation == "preprint_of":
        if is_preprint_pub(new_row) and is_journal_pub(present):
            from_id, to_id = new_id, present_id
        elif is_preprint_pub(present) and is_journal_pub(new_row):
            from_id, to_id = present_id, new_id
        elif missing_role == "from":
            from_id, to_id = new_id, present_id
        else:
            from_id, to_id = present_id, new_id
    else:
        from_id, to_id = present_id, new_id

    link = propose_version_link(
        conn,
        from_id,
        to_id,
        relation=relation,
        confidence=confidence,
        needs_confirm=1,
        evidence=ev,
    )
    return {
        "ok": True,
        "publication_id": new_id,
        "present_id": present_id,
        "missing_doi": doi,
        "created": created,
        "link": link,
        "metadata_only": True,
        "note": "Preview shell only — no PDF downloaded. Confirm link to merge.",
    }
