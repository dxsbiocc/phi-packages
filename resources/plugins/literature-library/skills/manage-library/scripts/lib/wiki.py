#!/usr/bin/env python3
"""Wiki vault helpers for the literature library.

Vault layout (one per collection, next to library.sqlite):

  wiki-vault/
    AGENTS.md          conventions for whoever edits the vault
    vault.json         schema label
    raw/publications/<safe_id>/{meta.json, body.txt|<file>}   snapshots, never rewritten
    wiki/
      index.md         table of contents; headings are the theme tree
      log.md           append-only activity log
      sources/<safe_id>.md    one page per paper (summary, key claims, methods)
      concepts/<slug>.md      one page per concept, linking the papers that mention it

The prose in source pages is written by the invoking agent (its own model). The scripts
only build the request, validate the agent's JSON, and write files. No API keys, no
embeddings, no keyword-list extraction.
"""
from __future__ import annotations

import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

VAULT_NAME = "wiki-vault"
SCHEMA_LABEL = "agent-wiki-v1"
COMPILE_SCHEMA = "wiki-compile-v1"
SUMMARY_MIN_CHARS = 40
MAX_CONCEPTS_PER_PAPER = 12

PLACEHOLDER_DEFINITION = "_Not written yet. Write it once two or more papers mention this concept._"

AGENTS_MD = """# Literature wiki vault

Conventions for the wiki beside `library.sqlite` (identity, deduplication and search
history live in SQLite; this vault is the readable knowledge layer).

## Layers

1. `raw/` holds snapshots of what each page was written from. Add, never rewrite.
2. `wiki/` is markdown you may edit. Prefer appending and linking over rewriting.
3. This file holds the conventions. Do not delete it.

## Pages

| Path | Role |
|--|--|
| `wiki/index.md` | Table of contents. Headings are the theme tree; a paper may sit under several. |
| `wiki/log.md` | Append-only activity log. |
| `wiki/sources/<id>.md` | One page per paper: Summary, Key claims, Methods. |
| `wiki/concepts/<slug>.md` | One page per concept, listing the papers that mention it. |

## Rules

- Ground every statement in the paper's own title, abstract or provided text. Do not
  invent identifiers, sample sizes, p-values or results. Say so when something is unclear.
- Describe findings; do not give clinical advice.
- Link with relative markdown links, for example `[title](../sources/<id>.md)`.
- Papers marked `needs_review` carry a banner. Do not merge their claims silently.
- A concept's Definition starts as a placeholder. Write it yourself when two or more
  papers mention the concept, using only what those papers say.
"""

INDEX_MD = """# Wiki index

Table of contents for this collection. Edit headings and links freely.

## Unfiled

<!-- New source pages land here until you move them under a theme. -->

## Themes

<!-- Group papers under your own headings. -->

## Concepts

<!-- Concept pages are linked here after compile. -->
"""

LOG_MD = """# Wiki log

Append-only. Newest entries at the bottom.

"""


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def safe_id(publication_id: str) -> str:
    """Filesystem-safe id derived from the canonical publication id."""
    s = publication_id.strip()
    s = s.replace(":", "_").replace("/", "_").replace("\\", "_")
    s = re.sub(r"[^\w.\-]+", "_", s, flags=re.UNICODE)
    s = re.sub(r"_+", "_", s).strip("._")
    return s[:180] or "unknown"


def slugify(text: str) -> str:
    s = text.strip().lower().replace("+", "-plus-")
    s = re.sub(r"[^\w\-]+", "-", s, flags=re.UNICODE)
    s = re.sub(r"-+", "-", s).strip("-")
    return s[:80] or "untitled"


def vault_root(col_dir: Path) -> Path:
    return col_dir / VAULT_NAME


def init_vault(col_dir: Path, *, collection_id: str, collection_name: str = "") -> dict[str, Any]:
    """Create the vault skeleton. Idempotent; never overwrites existing pages."""
    root = vault_root(col_dir)
    created: list[str] = []
    for sub in (
        root,
        root / "raw" / "publications",
        root / "wiki" / "sources",
        root / "wiki" / "concepts",
    ):
        if not sub.exists():
            sub.mkdir(parents=True, exist_ok=True)
            created.append(str(sub.relative_to(root)) if sub != root else ".")

    files = {
        root / "AGENTS.md": AGENTS_MD,
        root / "wiki" / "index.md": INDEX_MD.replace(
            "# Wiki index", f"# Wiki index — {collection_name or collection_id}", 1
        ),
        root / "wiki" / "log.md": LOG_MD,
        root / "vault.json": json.dumps(
            {
                "collection_id": collection_id,
                "collection_name": collection_name,
                "created_at": utc_now(),
                "schema": SCHEMA_LABEL,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
    }
    for path, body in files.items():
        if not path.exists():
            path.write_text(body, encoding="utf-8")
            created.append(str(path.relative_to(root)))
    return {
        "ok": True,
        "vault_root": str(root),
        "created": created,
        "already_present": not created,
        "schema": SCHEMA_LABEL,
    }


def append_log(vault: Path, line: str) -> None:
    log = vault / "wiki" / "log.md"
    log.parent.mkdir(parents=True, exist_ok=True)
    if not log.exists():
        log.write_text(LOG_MD, encoding="utf-8")
    with log.open("a", encoding="utf-8") as f:
        f.write(f"\n- `{utc_now()}` {line}\n")


def _md(text: Any) -> str:
    return str(text or "").replace("\r\n", "\n").replace("\r", "\n")


def parse_tags(raw_tags: Any) -> list[str]:
    if raw_tags is None:
        return []
    if isinstance(raw_tags, list):
        return [str(t).strip() for t in raw_tags if str(t).strip()]
    if isinstance(raw_tags, str):
        s = raw_tags.strip()
        if not s:
            return []
        try:
            data = json.loads(s)
        except json.JSONDecodeError:
            return [t.strip() for t in s.split(",") if t.strip()]
        return parse_tags(data) if isinstance(data, list) else []
    return []


# ---------------------------------------------------------------------------
# Compile payload (what the agent writes per paper)
# ---------------------------------------------------------------------------

RESPONSE_SCHEMA: dict[str, Any] = {
    "type": "object",
    "required": ["publication_id", "summary", "key_claims"],
    "properties": {
        "publication_id": {"type": "string", "description": "copy from the request item"},
        "summary": {"type": "string", "minLength": SUMMARY_MIN_CHARS,
                    "description": "2-4 sentences of grounded synthesis, not a copy of the abstract"},
        "key_claims": {"type": "array", "items": {"type": "string"}, "minItems": 1,
                       "description": "3-6 claims the text supports"},
        "methods": {"type": "array", "items": {"type": "string"}},
        "concepts": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["name"],
                "properties": {"name": {"type": "string"}, "evidence": {"type": "string"}},
            },
            "description": "topics, diseases, methods, genes or proteins the paper is about",
        },
        "tags": {"type": "array", "items": {"type": "string"}},
        "confidence_notes": {"type": "string"},
    },
}

INSTRUCTIONS = (
    "For each item, write one result object matching response_schema, using only that item's "
    "title, abstract and fulltext_excerpt. Do not invent identifiers, numbers or findings; leave "
    "a field out or say so in confidence_notes when the text does not support it. Keep "
    "publication_id exactly as given. Save {\"results\": [...]} to response_path, then call "
    "lit_wiki_apply."
)


def _need_str(v: Any, field: str) -> str:
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f"{field} must be a non-empty string")
    return v.strip()


def _str_list(v: Any, field: str, *, min_items: int = 0) -> list[str]:
    if v is None:
        v = []
    if not isinstance(v, list):
        raise ValueError(f"{field} must be a list")
    out = []
    for i, item in enumerate(v):
        if not isinstance(item, str) or not item.strip():
            raise ValueError(f"{field}[{i}] must be a non-empty string")
        out.append(item.strip())
    if len(out) < min_items:
        raise ValueError(f"{field} needs at least {min_items} item(s)")
    return out


def validate_compile_payload(obj: Any) -> dict[str, Any]:
    """Validate one agent result. Raises ValueError with a message the agent can act on."""
    if not isinstance(obj, dict):
        raise ValueError("result must be an object")
    summary = _need_str(obj.get("summary"), "summary")
    if len(summary) < SUMMARY_MIN_CHARS:
        raise ValueError(f"summary is shorter than {SUMMARY_MIN_CHARS} characters")
    concepts: list[dict[str, str]] = []
    # "entities" is accepted from older responses and treated as concepts.
    for field in ("concepts", "entities"):
        raw = obj.get(field) or []
        if not isinstance(raw, list):
            raise ValueError(f"{field} must be a list")
        for i, c in enumerate(raw):
            if isinstance(c, str):
                name, evidence = c.strip(), ""
            elif isinstance(c, dict):
                name = str(c.get("name") or c.get("symbol") or "").strip()
                evidence = str(c.get("evidence") or c.get("note") or "").strip()
            else:
                raise ValueError(f"{field}[{i}] must be a string or object")
            if not name:
                raise ValueError(f"{field}[{i}] has no name")
            if all(name.lower() != x["name"].lower() for x in concepts):
                concepts.append({"name": name, "evidence": evidence})
    notes = obj.get("confidence_notes")
    if notes is not None and not isinstance(notes, str):
        raise ValueError("confidence_notes must be a string")
    return {
        "publication_id": str(obj.get("publication_id") or "").strip(),
        "summary": summary,
        "key_claims": _str_list(obj.get("key_claims"), "key_claims", min_items=1),
        "methods": _str_list(obj.get("methods"), "methods"),
        "concepts": concepts[:MAX_CONCEPTS_PER_PAPER],
        "tags": _str_list(obj.get("tags"), "tags"),
        "confidence_notes": (notes or "").strip(),
    }


def request_item(pub: dict[str, Any], *, fulltext_excerpt: str = "") -> dict[str, Any]:
    """One paper's entry in a compile batch request."""
    title = (pub.get("title") or "").strip()
    abstract = (pub.get("abstract") or "").strip()
    item: dict[str, Any] = {
        "publication_id": str(pub.get("id") or ""),
        "title": title,
        "year": pub.get("year"),
        "venue": pub.get("venue") or "",
        "doi": pub.get("doi") or "",
        "pmid": pub.get("pmid") or "",
        "abstract": abstract,
        "basis": "fulltext" if fulltext_excerpt else "abstract",
    }
    if fulltext_excerpt:
        item["fulltext_excerpt"] = fulltext_excerpt
    if int(pub.get("needs_review") or 0):
        item["needs_review"] = True
    if not abstract and not fulltext_excerpt:
        item["warning"] = "no abstract: run lit_curate enrich first, or skip this paper"
    return item


# ---------------------------------------------------------------------------
# Writing pages
# ---------------------------------------------------------------------------


def source_page_markdown(
    pub: dict[str, Any],
    payload: dict[str, Any],
    *,
    sid: str,
    rel_raw: str,
    basis: str,
) -> str:
    title = pub.get("title") or pub.get("id") or sid
    needs = int(pub.get("needs_review") or 0)
    banner = (
        "\n> **needs_review**: identity is uncertain (soft match or conflicting IDs). "
        "Do not merge this paper's claims silently.\n"
        if needs
        else ""
    )
    ids = [
        f"- {label}: `{pub[key]}`"
        for label, key in (("DOI", "doi"), ("PMID", "pmid"), ("PMCID", "pmcid"), ("arXiv", "arxiv"))
        if pub.get(key)
    ]
    claims = "\n".join(f"- {_md(c)}" for c in payload["key_claims"])
    methods = ", ".join(f"`{m}`" for m in payload["methods"]) or "_None stated._"
    concepts = "\n".join(
        f"- [{c['name']}](../concepts/{slugify(c['name'])}.md)" for c in payload["concepts"]
    ) or "_None._"
    notes = f"\n## Confidence notes\n\n{_md(payload['confidence_notes'])}\n" if payload["confidence_notes"] else ""
    return f"""# {_md(title)}
{banner}
- publication_id: `{pub.get('id')}`
- year: {pub.get('year') or ''}
- venue: {pub.get('venue') or ''}
- needs_review: {needs}
- written from: {basis}
- raw: [`{rel_raw}`](../../raw/{rel_raw})

## Identifiers

{chr(10).join(ids) if ids else '_No strong IDs._'}

## Summary

{_md(payload['summary'])}

## Key claims

{claims}

## Methods

{methods}

## Abstract

{_md(pub.get('abstract')) or '_No abstract saved._'}

## Concepts

{concepts}
{notes}
## Related sources

<!-- Link other wiki/sources pages here. -->
"""


def ensure_section_link(index_path: Path, section: str, link_text: str, rel_target: str) -> bool:
    """Add a bullet under `## <section>` unless the target is already linked."""
    text = index_path.read_text(encoding="utf-8") if index_path.exists() else INDEX_MD
    if f"]({rel_target})" in text or f"](./{rel_target})" in text:
        return False
    bullet = f"- [{link_text}]({rel_target})"
    marker = f"## {section}"
    if marker not in text:
        index_path.write_text(text.rstrip() + f"\n\n{marker}\n\n{bullet}\n", encoding="utf-8")
        return True
    head, rest = text.split(marker, 1)
    lines = rest.splitlines(keepends=True)
    i = 0
    while i < len(lines) and (lines[i].strip() == "" or lines[i].lstrip().startswith("<!--")):
        i += 1
    new_rest = "".join(lines[:i])
    if not new_rest.endswith("\n"):
        new_rest += "\n"
    index_path.write_text(
        head + marker + new_rest + bullet + "\n" + "".join(lines[i:]), encoding="utf-8"
    )
    return True


def _upsert_concept(path: Path, *, name: str, slug: str, source_rel: str, source_title: str, evidence: str) -> bool:
    """Create the concept page or add this paper to its Sources. True when created."""
    bullet = f"- [{source_title}]({source_rel})" + (f" — {evidence[:240]}" if evidence else "")
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            f"# {name}\n\nSlug: `{slug}`\n\n## Definition\n\n{PLACEHOLDER_DEFINITION}\n\n"
            f"## Sources\n\n{bullet}\n\n## Related concepts\n\n",
            encoding="utf-8",
        )
        return True
    text = path.read_text(encoding="utf-8")
    if f"]({source_rel})" in text:
        return False
    marker = "## Sources\n"
    if marker in text:
        head, _, rest = text.partition(marker)
        path.write_text(head + marker + "\n" + bullet + "\n" + rest.lstrip("\n"), encoding="utf-8")
    else:
        path.write_text(text.rstrip() + f"\n\n## Sources\n\n{bullet}\n", encoding="utf-8")
    return False


def write_paper(
    vault: Path,
    pub: dict[str, Any],
    payload: dict[str, Any],
    *,
    raw_text: Optional[str] = None,
    basis: str = "abstract",
) -> dict[str, Any]:
    """Write raw snapshot, source page, concept pages, index links and a log line."""
    pub_id = str(pub.get("id") or "")
    if not pub_id:
        return {"ok": False, "error": "publication missing id"}
    sid = safe_id(pub_id)
    title = str(pub.get("title") or pub_id)

    raw_dir = vault / "raw" / "publications" / sid
    raw_dir.mkdir(parents=True, exist_ok=True)
    snapshot = raw_text
    if snapshot is None:
        snapshot = f"Title: {title}\n\nAbstract:\n{pub.get('abstract') or ''}\n"
    (raw_dir / "body.txt").write_text(snapshot, encoding="utf-8")
    (raw_dir / "meta.json").write_text(
        json.dumps(
            {
                "publication_id": pub_id,
                "safe_id": sid,
                "written_at": utc_now(),
                "title": pub.get("title"),
                "year": pub.get("year"),
                "doi": pub.get("doi"),
                "pmid": pub.get("pmid"),
                "needs_review": int(pub.get("needs_review") or 0),
                "basis": basis,
                "schema": COMPILE_SCHEMA,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    page = vault / "wiki" / "sources" / f"{sid}.md"
    page.write_text(
        source_page_markdown(pub, payload, sid=sid, rel_raw=f"publications/{sid}/body.txt", basis=basis),
        encoding="utf-8",
    )

    short = title if len(title) <= 80 else title[:77] + "…"
    created: list[str] = []
    updated: list[str] = []
    for c in payload["concepts"]:
        slug = slugify(c["name"])
        cpath = vault / "wiki" / "concepts" / f"{slug}.md"
        was_new = _upsert_concept(
            cpath, name=c["name"], slug=slug, source_rel=f"../sources/{sid}.md",
            source_title=short, evidence=c.get("evidence", ""),
        )
        (created if was_new else updated).append(str(cpath.relative_to(vault)))

    index = vault / "wiki" / "index.md"
    linked = ensure_section_link(index, "Unfiled", short, f"sources/{sid}.md")
    for c in payload["concepts"][:6]:
        ensure_section_link(index, "Concepts", c["name"], f"concepts/{slugify(c['name'])}.md")

    append_log(
        vault,
        f"wrote source `{sid}` (`{pub_id}`) from {basis}; claims={len(payload['key_claims'])} "
        f"concepts=+{len(created)}/~{len(updated)}",
    )
    return {
        "ok": True,
        "publication_id": pub_id,
        "source_page": str(page.relative_to(vault)),
        "concepts_created": created,
        "concepts_updated": updated,
        "index_linked": linked,
        "suggested_tags": payload["tags"],
    }


def concepts_needing_definition(vault: Path) -> list[dict[str, Any]]:
    """Concept pages that list two or more papers but still have the placeholder Definition."""
    out: list[dict[str, Any]] = []
    cdir = vault / "wiki" / "concepts"
    if not cdir.is_dir():
        return out
    for p in sorted(cdir.glob("*.md")):
        text = p.read_text(encoding="utf-8")
        if PLACEHOLDER_DEFINITION not in text:
            continue
        body = text.split("## Sources", 1)[-1].split("\n## ", 1)[0]
        n = sum(1 for ln in body.splitlines() if ln.strip().startswith("- ["))
        if n >= 2:
            out.append({"slug": p.stem, "path": str(p.relative_to(vault)), "n_sources": n})
    return out


def sources_written(vault: Path) -> set[str]:
    sdir = vault / "wiki" / "sources"
    return {p.stem for p in sdir.glob("*.md")} if sdir.is_dir() else set()


def lint_vault(vault: Path) -> dict[str, Any]:
    """Missing index entries, broken links, orphans, stale banners, raw without page."""
    issues: list[dict[str, Any]] = []
    wiki = vault / "wiki"
    index_path = wiki / "index.md"
    index_text = index_path.read_text(encoding="utf-8") if index_path.exists() else ""
    source_files = sorted((wiki / "sources").glob("*.md")) if (wiki / "sources").exists() else []
    source_names = {p.name for p in source_files}

    for p in source_files:
        rel = f"sources/{p.name}"
        if rel not in index_text and f"./{rel}" not in index_text:
            issues.append({"code": "missing_index_entry", "path": rel, "severity": "error"})

    inbound: set[str] = set()
    vault_resolved = vault.resolve()
    for md in (list(wiki.rglob("*.md")) if wiki.exists() else []):
        for link in extract_md_links(md.read_text(encoding="utf-8")):
            if link.startswith(("http://", "https://", "mailto:")):
                continue
            link_path = link.split("#", 1)[0]
            if not link_path:
                continue
            target = (md.parent / link_path).resolve()
            try:
                target.relative_to(vault_resolved)
            except ValueError:
                issues.append({"code": "link_outside_vault", "from": str(md.relative_to(vault)),
                               "link": link, "severity": "error"})
                continue
            if not target.exists():
                issues.append({"code": "broken_link", "from": str(md.relative_to(vault)),
                               "link": link, "severity": "error"})
            else:
                inbound.add(str(target.relative_to(vault_resolved)))

    for kind in ("concepts", "entities"):  # entities/ may exist in vaults made by older versions
        d = wiki / kind
        if d.exists():
            for p in d.glob("*.md"):
                rel = str(p.relative_to(vault))
                if rel not in inbound:
                    issues.append({"code": "orphan_page", "path": rel, "severity": "warn"})

    for p in source_files:
        meta_path = vault / "raw" / "publications" / p.stem / "meta.json"
        if not meta_path.exists():
            continue
        try:
            needs = int(json.loads(meta_path.read_text(encoding="utf-8")).get("needs_review") or 0)
        except (json.JSONDecodeError, ValueError):
            issues.append({"code": "bad_raw_meta", "path": str(meta_path.relative_to(vault)), "severity": "warn"})
            continue
        if needs and "needs_review" not in p.read_text(encoding="utf-8"):
            issues.append({"code": "stale_needs_review_banner", "path": str(p.relative_to(vault)),
                           "severity": "warn"})

    raw_pubs = vault / "raw" / "publications"
    if raw_pubs.exists():
        for d in raw_pubs.iterdir():
            if d.is_dir() and f"{d.name}.md" not in source_names:
                issues.append({"code": "raw_without_source_page", "path": str(d.relative_to(vault)),
                               "severity": "error"})

    for c in concepts_needing_definition(vault):
        issues.append({"code": "concept_needs_definition", "path": c["path"],
                       "n_sources": c["n_sources"], "severity": "info"})

    count = lambda sev: sum(1 for i in issues if i.get("severity") == sev)  # noqa: E731
    return {
        "ok": True,
        "n_issues": len(issues),
        "n_errors": count("error"),
        "n_warns": count("warn"),
        "n_info": count("info"),
        "n_sources": len(source_files),
        "issues": issues,
    }




def query_wiki_markdown(vault: Path, q: str, *, limit: int = 50) -> list[dict[str, Any]]:
    """Keyword search over wiki markdown (case-insensitive). No vectors."""
    q = (q or "").strip()
    if not q:
        return []
    limit = max(1, min(int(limit), 200))
    wiki = vault / "wiki"
    if not wiki.exists():
        return []
    terms = [t for t in re.split(r"\s+", q.lower()) if t]
    hits: list[dict[str, Any]] = []
    for md in sorted(wiki.rglob("*.md")):
        text = md.read_text(encoding="utf-8")
        low = text.lower()
        if not all(t in low for t in terms):
            continue
        # snippet: first line with a term
        snippet = ""
        for line in text.splitlines():
            if any(t in line.lower() for t in terms):
                snippet = line.strip()[:240]
                break
        if not snippet:
            snippet = text.strip().splitlines()[0][:240] if text.strip() else ""
        hits.append(
            {
                "path": str(md.relative_to(vault)),
                "title": next((ln[2:].strip() for ln in text.splitlines() if ln.startswith("# ")), md.stem),
                "snippet": snippet,
            }
        )
        if len(hits) >= limit:
            break
    return hits


def extract_md_links(text: str) -> list[str]:
    # Ignore markdown inside HTML comments (template placeholders)
    cleaned = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    return re.findall(r"\[[^\]]*\]\(([^)]+)\)", cleaned)


def _replace_publication_id_in_source_md(text: str, old_id: str, new_id: str, old_safe: str, new_safe: str) -> str:
    """Update publication_id bullet and safe_id path references inside a source page."""
    out = text
    out = out.replace(f"- publication_id: `{old_id}`", f"- publication_id: `{new_id}`")
    out = out.replace(f"publications/{old_safe}/", f"publications/{new_safe}/")
    out = out.replace(f"sources/{old_safe}.md", f"sources/{new_safe}.md")
    out = out.replace(f"`{old_safe}`", f"`{new_safe}`")
    return out


def _replace_links_in_text(text: str, old_safe: str, new_safe: str) -> str:
    if old_safe == new_safe:
        return text
    return (
        text.replace(f"sources/{old_safe}.md", f"sources/{new_safe}.md")
        .replace(f"]({old_safe}.md)", f"]({new_safe}.md)")
        .replace(f"/{old_safe}.md", f"/{new_safe}.md")
        .replace(f"publications/{old_safe}/", f"publications/{new_safe}/")
    )


def wiki_repoint_publication(
    project_root: Path,
    old_id: str,
    new_id: str,
    *,
    collection_id: Optional[str] = None,
    db_path: Optional[Path] = None,
    leave_stub: bool = True,
) -> dict[str, Any]:
    """Rename wiki source/raw paths when a publication canonical id changes.

    For each matching collection vault:
      - move wiki/sources/<old_safe>.md → <new_safe>.md (update publication_id)
      - optionally leave a stub redirect at the old path
      - rename raw/publications/<old_safe>/ when present
      - cheap string replace in index.md / log.md
    """
    root = Path(project_root).resolve()
    if db_path is None:
        db_path = root / "library.sqlite"
    old_safe = safe_id(old_id)
    new_safe = safe_id(new_id)
    if old_id == new_id or old_safe == new_safe:
        return {
            "ok": True,
            "skipped": True,
            "reason": "same_id",
            "old_id": old_id,
            "new_id": new_id,
            "vaults": [],
        }

    col_dirs: list[tuple[str, Path]] = []
    import sqlite3

    if db_path.exists():
        conn = sqlite3.connect(str(db_path))
        conn.row_factory = sqlite3.Row
        try:
            if collection_id:
                rows = conn.execute(
                    "SELECT id, path FROM collections WHERE id=?", (collection_id,)
                ).fetchall()
            else:
                rows = conn.execute("SELECT id, path FROM collections").fetchall()
            for r in rows:
                p = Path(r["path"]) if r["path"] else root / "collections" / r["id"]
                if not p.is_absolute():
                    p = root / p
                col_dirs.append((r["id"], p))
        finally:
            conn.close()
    elif collection_id:
        col_dirs.append((collection_id, root / "collections" / collection_id))

    vault_results: list[dict[str, Any]] = []
    for cid, col_dir in col_dirs:
        vault = vault_root(col_dir)
        if not vault.is_dir():
            continue
        sources_dir = vault / "wiki" / "sources"
        old_page = sources_dir / f"{old_safe}.md"
        new_page = sources_dir / f"{new_safe}.md"
        raw_old = vault / "raw" / "publications" / old_safe
        raw_new = vault / "raw" / "publications" / new_safe
        actions: list[str] = []

        if old_page.is_file():
            body = old_page.read_text(encoding="utf-8")
            body = _replace_publication_id_in_source_md(body, old_id, new_id, old_safe, new_safe)
            sources_dir.mkdir(parents=True, exist_ok=True)
            if new_page.exists() and new_page.resolve() != old_page.resolve():
                actions.append("new_source_exists")
            else:
                new_page.write_text(body, encoding="utf-8")
                actions.append("wrote_new_source")
            if leave_stub:
                stub = (
                    "# Redirect\n\n"
                    "This source moved after identity upgrade/merge.\n\n"
                    f"- old publication_id: `{old_id}`\n"
                    f"- new publication_id: `{new_id}`\n"
                    f"- see: [`{new_safe}.md`]({new_safe}.md)\n"
                )
                old_page.write_text(stub, encoding="utf-8")
                actions.append("left_stub")
            else:
                if old_page.resolve() != new_page.resolve():
                    old_page.unlink()
                    actions.append("removed_old_source")
        else:
            actions.append("no_old_source")

        if raw_old.is_dir():
            raw_new.parent.mkdir(parents=True, exist_ok=True)
            if raw_new.exists():
                for child in raw_old.iterdir():
                    dest = raw_new / child.name
                    if not dest.exists():
                        shutil.move(str(child), str(dest))
                shutil.rmtree(raw_old, ignore_errors=True)
                actions.append("merged_raw")
            else:
                shutil.move(str(raw_old), str(raw_new))
                actions.append("moved_raw")
            meta = raw_new / "meta.json"
            if meta.is_file():
                try:
                    data = json.loads(meta.read_text(encoding="utf-8"))
                    data["publication_id"] = new_id
                    data["safe_id"] = new_safe
                    data["repointed_from"] = old_id
                    meta.write_text(
                        json.dumps(data, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8",
                    )
                    actions.append("updated_raw_meta")
                except (json.JSONDecodeError, OSError):
                    actions.append("raw_meta_update_failed")

        for name in ("index.md", "log.md"):
            p = vault / "wiki" / name
            if not p.is_file():
                continue
            before = p.read_text(encoding="utf-8")
            after = _replace_links_in_text(before, old_safe, new_safe)
            if after != before:
                p.write_text(after, encoding="utf-8")
                actions.append(f"updated_{name}")

        if any(a for a in actions if a not in ("no_old_source",)):
            append_log(
                vault,
                f"wiki_repoint {old_id} → {new_id} ({old_safe} → {new_safe})",
            )
            actions.append("logged")

        vault_results.append(
            {
                "collection_id": cid,
                "vault": str(vault),
                "old_safe": old_safe,
                "new_safe": new_safe,
                "actions": actions,
            }
        )

    return {
        "ok": True,
        "old_id": old_id,
        "new_id": new_id,
        "old_safe": old_safe,
        "new_safe": new_safe,
        "n_vaults": len(vault_results),
        "vaults": vault_results,
    }


def wiki_repoint_remaps(
    project_root: Path,
    remaps: list[dict[str, str]],
    *,
    collection_id: Optional[str] = None,
    db_path: Optional[Path] = None,
    leave_stub: bool = True,
) -> dict[str, Any]:
    """Apply wiki_repoint_publication for a list of {old_id,new_id} remaps."""
    results = []
    for r in remaps or []:
        old_id = r.get("old_id") or r.get("from") or ""
        new_id = r.get("new_id") or r.get("to") or ""
        if not old_id or not new_id:
            continue
        results.append(
            wiki_repoint_publication(
                project_root,
                old_id,
                new_id,
                collection_id=collection_id,
                db_path=db_path,
                leave_stub=leave_stub,
            )
        )
    return {
        "ok": True,
        "n": len(results),
        "results": results,
    }
