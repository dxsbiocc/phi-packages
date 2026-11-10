"""Reference formats for lit_export: BibTeX, RIS, CSL-JSON and CSV."""
from __future__ import annotations

import csv
import io
import json
import re
import unicodedata
from typing import Any, Iterable

FORMATS = ("bibtex", "ris", "csl-json", "csv")

_BIBTEX_ESCAPES = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_"}


def parse_authors(raw: Any) -> list[dict[str, str]]:
    """Return [{"family","given"}] from authors_json (a JSON list of names, or a plain list)."""
    if isinstance(raw, str):
        try:
            raw = json.loads(raw)
        except json.JSONDecodeError:
            raw = [a for a in re.split(r";|\band\b", raw) if a.strip()]
    out: list[dict[str, str]] = []
    for a in raw or []:
        if isinstance(a, dict):
            family = str(a.get("family") or a.get("last") or a.get("name") or "").strip()
            given = str(a.get("given") or a.get("first") or "").strip()
        else:
            name = str(a).strip()
            if not name:
                continue
            if "," in name:
                family, given = [x.strip() for x in name.split(",", 1)]
            else:
                parts = name.split()
                family, given = (parts[-1], " ".join(parts[:-1])) if len(parts) > 1 else (name, "")
        if family:
            out.append({"family": family, "given": given})
    return out


def _ascii(text: str) -> str:
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")


def cite_keys(pubs: Iterable[dict[str, Any]]) -> dict[str, str]:
    """Stable, unique keys like smith2020tcell, keyed by publication id."""
    used: dict[str, int] = {}
    keys: dict[str, str] = {}
    for p in pubs:
        authors = parse_authors(p.get("authors_json"))
        family = re.sub(r"[^a-z]", "", _ascii(authors[0]["family"]).lower()) if authors else "anon"
        word = next(
            (w for w in re.findall(r"[A-Za-z]{4,}", _ascii(p.get("title") or "")) if w.lower() not in {"with", "from", "that", "this", "their"}),
            "",
        ).lower()
        base = f"{family or 'anon'}{p.get('year') or 'nd'}{word}"
        n = used.get(base, 0)
        used[base] = n + 1
        keys[p["id"]] = base if n == 0 else f"{base}{chr(ord('a') + n - 1)}"
    return keys


def _bib(text: Any) -> str:
    return "".join(_BIBTEX_ESCAPES.get(c, c) for c in str(text or ""))


def _entry_type(p: dict[str, Any]) -> str:
    t = (p.get("type") or "").lower()
    doi = (p.get("doi") or "").lower()
    if "preprint" in t or doi.startswith("10.1101/") or p.get("arxiv"):
        return "misc"
    if "article" in t or p.get("venue"):
        return "article"
    return "misc"


def to_bibtex(pubs: list[dict[str, Any]]) -> str:
    keys = cite_keys(pubs)
    chunks = []
    for p in pubs:
        etype = _entry_type(p)
        authors = " and ".join(
            f"{a['family']}, {a['given']}" if a["given"] else a["family"] for a in parse_authors(p.get("authors_json"))
        )
        fields = [
            ("author", authors),
            ("title", "{" + _bib(p.get("title")) + "}" if p.get("title") else ""),
            ("journal" if etype == "article" else "howpublished", _bib(p.get("venue"))),
            ("year", str(p.get("year") or "")),
            ("doi", p.get("doi") or ""),
            ("pmid", p.get("pmid") or ""),
            ("eprint", p.get("arxiv") or ""),
            ("abstract", _bib(p.get("abstract"))),
        ]
        body = ",\n".join(f"  {k} = {{{v}}}" for k, v in fields if v)
        chunks.append(f"@{etype}{{{keys[p['id']]},\n{body}\n}}")
    return "\n\n".join(chunks) + ("\n" if chunks else "")


def to_ris(pubs: list[dict[str, Any]]) -> str:
    out = []
    for p in pubs:
        lines = [f"TY  - {'JOUR' if _entry_type(p) == 'article' else 'GEN'}"]
        lines += [f"AU  - {a['family']}, {a['given']}" if a["given"] else f"AU  - {a['family']}" for a in parse_authors(p.get("authors_json"))]
        for tag, value in (
            ("TI", p.get("title")), ("T2", p.get("venue")), ("PY", p.get("year")),
            ("DO", p.get("doi")), ("AN", f"PMID:{p['pmid']}" if p.get("pmid") else ""), ("AB", p.get("abstract")),
        ):
            if value:
                lines.append(f"{tag}  - {str(value).replace(chr(10), ' ')}")
        lines.append("ER  - ")
        out.append("\n".join(lines))
    return "\n\n".join(out) + ("\n" if out else "")


def to_csl_json(pubs: list[dict[str, Any]]) -> str:
    keys = cite_keys(pubs)
    items = []
    for p in pubs:
        item: dict[str, Any] = {
            "id": keys[p["id"]],
            "type": "article-journal" if _entry_type(p) == "article" else "article",
            "title": p.get("title") or "",
            "author": parse_authors(p.get("authors_json")),
        }
        if p.get("year"):
            item["issued"] = {"date-parts": [[int(p["year"])]]}
        for key, value in (("DOI", p.get("doi")), ("PMID", p.get("pmid")), ("container-title", p.get("venue")), ("abstract", p.get("abstract"))):
            if value:
                item[key] = value
        items.append(item)
    return json.dumps(items, ensure_ascii=False, indent=2) + "\n"


def to_csv(pubs: list[dict[str, Any]]) -> str:
    keys = cite_keys(pubs)
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["cite_key", "publication_id", "year", "title", "authors", "venue", "doi", "pmid", "tags", "needs_review"])
    for p in pubs:
        authors = "; ".join(f"{a['family']}, {a['given']}".strip(", ") for a in parse_authors(p.get("authors_json")))
        tags = p.get("tags_json") or ""
        w.writerow([keys[p["id"]], p["id"], p.get("year") or "", p.get("title") or "", authors,
                    p.get("venue") or "", p.get("doi") or "", p.get("pmid") or "", tags, p.get("needs_review") or 0])
    return buf.getvalue()


RENDERERS = {"bibtex": to_bibtex, "ris": to_ris, "csl-json": to_csl_json, "csv": to_csv}
EXTENSIONS = {"bibtex": ".bib", "ris": ".ris", "csl-json": ".json", "csv": ".csv"}
