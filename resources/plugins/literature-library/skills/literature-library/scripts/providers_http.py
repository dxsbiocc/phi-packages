#!/usr/bin/env python3
"""HTTP providers for lit-search (draft CLI path; no MCP server in-plugin).

Public REST only. Optional env (never packaged keys):
  OPENALEX_MAILTO / CROSSREF_MAILTO / CONTACT_EMAIL — polite mailto for OpenAlex/Crossref
  S2_API_KEY / SEMANTIC_SCHOLAR_API_KEY — Semantic Scholar (soft-fail if missing/unauthorized)

Sources:
  pubmed, europepmc — live (existing)
  openalex, crossref, arxiv — live keyword search
  biorxiv / medrxiv — DOI detail, date browse, or Europe PMC JOURNAL keyword fallback
  pmc — NCBI E-utilities db=pmc (esearch + esummary)
  doaj — DOAJ public article search API
  datacite — DataCite REST dois?query=
  zenodo — Zenodo/Invenio records API
  semantic_scholar / S2 — live when key accepted; meta error otherwise

Expand helpers (lit-expand; mirror phi-literature related_publications over HTTP):
  pubmed elink neighbors, OpenAlex related/referenced/cites, OpenCitations DOI fallback.

No MCP server in-plugin; stdlib urllib/ssl/json/xml only.
"""
from __future__ import annotations

import json
import os
import re
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date, timedelta
from typing import Any, Optional
from xml.etree import ElementTree as ET

USER_AGENT = "phi-literature-library-draft/0.4.9 (local; research; +https://github.com/)"
TIMEOUT = 30

ATOM_NS = {
    "atom": "http://www.w3.org/2005/Atom",
    "arxiv": "http://arxiv.org/schemas/atom",
    "opensearch": "http://a9.com/-/spec/opensearch/1.1/",
}

S2_FIELDS = (
    "paperId,title,abstract,authors,externalIds,url,openAccessPdf,"
    "publicationDate,venue,journal,publicationTypes"
)


def _contact_email() -> Optional[str]:
    for key in ("OPENALEX_MAILTO", "CROSSREF_MAILTO", "CONTACT_EMAIL"):
        v = (os.environ.get(key) or "").strip()
        if v:
            return v
    return None


def _get(url: str, *, headers: Optional[dict[str, str]] = None) -> bytes:
    hdrs = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    if headers:
        hdrs.update(headers)
    req = urllib.request.Request(url, headers=hdrs)
    ctx = ssl.create_default_context()
    with urllib.request.urlopen(req, timeout=TIMEOUT, context=ctx) as resp:
        return resp.read()


def _get_json(url: str, *, headers: Optional[dict[str, str]] = None) -> Any:
    return json.loads(_get(url, headers=headers).decode("utf-8"))


def _strip_tags(s: Optional[str]) -> str:
    if not s:
        return ""
    return re.sub(r"<[^>]+>", " ", s).replace("\n", " ").strip()


def _doi_norm(raw: Optional[str]) -> Optional[str]:
    if not raw:
        return None
    s = str(raw).strip()
    s = re.sub(r"^https?://(dx\.)?doi\.org/", "", s, flags=re.I)
    m = re.search(r"(10\.\d{4,9}/\S+)", s, re.I)
    return m.group(1).rstrip(").,;") if m else None


def _year_from(value: Any) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, int):
        return str(value)
    s = str(value).strip()
    if not s:
        return None
    m = re.match(r"(\d{4})", s)
    return m.group(1) if m else None


def _inverted_abstract(idx: Any) -> Optional[str]:
    if not isinstance(idx, dict) or not idx:
        return None
    pairs: list[tuple[int, str]] = []
    for word, positions in idx.items():
        if not isinstance(positions, list):
            continue
        for pos in positions:
            if isinstance(pos, int):
                pairs.append((pos, str(word)))
    if not pairs:
        return None
    pairs.sort(key=lambda x: x[0])
    return " ".join(w for _, w in pairs)


def _map_type(value: Optional[str]) -> str:
    if not value:
        return "journal-article"
    token = value.lower().replace("_", "-")
    if "preprint" in token or token == "posted-content":
        return "preprint"
    if "review" in token:
        return "review"
    if "book-chapter" in token or token == "chapter":
        return "book-chapter"
    if "proceedings" in token:
        return "proceedings-article"
    if token == "article" or "journal" in token:
        return "journal-article"
    return "other"


# ---------------------------------------------------------------------------
# PubMed / Europe PMC (existing)
# ---------------------------------------------------------------------------


def search_pubmed(
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
) -> list[dict[str, Any]]:
    term = query
    if date_from or date_to:
        df = (date_from or "1900")[:4]
        dt = (date_to or "2099")[:4]
        term = f"({query}) AND ({df}:{dt}[pdat])"
    esearch = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode(
            {"db": "pubmed", "retmode": "json", "retmax": str(limit), "term": term, "sort": "pub+date"}
        )
    )
    data = json.loads(_get(esearch).decode("utf-8"))
    ids = data.get("esearchresult", {}).get("idlist") or []
    if not ids:
        return []
    efetch = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        + urllib.parse.urlencode({"db": "pubmed", "retmode": "xml", "id": ",".join(ids)})
    )
    return parse_pubmed_xml(_get(efetch))


def search_europepmc(
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
) -> list[dict[str, Any]]:
    q = query
    if date_from or date_to:
        df = date_from or "1900-01-01"
        dt = date_to or "2099-12-31"
        q = f"{query} AND (FIRST_PDATE:[{df} TO {dt}])"
    url = (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
        + urllib.parse.urlencode(
            {"query": q, "format": "json", "pageSize": str(limit), "resultType": "core"}
        )
    )
    data = json.loads(_get(url).decode("utf-8"))
    results = (data.get("resultList") or {}).get("result") or []
    hits: list[dict[str, Any]] = []
    for r in results:
        authors = []
        al = r.get("authorString") or ""
        if al:
            authors = [a.strip() for a in al.split(",") if a.strip()]
        hits.append(
            {
                "source": "europepmc",
                "pmid": r.get("pmid"),
                "pmcid": r.get("pmcid"),
                "doi": r.get("doi"),
                "title": r.get("title"),
                "abstract": r.get("abstractText"),
                "authors": authors,
                "year": r.get("pubYear"),
                "venue": r.get("journalTitle"),
                "type": r.get("pubType") or "journal-article",
                "url": (
                    f"https://europepmc.org/article/MED/{r['pmid']}"
                    if r.get("pmid")
                    else None
                ),
            }
        )
    return hits


# ---------------------------------------------------------------------------
# OpenAlex
# ---------------------------------------------------------------------------


def _openalex_url(path: str, params: dict[str, str]) -> str:
    url = f"https://api.openalex.org/{path}?" + urllib.parse.urlencode(params)
    email = _contact_email()
    if email:
        sep = "&" if "?" in url else "?"
        url = f"{url}{sep}mailto={urllib.parse.quote(email)}"
    return url


def _hit_from_openalex(work: dict[str, Any]) -> Optional[dict[str, Any]]:
    if work.get("error") or (isinstance(work.get("message"), str) and "not found" in work["message"]):
        return None
    ids = work.get("ids") or {}
    doi = _doi_norm(work.get("doi") or ids.get("doi"))
    pmid = None
    pmcid = None
    if ids.get("pmid"):
        m = re.search(r"(\d{1,9})$", str(ids["pmid"]))
        pmid = m.group(1) if m else None
    if ids.get("pmcid"):
        m = re.search(r"(PMC\d+)", str(ids["pmcid"]), re.I)
        pmcid = m.group(1).upper() if m else None
    authors: list[str] = []
    for row in work.get("authorships") or []:
        author = (row or {}).get("author") or {}
        name = author.get("display_name")
        if name:
            authors.append(name)
    primary = work.get("primary_location") or {}
    source = primary.get("source") or {}
    venue = source.get("display_name")
    year = work.get("publication_year")
    if year is not None:
        year = str(year)
    abstract = _inverted_abstract(work.get("abstract_inverted_index"))
    landing = primary.get("landing_page_url")
    oa_id = None
    for cand in (ids.get("openalex"), work.get("id")):
        if cand:
            m = re.search(r"(W\d+)", str(cand), re.I)
            if m:
                oa_id = m.group(1).upper()
                break
    return {
        "source": "openalex",
        "doi": doi,
        "pmid": pmid,
        "pmcid": pmcid,
        "openalex": oa_id,
        "title": work.get("display_name"),
        "abstract": abstract,
        "authors": authors,
        "year": year,
        "venue": venue,
        "type": _map_type(work.get("type")),
        "url": landing or (f"https://doi.org/{doi}" if doi else (f"https://openalex.org/{oa_id}" if oa_id else None)),
    }


def search_openalex(
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
) -> list[dict[str, Any]]:
    filters: list[str] = []
    if date_from:
        df = date_from if len(date_from) > 4 else f"{date_from}-01-01"
        filters.append(f"from_publication_date:{df}")
    if date_to:
        dt = date_to if len(date_to) > 4 else f"{date_to}-12-31"
        filters.append(f"to_publication_date:{dt}")
    params: dict[str, str] = {
        "search": query,
        "per-page": str(min(limit, 50)),
        "page": "1",
        "sort": "relevance_score:desc",
    }
    if filters:
        params["filter"] = ",".join(filters)
    url = _openalex_url("works", params)
    body = None
    last_err: Optional[Exception] = None
    for attempt in range(3):
        try:
            body = _get_json(url)
            break
        except urllib.error.HTTPError as e:
            last_err = e
            if e.code == 429 and attempt < 2:
                time.sleep(1.5 * (attempt + 1))
                continue
            raise
    if body is None:
        raise last_err or RuntimeError("openalex empty body")
    hits: list[dict[str, Any]] = []
    for row in body.get("results") or []:
        h = _hit_from_openalex(row if isinstance(row, dict) else {})
        if h and h.get("title"):
            hits.append(h)
        if len(hits) >= limit:
            break
    return hits


# ---------------------------------------------------------------------------
# Crossref
# ---------------------------------------------------------------------------


def _crossref_url(path: str, params: Optional[dict[str, str]] = None) -> str:
    base = f"https://api.crossref.org/{path}"
    if params:
        base = base + "?" + urllib.parse.urlencode(params)
    email = _contact_email()
    if email:
        sep = "&" if "?" in base else "?"
        base = f"{base}{sep}mailto={urllib.parse.quote(email)}"
    return base


def _parts_year(parts: Any) -> Optional[str]:
    if not isinstance(parts, list) or not parts:
        return None
    first = parts[0]
    if isinstance(first, list) and first and isinstance(first[0], int):
        return str(first[0])
    return None


def _hit_from_crossref(work: dict[str, Any]) -> Optional[dict[str, Any]]:
    doi = _doi_norm(work.get("DOI"))
    if not doi:
        return None
    titles = work.get("title") or []
    title = _strip_tags(titles[0] if titles else "")
    authors: list[str] = []
    for row in work.get("author") or []:
        name = row.get("name")
        if not name:
            given = row.get("given") or ""
            family = row.get("family") or ""
            name = f"{family}, {given}".strip(", ") if family else given
        if name:
            authors.append(name)
    containers = work.get("container-title") or []
    venue = containers[0] if containers else None
    year = _parts_year((work.get("issued") or {}).get("date-parts"))
    abstract = _strip_tags(work.get("abstract"))
    return {
        "source": "crossref",
        "doi": doi,
        "pmid": None,
        "pmcid": None,
        "title": title,
        "abstract": abstract or None,
        "authors": authors,
        "year": year,
        "venue": venue,
        "type": _map_type(work.get("type")),
        "url": f"https://doi.org/{doi}",
    }


def search_crossref(
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
) -> list[dict[str, Any]]:
    filters: list[str] = []
    if date_from:
        filters.append(f"from-pub-date:{date_from}")
    if date_to:
        filters.append(f"until-pub-date:{date_to}")
    params: dict[str, str] = {
        "query.bibliographic": query,
        "rows": str(min(limit, 50)),
        "offset": "0",
        "sort": "relevance",
        "order": "desc",
    }
    if filters:
        params["filter"] = ",".join(filters)
    body = _get_json(_crossref_url("works", params))
    message = body.get("message") or {}
    hits: list[dict[str, Any]] = []
    for row in message.get("items") or []:
        h = _hit_from_crossref(row if isinstance(row, dict) else {})
        if h and h.get("title"):
            hits.append(h)
        if len(hits) >= limit:
            break
    return hits


# ---------------------------------------------------------------------------
# arXiv
# ---------------------------------------------------------------------------


def _atom_text(el: Optional[ET.Element]) -> str:
    if el is None:
        return ""
    return "".join(el.itertext()).strip()


def _hit_from_arxiv_entry(entry: ET.Element) -> Optional[dict[str, Any]]:
    id_el = entry.find("atom:id", ATOM_NS)
    id_text = _atom_text(id_el)
    m = re.search(r"arxiv\.org/abs/([^v\s/#]+)(?:v\d+)?", id_text, re.I)
    arxiv_id = m.group(1) if m else None
    title = re.sub(r"\s+", " ", _atom_text(entry.find("atom:title", ATOM_NS)))
    summary = re.sub(r"\s+", " ", _atom_text(entry.find("atom:summary", ATOM_NS)))
    authors: list[str] = []
    for au in entry.findall("atom:author", ATOM_NS):
        name = _atom_text(au.find("atom:name", ATOM_NS))
        if name:
            authors.append(name)
    published = _atom_text(entry.find("atom:published", ATOM_NS))
    year = _year_from(published)
    doi_el = entry.find("arxiv:doi", ATOM_NS)
    doi = _doi_norm(_atom_text(doi_el)) if doi_el is not None else None
    # Prefer abs URL
    url = id_text if id_text.startswith("http") else (f"https://arxiv.org/abs/{arxiv_id}" if arxiv_id else None)
    if not title:
        return None
    return {
        "source": "arxiv",
        "doi": doi,
        "pmid": None,
        "pmcid": None,
        "arxiv": arxiv_id,
        "title": title,
        "abstract": summary or None,
        "authors": authors,
        "year": year,
        "venue": "arXiv",
        "type": "preprint",
        "url": url,
    }


def search_arxiv(
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,  # noqa: ARG001 — API has no simple date filter
    date_to: Optional[str] = None,  # noqa: ARG001
) -> list[dict[str, Any]]:
    q = (query or "").strip()
    if not q:
        search = "cat:q-bio*"
    elif re.search(r"\bcat:", q, re.I):
        search = q
    else:
        # Prefer q-bio for life-science draft default (mirrors MCP)
        search = f"cat:q-bio* AND all:{q}"
    url = (
        "https://export.arxiv.org/api/query?"
        + urllib.parse.urlencode(
            {
                "search_query": search,
                "start": "0",
                "max_results": str(min(limit, 50)),
                "sortBy": "relevance",
                "sortOrder": "descending",
            }
        )
    )
    root = ET.fromstring(_get(url))
    hits: list[dict[str, Any]] = []
    for entry in root.findall("atom:entry", ATOM_NS):
        h = _hit_from_arxiv_entry(entry)
        if h:
            hits.append(h)
        if len(hits) >= limit:
            break
    # Fallback: if q-bio filter yielded nothing, retry all:query
    if not hits and q and not re.search(r"\bcat:", q, re.I):
        url2 = (
            "https://export.arxiv.org/api/query?"
            + urllib.parse.urlencode(
                {
                    "search_query": f"all:{q}",
                    "start": "0",
                    "max_results": str(min(limit, 50)),
                    "sortBy": "relevance",
                    "sortOrder": "descending",
                }
            )
        )
        root2 = ET.fromstring(_get(url2))
        for entry in root2.findall("atom:entry", ATOM_NS):
            h = _hit_from_arxiv_entry(entry)
            if h:
                hits.append(h)
            if len(hits) >= limit:
                break
    return hits


# ---------------------------------------------------------------------------
# bioRxiv / medRxiv
# ---------------------------------------------------------------------------


def _looks_like_doi(q: str) -> Optional[str]:
    return _doi_norm(q.strip()) if re.search(r"10\.\d{4,9}/", q) else None


def _hit_from_biorxiv(row: dict[str, Any], server: str) -> Optional[dict[str, Any]]:
    doi = _doi_norm(row.get("doi"))
    title = (row.get("title") or "").strip()
    if not title and not doi:
        return None
    authors_raw = row.get("authors") or ""
    authors = [a.strip() for a in re.split(r";\s*", authors_raw) if a.strip()]
    year = _year_from(row.get("date"))
    version = row.get("version")
    url = f"https://www.{server}.org/content/{doi}" if doi else None
    return {
        "source": server,
        "doi": doi,
        "pmid": None,
        "pmcid": None,
        "title": title,
        "abstract": row.get("abstract"),
        "authors": authors,
        "year": year,
        "venue": server,
        "type": "preprint",
        "url": url,
        "version": version,
    }


def search_biorxiv_family(
    server: str,
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """biorxiv/medrxiv: DOI detail, date browse, or Europe PMC JOURNAL keyword.

    api.biorxiv.org has no keyword search. Strategies:
      1) query looks like DOI → details/{server}/{doi}/na/json
      2) date_from/date_to set (or default last 7d if query empty) → date window browse
      3) otherwise → Europe PMC JOURNAL:\"bioRxiv\"|\"medRxiv\" keyword (real hits)
    """
    server = server.lower().strip()
    assert server in ("biorxiv", "medrxiv")
    journal = "bioRxiv" if server == "biorxiv" else "medRxiv"
    q = (query or "").strip()

    # 1) DOI detail
    doi = _looks_like_doi(q) if q else None
    if doi:
        url = f"https://api.biorxiv.org/details/{server}/{doi}/na/json"
        body = _get_json(url)
        rows = body.get("collection") or []
        if rows:
            latest = rows[-1] if isinstance(rows[-1], dict) else {}
            h = _hit_from_biorxiv(latest, server)
            hits = [h] if h else []
            return hits, {
                "provider": "http",
                "source": server,
                "via": "doi_detail",
                "note": "api.biorxiv.org has no keyword search; DOI detail used",
            }
        # try other server once
        other = "medrxiv" if server == "biorxiv" else "biorxiv"
        body2 = _get_json(f"https://api.biorxiv.org/details/{other}/{doi}/na/json")
        rows2 = body2.get("collection") or []
        if rows2:
            latest = rows2[-1]
            h = _hit_from_biorxiv(latest if isinstance(latest, dict) else {}, other)
            return ([h] if h else []), {
                "provider": "http",
                "source": server,
                "via": f"doi_detail_{other}",
            }

    # 2) Date browse when dates given, or empty query
    if date_from or date_to or not q:
        dt = date_to or date.today().isoformat()
        df = date_from or (date.today() - timedelta(days=7)).isoformat()
        # normalize YYYY only
        if len(df) == 4:
            df = f"{df}-01-01"
        if len(dt) == 4:
            dt = f"{dt}-12-31"
        body = _get_json(f"https://api.biorxiv.org/details/{server}/{df}/{dt}/0/json")
        rows = body.get("collection") or []
        hits: list[dict[str, Any]] = []
        # newest first: reverse
        for row in reversed(rows):
            if not isinstance(row, dict):
                continue
            if q:
                blob = f"{row.get('title','')} {row.get('abstract','')}".lower()
                if q.lower() not in blob:
                    continue
            h = _hit_from_biorxiv(row, server)
            if h:
                hits.append(h)
            if len(hits) >= limit:
                break
        return hits, {
            "provider": "http",
            "source": server,
            "via": "date_browse",
            "interval": f"{df}/{dt}",
            "note": (
                "api.biorxiv.org is date-window browse; optional substring filter on title/abstract. "
                "For topic search prefer europepmc JOURNAL filter (auto when no dates)."
            ),
        }

    # 3) Keyword via Europe PMC
    epmc_q = f'({q}) AND JOURNAL:"{journal}"'
    epmc_hits = search_europepmc(epmc_q, limit=limit, date_from=date_from, date_to=date_to)
    hits = []
    for h in epmc_hits:
        hits.append(
            {
                **h,
                "source": server,
                "type": "preprint",
                "venue": h.get("venue") or journal,
            }
        )
    return hits, {
        "provider": "http",
        "source": server,
        "via": "europepmc_journal",
        "note": (
            f"api.biorxiv.org has no keyword API; keyword search delegated to Europe PMC "
            f'JOURNAL:"{journal}"'
        ),
    }


# ---------------------------------------------------------------------------
# Semantic Scholar
# ---------------------------------------------------------------------------


def _s2_headers() -> dict[str, str]:
    key = (os.environ.get("S2_API_KEY") or os.environ.get("SEMANTIC_SCHOLAR_API_KEY") or "").strip()
    hdrs: dict[str, str] = {"Accept": "application/json"}
    if key:
        hdrs["x-api-key"] = key
    return hdrs


def _hit_from_s2(paper: dict[str, Any]) -> Optional[dict[str, Any]]:
    external = paper.get("externalIds") or {}
    title = paper.get("title")
    if not title:
        return None
    authors = []
    for a in paper.get("authors") or []:
        name = (a or {}).get("name")
        if name:
            authors.append(name)
    journal = paper.get("journal") or {}
    venue = journal.get("name") or paper.get("venue")
    year = _year_from(paper.get("publicationDate"))
    types = paper.get("publicationTypes") or []
    ptype = "journal-article"
    for t in types:
        if t and re.search(r"review", str(t), re.I):
            ptype = "review"
            break
        if t and re.search(r"preprint", str(t), re.I):
            ptype = "preprint"
            break
    return {
        "source": "semantic_scholar",
        "doi": _doi_norm(external.get("DOI")),
        "pmid": external.get("PubMed"),
        "pmcid": external.get("PubMedCentral"),
        "arxiv": external.get("ArXiv"),
        "s2": paper.get("paperId"),
        "title": title,
        "abstract": paper.get("abstract"),
        "authors": authors,
        "year": year,
        "venue": venue,
        "type": ptype,
        "url": paper.get("url"),
    }


def search_semantic_scholar(
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    params: dict[str, str] = {
        "query": query,
        "offset": "0",
        "limit": str(min(limit, 50)),
        "fields": S2_FIELDS,
    }
    if date_from or date_to:
        y0 = (date_from or "")[:4]
        y1 = (date_to or "")[:4]
        params["year"] = f"{y0}-{y1}"
    url = "https://api.semanticscholar.org/graph/v1/paper/search?" + urllib.parse.urlencode(params)
    try:
        body = _get_json(url, headers=_s2_headers())
    except urllib.error.HTTPError as e:
        if e.code in (401, 403, 429):
            return [], {
                "provider": "error",
                "source": "semantic_scholar",
                "message": (
                    f"Semantic Scholar HTTP {e.code}; set S2_API_KEY / SEMANTIC_SCHOLAR_API_KEY "
                    "in the environment (never store keys in the package). Soft-fail."
                ),
                "http_status": e.code,
            }
        raise
    hits: list[dict[str, Any]] = []
    for row in body.get("data") or []:
        h = _hit_from_s2(row if isinstance(row, dict) else {})
        if h:
            hits.append(h)
        if len(hits) >= limit:
            break
    meta: dict[str, Any] = {"provider": "http", "source": "semantic_scholar"}
    if not (os.environ.get("S2_API_KEY") or os.environ.get("SEMANTIC_SCHOLAR_API_KEY")):
        meta["note"] = "no S2_API_KEY; using unauthenticated quota"
    return hits, meta



# ---------------------------------------------------------------------------
# PMC / DOAJ / DataCite / Zenodo (0.4.9 open HTTP additions)
# ---------------------------------------------------------------------------


def _pmc_id_display(uid: str) -> str:
    u = (uid or "").strip()
    if not u:
        return ""
    u = u.upper()
    if u.startswith("PMC"):
        return u
    return f"PMC{u}"


def _hit_from_pmc_esummary(uid: str, row: dict[str, Any]) -> Optional[dict[str, Any]]:
    if not isinstance(row, dict):
        return None
    title = (row.get("title") or "").strip()
    if not title:
        return None
    pmid = None
    pmcid = _pmc_id_display(str(row.get("pmcid") or uid))
    doi = None
    for aid in row.get("articleids") or []:
        if not isinstance(aid, dict):
            continue
        idtype = (aid.get("idtype") or "").lower()
        val = aid.get("value")
        if not val:
            continue
        if idtype == "pmid":
            pmid = str(val)
        elif idtype in ("pmcid", "pmc"):
            pmcid = _pmc_id_display(str(val))
        elif idtype == "doi":
            doi = _doi_norm(str(val))
    authors = []
    for a in row.get("authors") or []:
        name = (a or {}).get("name")
        if name:
            authors.append(name)
    year = _year_from(row.get("pubdate") or row.get("epubdate") or row.get("sortdate"))
    venue = row.get("fulljournalname") or row.get("source")
    return {
        "source": "pmc",
        "doi": doi,
        "pmid": pmid,
        "pmcid": pmcid or None,
        "title": title,
        "abstract": None,  # esummary has no abstract; enrich later
        "authors": authors,
        "year": year,
        "venue": venue,
        "type": "journal-article",
        "url": f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/" if pmcid else None,
    }


def search_pmc(
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
) -> list[dict[str, Any]]:
    """NCBI E-utilities db=pmc: esearch → esummary (PMCIDs + DOI/PMID when present)."""
    term = query
    if date_from or date_to:
        df = (date_from or "1900")[:4]
        dt = (date_to or "2099")[:4]
        term = f"({query}) AND ({df}:{dt}[pubdate])"
    esearch = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urllib.parse.urlencode(
            {"db": "pmc", "retmode": "json", "retmax": str(limit), "term": term, "sort": "pub+date"}
        )
    )
    data = json.loads(_get(esearch).decode("utf-8"))
    ids = data.get("esearchresult", {}).get("idlist") or []
    if not ids:
        return []
    esummary = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode({"db": "pmc", "retmode": "json", "id": ",".join(ids)})
    )
    body = json.loads(_get(esummary).decode("utf-8"))
    result = body.get("result") or {}
    hits: list[dict[str, Any]] = []
    for uid in ids:
        row = result.get(str(uid))
        if not isinstance(row, dict):
            continue
        h = _hit_from_pmc_esummary(str(uid), row)
        if h:
            hits.append(h)
        if len(hits) >= limit:
            break
    return hits


def _hit_from_doaj(row: dict[str, Any]) -> Optional[dict[str, Any]]:
    if not isinstance(row, dict):
        return None
    bj = row.get("bibjson") or {}
    if not isinstance(bj, dict):
        return None
    title = (bj.get("title") or "").strip()
    if not title:
        return None
    doi = None
    for ident in bj.get("identifier") or []:
        if not isinstance(ident, dict):
            continue
        if (ident.get("type") or "").lower() == "doi":
            doi = _doi_norm(ident.get("id"))
            if doi:
                break
    authors = []
    for a in bj.get("author") or []:
        if isinstance(a, dict) and a.get("name"):
            authors.append(a["name"])
        elif isinstance(a, str) and a.strip():
            authors.append(a.strip())
    journal = bj.get("journal") or {}
    venue = journal.get("title") if isinstance(journal, dict) else None
    url = None
    for link in bj.get("link") or []:
        if isinstance(link, dict) and link.get("url"):
            url = link["url"]
            break
    if not url and doi:
        url = f"https://doi.org/{doi}"
    return {
        "source": "doaj",
        "doi": doi,
        "pmid": None,
        "pmcid": None,
        "title": title,
        "abstract": bj.get("abstract"),
        "authors": authors,
        "year": _year_from(bj.get("year")),
        "venue": venue,
        "type": "journal-article",
        "url": url,
        "doaj_id": row.get("id"),
    }


def search_doaj(
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
) -> list[dict[str, Any]]:
    """DOAJ public article search (v4 path via /api/search/articles/<query>)."""
    q = query
    # Optional year filter via Lucene-ish field (soft: ignore if API rejects — outer soft-fail)
    year_bits = []
    if date_from:
        year_bits.append(f"year:[{(date_from or '')[:4]} TO *]")
    if date_to and not date_from:
        year_bits.append(f"year:[* TO {(date_to or '')[:4]}]")
    elif date_from and date_to:
        year_bits = [f"year:[{(date_from or '')[:4]} TO {(date_to or '')[:4]}]"]
    if year_bits:
        q = f"({query}) AND {year_bits[0]}"
    # DOAJ expects the query in the path; encode carefully
    path_q = urllib.parse.quote(q, safe="")
    url = (
        f"https://doaj.org/api/search/articles/{path_q}?"
        + urllib.parse.urlencode({"pageSize": str(min(limit, 100)), "page": "1"})
    )
    body = _get_json(url)
    hits: list[dict[str, Any]] = []
    for row in body.get("results") or []:
        h = _hit_from_doaj(row if isinstance(row, dict) else {})
        if h:
            hits.append(h)
        if len(hits) >= limit:
            break
    return hits


def _hit_from_datacite(item: dict[str, Any]) -> Optional[dict[str, Any]]:
    if not isinstance(item, dict):
        return None
    attrs = item.get("attributes") or {}
    if not isinstance(attrs, dict):
        attrs = {}
    doi = _doi_norm(attrs.get("doi") or item.get("id"))
    titles = attrs.get("titles") or []
    title = None
    if titles and isinstance(titles[0], dict):
        title = (titles[0].get("title") or "").strip()
    elif titles and isinstance(titles[0], str):
        title = titles[0].strip()
    if not title:
        return None
    authors = []
    for c in attrs.get("creators") or []:
        if not isinstance(c, dict):
            continue
        name = c.get("name") or " ".join(
            x for x in [(c.get("givenName") or ""), (c.get("familyName") or "")] if x
        ).strip()
        if name:
            authors.append(name)
    abstract = None
    for d in attrs.get("descriptions") or []:
        if isinstance(d, dict) and (d.get("descriptionType") or "").lower() in ("abstract", "description", ""):
            abstract = d.get("description")
            if abstract:
                break
    types = attrs.get("types") or {}
    rtype = ""
    if isinstance(types, dict):
        rtype = types.get("resourceTypeGeneral") or types.get("resourceType") or ""
    ptype = _map_type(str(rtype))
    if "dataset" in str(rtype).lower():
        ptype = "other"
    publisher = attrs.get("publisher")
    if isinstance(publisher, dict):
        publisher = publisher.get("name")
    url = attrs.get("url") or (f"https://doi.org/{doi}" if doi else None)
    return {
        "source": "datacite",
        "doi": doi,
        "pmid": None,
        "pmcid": None,
        "title": title,
        "abstract": abstract,
        "authors": authors,
        "year": _year_from(attrs.get("publicationYear")),
        "venue": publisher,
        "type": ptype,
        "url": url,
    }


def search_datacite(
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
) -> list[dict[str, Any]]:
    """DataCite REST: GET /dois?query=… (datasets/preprints/articles with DOI)."""
    params: dict[str, str] = {
        "query": query,
        "page[size]": str(min(limit, 100)),
        "page[number]": "1",
    }
    # publicationYear filter when dates given
    if date_from or date_to:
        y0 = (date_from or "1900")[:4]
        y1 = (date_to or "2099")[:4]
        params["query"] = f"({query}) AND publicationYear:[{y0} TO {y1}]"
    url = "https://api.datacite.org/dois?" + urllib.parse.urlencode(params)
    headers = {"Accept": "application/vnd.api+json"}
    email = _contact_email()
    if email:
        headers["User-Agent"] = f"{USER_AGENT} mailto:{email}"
    body = _get_json(url, headers=headers)
    hits: list[dict[str, Any]] = []
    for item in body.get("data") or []:
        h = _hit_from_datacite(item if isinstance(item, dict) else {})
        if h:
            hits.append(h)
        if len(hits) >= limit:
            break
    return hits


def _hit_from_zenodo(row: dict[str, Any]) -> Optional[dict[str, Any]]:
    if not isinstance(row, dict):
        return None
    md = row.get("metadata") or {}
    if not isinstance(md, dict):
        md = {}
    title = (row.get("title") or md.get("title") or "").strip()
    if not title:
        return None
    doi = _doi_norm(row.get("doi") or md.get("doi"))
    authors = []
    for c in md.get("creators") or []:
        if isinstance(c, dict) and c.get("name"):
            authors.append(c["name"])
    rtype = md.get("resource_type") or {}
    type_str = ""
    if isinstance(rtype, dict):
        type_str = rtype.get("type") or rtype.get("title") or ""
    ptype = _map_type(str(type_str))
    if str(type_str).lower() in ("dataset", "software", "poster", "image"):
        ptype = "other"
    elif str(type_str).lower() == "publication":
        subtype = (rtype.get("subtype") or "") if isinstance(rtype, dict) else ""
        if "preprint" in str(subtype).lower():
            ptype = "preprint"
        else:
            ptype = "journal-article"
    journal = md.get("journal") or {}
    venue = None
    if isinstance(journal, dict):
        venue = journal.get("title")
    venue = venue or md.get("journal_title")
    links = row.get("links") or {}
    url = None
    if isinstance(links, dict):
        url = links.get("self_html") or links.get("html") or links.get("doi")
    if not url and doi:
        url = f"https://doi.org/{doi}"
    abstract = md.get("description")
    if abstract:
        abstract = _strip_tags(abstract)
    return {
        "source": "zenodo",
        "doi": doi,
        "pmid": None,
        "pmcid": None,
        "title": title,
        "abstract": abstract,
        "authors": authors,
        "year": _year_from(md.get("publication_date") or row.get("created")),
        "venue": venue or "Zenodo",
        "type": ptype,
        "url": url,
        "zenodo_id": row.get("id") or row.get("recid"),
    }


def search_zenodo(
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
) -> list[dict[str, Any]]:
    """Zenodo / Invenio records search (publications, preprints, datasets)."""
    params: dict[str, str] = {
        "q": query,
        "size": str(min(limit, 100)),
        "page": "1",
        "sort": "bestmatch",
    }
    # Optional publication date range (Invenio syntax)
    if date_from or date_to:
        df = date_from or "1900-01-01"
        dt = date_to or "2099-12-31"
        # keep query simple; append lucene range on publication_date when both look like dates/years
        params["q"] = f"({query}) AND publication_date:[{df[:10]} TO {dt[:10]}]"
    url = "https://zenodo.org/api/records?" + urllib.parse.urlencode(params)
    body = _get_json(url)
    hits: list[dict[str, Any]] = []
    for row in (body.get("hits") or {}).get("hits") or []:
        h = _hit_from_zenodo(row if isinstance(row, dict) else {})
        if h:
            hits.append(h)
        if len(hits) >= limit:
            break
    return hits


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


SOURCE_ALIASES = {
    "s2": "semantic_scholar",
    "semantic-scholar": "semantic_scholar",
    "europe_pmc": "europepmc",
    "epmc": "europepmc",
    "pubmed_central": "pmc",
    "pubmed-central": "pmc",
    "doaj_articles": "doaj",
    "data_cite": "datacite",
    "data-cite": "datacite",
}


def normalize_source(source: str) -> str:
    s = source.lower().strip()
    return SOURCE_ALIASES.get(s, s)


def search_source(
    source: str,
    query: str,
    *,
    limit: int = 20,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Return (hits, meta). meta may include error flags; never raises for soft-fail sources."""
    source = normalize_source(source)
    try:
        if source == "pubmed":
            return search_pubmed(query, limit=limit, date_from=date_from, date_to=date_to), {
                "provider": "http",
                "source": source,
            }
        if source == "europepmc":
            return search_europepmc(query, limit=limit, date_from=date_from, date_to=date_to), {
                "provider": "http",
                "source": source,
            }
        if source == "openalex":
            return search_openalex(query, limit=limit, date_from=date_from, date_to=date_to), {
                "provider": "http",
                "source": source,
            }
        if source == "crossref":
            return search_crossref(query, limit=limit, date_from=date_from, date_to=date_to), {
                "provider": "http",
                "source": source,
            }
        if source == "arxiv":
            return search_arxiv(query, limit=limit, date_from=date_from, date_to=date_to), {
                "provider": "http",
                "source": source,
            }
        if source in ("biorxiv", "medrxiv"):
            return search_biorxiv_family(
                source, query, limit=limit, date_from=date_from, date_to=date_to
            )
        if source == "pmc":
            return search_pmc(query, limit=limit, date_from=date_from, date_to=date_to), {
                "provider": "http",
                "source": source,
            }
        if source == "doaj":
            return search_doaj(query, limit=limit, date_from=date_from, date_to=date_to), {
                "provider": "http",
                "source": source,
            }
        if source == "datacite":
            return search_datacite(query, limit=limit, date_from=date_from, date_to=date_to), {
                "provider": "http",
                "source": source,
            }
        if source == "zenodo":
            return search_zenodo(query, limit=limit, date_from=date_from, date_to=date_to), {
                "provider": "http",
                "source": source,
            }
        if source == "semantic_scholar":
            return search_semantic_scholar(
                query, limit=limit, date_from=date_from, date_to=date_to
            )
        return [], {"provider": "error", "source": source, "message": f"unknown source {source}"}
    except urllib.error.HTTPError as e:
        return [], {
            "provider": "error",
            "source": source,
            "message": f"HTTP {e.code}: {e.reason}",
            "http_status": e.code,
        }
    except (urllib.error.URLError, TimeoutError, ET.ParseError, json.JSONDecodeError, OSError, ValueError) as e:
        return [], {"provider": "error", "source": source, "message": str(e)}


def parse_pubmed_xml(xml_bytes: bytes) -> list[dict[str, Any]]:
    """Parse PubmedArticle XML into hit dicts (shared by search + enrich)."""
    root = ET.fromstring(xml_bytes)
    hits: list[dict[str, Any]] = []
    for art in root.findall(".//PubmedArticle"):
        medline = art.find("MedlineCitation")
        if medline is None:
            continue
        pmid_el = medline.find("PMID")
        pmid = pmid_el.text if pmid_el is not None else None
        article = medline.find("Article")
        title = ""
        abstract = ""
        journal = ""
        year = None
        authors: list[str] = []
        doi = None
        if article is not None:
            te = article.find("ArticleTitle")
            title = "".join(te.itertext()) if te is not None else ""
            abs_el = article.find("Abstract")
            if abs_el is not None:
                parts = []
                for at in abs_el.findall("AbstractText"):
                    parts.append("".join(at.itertext()))
                abstract = " ".join(parts)
            j = article.find("Journal")
            if j is not None:
                jt = j.find("Title")
                journal = jt.text if jt is not None and jt.text else ""
                y = j.find("JournalIssue/PubDate/Year")
                if y is not None and y.text:
                    year = y.text
                else:
                    md = j.find("JournalIssue/PubDate/MedlineDate")
                    if md is not None and md.text:
                        year = md.text[:4]
            for au in article.findall("AuthorList/Author"):
                last = au.findtext("LastName") or ""
                fore = au.findtext("ForeName") or ""
                collective = au.findtext("CollectiveName")
                if collective:
                    authors.append(collective)
                elif last:
                    authors.append(f"{last}, {fore}".strip(", "))
            for id_el in article.findall("ELocationID"):
                if id_el.get("EIdType") == "doi" and id_el.text:
                    doi = id_el.text.strip()
        if not doi:
            for aid in art.findall(".//ArticleId"):
                if aid.get("IdType") == "doi" and aid.text:
                    doi = aid.text.strip()
        hits.append(
            {
                "source": "pubmed",
                "pmid": pmid,
                "doi": doi,
                "title": title,
                "abstract": abstract,
                "authors": authors,
                "year": year,
                "venue": journal,
                "type": "journal-article",
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None,
            }
        )
    return hits


def fetch_by_pmid(pmid: str) -> tuple[Optional[dict[str, Any]], dict[str, Any]]:
    """Efetch one PMID. Returns (hit_or_None, meta)."""
    pmid = (pmid or "").strip()
    if not pmid:
        return None, {"provider": "error", "message": "empty pmid"}
    try:
        efetch = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            + urllib.parse.urlencode({"db": "pubmed", "retmode": "xml", "id": pmid})
        )
        hits = parse_pubmed_xml(_get(efetch))
        if not hits:
            return None, {"provider": "http", "source": "pubmed", "message": "no hit"}
        return hits[0], {"provider": "http", "source": "pubmed", "via": "pmid"}
    except (urllib.error.URLError, TimeoutError, ET.ParseError, OSError) as e:
        return None, {"provider": "error", "source": "pubmed", "message": str(e)}


def fetch_by_doi(doi: str) -> tuple[Optional[dict[str, Any]], dict[str, Any]]:
    """Resolve DOI via PubMed doi term, then Europe PMC DOI query."""
    doi = (doi or "").strip()
    if not doi:
        return None, {"provider": "error", "message": "empty doi"}
    pubmed_err = None
    try:
        hits, meta = search_source("pubmed", f"{doi}[doi]", limit=3)
        for h in hits:
            if (h.get("doi") or "").lower() == doi.lower() or not h.get("doi"):
                if h.get("title") and (h.get("abstract") or h.get("pmid")):
                    meta = {**meta, "via": "doi_pubmed"}
                    return h, meta
        if hits:
            return hits[0], {**meta, "via": "doi_pubmed_first"}
        pubmed_err = (meta or {}).get("message")
    except Exception as e:  # noqa: BLE001
        pubmed_err = str(e)

    hits, meta = search_source("europepmc", f'DOI:"{doi}"', limit=3)
    for h in hits:
        hd = (h.get("doi") or "").lower()
        if hd == doi.lower() or not hd:
            if h.get("title"):
                return h, {**meta, "via": "doi_europepmc"}
    if hits:
        return hits[0], {**meta, "via": "doi_europepmc_first"}
    return None, {
        "provider": "error",
        "source": "doi",
        "message": pubmed_err or meta.get("message") or "no hit for doi",
    }


# ---------------------------------------------------------------------------
# Expand helpers (related / citations / references) — HTTP only
# Mirror phi-literature related_publications semantics without Node/MCP.
# ---------------------------------------------------------------------------

OPENALEX_SEARCH_SELECT = (
    "id,doi,display_name,publication_date,publication_year,type,ids,authorships,"
    "primary_location,best_oa_location,locations,open_access,abstract_inverted_index,"
    "cited_by_count"
)
OPENALEX_DETAIL_SELECT = OPENALEX_SEARCH_SELECT + ",referenced_works,related_works"


def _openalex_short_id(value: Any) -> Optional[str]:
    if not value:
        return None
    m = re.search(r"(W\d+)", str(value), re.I)
    return m.group(1).upper() if m else None


def dedup_hits(hits: list[dict[str, Any]], *, limit: Optional[int] = None) -> list[dict[str, Any]]:
    """Dedup by DOI > PMID > PMCID > arxiv > openalex; preserve order."""
    seen: set[str] = set()
    out: list[dict[str, Any]] = []
    for h in hits:
        keys: list[str] = []
        doi = _doi_norm(h.get("doi"))
        if doi:
            keys.append(f"doi:{doi.lower()}")
        pmid = str(h.get("pmid") or "").strip()
        if pmid:
            keys.append(f"pmid:{pmid}")
        pmcid = str(h.get("pmcid") or "").strip().upper()
        if pmcid:
            keys.append(f"pmcid:{pmcid}")
        arxiv = str(h.get("arxiv") or "").strip().lower()
        if arxiv:
            keys.append(f"arxiv:{arxiv}")
        oa = _openalex_short_id(h.get("openalex") or h.get("openalex_id"))
        if oa:
            keys.append(f"openalex:{oa}")
        if keys and any(k in seen for k in keys):
            continue
        for k in keys:
            seen.add(k)
        if not keys:
            # untitled / no id — keep once by title
            title_key = f"title:{(h.get('title') or '').strip().lower()}"
            if title_key == "title:" or title_key in seen:
                continue
            seen.add(title_key)
        out.append(h)
        if limit is not None and len(out) >= limit:
            break
    return out


def pubmed_neighbors(pmid: str) -> tuple[list[str], dict[str, Any]]:
    """PubMed elink neighbor_score → related PMIDs (excludes seed)."""
    pmid = re.sub(r"\D", "", str(pmid or ""))
    if not pmid:
        return [], {"provider": "error", "source": "pubmed", "message": "empty pmid"}
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?"
        + urllib.parse.urlencode(
            {
                "dbfrom": "pubmed",
                "db": "pubmed",
                "id": pmid,
                "cmd": "neighbor_score",
                "retmode": "json",
            }
        )
    )
    try:
        body = _get_json(url)
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
        return [], {"provider": "error", "source": "pubmed", "message": str(e)}
    linksets = body.get("linksets") or body.get("linkSets") or []
    if not linksets:
        return [], {"provider": "http", "source": "pubmed", "via": "elink", "n": 0}
    linkset = linksets[0] if isinstance(linksets[0], dict) else {}
    dbs = linkset.get("linksetdbs") or linkset.get("linksetDbs") or []
    chosen = None
    for row in dbs:
        if not isinstance(row, dict):
            continue
        name = str(row.get("linkname") or row.get("linkName") or "")
        if "pubmed" in name.lower():
            chosen = row
            break
    if chosen is None and dbs:
        chosen = dbs[0] if isinstance(dbs[0], dict) else {}
    links = (chosen or {}).get("links") or []
    ids: list[str] = []
    for link in links:
        if isinstance(link, dict):
            lid = str(link.get("id") or "").strip()
        else:
            lid = str(link).strip()
        if lid and lid != pmid and lid not in ids:
            ids.append(lid)
    return ids, {"provider": "http", "source": "pubmed", "via": "elink", "n": len(ids)}


def fetch_pubmed_by_ids(ids: list[str]) -> list[dict[str, Any]]:
    """Efetch a list of PMIDs → hit dicts."""
    clean = [re.sub(r"\D", "", str(i)) for i in ids if i]
    clean = [c for c in clean if c]
    if not clean:
        return []
    efetch = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        + urllib.parse.urlencode({"db": "pubmed", "retmode": "xml", "id": ",".join(clean)})
    )
    return parse_pubmed_xml(_get(efetch))


def get_openalex_work(
    *,
    openalex_id: Optional[str] = None,
    doi: Optional[str] = None,
) -> tuple[Optional[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    """Fetch one OpenAlex work with related_works / referenced_works.

    Returns (hit_or_None, graph, meta).
    graph keys: openalex_id, related_works, referenced_works (short W… ids).
    """
    path: Optional[str] = None
    if openalex_id:
        short = _openalex_short_id(openalex_id)
        if short:
            path = f"works/{short}"
    if path is None and doi:
        d = _doi_norm(doi)
        if d:
            path = f"works/doi:{d}"
    if not path:
        return None, {}, {"provider": "error", "source": "openalex", "message": "need openalex id or doi"}
    url = _openalex_url(path, {"select": OPENALEX_DETAIL_SELECT})
    try:
        work = _get_json(url)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None, {}, {"provider": "http", "source": "openalex", "message": "not found", "http_status": 404}
        return None, {}, {
            "provider": "error",
            "source": "openalex",
            "message": f"HTTP {e.code}: {e.reason}",
            "http_status": e.code,
        }
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
        return None, {}, {"provider": "error", "source": "openalex", "message": str(e)}
    if not isinstance(work, dict) or work.get("error"):
        return None, {}, {"provider": "http", "source": "openalex", "message": "empty/error body"}
    hit = _hit_from_openalex(work)
    oa = _openalex_short_id(work.get("id") or (work.get("ids") or {}).get("openalex"))
    related = [_openalex_short_id(x) for x in (work.get("related_works") or [])]
    referenced = [_openalex_short_id(x) for x in (work.get("referenced_works") or [])]
    graph = {
        "openalex_id": oa,
        "related_works": [x for x in related if x],
        "referenced_works": [x for x in referenced if x],
        "cited_by_count": work.get("cited_by_count"),
    }
    return hit, graph, {"provider": "http", "source": "openalex", "via": "work_detail", "openalex_id": oa}


def hydrate_openalex_ids(openalex_ids: list[str], *, limit: int = 20) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Batch-hydrate OpenAlex W… ids via filter=openalex:W1|W2|…"""
    shorts: list[str] = []
    for raw in openalex_ids:
        s = _openalex_short_id(raw)
        if s and s not in shorts:
            shorts.append(s)
        if len(shorts) >= limit:
            break
    if not shorts:
        return [], {"provider": "http", "source": "openalex", "via": "hydrate", "n": 0}
    params = {
        "filter": "openalex:" + "|".join(shorts),
        "per-page": str(len(shorts)),
        "select": OPENALEX_SEARCH_SELECT,
    }
    try:
        body = _get_json(_openalex_url("works", params))
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError, urllib.error.HTTPError) as e:
        return [], {"provider": "error", "source": "openalex", "message": str(e), "via": "hydrate"}
    by_id: dict[str, dict[str, Any]] = {}
    for row in body.get("results") or []:
        if not isinstance(row, dict):
            continue
        h = _hit_from_openalex(row)
        if h and h.get("openalex"):
            by_id[str(h["openalex"]).upper()] = h
    ordered = [by_id[s] for s in shorts if s in by_id]
    return ordered, {"provider": "http", "source": "openalex", "via": "hydrate", "n": len(ordered)}


def openalex_by_filter(
    filter_str: str,
    *,
    limit: int = 20,
    offset: int = 0,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """OpenAlex works filter (e.g. cites:W… or doi:…). Soft-fail on HTTP errors."""
    # OpenAlex pages are 1-indexed; keep offset small for draft expand
    per_page = min(max(limit, 1), 50)
    page = (offset // per_page) + 1
    skip = offset % per_page
    params = {
        "filter": filter_str,
        "per-page": str(per_page),
        "page": str(page),
        "select": OPENALEX_SEARCH_SELECT,
        "sort": "publication_date:desc",
    }
    try:
        body = _get_json(_openalex_url("works", params))
    except urllib.error.HTTPError as e:
        return [], {
            "provider": "error",
            "source": "openalex",
            "message": f"HTTP {e.code}: {e.reason}",
            "http_status": e.code,
            "via": "filter",
            "filter": filter_str,
        }
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
        return [], {
            "provider": "error",
            "source": "openalex",
            "message": str(e),
            "via": "filter",
            "filter": filter_str,
        }
    rows = body.get("results") or []
    hits: list[dict[str, Any]] = []
    for row in rows[skip : skip + limit]:
        if not isinstance(row, dict):
            continue
        h = _hit_from_openalex(row)
        if h and h.get("title"):
            hits.append(h)
    meta_count = (body.get("meta") or {}).get("count")
    return hits, {
        "provider": "http",
        "source": "openalex",
        "via": "filter",
        "filter": filter_str,
        "n": len(hits),
        "total": meta_count,
    }


def opencitations_dois(doi: str, kind: str) -> tuple[list[str], dict[str, Any]]:
    """OpenCitations index API → list of DOIs (citing or cited). Soft-fail."""
    d = _doi_norm(doi)
    if not d:
        return [], {"provider": "error", "source": "opencitations", "message": "empty doi"}
    if kind not in ("citations", "references"):
        return [], {"provider": "error", "source": "opencitations", "message": f"bad kind {kind}"}
    path = "citations" if kind == "citations" else "references"
    field = "citing" if kind == "citations" else "cited"
    url = f"https://opencitations.net/index/api/v1/{path}/doi:{urllib.parse.quote(d)}"
    try:
        body = _get_json(url)
    except urllib.error.HTTPError as e:
        return [], {
            "provider": "error",
            "source": "opencitations",
            "message": f"HTTP {e.code}: {e.reason}",
            "http_status": e.code,
        }
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
        return [], {"provider": "error", "source": "opencitations", "message": str(e)}
    if not isinstance(body, list):
        return [], {"provider": "http", "source": "opencitations", "message": "unexpected body", "n": 0}
    dois: list[str] = []
    for row in body:
        if not isinstance(row, dict):
            continue
        raw = row.get(field)
        nd = _doi_norm(raw)
        if nd and nd.lower() not in {x.lower() for x in dois}:
            dois.append(nd)
    return dois, {"provider": "http", "source": "opencitations", "via": path, "n": len(dois)}


def expand_related(
    *,
    relation: str = "similar",
    limit: int = 10,
    pmid: Optional[str] = None,
    doi: Optional[str] = None,
    openalex: Optional[str] = None,
    pmcid: Optional[str] = None,  # noqa: ARG001 — reserved for future PMC neighbors
    arxiv: Optional[str] = None,  # noqa: ARG001
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Expand from seed IDs. Returns (hits, meta). Never writes DB.

    relation:
      similar    — PubMed elink when PMID known; else OpenAlex related_works
      references — OpenAlex referenced_works; OpenCitations DOI fallback
      citations  — OpenAlex filter cites:W…; OpenCitations DOI fallback
    """
    relation = (relation or "similar").strip().lower()
    limit = max(1, min(int(limit or 10), 50))
    pmid = re.sub(r"\D", "", str(pmid)) if pmid else None
    doi = _doi_norm(doi) if doi else None
    openalex = _openalex_short_id(openalex) if openalex else None
    providers: list[dict[str, Any]] = []
    hits: list[dict[str, Any]] = []

    if relation not in ("similar", "citations", "references"):
        return [], {
            "ok": False,
            "provider": "error",
            "message": f"unsupported relation {relation}",
            "relation": relation,
            "providers": providers,
        }

    # similar + PMID → PubMed neighbors
    if relation == "similar" and pmid:
        neighbor_ids, meta = pubmed_neighbors(pmid)
        providers.append(meta)
        if neighbor_ids:
            try:
                batch = fetch_pubmed_by_ids(neighbor_ids[:limit])
                hits = dedup_hits(batch, limit=limit)
                return hits, {
                    "ok": True,
                    "relation": relation,
                    "provider": "pubmed",
                    "via": "elink",
                    "seed": {"pmid": pmid, "doi": doi, "openalex": openalex},
                    "n_candidates": len(neighbor_ids),
                    "providers": providers,
                }
            except (urllib.error.URLError, TimeoutError, ET.ParseError, OSError) as e:
                providers.append({"provider": "error", "source": "pubmed", "message": str(e), "via": "efetch"})

    # Ensure OpenAlex graph when needed
    graph: dict[str, Any] = {}
    if openalex or doi:
        _hit, graph, gmeta = get_openalex_work(openalex_id=openalex, doi=doi)
        providers.append(gmeta)
        if not openalex and graph.get("openalex_id"):
            openalex = graph["openalex_id"]
        if _hit and not doi:
            doi = _doi_norm(_hit.get("doi"))
        if _hit and not pmid:
            pmid = str(_hit.get("pmid") or "") or pmid

    if relation == "similar":
        related = graph.get("related_works") or []
        if related:
            hydrated, hmeta = hydrate_openalex_ids(related, limit=limit)
            providers.append(hmeta)
            hits = dedup_hits(hydrated, limit=limit)
            if hits:
                return hits, {
                    "ok": True,
                    "relation": relation,
                    "provider": "openalex",
                    "via": "related_works",
                    "seed": {"pmid": pmid, "doi": doi, "openalex": openalex},
                    "n_candidates": len(related),
                    "providers": providers,
                }
        # last chance: if we now have pmid from OpenAlex, try elink
        if pmid and not any(p.get("via") == "elink" for p in providers):
            neighbor_ids, meta = pubmed_neighbors(pmid)
            providers.append(meta)
            if neighbor_ids:
                try:
                    batch = fetch_pubmed_by_ids(neighbor_ids[:limit])
                    hits = dedup_hits(batch, limit=limit)
                    if hits:
                        return hits, {
                            "ok": True,
                            "relation": relation,
                            "provider": "pubmed",
                            "via": "elink_after_openalex",
                            "seed": {"pmid": pmid, "doi": doi, "openalex": openalex},
                            "n_candidates": len(neighbor_ids),
                            "providers": providers,
                        }
                except (urllib.error.URLError, TimeoutError, ET.ParseError, OSError) as e:
                    providers.append({"provider": "error", "source": "pubmed", "message": str(e)})

    if relation == "references":
        referenced = graph.get("referenced_works") or []
        if referenced:
            hydrated, hmeta = hydrate_openalex_ids(referenced, limit=limit)
            providers.append(hmeta)
            hits = dedup_hits(hydrated, limit=limit)
            if hits:
                return hits, {
                    "ok": True,
                    "relation": relation,
                    "provider": "openalex",
                    "via": "referenced_works",
                    "seed": {"pmid": pmid, "doi": doi, "openalex": openalex},
                    "n_candidates": len(referenced),
                    "providers": providers,
                }

    if relation == "citations" and openalex:
        cited, cmeta = openalex_by_filter(f"cites:{openalex}", limit=limit)
        providers.append(cmeta)
        if cmeta.get("provider") == "error":
            # soft-fail note kept; try OpenCitations below
            pass
        else:
            hits = dedup_hits(cited, limit=limit)
            if hits:
                return hits, {
                    "ok": True,
                    "relation": relation,
                    "provider": "openalex",
                    "via": "cites_filter",
                    "seed": {"pmid": pmid, "doi": doi, "openalex": openalex},
                    "n_candidates": cmeta.get("total") or len(cited),
                    "providers": providers,
                }

    # OpenCitations fallback for citations/references when DOI known
    if doi and relation in ("citations", "references"):
        dois, ocmeta = opencitations_dois(doi, relation)
        providers.append(ocmeta)
        if dois:
            slice_dois = dois[:limit]
            # hydrate via OpenAlex doi: filter
            filt = "doi:" + "|".join(slice_dois)
            hydrated, hmeta = openalex_by_filter(filt, limit=len(slice_dois))
            providers.append(hmeta)
            hits = dedup_hits(hydrated, limit=limit)
            unresolved = [
                d for d in slice_dois if d.lower() not in {(h.get("doi") or "").lower() for h in hits}
            ]
            return hits, {
                "ok": True if hits else False,
                "relation": relation,
                "provider": "openalex" if hits else "opencitations",
                "via": "opencitations",
                "seed": {"pmid": pmid, "doi": doi, "openalex": openalex},
                "n_candidates": len(dois),
                "unresolved_dois": unresolved,
                "providers": providers,
                "message": None if hits else f"No {relation} hydrated for {doi}",
            }

    return [], {
        "ok": False,
        "relation": relation,
        "provider": "none",
        "seed": {"pmid": pmid, "doi": doi, "openalex": openalex},
        "providers": providers,
        "message": f"No {relation} were found for seed",
    }


# ---------------------------------------------------------------------------
# Version ties (preprint↔journal) via OpenAlex / Crossref — soft-fail HTTP
# ---------------------------------------------------------------------------

CROSSREF_PREPRINT_RELATIONS = {
    "is-preprint-of": ("preprint_of", 0.95),
    "has-preprint": ("published_as", 0.95),  # journal has-preprint → journal is to_id side when oriented
    "is-version-of": ("same_work", 0.9),
    "has-version": ("same_work", 0.85),
    "is-identical-to": ("same_work", 0.9),
}


def get_crossref_work(doi: str) -> tuple[Optional[dict[str, Any]], dict[str, Any]]:
    """Fetch one Crossref work by DOI. Returns (message_or_None, meta). Soft-fail."""
    d = _doi_norm(doi)
    if not d:
        return None, {"provider": "error", "source": "crossref", "message": "need doi"}
    url = _crossref_url(f"works/{urllib.parse.quote(d)}")
    try:
        body = _get_json(url)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None, {
                "provider": "http",
                "source": "crossref",
                "message": "not found",
                "http_status": 404,
            }
        return None, {
            "provider": "error",
            "source": "crossref",
            "message": f"HTTP {e.code}: {e.reason}",
            "http_status": e.code,
        }
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
        return None, {"provider": "error", "source": "crossref", "message": str(e)}
    message = body.get("message") if isinstance(body, dict) else None
    if not isinstance(message, dict):
        return None, {"provider": "http", "source": "crossref", "message": "empty/error body"}
    return message, {"provider": "http", "source": "crossref", "via": "work_detail", "doi": d}


def _dois_from_crossref_relation_entries(entries: Any) -> list[str]:
    out: list[str] = []
    if not isinstance(entries, list):
        return out
    for ent in entries:
        if not isinstance(ent, dict):
            continue
        id_type = str(ent.get("id-type") or ent.get("id_type") or "").lower()
        raw = ent.get("id") or ent.get("DOI") or ent.get("doi")
        if id_type and id_type not in ("doi", "doi-url", ""):
            # still try if value looks like a DOI
            pass
        d = _doi_norm(raw)
        if d and d not in out:
            out.append(d)
    return out


def extract_crossref_version_ties(work: dict[str, Any]) -> list[dict[str, Any]]:
    """Parse Crossref message.relation into preprint↔journal DOI ties."""
    seed = _doi_norm(work.get("DOI") or work.get("doi"))
    if not seed:
        return []
    relation = work.get("relation") or {}
    if not isinstance(relation, dict):
        return []
    ties: list[dict[str, Any]] = []
    for rel_key, entries in relation.items():
        mapped = CROSSREF_PREPRINT_RELATIONS.get(str(rel_key).lower())
        if not mapped:
            continue
        rel_name, confidence = mapped
        for other in _dois_from_crossref_relation_entries(entries):
            if other == seed:
                continue
            # Orient: preprint DOI → journal DOI when possible
            seed_pre = seed.startswith("10.1101/")
            other_pre = other.startswith("10.1101/")
            if rel_key.lower() == "has-preprint":
                # seed is journal (or version that has a preprint)
                from_doi, to_doi = other, seed
                oriented = "preprint_of"
            elif rel_key.lower() == "is-preprint-of":
                from_doi, to_doi = seed, other
                oriented = "preprint_of"
            elif seed_pre and not other_pre:
                from_doi, to_doi = seed, other
                oriented = "preprint_of"
            elif other_pre and not seed_pre:
                from_doi, to_doi = other, seed
                oriented = "preprint_of"
            else:
                from_doi, to_doi = seed, other
                oriented = rel_name
            ties.append(
                {
                    "from_doi": from_doi,
                    "to_doi": to_doi,
                    "relation": oriented,
                    "confidence": confidence,
                    "high_confidence": confidence >= 0.85,
                    "provider": "crossref",
                    "evidence": {
                        "crossref_relation": rel_key,
                        "seed_doi": seed,
                        "linked_doi": other,
                    },
                }
            )
    return ties


def _location_dois(work: dict[str, Any]) -> list[dict[str, Any]]:
    """Collect DOIs + version hints from OpenAlex locations."""
    found: list[dict[str, Any]] = []
    seen: set[str] = set()
    for loc in work.get("locations") or []:
        if not isinstance(loc, dict):
            continue
        landing = loc.get("landing_page_url") or ""
        d = _doi_norm(landing)
        if not d:
            # some locations nest doi
            d = _doi_norm((loc.get("source") or {}).get("host_organization_lineage_names"))
        if not d:
            continue
        if d in seen:
            continue
        seen.add(d)
        found.append(
            {
                "doi": d,
                "version": loc.get("version"),
                "is_accepted": loc.get("is_accepted"),
                "is_published": loc.get("is_published"),
                "source_name": ((loc.get("source") or {}).get("display_name")),
            }
        )
    return found


def extract_openalex_version_ties(work: dict[str, Any]) -> list[dict[str, Any]]:
    """Infer preprint↔journal DOI ties from one OpenAlex work (ids + locations).

    High confidence when locations expose a non-preprint DOI beside a preprint
    seed DOI (or vice versa). related_works alone are not enough for high confidence.
    """
    ids = work.get("ids") or {}
    seed = _doi_norm(work.get("doi") or ids.get("doi"))
    wtype = _map_type(work.get("type"))
    locs = _location_dois(work)
    alt_dois = [x["doi"] for x in locs if x.get("doi")]
    # Also check ids for a second DOI-like value (rare)
    for key in ("doi",):
        d = _doi_norm(ids.get(key))
        if d and d not in alt_dois:
            alt_dois.append(d)

    ties: list[dict[str, Any]] = []
    if not seed:
        return ties

    for loc in locs:
        other = loc.get("doi")
        if not other or other == seed:
            continue
        seed_pre = seed.startswith("10.1101/")
        other_pre = other.startswith("10.1101/")
        ver = str(loc.get("version") or "").lower()
        publishedish = (
            ver in ("publishedversion", "published_version", "versionofrecord")
            or bool(loc.get("is_published"))
            or (other_pre is False and seed_pre)
        )
        if seed_pre and not other_pre and publishedish:
            confidence = 0.92
            from_doi, to_doi = seed, other
        elif other_pre and not seed_pre:
            confidence = 0.9
            from_doi, to_doi = other, seed
        elif seed_pre and not other_pre:
            confidence = 0.85
            from_doi, to_doi = seed, other
        elif wtype == "preprint" and not other_pre:
            confidence = 0.8
            from_doi, to_doi = seed, other
        else:
            continue
        ties.append(
            {
                "from_doi": from_doi,
                "to_doi": to_doi,
                "relation": "preprint_of",
                "confidence": confidence,
                "high_confidence": confidence >= 0.85,
                "provider": "openalex",
                "evidence": {
                    "seed_doi": seed,
                    "location_doi": other,
                    "location_version": loc.get("version"),
                    "work_type": work.get("type"),
                    "openalex": _openalex_short_id(work.get("id") or ids.get("openalex")),
                },
            }
        )
    return ties


def discover_version_ties(
    *,
    doi: Optional[str] = None,
    openalex_id: Optional[str] = None,
    providers: Optional[list[str]] = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Query OpenAlex and/or Crossref for version ties. Soft-fail per provider.

    Returns (ties, provider_metas). Never raises on network errors.
    """
    wanted = [p.lower() for p in (providers or ["openalex", "crossref"])]
    ties: list[dict[str, Any]] = []
    metas: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()

    def _add(t: dict[str, Any]) -> None:
        key = (t["from_doi"], t["to_doi"], t["relation"])
        if key in seen:
            return
        seen.add(key)
        ties.append(t)

    if "openalex" in wanted and (doi or openalex_id):
        try:
            path = None
            if openalex_id:
                short = _openalex_short_id(openalex_id)
                if short:
                    path = f"works/{short}"
            if path is None and doi:
                d = _doi_norm(doi)
                if d:
                    path = f"works/doi:{d}"
            if not path:
                metas.append(
                    {
                        "provider": "error",
                        "source": "openalex",
                        "message": "need openalex id or doi",
                    }
                )
            else:
                try:
                    work = _get_json(_openalex_url(path, {"select": OPENALEX_DETAIL_SELECT}))
                except urllib.error.HTTPError as e:
                    metas.append(
                        {
                            "provider": "http" if e.code == 404 else "error",
                            "source": "openalex",
                            "message": "not found" if e.code == 404 else f"HTTP {e.code}",
                            "http_status": e.code,
                        }
                    )
                    work = None
                except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
                    metas.append({"provider": "error", "source": "openalex", "message": str(e)})
                    work = None
                if isinstance(work, dict) and not work.get("error"):
                    metas.append(
                        {
                            "provider": "http",
                            "source": "openalex",
                            "via": "version_ties",
                            "openalex_id": _openalex_short_id(
                                work.get("id") or (work.get("ids") or {}).get("openalex")
                            ),
                        }
                    )
                    for t in extract_openalex_version_ties(work):
                        _add(t)
                elif work is not None:
                    metas.append(
                        {
                            "provider": "http",
                            "source": "openalex",
                            "message": "empty/error body",
                        }
                    )
        except Exception as e:  # noqa: BLE001
            metas.append({"provider": "error", "source": "openalex", "message": str(e)})

    if "crossref" in wanted and doi:
        try:
            work, meta = get_crossref_work(doi)
            metas.append(meta)
            if work:
                for t in extract_crossref_version_ties(work):
                    _add(t)
        except Exception as e:  # noqa: BLE001
            metas.append({"provider": "error", "source": "crossref", "message": str(e)})

    return ties, metas
