#!/usr/bin/env python3
"""HTTP providers for lit_search, lit_expand and metadata lookups.

Public REST only, standard library only. Optional environment (nothing is packaged):
  OPENALEX_MAILTO / CROSSREF_MAILTO / CONTACT_EMAIL  polite contact address for OpenAlex and Crossref

Sources (SUPPORTED_SOURCES):
  pubmed, europepmc            keyword search
  openalex, crossref, arxiv    keyword search
  biorxiv / medrxiv            DOI detail, date browse, or a Europe PMC journal-filter fallback

Related-paper helpers: PubMed elink neighbours, OpenAlex related / referenced / cited-by,
OpenCitations as a DOI fallback.
"""
from __future__ import annotations

import gzip
import http.client
import json
import os
import re
import ssl
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from datetime import date, timedelta
from difflib import SequenceMatcher
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


NCBI_HOST = "eutils.ncbi.nlm.nih.gov"
_last_ncbi_call = 0.0


def _ncbi_url(url: str) -> str:
    """Add the identification NCBI asks for (tool, email) and the optional NCBI_API_KEY."""
    extra = {"tool": "phi-literature-library"}
    email = _contact_email()
    if email:
        extra["email"] = email
    key = (os.environ.get("NCBI_API_KEY") or "").strip()
    if key:
        extra["api_key"] = key
    parts = urllib.parse.urlsplit(url)
    have = urllib.parse.parse_qs(parts.query)
    add = urllib.parse.urlencode({k: v for k, v in extra.items() if k not in have})
    return urllib.parse.urlunsplit(parts._replace(query=f"{parts.query}&{add}" if add else parts.query))


def _throttle_ncbi() -> None:
    """NCBI allows 3 requests/s without a key and 10/s with one."""
    global _last_ncbi_call
    gap = 0.11 if os.environ.get("NCBI_API_KEY") else 0.35
    wait = _last_ncbi_call + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_ncbi_call = time.monotonic()


def _get(url: str, *, headers: Optional[dict[str, str]] = None) -> bytes:
    hdrs = {"User-Agent": USER_AGENT, "Accept": "*/*", "Accept-Encoding": "gzip"}
    if headers:
        hdrs.update(headers)
    is_ncbi = urllib.parse.urlsplit(url).netloc == NCBI_HOST
    if is_ncbi:
        url = _ncbi_url(url)
    req = urllib.request.Request(url, headers=hdrs)
    ctx = ssl.create_default_context()
    last: Exception | None = None
    for attempt in range(3):
        if is_ncbi:
            _throttle_ncbi()
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT, context=ctx) as resp:
                body = resp.read()
                if resp.headers.get("Content-Encoding", "").lower() == "gzip":
                    body = gzip.decompress(body)
                return body
        except urllib.error.HTTPError as e:
            if e.code in (429, 503) and attempt < 2:
                time.sleep(1.5 * (attempt + 1))
                continue
            raise
        except (http.client.HTTPException, EOFError, gzip.BadGzipFile) as e:
            # Large responses are sometimes cut off mid-body; retry up to three times, then report as OSError
            # so callers' network soft-fail handling applies.
            last = e
            time.sleep(0.5)
    raise OSError(f"incomplete response from {urllib.parse.urlsplit(url).netloc}: {last}")


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
            {"db": "pubmed", "retmode": "json", "retmax": str(limit), "term": term}
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


def openalex_oa_urls(doi: str) -> list[str]:
    """Open-access PDF/landing URLs OpenAlex knows for a DOI (best location first). [] on any failure."""
    d = _doi_norm(doi)
    if not d:
        return []
    try:
        work = _get_json(_openalex_url(f"works/doi:{urllib.parse.quote(d)}", {"select": "open_access,best_oa_location,locations"}))
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError):
        return []
    urls: list[str] = []

    def add(loc: Any) -> None:
        if isinstance(loc, dict):
            for key in ("pdf_url", "landing_page_url"):
                u = loc.get(key)
                if isinstance(u, str) and u.startswith("http") and u not in urls:
                    urls.append(u)

    add(work.get("best_oa_location"))
    oa = work.get("open_access") or {}
    u = oa.get("oa_url")
    if isinstance(u, str) and u.startswith("http") and u not in urls:
        urls.append(u)
    for loc in work.get("locations") or []:
        if isinstance(loc, dict) and loc.get("is_oa"):
            add(loc)
    return urls[:5]


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









# ---------------------------------------------------------------------------
# PMC / DOAJ / DataCite / Zenodo (0.4.9 open HTTP additions)
# ---------------------------------------------------------------------------




















# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


SOURCE_ALIASES = {
    "europe_pmc": "europepmc",
    "epmc": "europepmc",
}

SUPPORTED_SOURCES = ("pubmed", "europepmc", "openalex", "crossref", "arxiv", "biorxiv", "medrxiv")


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
    # Journals outside PubMed / Europe PMC (e.g. many Elsevier titles): Crossref knows the DOI.
    work, cmeta = get_crossref_work(doi)
    hit = _hit_from_crossref(work) if work else None
    if hit and hit.get("title"):
        return hit, {**cmeta, "via": "doi_crossref"}
    return None, {
        "provider": "error",
        "source": "doi",
        "message": pubmed_err or meta.get("message") or "no hit for doi",
    }


def fetch_by_arxiv(arxiv_id: str) -> tuple[Optional[dict[str, Any]], dict[str, Any]]:
    """One arXiv entry by id (e.g. 2205.12454). Returns (hit_or_None, meta)."""
    aid = re.sub(r"^arxiv:", "", (arxiv_id or "").strip(), flags=re.I)
    if not aid:
        return None, {"provider": "error", "source": "arxiv", "message": "empty arxiv id"}
    url = "https://export.arxiv.org/api/query?" + urllib.parse.urlencode({"id_list": aid, "max_results": "1"})
    try:
        root = ET.fromstring(_get(url))
    except (urllib.error.URLError, TimeoutError, ET.ParseError, OSError) as e:
        return None, {"provider": "error", "source": "arxiv", "message": str(e)}
    for entry in root.findall("atom:entry", ATOM_NS):
        h = _hit_from_arxiv_entry(entry)
        if h and h.get("title"):
            return h, {"provider": "http", "source": "arxiv", "via": "arxiv_id"}
    return None, {"provider": "http", "source": "arxiv", "message": "no hit"}


def _norm_title(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", unicodedata.normalize("NFKD", text or "").lower()).strip()


def title_similarity(a: str, b: str) -> float:
    return SequenceMatcher(None, _norm_title(a), _norm_title(b)).ratio()


TITLE_MATCH_MIN = 0.95


def find_by_title(title: str) -> tuple[Optional[dict[str, Any]], dict[str, Any]]:
    """Find a paper by its exact title in Crossref. Only a near-identical title is accepted."""
    t = re.sub(r"\s+", " ", title or "").strip()
    if len(t.split()) < 4:
        return None, {"provider": "error", "source": "crossref", "message": "title too short to match"}
    try:
        body = _get_json(_crossref_url("works", {"query.bibliographic": t, "rows": "10"}))
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
        return None, {"provider": "error", "source": "crossref", "message": str(e)}
    best, best_score, best_cited = None, 0.0, -1
    for row in ((body.get("message") or {}).get("items") or []):
        row = row if isinstance(row, dict) else {}
        h = _hit_from_crossref(row)
        if not h or not h.get("title"):
            continue
        score = title_similarity(t, h["title"])
        cited = int(row.get("is-referenced-by-count") or 0)
        # Conference abstracts often reuse a paper's exact title: on a (near) tie, the
        # version that others cite is the article.
        if score > best_score + 0.01 or (abs(score - best_score) <= 0.01 and cited > best_cited):
            best, best_score, best_cited = h, max(score, best_score), cited
    if best and best_score >= TITLE_MATCH_MIN:
        return best, {"provider": "http", "source": "crossref", "via": "title_crossref", "title_similarity": round(best_score, 3)}
    return None, {"provider": "http", "source": "crossref", "message": "no near-identical title", "best_similarity": round(best_score, 3)}


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
    # The index API moved to api.opencitations.net and takes a bare DOI (no "doi:" prefix).
    url = f"https://api.opencitations.net/index/v1/{path}/{urllib.parse.quote(d)}"
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
