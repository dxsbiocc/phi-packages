"""Optional PubChem resolution. Off unless --resolve-compounds is passed."""

from __future__ import annotations

import json
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import Request, urlopen

from pgx.download import USER_AGENT, DownloadError
from pgx.ids import CHEMBL_RE, CompoundRecord, apply_pubchem

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def _get(url: str) -> dict | None:
    request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    try:
        with urlopen(request, timeout=60) as response:
            return json.loads(response.read().decode())
    except HTTPError as exc:
        if exc.code == 404:
            return None
        raise DownloadError(f"PubChem {url} returned HTTP {exc.code}") from exc


def _cids(payload: dict | None) -> list[str]:
    if not payload:
        return []
    values = payload.get("IdentifierList", {}).get("CID") or []
    return [str(value) for value in values]


def lookup_name(name: str) -> list[str]:
    return _cids(_get(f"{PUG}/compound/name/{quote(name)}/cids/JSON"))


def properties(cid: str) -> tuple[str | None, str | None, str | None]:
    prop = _get(f"{PUG}/compound/cid/{cid}/property/InChIKey,Title/JSON") or {}
    rows = prop.get("PropertyTable", {}).get("Properties") or [{}]
    inchikey = rows[0].get("InChIKey")
    synonyms = _get(f"{PUG}/compound/cid/{cid}/synonyms/JSON") or {}
    info = (synonyms.get("InformationList", {}).get("Information") or [{}])[0]
    chembl = [item for item in info.get("Synonym") or [] if CHEMBL_RE.match(item)]
    chembl_id = chembl[0] if len(chembl) == 1 else None
    return inchikey, chembl_id, rows[0].get("Title")


def resolve_records(records: list[CompoundRecord], cache_dir: Path, pause_s: float = 0.25, limit: int | None = None) -> int:
    cache_dir.mkdir(parents=True, exist_ok=True)
    resolved = 0
    for record in records:
        if limit is not None and resolved >= limit:
            break
        if record.status == "mapped" or record.status == "ambiguous":
            continue
        query = record.name
        if not query:
            continue
        cache = cache_dir / f"name-{quote(query, safe='')}.json"
        if cache.exists():
            cached = json.loads(cache.read_text())
        else:
            cids = lookup_name(query)
            cached = {"cids": cids}
            if len(cids) == 1:
                inchikey, chembl_id, _title = properties(cids[0])
                cached.update({"inchikey": inchikey, "chembl_id": chembl_id})
            cache.write_text(json.dumps(cached))
            time.sleep(pause_s)
        cids = cached.get("cids") or []
        if len(cids) == 1:
            apply_pubchem(record, cids[0], cached.get("inchikey"), cached.get("chembl_id"), "pubchem_name_unique", f"PubChem name {query}")
            resolved += 1
        elif len(cids) > 1:
            apply_pubchem(record, None, None, None, "pubchem_name_ambiguous", f"{query} matched {len(cids)} CIDs")
        else:
            apply_pubchem(record, None, None, None, "pubchem_name_miss", f"PubChem has no CID for {query}")
    return resolved
