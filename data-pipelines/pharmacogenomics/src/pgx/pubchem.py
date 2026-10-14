"""PubChem PUG-REST and ChEMBL lookups with a resumable cache.

Name lists are not posted together: a multi-name PUG body returns one CID list
and drops the name-to-CID pairing. Each query is cached on its own. Requests
stay under PubChem's published rate.
"""

from __future__ import annotations

import hashlib
import json
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

from pgx.download import USER_AGENT, DownloadError
from pgx.ids import CompoundRecord

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"


def _cache_path(cache_dir: Path, kind: str, query: str) -> Path:
    digest = hashlib.sha256(f"{kind}\n{query}".encode()).hexdigest()
    return cache_dir / kind / f"{digest}.json"


class PubChemClient:
    def __init__(self, cache_dir: Path, pause_s: float = 0.25) -> None:
        self.cache_dir = cache_dir
        self.pause_s = pause_s
        self._next = 0.0
        self._lock = threading.Lock()
        cache_dir.mkdir(parents=True, exist_ok=True)

    def _read(self, kind: str, query: str) -> dict | None:
        path = _cache_path(self.cache_dir, kind, query)
        if not path.exists():
            return None
        try:
            return json.loads(path.read_text())
        except json.JSONDecodeError:
            return None

    def _write(self, kind: str, query: str, payload: dict) -> None:
        path = _cache_path(self.cache_dir, kind, query)
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(payload))
        temporary.replace(path)

    def _reserve(self) -> None:
        if self.pause_s <= 0:
            return
        with self._lock:
            now = time.monotonic()
            delay = max(0.0, self._next - now)
            self._next = max(now, self._next) + self.pause_s
        if delay:
            time.sleep(delay)

    def _request(self, url: str, data: bytes | None = None, accept: str = "application/json") -> dict | None:
        last_error: Exception | None = None
        for attempt in range(5):
            self._reserve()
            request = Request(url, data=data, headers={"User-Agent": USER_AGENT, "Accept": accept})
            if data is not None:
                request.add_header("Content-Type", "application/x-www-form-urlencoded")
            try:
                with urlopen(request, timeout=90) as response:
                    return json.loads(response.read().decode())
            except HTTPError as exc:
                if exc.code in {400, 404}:
                    return None
                last_error = exc
                if exc.code not in {429, 500, 503}:
                    raise DownloadError(f"PubChem/ChEMBL {url} returned HTTP {exc.code}") from exc
            except URLError as exc:
                last_error = exc
            time.sleep(min(60, 2 ** attempt))
        raise DownloadError(f"PubChem/ChEMBL {url} failed after retries: {last_error}")

    def _cids_from(self, payload: dict | None) -> list[str]:
        if not payload:
            return []
        values = payload.get("IdentifierList", {}).get("CID") or []
        return [str(value) for value in values]

    def _fetch_cids(self, kind: str, query: str) -> list[str]:
        if kind == "inchikey":
            return self._cids_from(self._request(f"{PUG}/compound/inchikey/{quote(query)}/cids/JSON"))
        if kind == "smiles":
            body = urlencode({"smiles": query}).encode()
            return self._cids_from(self._request(f"{PUG}/compound/smiles/cids/JSON", data=body))
        return self._cids_from(self._request(f"{PUG}/compound/name/{quote(query, safe='')}/cids/JSON"))

    def _fetch_and_store(self, kind: str, query: str) -> list[str]:
        existing = self._read(kind, query)
        if existing is not None and "cids" in existing:
            return [str(cid) for cid in existing["cids"]]
        cids = self._fetch_cids(kind, query)
        self._write(kind, query, {"kind": kind, "query": query, "cids": cids})
        return cids

    def cids(self, kind: str, query: str) -> list[str]:
        existing = self._read(kind, query)
        if existing is not None and "cids" in existing:
            return [str(cid) for cid in existing["cids"]]
        return self._fetch_and_store(kind, query)

    def prefetch(self, queries: list[tuple[str, str]]) -> None:
        missing: list[tuple[str, str]] = []
        seen: set[tuple[str, str]] = set()
        for kind, query in queries:
            if (kind, query) in seen:
                continue
            seen.add((kind, query))
            existing = self._read(kind, query)
            if existing is None or "cids" not in existing:
                missing.append((kind, query))
        if not missing:
            return
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = [pool.submit(self._fetch_and_store, kind, query) for kind, query in missing]
            total = len(futures)
            for done, future in enumerate(as_completed(futures), start=1):
                future.result()
                if done == total or done % 200 == 0:
                    print(f"pubchem {done}/{total}", file=sys.stderr, flush=True)

    def properties_for(self, cids: list[str]) -> dict[str, dict[str, str | None]]:
        wanted = list(dict.fromkeys(cids))
        found: dict[str, dict[str, str | None]] = {}
        missing: list[str] = []
        for cid in wanted:
            cached = self._read("cid", cid)
            if cached is not None and "inchikey" in cached:
                found[cid] = {"inchikey": cached.get("inchikey"), "title": cached.get("title")}
            else:
                missing.append(cid)
        for start in range(0, len(missing), 80):
            batch = missing[start:start + 80]
            payload = self._request(f"{PUG}/compound/cid/{','.join(batch)}/property/InChIKey,Title/JSON") or {}
            rows = payload.get("PropertyTable", {}).get("Properties") or []
            returned: set[str] = set()
            for row in rows:
                cid = str(row.get("CID"))
                returned.add(cid)
                stored = {"inchikey": row.get("InChIKey"), "title": row.get("Title")}
                self._write("cid", cid, stored)
                found[cid] = stored
            for cid in batch:
                if cid not in returned:
                    stored = {"inchikey": None, "title": None}
                    self._write("cid", cid, stored)
                    found[cid] = stored
        return found

    def chembl_for(self, inchikeys: list[str]) -> dict[str, list[str]]:
        wanted = list(dict.fromkeys(inchikeys))
        found: dict[str, list[str]] = {}
        missing: list[str] = []
        for key in wanted:
            cached = self._read("chembl", key)
            if cached is not None and "chembl_ids" in cached:
                found[key] = [str(item) for item in cached["chembl_ids"]]
            else:
                missing.append(key)
        for start in range(0, len(missing), 40):
            batch = missing[start:start + 40]
            query = ",".join(batch)
            url = f"{CHEMBL}?molecule_structures__standard_inchi_key__in={quote(query, safe=',')}&limit=200"
            payload = self._request(url) or {}
            grouped: dict[str, list[str]] = {key: [] for key in batch}
            for molecule in payload.get("molecules") or []:
                structures = molecule.get("molecule_structures") or {}
                key = structures.get("standard_inchi_key")
                chembl_id = molecule.get("molecule_chembl_id")
                if key in grouped and chembl_id and chembl_id not in grouped[key]:
                    grouped[key].append(chembl_id)
            for key, chembl_ids in grouped.items():
                self._write("chembl", key, {"chembl_ids": chembl_ids})
                found[key] = chembl_ids
        return found


def resolve_records(
    records: list[CompoundRecord],
    cache_dir: Path,
    pause_s: float = 0.25,
    limit: int | None = None,
) -> int:
    from pgx.normalize.compounds import resolve_compounds

    client = PubChemClient(cache_dir, pause_s=pause_s)
    return resolve_compounds(records, client, limit=limit)
