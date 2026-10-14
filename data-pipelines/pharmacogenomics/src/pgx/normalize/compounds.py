"""Compound name normalization and PubChem / ChEMBL resolution.

Salt forms and stereo parentheticals are removed only from the lookup string.
The source name stored on the record is unchanged. A lookup that returns more
than one CID is ambiguous and stops; a looser name is not tried after that.
"""

from __future__ import annotations

import re
import unicodedata
from collections.abc import Iterable
from typing import Protocol

from pgx.ids import CompoundRecord, apply_pubchem

SALT_WORDS = (
    "dihydrochloride",
    "dihydrobromide",
    "hydrochloride",
    "hydrobromide",
    "methanesulfonate",
    "methochloride",
    "monohydrate",
    "hemihydrate",
    "dihydrate",
    "disodium",
    "mesylate",
    "mesilate",
    "maleate",
    "sulfate",
    "sulphate",
    "acetate",
    "citrate",
    "phosphate",
    "tartrate",
    "fumarate",
    "besylate",
    "besilate",
    "tosylate",
    "camsylate",
    "succinate",
    "benzoate",
    "oxalate",
    "lactate",
    "nitrate",
    "bromide",
    "chloride",
    "iodide",
    "sodium",
    "potassium",
    "calcium",
    "hydrate",
    "2hcl",
    "2hbr",
    "hcl",
    "hbr",
    "salt",
)
_SALT = "(?:" + "|".join(SALT_WORDS) + ")"
_SALT_TRAIL = re.compile(rf"(?:[\s\-]+|\s*\(){_SALT}\)?\s*$", re.IGNORECASE)
_STEREO = re.compile(
    r"(?:\s*-\s*)?\(\s*(?:[+-]|±|\+/-|R|S|RS|rac|E|Z|cis|trans)\s*\)",
    re.IGNORECASE,
)
_BRACKET = re.compile(r"\s*\[[^\]]*\]\s*$")


class CompoundClient(Protocol):
    def cids(self, kind: str, query: str) -> list[str]: ...

    def properties_for(self, cids: list[str]) -> dict[str, dict[str, str | None]]: ...

    def chembl_for(self, inchikeys: list[str]) -> dict[str, list[str]]: ...


def normalize_compound_name(name: str) -> str:
    text = unicodedata.normalize("NFKC", name)
    text = text.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-").replace("\u00a0", " ")
    return re.sub(r"\s+", " ", text).strip()


def strip_stereo(name: str) -> str:
    text = name
    previous = None
    while previous != text:
        previous = text
        text = _STEREO.sub(" ", text)
        text = re.sub(r"\s+", " ", text).strip()
    bracket = _BRACKET.search(text)
    if bracket and re.search(r"[RS]", bracket.group(0)):
        text = text[: bracket.start()]
    return text.strip(" -")


def strip_salts(name: str) -> str:
    text = name.strip()
    while True:
        match = _SALT_TRAIL.search(text)
        if not match:
            return text
        nxt = text[: match.start()].strip(" -")
        if len(re.sub(r"[^A-Za-z0-9]", "", nxt)) < 3:
            return text
        text = nxt


def name_variants(original: str) -> list[tuple[str, str, str]]:
    """Lookup stages for one source string: `(stage, query, confidence)`.

    Stages are `exact`, `stereo`, and `salt`. The original string is not modified.
    """
    base = normalize_compound_name(original)
    if not base:
        return []
    stereo = strip_stereo(base)
    salted = strip_salts(stereo)
    variants = [("exact", base, "medium")]
    if stereo and stereo.casefold() != base.casefold():
        variants.append(("stereo", stereo, "low"))
    if salted and salted.casefold() not in {base.casefold(), stereo.casefold()}:
        variants.append(("salt", salted, "low"))
    return variants


def iter_queries(record: CompoundRecord) -> list[tuple[str, str, str, str]]:
    """Ordered `(kind, query, method, confidence)` lookups for one compound."""
    queries: list[tuple[str, str, str, str]] = []
    seen: set[tuple[str, str]] = set()

    def add(kind: str, query: str, method: str, confidence: str) -> None:
        key = (kind, query.casefold())
        if not query or key in seen:
            return
        seen.add(key)
        queries.append((kind, query, method, confidence))

    if record.inchikey:
        add("inchikey", record.inchikey, "pubchem_inchikey", "high")
    if record.smiles:
        add("smiles", record.smiles, "pubchem_smiles", "high")
    if record.name:
        _add_name(add, record.name, "pubchem_name", "pubchem_name_stereo", "pubchem_name_salt")
    for synonym in sorted(record.synonyms):
        _add_name(add, synonym, "pubchem_synonym", "pubchem_synonym_stereo", "pubchem_synonym_salt")
    if record.brd_id:
        add("name", record.brd_id, "pubchem_brd", "medium")
    return queries


def _add_name(add, original: str, exact: str, stereo: str, salt: str) -> None:
    methods = {"exact": exact, "stereo": stereo, "salt": salt}
    for stage, query, confidence in name_variants(original):
        add("name", query, methods[stage], confidence)


def resolve_compounds(
    records: list[CompoundRecord],
    client: CompoundClient,
    limit: int | None = None,
) -> int:
    """Resolve unmatched compounds. Returns how many received a single CID."""
    pending: list[tuple[CompoundRecord, list[tuple[str, str, str, str]]]] = []
    for record in records:
        if record.status in {"mapped", "ambiguous"}:
            continue
        queries = iter_queries(record)
        if queries:
            pending.append((record, queries))
        if limit is not None and len(pending) >= limit:
            break
    active = [(record, queries, 0) for record, queries in pending]
    mapped = 0
    while active:
        if hasattr(client, "prefetch"):
            client.prefetch([(queries[index][0], queries[index][1]) for _record, queries, index in active])
        nxt: list[tuple[CompoundRecord, list[tuple[str, str, str, str]], int]] = []
        for record, queries, index in active:
            kind, query, method, confidence = queries[index]
            cids = client.cids(kind, query)
            if len(cids) > 1:
                apply_pubchem(record, None, None, None, method, f"{query} matched {len(cids)} CIDs")
                record.status = "ambiguous"
                record.confidence = None
            elif len(cids) == 1:
                apply_pubchem(record, cids[0], record.inchikey, record.chembl_id, method, f"{method} {query}", confidence)
                if record.status == "mapped":
                    mapped += 1
            elif index + 1 < len(queries):
                nxt.append((record, queries, index + 1))
            else:
                apply_pubchem(record, None, None, None, "pubchem_miss", "no unique PubChem CID")
        active = nxt
    _enrich(records, client)
    return mapped


def _enrich(records: Iterable[CompoundRecord], client: CompoundClient) -> None:
    wanted = [record.pubchem_cid for record in records if record.status == "mapped" and record.pubchem_cid]
    props = client.properties_for([cid for cid in wanted if cid]) if wanted else {}
    for record in records:
        if record.status != "mapped" or not record.pubchem_cid:
            continue
        row = props.get(record.pubchem_cid) or {}
        if row.get("inchikey") and not record.inchikey:
            record.inchikey = row["inchikey"]
    keys = [record.inchikey for record in records if record.status == "mapped" and record.inchikey and not record.chembl_id]
    chembl = client.chembl_for([key for key in keys if key]) if keys else {}
    for record in records:
        if record.status != "mapped" or not record.inchikey or record.chembl_id:
            continue
        found = chembl.get(record.inchikey) or []
        if len(found) == 1:
            record.chembl_id = found[0]
