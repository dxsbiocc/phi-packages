"""Canonical cell-line names and Cellosaurus accession matching.

`NCI-H460` and `NCIH460` share a canonical key. `H460` does not: the NCI
prefix is not stripped. A short name matches only when Cellosaurus lists it
as a name or synonym. DepMap `CCLE_Name` tissue suffixes such as `_LUNG` are
removed for a second lookup key and the original string is kept.
"""

from __future__ import annotations

import re
import unicodedata
from dataclasses import dataclass, field
from pathlib import Path

CVCL_RE = re.compile(r"^CVCL_[A-Z0-9]{4}$")

# DepMap CCLE_Name appends an Oncotree lineage. Longest suffix wins.
CCLE_LINEAGES = (
    "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
    "UPPER_AERODIGESTIVE_TRACT",
    "CENTRAL_NERVOUS_SYSTEM",
    "AUTONOMIC_GANGLIA",
    "LARGE_INTESTINE",
    "SMALL_INTESTINE",
    "SALIVARY_GLAND",
    "ADRENAL_CORTEX",
    "BILIARY_TRACT",
    "URINARY_TRACT",
    "SOFT_TISSUE",
    "HEAD_AND_NECK",
    "ADRENAL_GLAND",
    "PLASMA_CELL",
    "ENDOMETRIUM",
    "FIBROBLAST",
    "ENGINEERED",
    "LYMPHOCYTE",
    "ESOPHAGUS",
    "PANCREAS",
    "STOMACH",
    "THYROID",
    "PROSTATE",
    "BREAST",
    "OVARY",
    "KIDNEY",
    "LIVER",
    "CERVIX",
    "PLEURA",
    "VULVA",
    "COLON",
    "BLOOD",
    "BONE",
    "SKIN",
    "LUNG",
    "EYE",
    "UNKNOWN",
    "LYMPHOID",
    "MYELOID",
    "EMBRYONIC",
    "TESTIS",
    "BLADDER",
)

IDENTIFIER_METHODS = {
    "depmap_rrid",
    "cmp_broad_id",
    "cmp_sidm",
    "cmp_cosmic",
    "cmp_rrid",
    "depmap_ach",
}

_LINEAGE_TAILS = tuple("_" + suffix for suffix in sorted(CCLE_LINEAGES, key=len, reverse=True))


def canonical_name(value: str) -> str:
    text = unicodedata.normalize("NFKC", value).casefold()
    return re.sub(r"[^a-z0-9]", "", text)


def ccle_lookup_keys(name: str) -> list[tuple[str, str]]:
    """Return `(canonical key, kind)` pairs. Kind is `full` or `stem`."""
    full = canonical_name(name)
    keys: list[tuple[str, str]] = []
    if full:
        keys.append((full, "full"))
    upper = name.strip().upper()
    for tail in _LINEAGE_TAILS:
        if upper.endswith(tail) and len(upper) > len(tail):
            stem = canonical_name(upper[: -len(tail)])
            if stem and stem != full:
                keys.append((stem, "stem"))
            break
    return keys


@dataclass
class CellLineRecord:
    accession: str
    name: str
    problem: str | None = None


@dataclass
class CellosaurusIndex:
    by_ac: dict[str, CellLineRecord] = field(default_factory=dict)
    alias_ac: dict[str, str] = field(default_factory=dict)
    by_key: dict[str, set[str]] = field(default_factory=dict)
    official: dict[str, set[str]] = field(default_factory=dict)

    def _add_key(self, name: str, accession: str, official: bool) -> None:
        key = canonical_name(name)
        if len(key) < 2:
            return
        self.by_key.setdefault(key, set()).add(accession)
        if official:
            self.official.setdefault(key, set()).add(accession)

    def primary(self, accession: str) -> str | None:
        if accession in self.by_ac:
            return accession
        return self.alias_ac.get(accession)

    def flag(self, accession: str) -> str | None:
        primary = self.primary(accession)
        if not primary:
            return None
        return self.by_ac[primary].problem

    @classmethod
    def load(cls, path: Path) -> CellosaurusIndex:
        index = cls()
        current: dict[str, list[str]] = {}

        def finish() -> None:
            accession = next((item for item in current.get("AC", []) if CVCL_RE.match(item)), "")
            name = current.get("ID", [""])[0].strip()
            if not accession or not name:
                current.clear()
                return
            problem = _problem_flag(current.get("CC", []))
            index.by_ac[accession] = CellLineRecord(accession, name, problem)
            index._add_key(name, accession, official=True)
            for synonym in current.get("SY", []):
                for part in synonym.split(";"):
                    token = part.strip()
                    if token:
                        index._add_key(token, accession, official=False)
            for secondary in current.get("AS", []):
                if CVCL_RE.match(secondary) and secondary != accession:
                    index.alias_ac[secondary] = accession
            current.clear()

        with path.open(errors="replace") as handle:
            for line in handle:
                if line.startswith("//"):
                    finish()
                    continue
                if len(line) < 5 or line[2:5] != "   ":
                    continue
                tag = line[:2]
                value = line[5:].rstrip("\n")
                if tag in {"ID", "AC", "SY", "AS", "CC"}:
                    current.setdefault(tag, []).append(value)
            if current:
                finish()
        return index

    def decide_names(self, names: list[object]) -> object:
        from pgx.ids import Decision, clean

        found: list[tuple[str, str]] = []
        ambiguous = False
        for raw in names:
            token = clean(raw)
            if not token:
                continue
            for key, kind in ccle_lookup_keys(token):
                accessions = self.by_key.get(key)
                if not accessions:
                    continue
                if len(accessions) != 1:
                    ambiguous = True
                    continue
                accession = next(iter(accessions))
                found.append((accession, self._method(key, accession, kind)))
        distinct = {accession for accession, _method in found}
        if ambiguous or len(distinct) > 1:
            return Decision(None, "cellosaurus_name_ambiguous", "ambiguous", "name matched more than one Cellosaurus accession")
        if len(distinct) == 1:
            accession = next(iter(distinct))
            methods = {method for cvcl, method in found if cvcl == accession}
            if "cellosaurus_name" in methods:
                method = "cellosaurus_name"
            elif "cellosaurus_synonym" in methods:
                method = "cellosaurus_synonym"
            else:
                method = "cellosaurus_ccle_stem"
            return Decision(accession, method, "mapped", f"Cellosaurus {accession}", self.flag(accession))
        return Decision(None, "unmatched", "unmatched", "no Cellosaurus name or synonym")

    def _method(self, key: str, accession: str, kind: str) -> str:
        if kind == "stem":
            return "cellosaurus_ccle_stem"
        if accession in self.official.get(key, set()):
            return "cellosaurus_name"
        return "cellosaurus_synonym"

    def harmonize(self, decision: object, names: list[object], identifier_names: list[object] | None = None) -> object:
        """Apply Cellosaurus after an identifier decision. Never auto-resolve ambiguity."""
        from pgx.ids import Decision

        assert isinstance(decision, Decision)
        if decision.match_method in {"conflict_identifier", "unmatched_identifier"}:
            return decision
        check_names = names if identifier_names is None else identifier_names
        if decision.match_method in IDENTIFIER_METHODS and decision.status == "mapped" and decision.cvcl:
            primary = self.primary(decision.cvcl)
            if primary is None:
                return Decision(None, "unmatched_identifier", "unmatched", f"{decision.cvcl} is not in Cellosaurus")
            name_decision = self.decide_names(check_names)
            assert isinstance(name_decision, Decision)
            if name_decision.status == "mapped" and name_decision.cvcl != primary:
                return Decision(
                    None,
                    "conflict_identifier",
                    "ambiguous",
                    f"{decision.cvcl} disagrees with Cellosaurus name {name_decision.cvcl}",
                )
            evidence = decision.evidence
            if primary != decision.cvcl:
                evidence = f"{decision.evidence}; secondary {decision.cvcl} -> {primary}"
            return Decision(primary, decision.match_method, "mapped", evidence, self.flag(primary))
        return self.decide_names(names)


def _problem_flag(comments: list[str]) -> str | None:
    text = " ".join(comments)
    if "Problematic cell line" not in text and "Misidentified" not in text and "Contaminated" not in text:
        return None
    if re.search(r"Misidentified", text, re.IGNORECASE):
        return "misidentified"
    if re.search(r"Contaminated", text, re.IGNORECASE):
        return "contaminated"
    if "Problematic cell line" in text:
        return "problematic"
    return None
