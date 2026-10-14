"""Identifier harmonization.

Cell lines are keyed by a Cellosaurus CVCL accession. Compounds are keyed by
a PubChem CID. A crosswalk is kept only when one identifier says so. Conflicts
and unknown aliases stay unresolved.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

CVCL_RE = re.compile(r"^CVCL_[A-Z0-9]{4}$")
GENE_RE = re.compile(r"^(?P<symbol>.+?) \((?P<entrez>\d+)\)$")
SIDM_RE = re.compile(r"^SIDM\d+$", re.IGNORECASE)
ACH_RE = re.compile(r"^ACH-\d+$", re.IGNORECASE)
CHEMBL_RE = re.compile(r"^CHEMBL\d+$")
EMPTY = {"", "NA", "N/A", "NONE", "NULL", "NAN", "-"}

SILENT_CONSEQUENCES = {
    "synonymous_variant",
    "intron_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "5_prime_utr_variant",
    "3_prime_utr_variant",
    "non_coding_transcript_exon_variant",
    "intergenic_variant",
}


def clean(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if text.upper() in EMPTY:
        return None
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    return text


def norm_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", value).upper()


def parse_gene_label(label: str) -> tuple[str, str] | None:
    match = GENE_RE.match(label.strip())
    if not match:
        return None
    return match.group("symbol"), match.group("entrez")


def cvcl_values(value: object) -> list[str]:
    text = clean(value)
    if not text:
        return []
    return [part.strip() for part in re.split(r"[;|]", text) if CVCL_RE.match(part.strip())]


def parse_cvcl(value: object) -> str | None:
    found = cvcl_values(value)
    if len(found) == 1:
        return found[0]
    return None


def parse_sidm(value: object) -> str | None:
    text = clean(value)
    if text and SIDM_RE.match(text):
        return text.upper()
    return None


def parse_ach(value: object) -> str | None:
    text = clean(value)
    if text and ACH_RE.match(text):
        return text.upper()
    return None


def parse_cosmic(value: object) -> str | None:
    text = clean(value)
    if text and re.fullmatch(r"\d+", text):
        return text
    return None


def parse_brd(value: object) -> str | None:
    text = clean(value)
    if not text:
        return None
    text = text.removeprefix("BRD:")
    text = re.sub(r"\s+-\s+QC Failure$", "", text, flags=re.IGNORECASE).strip()
    if text.upper().startswith("BRD-"):
        return text
    return None


@dataclass
class Decision:
    cvcl: str | None
    match_method: str
    status: str
    evidence: str


@dataclass
class CmpIndex:
    by_sidm: dict[str, set[str]] = field(default_factory=dict)
    by_ach: dict[str, set[str]] = field(default_factory=dict)
    by_cosmic: dict[str, set[str]] = field(default_factory=dict)
    by_name: dict[str, set[str]] = field(default_factory=dict)
    display: dict[str, dict[str, str | None]] = field(default_factory=dict)

    def add(self, bucket: dict[str, set[str]], key: str | None, cvcl: str) -> None:
        if key:
            bucket.setdefault(key, set()).add(cvcl)


def _lookup(bucket: dict[str, set[str]], key: str | None) -> tuple[str | None, str]:
    if not key:
        return None, "missing"
    found = bucket.get(key)
    if not found:
        return None, "missing"
    if len(found) == 1:
        return next(iter(found)), "unique"
    return None, "ambiguous"


def index_cmp(rows: list[dict[str, str]]) -> CmpIndex:
    index = CmpIndex()
    for row in rows:
        found = cvcl_values(row.get("RRID"))
        sidm = parse_sidm(row.get("model_id"))
        if not sidm or not found:
            continue
        for cvcl in found:
            index.add(index.by_sidm, sidm, cvcl)
        if len(found) != 1:
            continue
        cvcl = found[0]
        index.add(index.by_ach, parse_ach(row.get("BROAD_ID")), cvcl)
        index.add(index.by_cosmic, parse_cosmic(row.get("COSMIC_ID")), cvcl)
        names = [row.get("model_name") or ""]
        synonym = clean(row.get("synonyms"))
        if synonym:
            names.extend(part.strip() for part in synonym.split(";"))
        for name in names:
            token = clean(name)
            if token and len(norm_name(token)) >= 3:
                index.add(index.by_name, norm_name(token), cvcl)
        species = clean(row.get("species"))
        index.display[cvcl] = {
            "display_name": clean(row.get("model_name")),
            "tissue": clean(row.get("tissue")),
            "disease": clean(row.get("cancer_type")),
            "species": species,
        }
    return index


def _identifier_decision(
    rrid: str | None,
    ach: tuple[str | None, str],
    sidm: tuple[str | None, str],
    cosmic: tuple[str | None, str],
) -> Decision | None:
    known: list[tuple[str, str]] = []
    if rrid:
        known.append((rrid, "depmap_rrid"))
    if ach[1] == "unique" and ach[0]:
        known.append((ach[0], "cmp_broad_id"))
    if sidm[1] == "unique" and sidm[0]:
        known.append((sidm[0], "cmp_sidm"))
    if cosmic[1] == "unique" and cosmic[0]:
        known.append((cosmic[0], "cmp_cosmic"))
    ambiguous = [label for label, state in (("ACH", ach[1]), ("SIDM", sidm[1]), ("COSMIC", cosmic[1])) if state == "ambiguous"]
    distinct = {cvcl for cvcl, _method in known}
    if ambiguous or len(distinct) > 1:
        bits = []
        if distinct:
            bits.append("values " + ",".join(sorted(distinct)))
        if ambiguous:
            bits.append("ambiguous " + ",".join(ambiguous))
        evidence = "conflicting identifiers " + "; ".join(bits)
        return Decision(None, "conflict_identifier", "ambiguous", evidence)
    if len(distinct) == 1:
        cvcl, method = known[0]
        return Decision(cvcl, method, "mapped", method)
    return None


def classify_depmap(row: dict[str, str], cmp: CmpIndex) -> Decision:
    rrid = parse_cvcl(row.get("RRID"))
    sidm = parse_sidm(row.get("SangerModelID"))
    cosmic = parse_cosmic(row.get("COSMICID"))
    decision = _identifier_decision(
        rrid,
        _lookup(cmp.by_ach, parse_ach(row.get("ModelID"))),
        _lookup(cmp.by_sidm, sidm),
        _lookup(cmp.by_cosmic, cosmic),
    )
    if decision:
        return decision
    if rrid or sidm or cosmic:
        return Decision(None, "unmatched_identifier", "unmatched", "stable id did not join to one CVCL")
    return _name_decision(row.get("CellLineName"), row.get("StrippedCellLineName"), cmp)


def _name_decision(primary: object, secondary: object, cmp: CmpIndex) -> Decision:
    hits: set[str] = set()
    saw_ambiguous = False
    for raw in (primary, secondary):
        token = clean(raw)
        if not token or len(norm_name(token)) < 3:
            continue
        cvcl, state = _lookup(cmp.by_name, norm_name(token))
        if state == "ambiguous":
            saw_ambiguous = True
        elif state == "unique" and cvcl:
            hits.add(cvcl)
    if saw_ambiguous or len(hits) > 1:
        return Decision(None, "name_ambiguous", "ambiguous", "name matched more than one CVCL")
    if len(hits) == 1:
        return Decision(next(iter(hits)), "name_unique", "mapped", "unique normalized name")
    return Decision(None, "unmatched", "unmatched", "no identifier or unique name")


def classify_gdsc(sidm: str | None, cosmic: str | None, name: str | None, cmp: CmpIndex) -> Decision:
    decision = _identifier_decision(
        None,
        (None, "missing"),
        _lookup(cmp.by_sidm, sidm),
        _lookup(cmp.by_cosmic, cosmic),
    )
    if decision:
        return decision
    if sidm or cosmic:
        return Decision(None, "unmatched_identifier", "unmatched", "GDSC id did not join to one CVCL")
    return _name_decision(name, None, cmp)


def classify_ach(ach: str, depmap: dict[str, Decision], cmp: CmpIndex) -> Decision:
    if ach in depmap:
        prior = depmap[ach]
        if prior.status == "mapped" and prior.cvcl:
            return Decision(prior.cvcl, "depmap_ach", "mapped", f"DepMap {ach}")
        return Decision(None, prior.match_method, prior.status, prior.evidence)
    cvcl, state = _lookup(cmp.by_ach, ach)
    if state == "unique" and cvcl:
        return Decision(cvcl, "cmp_broad_id", "mapped", f"CMP BROAD_ID {ach}")
    if state == "ambiguous":
        return Decision(None, "conflict_identifier", "ambiguous", f"CMP BROAD_ID {ach} is ambiguous")
    return Decision(None, "unmatched", "unmatched", f"{ach} is absent from DepMap and CMP")


def xref_rows(
    source: str,
    source_model_id: str,
    decision: Decision,
    aliases: list[tuple[str, str]],
    display_name: str | None = None,
    tissue: str | None = None,
    disease: str | None = None,
    species: str | None = None,
) -> list[dict[str, str | None]]:
    rows: list[dict[str, str | None]] = []
    seen: set[tuple[str, str]] = set()
    items = [("primary", source_model_id), *aliases]
    for alias_type, alias in items:
        token = clean(alias)
        if not token:
            continue
        key = (alias_type, token)
        if key in seen:
            continue
        seen.add(key)
        rows.append({
            "source": source,
            "source_model_id": source_model_id,
            "alias": token,
            "alias_type": alias_type,
            "cvcl": decision.cvcl,
            "match_method": decision.match_method,
            "status": decision.status,
            "evidence": decision.evidence,
            "display_name": display_name,
            "tissue": tissue,
            "disease": disease,
            "species": species,
        })
    return rows


def build_model_xref(
    cmp_rows: list[dict[str, str]],
    depmap_rows: list[dict[str, str]],
    gdsc_rows: list[dict[str, str]] | None = None,
    prism_rows: list[dict[str, str]] | None = None,
) -> list[dict[str, str | None]]:
    cmp = index_cmp(cmp_rows)
    out: list[dict[str, str | None]] = []
    depmap_decisions: dict[str, Decision] = {}
    for row in depmap_rows:
        ach = parse_ach(row.get("ModelID"))
        if not ach:
            continue
        decision = classify_depmap(row, cmp)
        depmap_decisions[ach] = decision
        display = cmp.display.get(decision.cvcl or "", {})
        out.extend(xref_rows(
            "depmap",
            ach,
            decision,
            [
                ("rrid", row.get("RRID") or ""),
                ("sidm", row.get("SangerModelID") or ""),
                ("cosmic", row.get("COSMICID") or ""),
                ("cell_line_name", row.get("CellLineName") or ""),
                ("stripped_name", row.get("StrippedCellLineName") or ""),
                ("ccle_name", row.get("CCLEName") or ""),
            ],
            display_name=clean(row.get("CellLineName")) or display.get("display_name"),
            tissue=clean(row.get("OncotreeLineage")) or display.get("tissue"),
            disease=clean(row.get("OncotreePrimaryDisease")) or display.get("disease"),
            species=display.get("species") or "Homo sapiens",
        ))
    for row in cmp_rows:
        sidm = parse_sidm(row.get("model_id"))
        cvcl = parse_cvcl(row.get("RRID"))
        if not sidm:
            continue
        _found, state = _lookup(cmp.by_sidm, sidm)
        if state == "ambiguous":
            decision = Decision(None, "conflict_identifier", "ambiguous", f"{sidm} has more than one RRID")
        elif cvcl:
            decision = Decision(cvcl, "cmp_rrid", "mapped", "Cell Model Passports RRID")
        else:
            decision = Decision(None, "unmatched", "unmatched", "CMP model has no RRID")
        out.extend(xref_rows(
            "sanger_cmp",
            sidm,
            decision,
            [
                ("rrid", row.get("RRID") or ""),
                ("ach", row.get("BROAD_ID") or ""),
                ("cosmic", row.get("COSMIC_ID") or ""),
                ("cell_line_name", row.get("model_name") or ""),
                ("ccle_name", row.get("CCLE_ID") or ""),
            ],
            display_name=clean(row.get("model_name")),
            tissue=clean(row.get("tissue")),
            disease=clean(row.get("cancer_type")),
            species=clean(row.get("species")),
        ))
    for row in gdsc_rows or []:
        sidm = parse_sidm(row.get("SANGER_MODEL_ID"))
        cosmic = parse_cosmic(row.get("COSMIC_ID"))
        name = clean(row.get("CELL_LINE_NAME"))
        source_id = sidm or (f"COSMIC:{cosmic}" if cosmic else name)
        if not source_id:
            continue
        decision = classify_gdsc(sidm, cosmic, name, cmp)
        display = cmp.display.get(decision.cvcl or "", {})
        out.extend(xref_rows(
            "gdsc",
            source_id,
            decision,
            [("sidm", sidm or ""), ("cosmic", cosmic or ""), ("cell_line_name", name or "")],
            display_name=name or display.get("display_name"),
            tissue=display.get("tissue"),
            disease=display.get("disease"),
            species=display.get("species"),
        ))
    seen_prism: set[str] = set()
    for row in prism_rows or []:
        ach = parse_ach(row.get("depmap_id"))
        if not ach or ach in seen_prism:
            continue
        seen_prism.add(ach)
        decision = classify_ach(ach, depmap_decisions, cmp)
        display = cmp.display.get(decision.cvcl or "", {})
        out.extend(xref_rows(
            "prism",
            ach,
            decision,
            [("ach", ach), ("ccle_name", row.get("ccle_name") or "")],
            display_name=display.get("display_name"),
            tissue=clean(row.get("primary_tissue")) or display.get("tissue"),
            disease=display.get("disease"),
            species=display.get("species") or "Homo sapiens",
        ))
    return out


@dataclass
class CompoundRecord:
    source: str
    source_compound_id: str
    name: str | None = None
    synonyms: set[str] = field(default_factory=set)
    brd_id: str | None = None
    gdsc_drug_id: str | None = None
    pubchem_cid: str | None = None
    chembl_id: str | None = None
    inchikey: str | None = None
    smiles: str | None = None
    match_method: str = "unmatched"
    status: str = "unmatched"
    evidence: str = "no PubChem CID in the source file"


def _cid(value: object) -> str | None:
    text = clean(value)
    if text and re.fullmatch(r"\d+", text):
        return text
    return None


def add_compound(records: dict[tuple[str, str], CompoundRecord], record: CompoundRecord) -> None:
    key = (record.source, record.source_compound_id)
    current = records.get(key)
    if current is None:
        records[key] = record
        return
    if record.name and not current.name:
        current.name = record.name
    current.synonyms.update(record.synonyms)
    if record.name and current.name and record.name.lower() != current.name.lower():
        current.synonyms.add(record.name)
    for attr in ("brd_id", "gdsc_drug_id", "smiles", "inchikey", "chembl_id"):
        incoming = getattr(record, attr)
        if incoming and not getattr(current, attr):
            setattr(current, attr, incoming)
    if record.pubchem_cid:
        if current.pubchem_cid and current.pubchem_cid != record.pubchem_cid:
            current.pubchem_cid = None
            current.status = "ambiguous"
            current.match_method = "conflict_pubchem"
            current.evidence = "source files disagree on PubChem CID"
        elif current.status != "ambiguous":
            current.pubchem_cid = record.pubchem_cid
            current.status = "mapped"
            current.match_method = record.match_method
            current.evidence = record.evidence


def compound_rows(records: dict[tuple[str, str], CompoundRecord]) -> list[dict[str, str | None]]:
    rows: list[dict[str, str | None]] = []
    for record in records.values():
        rows.append({
            "source": record.source,
            "source_compound_id": record.source_compound_id,
            "name": record.name,
            "synonyms": ";".join(sorted(record.synonyms)),
            "brd_id": record.brd_id,
            "gdsc_drug_id": record.gdsc_drug_id,
            "pubchem_cid": record.pubchem_cid,
            "chembl_id": record.chembl_id,
            "inchikey": record.inchikey,
            "smiles": record.smiles,
            "match_method": record.match_method,
            "status": record.status,
            "evidence": record.evidence,
        })
    return rows


def compounds_from_gdsc(rows: list[dict[str, str]]) -> list[CompoundRecord]:
    out: list[CompoundRecord] = []
    for row in rows:
        drug_id = clean(row.get("DRUG_ID"))
        if not drug_id:
            continue
        cid = _cid(row.get("PUBCHEM") or row.get("pubchem_cid"))
        chembl = clean(row.get("CHEMBL") or row.get("chembl_id"))
        if chembl and not CHEMBL_RE.match(chembl):
            chembl = None
        record = CompoundRecord(
            source="gdsc",
            source_compound_id=drug_id,
            name=clean(row.get("DRUG_NAME")),
            gdsc_drug_id=drug_id,
            pubchem_cid=cid,
            chembl_id=chembl,
            inchikey=clean(row.get("INCHIKEY") or row.get("inchikey")),
        )
        synonym = clean(row.get("SYNONYMS"))
        if synonym:
            record.synonyms.update(part.strip() for part in synonym.split(",") if clean(part))
        if cid:
            record.status = "mapped"
            record.match_method = "source_pubchem"
            record.evidence = "PubChem CID column"
        out.append(record)
    return out


def compounds_from_prism_names(rows: list[dict[str, str]], source: str, id_field: str) -> list[CompoundRecord]:
    out: list[CompoundRecord] = []
    for row in rows:
        brd = parse_brd(row.get(id_field) or row.get("broad_id") or row.get("IDs"))
        if not brd:
            continue
        name = clean(row.get("name") or row.get("Drug.Name") or row.get("DRUG_NAME"))
        record = CompoundRecord(source=source, source_compound_id=brd, name=name, brd_id=brd, smiles=clean(row.get("smiles")))
        cid = _cid(row.get("pubchem_cid") or row.get("PUBCHEM"))
        if cid:
            record.pubchem_cid = cid
            record.status = "mapped"
            record.match_method = "source_pubchem"
            record.evidence = "PubChem CID column"
        out.append(record)
    return out


def apply_pubchem(record: CompoundRecord, cid: str | None, inchikey: str | None, chembl_id: str | None, method: str, evidence: str) -> None:
    if record.status == "ambiguous":
        return
    if record.pubchem_cid and cid and record.pubchem_cid != cid:
        record.status = "ambiguous"
        record.match_method = "conflict_pubchem"
        record.evidence = f"existing CID {record.pubchem_cid} disagrees with {cid}"
        record.pubchem_cid = None
        return
    if not cid:
        record.status = "unmatched"
        record.match_method = method
        record.evidence = evidence
        return
    record.pubchem_cid = cid
    record.inchikey = inchikey or record.inchikey
    record.chembl_id = chembl_id or record.chembl_id
    record.status = "mapped"
    record.match_method = method
    record.evidence = evidence
