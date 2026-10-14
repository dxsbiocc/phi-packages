"""Parse DepMap gene headers and map them onto the HGNC approved set."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9@._-]*$")
ENSG_RE = re.compile(r"^ENSG\d{11}$", re.IGNORECASE)
HEADER_RE = re.compile(r"^(?P<symbol>.+?)\s*\((?P<ids>[^()]*)\)\s*$")
LEGACY_LABEL = re.compile(r"^(.*) \(([0-9]+)\)$")


@dataclass
class ParsedGene:
    source_label: str
    source_symbol: str | None
    entrez_ids: list[str]
    ensembl_ids: list[str]
    odd_tokens: list[str]
    parse_status: str


@dataclass(frozen=True)
class HgncGene:
    hgnc_id: str
    symbol: str
    entrez_id: str | None
    ensembl_gene_id: str | None


@dataclass(frozen=True)
class SymbolHit:
    gene: HgncGene
    method: str


@dataclass
class SymbolMatch:
    gene: HgncGene | None
    method: str
    status: str
    evidence: str


@dataclass
class GeneDecision:
    source_label: str
    source_symbol: str | None
    symbol: str | None
    ncbigene: str | None
    ensembl_gene_id: str | None
    match_method: str
    status: str
    evidence: str

    def row(self) -> dict[str, str | None]:
        return {
            "source_label": self.source_label,
            "source_symbol": self.source_symbol,
            "symbol": self.symbol,
            "ncbigene": self.ncbigene,
            "ensembl_gene_id": self.ensembl_gene_id,
            "match_method": self.match_method,
            "status": self.status,
            "evidence": self.evidence,
        }


def parse_gene_header(label: str) -> ParsedGene:
    """Split a DepMap-style header into a symbol and any parenthetical ids.

    `TP53 (7157)` yields one Entrez id. A header with no parentheses is a
    symbol-only token when it looks like a gene symbol. Several Entrez ids,
    Ensembl ids (`ENSG...`), or other parenthetical tokens are kept separate.
    """
    text = (label or "").strip()
    if not text:
        return ParsedGene(text, None, [], [], [], "empty")
    match = HEADER_RE.match(text)
    if not match:
        if SYMBOL_RE.match(text):
            return ParsedGene(text, text, [], [], [], "no_id")
        return ParsedGene(text, None, [], [], [text], "odd_symbol")
    symbol = match.group("symbol").strip()
    entrez: list[str] = []
    ensembl: list[str] = []
    odd: list[str] = []
    for token in (part.strip() for part in re.split(r"[;,|/]", match.group("ids"))):
        if not token:
            continue
        if re.fullmatch(r"\d+", token):
            entrez.append(str(int(token)))
        elif ENSG_RE.match(token):
            ensembl.append(token.upper())
        else:
            odd.append(token)
    if len(set(entrez)) > 1:
        status = "multiple_entrez"
    elif not SYMBOL_RE.match(symbol):
        status = "odd_symbol"
    elif not entrez and not ensembl and not odd:
        status = "no_id"
    else:
        status = "ok"
    source_symbol = symbol or None
    return ParsedGene(text, source_symbol, entrez, ensembl, odd, status)


def _pipe(value: str | None) -> list[str]:
    if not value:
        return []
    return [part.strip() for part in value.split("|") if part.strip()]


class HgncIndex:
    def __init__(self) -> None:
        self.by_entrez: dict[str, HgncGene] = {}
        self.ambiguous_entrez: set[str] = set()
        self.exact: dict[str, list[SymbolHit]] = {}
        self.folded: dict[str, list[SymbolHit]] = {}
        self.genes = 0

    def _add_symbol(self, symbol: str, gene: HgncGene, method: str) -> None:
        for bucket, key in ((self.exact, symbol), (self.folded, symbol.casefold())):
            hits = bucket.setdefault(key, [])
            if any(hit.gene.hgnc_id == gene.hgnc_id and hit.method == method for hit in hits):
                continue
            hits.append(SymbolHit(gene, method))

    def _add_entrez(self, entrez: str, gene: HgncGene) -> None:
        if entrez in self.ambiguous_entrez:
            return
        current = self.by_entrez.get(entrez)
        if current and current.hgnc_id != gene.hgnc_id:
            del self.by_entrez[entrez]
            self.ambiguous_entrez.add(entrez)
            return
        self.by_entrez[entrez] = gene

    def add_gene(self, gene: HgncGene, previous: list[str], aliases: list[str]) -> None:
        self.genes += 1
        if gene.entrez_id:
            self._add_entrez(gene.entrez_id, gene)
        self._add_symbol(gene.symbol, gene, "hgnc_symbol")
        for symbol in previous:
            self._add_symbol(symbol, gene, "hgnc_previous")
        for symbol in aliases:
            self._add_symbol(symbol, gene, "hgnc_alias")

    @classmethod
    def load(cls, path: Path) -> HgncIndex:
        index = cls()
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if (row.get("status") or "").strip() != "Approved":
                    continue
                symbol = (row.get("symbol") or "").strip()
                hgnc_id = (row.get("hgnc_id") or "").strip()
                if not symbol or not hgnc_id:
                    continue
                entrez = (row.get("entrez_id") or "").strip() or None
                ensembl = (row.get("ensembl_gene_id") or "").strip() or None
                gene = HgncGene(hgnc_id, symbol, entrez, ensembl.upper() if ensembl else None)
                index.add_gene(gene, _pipe(row.get("prev_symbol")), _pipe(row.get("alias_symbol")))
        return index

    def lookup_symbol(self, symbol: str) -> SymbolMatch:
        exact = self.exact.get(symbol, [])
        if exact:
            return _collapse(exact, f"exact symbol {symbol}")
        folded = self.folded.get(symbol.casefold(), [])
        if folded:
            return _collapse(folded, f"case-folded symbol {symbol}")
        return SymbolMatch(None, "unmatched", "unmatched", f"{symbol} is not an HGNC symbol")


def _collapse(hits: list[SymbolHit], how: str) -> SymbolMatch:
    approved: dict[str, HgncGene] = {}
    for hit in hits:
        if hit.method == "hgnc_symbol":
            approved[hit.gene.hgnc_id] = hit.gene
    if len(approved) == 1:
        return SymbolMatch(next(iter(approved.values())), "hgnc_symbol", "mapped", how)
    if len(approved) > 1:
        return SymbolMatch(None, "hgnc_symbol", "ambiguous", f"{how} matches {len(approved)} approved symbols")
    by_gene: dict[str, HgncGene] = {}
    methods: dict[str, set[str]] = {}
    for hit in hits:
        by_gene[hit.gene.hgnc_id] = hit.gene
        methods.setdefault(hit.gene.hgnc_id, set()).add(hit.method)
    if len(by_gene) != 1:
        return SymbolMatch(None, "hgnc_alias", "ambiguous", f"{how} matches {len(by_gene)} HGNC genes")
    gene = next(iter(by_gene.values()))
    method = "hgnc_previous" if "hgnc_previous" in methods[gene.hgnc_id] else "hgnc_alias"
    return SymbolMatch(gene, method, "mapped", how)


def _ensembl(parsed: ParsedGene, gene: HgncGene) -> tuple[str | None, str | None]:
    """Return (ensembl, error). A disagreeing header Ensembl id is an error."""
    approved = gene.ensembl_gene_id
    if len(parsed.ensembl_ids) > 1:
        return None, "header lists multiple Ensembl ids"
    if not parsed.ensembl_ids:
        return approved, None
    header = parsed.ensembl_ids[0]
    if approved and header != approved:
        return None, f"header {header} disagrees with HGNC {approved}"
    return approved or header, None


def _from_symbol(parsed: ParsedGene, match: SymbolMatch, extra: str = "") -> GeneDecision:
    gene = match.gene
    assert gene is not None
    ensembl, error = _ensembl(parsed, gene)
    if error:
        return GeneDecision(
            parsed.source_label, parsed.source_symbol, gene.symbol, gene.entrez_id, None,
            "conflict_ensembl", "ambiguous", error,
        )
    evidence = f"HGNC {gene.hgnc_id} via {match.method}"
    if parsed.source_symbol and parsed.source_symbol != gene.symbol:
        evidence += f"; source symbol {parsed.source_symbol}"
    if extra:
        evidence += f"; {extra}"
    if parsed.odd_tokens:
        evidence += "; odd tokens " + ",".join(parsed.odd_tokens)
    return GeneDecision(
        parsed.source_label,
        parsed.source_symbol,
        gene.symbol,
        gene.entrez_id,
        ensembl,
        match.method,
        "mapped",
        evidence,
    )


def resolve_gene(parsed: ParsedGene, index: HgncIndex) -> GeneDecision:
    """Map one parsed header onto a single approved HGNC gene, or send it to review."""
    label = parsed.source_label
    source = parsed.source_symbol

    def finish(status: str, method: str, evidence: str, symbol: str | None = None, ncbigene: str | None = None, ensembl: str | None = None) -> GeneDecision:
        return GeneDecision(label, source, symbol, ncbigene, ensembl, method, status, evidence)

    if parsed.parse_status == "empty":
        return finish("unmatched", "empty", "empty label")
    if parsed.parse_status == "multiple_entrez":
        return finish("ambiguous", "multiple_entrez", "header lists multiple Entrez ids: " + ",".join(parsed.entrez_ids))
    if len(parsed.ensembl_ids) > 1:
        return finish("ambiguous", "multiple_ensembl", "header lists multiple Ensembl ids")

    if len(parsed.entrez_ids) == 1:
        entrez = parsed.entrez_ids[0]
        if entrez in index.ambiguous_entrez:
            return finish("ambiguous", "hgnc_entrez", f"Entrez {entrez} matches more than one HGNC gene", ncbigene=entrez)
        gene = index.by_entrez.get(entrez)
        if gene is not None:
            if source and SYMBOL_RE.match(source):
                symbol_match = index.lookup_symbol(source)
                if symbol_match.status == "mapped" and symbol_match.gene and symbol_match.gene.hgnc_id != gene.hgnc_id:
                    return finish(
                        "ambiguous",
                        "conflict_symbol",
                        f"Entrez {entrez} is {gene.symbol} but {source} is {symbol_match.gene.symbol}",
                        symbol=gene.symbol,
                        ncbigene=entrez,
                    )
            ensembl, error = _ensembl(parsed, gene)
            if error:
                return finish("ambiguous", "conflict_ensembl", error, symbol=gene.symbol, ncbigene=entrez)
            evidence = f"HGNC {gene.hgnc_id} entrez {entrez}"
            if source and source != gene.symbol:
                evidence += f"; source symbol {source}"
            if parsed.odd_tokens:
                evidence += "; odd tokens " + ",".join(parsed.odd_tokens)
            return finish("mapped", "hgnc_entrez", evidence, gene.symbol, entrez, ensembl)
        if source and SYMBOL_RE.match(source):
            symbol_match = index.lookup_symbol(source)
            if symbol_match.status == "mapped" and symbol_match.gene and symbol_match.gene.entrez_id and symbol_match.gene.entrez_id != entrez:
                return finish(
                    "ambiguous",
                    "conflict_symbol",
                    f"header Entrez {entrez} disagrees with {symbol_match.gene.symbol} entrez {symbol_match.gene.entrez_id}",
                )
            if symbol_match.status == "ambiguous":
                return finish("ambiguous", symbol_match.method, symbol_match.evidence)
        return finish("unmatched", "hgnc_entrez_miss", f"Entrez {entrez} is not in the HGNC approved set", symbol=source, ncbigene=entrez)

    if parsed.parse_status == "odd_symbol" or not source:
        return finish("unmatched", "odd_symbol", f"symbol is not a gene-symbol token: {label}")
    symbol_match = index.lookup_symbol(source)
    if symbol_match.status == "ambiguous":
        return finish("ambiguous", symbol_match.method, symbol_match.evidence)
    if symbol_match.status == "mapped":
        return _from_symbol(parsed, symbol_match)
    return finish("unmatched", "unmatched", f"{source} is not an approved, previous, or alias symbol")


def resolve_labels(labels: list[str], index: HgncIndex) -> list[GeneDecision]:
    seen: set[str] = set()
    decisions: list[GeneDecision] = []
    for label in labels:
        if label in seen:
            continue
        seen.add(label)
        decisions.append(resolve_gene(parse_gene_header(label), index))
    return decisions


def legacy_parse_counts(labels: list[str]) -> dict[str, int]:
    """Coverage of the pre-normalization regex, which keeps only `SYMBOL (digits)`."""
    kept = 0
    for label in labels:
        if LEGACY_LABEL.match(label.strip()):
            kept += 1
    return {"labels": len(labels), "regex_parsed": kept, "regex_dropped": len(labels) - kept}


def decision_counts(decisions: list[GeneDecision]) -> dict[str, int]:
    counts: dict[str, int] = {"labels": len(decisions), "mapped": 0, "unmatched": 0, "ambiguous": 0, "symbol_updated": 0}
    for decision in decisions:
        counts[decision.status] = counts.get(decision.status, 0) + 1
        if decision.status == "mapped" and decision.source_symbol and decision.symbol and decision.source_symbol != decision.symbol:
            counts["symbol_updated"] += 1
        method_key = f"method:{decision.match_method}"
        counts[method_key] = counts.get(method_key, 0) + 1
    return counts


def header_labels(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        header = handle.readline()
    if not header:
        return []
    return [column.strip() for column in next(csv.reader([header])) if column.strip()]
