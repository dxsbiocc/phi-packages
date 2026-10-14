from pathlib import Path

from pgx.ids import CompoundRecord, build_model_xref
from pgx.normalize.cell_lines import CellosaurusIndex, canonical_name
from pgx.normalize.compounds import iter_queries, name_variants, resolve_compounds
from pgx.normalize.genes import HgncIndex, parse_gene_header, resolve_gene, resolve_labels


HGNC = """hgnc_id\tsymbol\tstatus\talias_symbol\tprev_symbol\tentrez_id\tensembl_gene_id
HGNC:11998\tTP53\tApproved\tp53|LFS1\tTP53L\t7157\tENSG00000141510
HGNC:1\tA1BG\tApproved\t\tA1B\t1\tENSG00000121410
HGNC:2\tONE\tApproved\tSHARED\t\t10\tENSG00000000010
HGNC:3\tTWO\tApproved\tSHARED\t\t11\tENSG00000000011
"""

CELLOSAURUS = """ID   NCI-H460
AC   CVCL_0459
SY   NCIH460; NCI-H460
OX   NCBI_TaxID=9606; ! Homo sapiens
//
ID   OtherLine
AC   CVCL_0001
SY   OnlyOther
//
ID   Alpha
AC   CVCL_0002
SY   Shared
//
ID   Beta
AC   CVCL_0003
SY   Shared
//
ID   BadLine
AC   CVCL_0008
SY   BadLine
CC   Problematic cell line: Misidentified. The patient origin is wrong.
//
ID   OldName
AC   CVCL_0100
AS   CVCL_9990
SY   OldName
//
"""


def cmp_row(**kwargs):
    base = {
        "model_id": "", "RRID": "", "BROAD_ID": "", "COSMIC_ID": "", "model_name": "",
        "synonyms": "", "CCLE_ID": "", "species": "Homo sapiens", "tissue": "", "cancer_type": "",
    }
    base.update(kwargs)
    return base


def depmap_row(**kwargs):
    base = {
        "ModelID": "", "RRID": "", "SangerModelID": "", "COSMICID": "", "CellLineName": "",
        "StrippedCellLineName": "", "CCLEName": "", "OncotreeLineage": "", "OncotreePrimaryDisease": "",
    }
    base.update(kwargs)
    return base


def write_hgnc(tmp_path: Path) -> HgncIndex:
    path = tmp_path / "hgnc.tsv"
    path.write_text(HGNC)
    return HgncIndex.load(path)


def write_cellosaurus(tmp_path: Path) -> CellosaurusIndex:
    path = tmp_path / "cellosaurus.txt"
    path.write_text(CELLOSAURUS)
    return CellosaurusIndex.load(path)


def test_tp53_header_splits_symbol_and_entrez(tmp_path):
    parsed = parse_gene_header("TP53 (7157)")
    assert parsed.source_symbol == "TP53"
    assert parsed.entrez_ids == ["7157"]
    assert parsed.parse_status == "ok"
    decision = resolve_gene(parsed, write_hgnc(tmp_path))
    assert decision.status == "mapped"
    assert decision.match_method == "hgnc_entrez"
    assert decision.symbol == "TP53"
    assert decision.ncbigene == "7157"
    assert decision.ensembl_gene_id == "ENSG00000141510"
    assert decision.source_symbol == "TP53"


def test_previous_symbol_maps_to_approved_symbol(tmp_path):
    index = write_hgnc(tmp_path)
    decision = resolve_gene(parse_gene_header("TP53L"), index)
    assert decision.status == "mapped"
    assert decision.match_method == "hgnc_previous"
    assert decision.symbol == "TP53"
    assert decision.source_symbol == "TP53L"
    assert decision.ncbigene == "7157"
    alias = resolve_gene(parse_gene_header("p53"), index)
    assert alias.match_method == "hgnc_alias"
    assert alias.symbol == "TP53"
    assert alias.source_symbol == "p53"


def test_gene_headers_without_one_entrez_id_go_to_review(tmp_path):
    index = write_hgnc(tmp_path)
    multi = resolve_gene(parse_gene_header("TP53 (7157;1)"), index)
    assert multi.status == "ambiguous"
    assert multi.match_method == "multiple_entrez"
    odd = resolve_gene(parse_gene_header("not a gene"), index)
    assert odd.status == "unmatched"
    assert odd.match_method == "odd_symbol"
    shared = resolve_gene(parse_gene_header("SHARED"), index)
    assert shared.status == "ambiguous"
    bare = resolve_gene(parse_gene_header("TP53"), index)
    assert bare.status == "mapped"
    assert bare.match_method == "hgnc_symbol"
    conflict = resolve_gene(parse_gene_header("TP53 (ENSG00000000000)"), index)
    assert conflict.status == "ambiguous"
    assert conflict.match_method == "conflict_ensembl"
    labels = resolve_labels(["TP53 (7157)", "TP53 (7157)", "A1B"], index)
    assert [item.source_label for item in labels] == ["TP53 (7157)", "A1B"]
    assert labels[1].symbol == "A1BG"
    assert labels[1].match_method == "hgnc_previous"


def test_nci_h460_canonical_form_does_not_drop_the_prefix(tmp_path):
    assert canonical_name("NCI-H460") == canonical_name("NCIH460") == "ncih460"
    assert canonical_name("H460") == "h460"
    assert canonical_name("NCI-H460") != canonical_name("H460")
    index = write_cellosaurus(tmp_path)
    by_name = {row["source_model_id"]: row for row in build_model_xref(
        [cmp_row(model_id="SIDM00001", RRID="CVCL_9999", model_name="Unrelated")],
        [
            depmap_row(ModelID="ACH-000460", CellLineName="NCI-H460", CCLEName="NCIH460_LUNG"),
            depmap_row(ModelID="ACH-000459", CCLEName="NCIH460_LUNG"),
            depmap_row(ModelID="ACH-000461", CellLineName="H460"),
            depmap_row(ModelID="ACH-000462", CellLineName="Shared"),
            depmap_row(ModelID="ACH-000008", RRID="CVCL_0008", CellLineName="BadLine"),
            depmap_row(ModelID="ACH-000100", RRID="CVCL_9990", CellLineName="OldName"),
        ],
        cellosaurus=index,
    ) if row["source"] == "depmap" and row["alias_type"] == "primary"}
    assert by_name["ACH-000460"]["status"] == "mapped"
    assert by_name["ACH-000460"]["cvcl"] == "CVCL_0459"
    assert by_name["ACH-000460"]["match_method"] == "cellosaurus_name"
    assert by_name["ACH-000459"]["status"] == "mapped"
    assert by_name["ACH-000459"]["cvcl"] == "CVCL_0459"
    assert by_name["ACH-000459"]["match_method"] == "cellosaurus_ccle_stem"
    assert by_name["ACH-000461"]["status"] == "unmatched"
    assert by_name["ACH-000461"]["cvcl"] is None
    assert by_name["ACH-000462"]["status"] == "ambiguous"
    assert by_name["ACH-000462"]["cvcl"] is None
    assert by_name["ACH-000008"]["status"] == "mapped"
    assert by_name["ACH-000008"]["cellosaurus_flag"] == "misidentified"
    assert by_name["ACH-000100"]["cvcl"] == "CVCL_0100"
    assert "secondary" in by_name["ACH-000100"]["evidence"]


def test_h460_matches_only_when_cellosaurus_lists_the_synonym(tmp_path):
    text = CELLOSAURUS.replace("SY   NCIH460; NCI-H460", "SY   NCIH460; NCI-H460; H460")
    path = tmp_path / "cellosaurus.txt"
    path.write_text(text)
    index = CellosaurusIndex.load(path)
    decision = index.decide_names(["H460"])
    assert decision.status == "mapped"
    assert decision.cvcl == "CVCL_0459"
    assert decision.match_method == "cellosaurus_synonym"
    assert canonical_name("H460") != canonical_name("NCI-H460")


class FakeClient:
    def __init__(self, table):
        self.table = table
        self.calls: list[tuple[str, str]] = []

    def cids(self, kind, query):
        self.calls.append((kind, query))
        return list(self.table.get((kind, query), []))

    def properties_for(self, cids):
        return {cid: {"inchikey": "BDUHCSBCVGXTJM-WUFINQPMSA-N", "title": "Nutlin 3"} for cid in cids}

    def chembl_for(self, inchikeys):
        return {key: ["CHEMBL191334"] for key in inchikeys}


def test_nutlin_parenthetical_is_lookup_only():
    variants = name_variants("Nutlin-3a (-)")
    assert variants[0] == ("exact", "Nutlin-3a (-)", "medium")
    assert ("stereo", "Nutlin-3a", "low") in variants
    record = CompoundRecord(source="gdsc", source_compound_id="1047", name="Nutlin-3a (-)", gdsc_drug_id="1047")
    assert record.name == "Nutlin-3a (-)"
    assert [item[1] for item in iter_queries(record)][:2] == ["Nutlin-3a (-)", "Nutlin-3a"]
    client = FakeClient({
        ("name", "Nutlin-3a (-)"): [],
        ("name", "Nutlin-3a"): ["11433190"],
    })
    mapped = resolve_compounds([record], client)
    assert mapped == 1
    assert record.name == "Nutlin-3a (-)"
    assert record.pubchem_cid == "11433190"
    assert record.match_method == "pubchem_name_stereo"
    assert record.confidence == "low"
    assert record.status == "mapped"
    assert record.chembl_id == "CHEMBL191334"
    assert record.inchikey == "BDUHCSBCVGXTJM-WUFINQPMSA-N"
    salted = name_variants("Erlotinib hydrochloride")
    assert ("salt", "Erlotinib", "low") in salted
    assert ("salt", "GSK-LSD1", "low") in name_variants("GSK-LSD1-2HCl")


def test_multiple_cids_are_ambiguous_and_stop():
    record = CompoundRecord(source="gdsc", source_compound_id="9", name="Ambiguous (-)")
    client = FakeClient({
        ("name", "Ambiguous (-)"): ["1", "2"],
        ("name", "Ambiguous"): ["9"],
    })
    resolve_compounds([record], client)
    assert record.status == "ambiguous"
    assert record.pubchem_cid is None
    assert client.calls == [("name", "Ambiguous (-)")]
