from pgx.ids import CompoundRecord, compounds_from_gdsc
from pgx.normalize.combinations import (
    combination_components,
    is_annotation_combination,
    is_internal_code,
    parse_annotation_combination,
    parse_screen_combination,
)
from pgx.normalize.compounds import iter_queries, resolve_compounds
from pgx.process import _write_rows, connect, load_combinations


def test_annotation_combo_is_anchor_then_library_not_one_compound():
    assert is_annotation_combination("BAY-MPS-combo-1 (paclitaxel 5 uM)")
    parsed = parse_annotation_combination("1369", "BAY-MPS-combo-1 (paclitaxel 5 uM)")
    assert parsed is not None
    assert parsed["combination_id"] == "gdsc.combination:annotation:1369"
    components = combination_components(parsed)
    assert [item["role"] for item in components] == ["anchor", "library"]
    assert components[0]["name"] == "BAY-MPS"
    assert components[1]["name"] == "paclitaxel"
    assert components[1]["concentrations"] == "5"
    assert components[1]["unit"] == "uM"
    spaced = parse_annotation_combination("1370", "BAY-MPS-combo 2 (paclitaxel 1 uM)")
    assert spaced is not None
    assert combination_components(spaced)[0]["name"] == "BAY-MPS"
    assert combination_components(spaced)[1]["concentrations"] == "1"
    rows = compounds_from_gdsc([{"DRUG_ID": "1369", "DRUG_NAME": "BAY-MPS-combo-1 (paclitaxel 5 uM)", "SYNONYMS": ""}])
    assert rows[0].status == "combination"
    assert rows[0].pubchem_cid is None
    assert rows[0].review_reason == "combination"


def test_screen_combi_id_keeps_anchor_before_library():
    parsed = parse_screen_combination(
        "1012", "Vorinostat", "1003", "Camptothecin",
        anchor_concentrations=[2.0, 0.5],
        library_concentrations=[0.1],
        screen="gdsc_combination_validation",
    )
    assert parsed is not None
    assert parsed["combination_id"] == "gdsc.combination:1012:1003"
    assert parsed["name"] == "Vorinostat + Camptothecin"
    components = combination_components(parsed)
    assert [item["source_compound_id"] for item in components] == ["1012", "1003"]
    assert components[0]["concentrations"] == "0.5;2"
    assert components[1]["concentrations"] == "0.1"
    assert parse_screen_combination("Vorinostat", "Vorinostat", "1003", "Camptothecin") is None
    triple = parse_screen_combination("1032|1372", "Afatinib|Trametinib", "1003", "Camptothecin", [1.0], [0.1])
    assert triple is not None
    triple_parts = combination_components(triple)
    assert [item["role"] for item in triple_parts] == ["anchor", "anchor", "library"]
    assert [item["name"] for item in triple_parts] == ["Afatinib", "Trametinib", "Camptothecin"]
    assert [item["source_compound_id"] for item in triple_parts] == ["1032", "1372", "1003"]
    assert triple["combination_id"] == "gdsc.combination:1032|1372:1003"


def test_counterion_parentheses_and_locant_commas_are_not_extra_drugs():
    from pgx.normalize.compounds import extra_lookup_names

    assert extra_lookup_names("POLYMYXIN B NONAPEPTIDE (TFA)") == ["POLYMYXIN B NONAPEPTIDE"]
    assert "vitamin C" in extra_lookup_names("ascorbate (vitamin C)")
    assert extra_lookup_names("MEBHYDROLINE-1,5-NAPHTALENEDISULFONATE") == []
    assert "Brivanib" in extra_lookup_names("Brivanib, BMS-540215")


def test_single_agents_are_not_combinations():
    for name in (
        "Nutlin-3a (-)",
        "Bleomycin (50 uM)",
        "VNLG/124",
        "Akti-1/2",
        "Brivanib, BMS-540215",
        "ascorbate (vitamin C)",
    ):
        assert is_annotation_combination(name) is False
        assert parse_annotation_combination("1", name) is None


def test_internal_codes_skip_lookup_and_public_names_do_not():
    assert is_internal_code("CDK9_5576", [])
    assert is_internal_code("KIN001-236", [])
    assert is_internal_code("123829", [])
    assert is_internal_code("THR-101", ["WIMM synthesis"])
    assert is_internal_code("Genentech Cpd 10", [])
    assert is_internal_code("BAY-HDAC11_1", [])
    assert is_internal_code("GSK2256098C", []) is False
    assert is_internal_code("CI-1033", []) is False
    assert is_internal_code("FTY-720", []) is False
    assert is_internal_code("BAY-MPS1", []) is False
    assert is_internal_code("AZD7969", []) is False
    client = _Client({})
    internal = CompoundRecord(source="gdsc", source_compound_id="286", name="KIN001-236", gdsc_drug_id="286")
    combo = CompoundRecord(source="gdsc", source_compound_id="1369", name="BAY-MPS-combo-1 (paclitaxel 5 uM)")
    resolve_compounds([internal, combo], client)
    assert client.calls == []
    assert internal.review_reason == "internal_code"
    assert internal.status == "unmatched"
    assert combo.status == "combination"
    assert combo.pubchem_cid is None
    assert iter_queries(combo) == []


class _Client:
    def __init__(self, table, chembl=None):
        self.table = table
        self.chembl = chembl or {}
        self.calls = []

    def cids(self, kind, query):
        self.calls.append((kind, query))
        return list(self.table.get((kind, query), []))

    def properties_for(self, cids):
        return {cid: {"inchikey": "RFTSSZJZXOSICM-GRSHGNNSSA-N", "title": "Obatoclax"} for cid in cids}

    def chembl_for(self, inchikeys):
        return {key: ["CHEMBL408194"] for key in inchikeys}

    def chembl_by_name(self, query):
        self.calls.append(("chembl_name", query))
        return list(self.chembl.get(query.casefold(), []))


def test_chembl_synonym_disambiguates_and_records_missing_structure():
    ambiguous = CompoundRecord(
        source="gdsc", source_compound_id="182", name="Obatoclax Mesylate",
        synonyms={"Obatoclax", "GX15-070"}, status="ambiguous", match_method="pubchem_name",
    )
    bleomycin = CompoundRecord(source="gdsc", source_compound_id="1378", name="Bleomycin (50 uM)", gdsc_drug_id="1378")
    client = _Client(
        {("inchikey", "RFTSSZJZXOSICM-GRSHGNNSSA-N"): ["49848060"]},
        {
            "obatoclax": [{"chembl_id": "CHEMBL408194", "inchikey": "RFTSSZJZXOSICM-GRSHGNNSSA-N", "pref_name": "OBATOCLAX"}],
            "bleomycin": [{"chembl_id": "CHEMBL6067492", "inchikey": None, "pref_name": "BLEOMYCIN"}],
        },
    )
    resolve_compounds([ambiguous, bleomycin], client)
    assert ambiguous.status == "mapped"
    assert ambiguous.pubchem_cid == "49848060"
    assert ambiguous.chembl_id == "CHEMBL408194"
    assert ambiguous.match_method == "chembl_synonym"
    assert ambiguous.confidence == "medium"
    assert bleomycin.status == "unmatched"
    assert bleomycin.review_reason == "no_public_structure"
    assert bleomycin.chembl_id == "CHEMBL6067492"
    assert bleomycin.pubchem_cid is None
    assert ("name", "Bleomycin") in client.calls
    assert ("name", "Bleomycin (50 uM)") in client.calls


def test_chembl_salt_without_a_pubchem_cid_tries_the_parent_name():
    record = CompoundRecord(
        source="gdsc", source_compound_id="182", name="Obatoclax Mesylate",
        synonyms={"Obatoclax"}, status="ambiguous",
    )
    client = _Client(
        {("inchikey", "RFTSSZJZXOSICM-GRSHGNNSSA-N"): ["49848060"]},
        {
            "obatoclax mesylate": [{"chembl_id": "CHEMBL2107358", "inchikey": "ZFKXDVMHNXPEIY-PEZBNFGJSA-N", "pref_name": "salt"}],
            "obatoclax": [{"chembl_id": "CHEMBL408194", "inchikey": "RFTSSZJZXOSICM-GRSHGNNSSA-N", "pref_name": "OBATOCLAX"}],
        },
    )
    resolve_compounds([record], client)
    assert record.status == "mapped"
    assert record.pubchem_cid == "49848060"
    assert record.chembl_id == "CHEMBL408194"
    assert record.review_reason is None


def test_ambiguous_pubchem_name_stays_ambiguous_when_chembl_has_no_cid():
    record = CompoundRecord(
        source="gdsc", source_compound_id="182", name="Obatoclax Mesylate",
        synonyms={"Obatoclax"}, status="ambiguous", evidence="Obatoclax matched 7 CIDs",
    )
    client = _Client(
        {},
        {"obatoclax": [{"chembl_id": "CHEMBL408194", "inchikey": "RFTSSZJZXOSICM-GRSHGNNSSA-N", "pref_name": "OBATOCLAX"}]},
    )
    resolve_compounds([record], client)
    assert record.status == "ambiguous"
    assert record.pubchem_cid is None
    assert record.review_reason is None
    assert "7 CIDs" in record.evidence


def test_combination_screen_rows_stay_out_of_single_agent_responses(tmp_path):
    raw = tmp_path / "raw"
    compounds = raw / "gdsc/8.4/screened_compounds_rel_8.4.csv"
    compounds.parent.mkdir(parents=True)
    compounds.write_text(
        "DRUG_ID,DRUG_NAME,SYNONYMS\n"
        "1369,BAY-MPS-combo-1 (paclitaxel 5 uM),\n"
        "11,Paclitaxel,\n"
        "1012,Vorinostat,\n"
        "1003,Camptothecin,\n"
    )
    screen = raw / "gdsc_combinations/2022/Validation_screen_All_tissues_fitted_data.csv"
    screen.parent.mkdir(parents=True)
    screen.write_text(
        "SIDM,COSMIC_ID,CELL_LINE_NAME,ANCHOR_ID,ANCHOR_NAME,ANCHOR_CONC,LIBRARY_ID,LIBRARY_NAME,LIBRARY_CONC,"
        "SYNERGY_DELTA_AUC,SYNERGY_DELTA_EMAX,SYNERGY_DELTA_XMID,SYNERGY_DELTA_XMID_uM,SYNERGY_DELTA_fAUC\n"
        "SIDM01042,100,EVSA-T,1012,Vorinostat,2,1003,Camptothecin,0.1,0.005,0.1,0.2,0.3,0.4\n"
        "SIDM01042,100,EVSA-T,1012,Vorinostat,0.5,1003,Camptothecin,0.1,0.008,0.2,0.3,0.4,0.5\n"
    )
    out = tmp_path / "processed"
    con = connect()
    _write_rows(con, [{
        "source": "sanger_cmp", "source_model_id": "SIDM01042", "alias_type": "primary",
        "status": "mapped", "cvcl": "CVCL_1202",
    }], out / "model_xref", ["source", "source_model_id", "alias_type", "status", "cvcl"])
    _write_rows(con, [
        {"source": "gdsc", "source_compound_id": "11", "name": "Paclitaxel", "pubchem_cid": "36314", "chembl_id": "CHEMBL428647", "status": "mapped"},
        {"source": "gdsc", "source_compound_id": "1012", "name": "Vorinostat", "pubchem_cid": "5311", "chembl_id": "CHEMBL98", "status": "mapped"},
        {"source": "gdsc", "source_compound_id": "1003", "name": "Camptothecin", "pubchem_cid": "24360", "chembl_id": "CHEMBL65", "status": "mapped"},
    ], out / "compound_xref", ["source", "source_compound_id", "name", "pubchem_cid", "chembl_id", "status"])
    document = {"sources": [
        {"id": "gdsc", "version": "8.4", "enabled": True, "files": [
            {"name": "screened_compounds_rel_8.4.csv", "role": "gdsc_compounds"},
        ]},
        {"id": "gdsc_combinations", "version": "2022", "enabled": True, "files": [
            {"name": "Validation_screen_All_tissues_fitted_data.csv", "role": "gdsc_combination_screen", "dataset": "gdsc_combination_validation"},
        ]},
    ]}
    counts = load_combinations(document, raw, out)
    assert counts["combinations"] == 2
    assert counts["combination_responses"] == 10
    combo = {
        row[0]: row
        for row in con.execute(
            "SELECT combination_id, anchor_name, library_name, anchor_pubchem_cid, library_pubchem_cid, anchor_concentrations, library_concentrations, status FROM read_parquet(?)",
            [str(out / "combination_xref" / "part.parquet")],
        ).fetchall()
    }
    screen_row = combo["gdsc.combination:1012:1003"]
    assert screen_row[1:5] == ("Vorinostat", "Camptothecin", "5311", "24360")
    assert screen_row[5] == "0.5;2"
    assert screen_row[7] == "mapped"
    annotation = combo["gdsc.combination:annotation:1369"]
    assert annotation[1:3] == ("BAY-MPS", "paclitaxel")
    assert annotation[4] == "36314"
    assert annotation[7] == "partial"
    response = con.execute(
        """
        SELECT response_class, combination_id, pubchem_cid, compound_name, metric, anchor_conc, library_conc, model_cvcl
        FROM read_parquet(?, hive_partitioning=true)
        WHERE metric = 'SYNERGY_DELTA_AUC'
        ORDER BY anchor_conc
        """,
        [str(out / "drug_response" / "**" / "*.parquet")],
    ).fetchall()
    assert response[0][0] == "combination"
    assert response[0][1] == "gdsc.combination:1012:1003"
    assert response[0][2] is None
    assert response[0][3] == "Vorinostat + Camptothecin"
    assert response[0][7] == "CVCL_1202"
    assert [row[5] for row in response] == [0.5, 2.0]
