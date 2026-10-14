from pgx.ids import build_model_xref, classify_depmap, index_cmp, parse_gene_label


def cmp_row(**kwargs):
    base = {"model_id": "", "RRID": "", "BROAD_ID": "", "COSMIC_ID": "", "model_name": "", "synonyms": "", "CCLE_ID": "", "species": "Homo sapiens", "tissue": "", "cancer_type": ""}
    base.update(kwargs)
    return base


def depmap_row(**kwargs):
    base = {"ModelID": "", "RRID": "", "SangerModelID": "", "COSMICID": "", "CellLineName": "", "StrippedCellLineName": "", "CCLEName": "", "OncotreeLineage": "", "OncotreePrimaryDisease": ""}
    base.update(kwargs)
    return base


def test_gene_label_parses_symbol_and_entrez():
    assert parse_gene_label("TP53 (7157)") == ("TP53", "7157")
    assert parse_gene_label("A1BG (1)") == ("A1BG", "1")
    assert parse_gene_label("not a gene") is None


def test_rrid_disagreement_is_not_guessed():
    rows = build_model_xref(
        [cmp_row(model_id="SIDM00001", RRID="CVCL_0002", BROAD_ID="ACH-000001", model_name="Alpha")],
        [depmap_row(ModelID="ACH-000001", RRID="CVCL_0001", SangerModelID="SIDM00001", CellLineName="Alpha")],
    )
    primary = next(row for row in rows if row["source"] == "depmap" and row["alias_type"] == "primary")
    assert primary["status"] == "ambiguous"
    assert primary["cvcl"] is None
    assert primary["match_method"] == "conflict_identifier"


def test_unique_name_maps_and_shared_name_does_not():
    cmp = [
        cmp_row(model_id="SIDM00010", RRID="CVCL_0010", model_name="OnlyOne"),
        cmp_row(model_id="SIDM00011", RRID="CVCL_0011", model_name="Shared"),
        cmp_row(model_id="SIDM00012", RRID="CVCL_0012", model_name="Shared"),
    ]
    rows = build_model_xref(cmp, [
        depmap_row(ModelID="ACH-000010", CellLineName="OnlyOne"),
        depmap_row(ModelID="ACH-000011", CellLineName="Shared"),
        depmap_row(ModelID="ACH-000012", SangerModelID="SIDM00999", CellLineName="OnlyOne"),
    ])
    by_id = {row["source_model_id"]: row for row in rows if row["source"] == "depmap" and row["alias_type"] == "primary"}
    assert by_id["ACH-000010"]["status"] == "mapped"
    assert by_id["ACH-000010"]["cvcl"] == "CVCL_0010"
    assert by_id["ACH-000010"]["match_method"] == "name_unique"
    assert by_id["ACH-000011"]["status"] == "ambiguous"
    assert by_id["ACH-000011"]["cvcl"] is None
    assert by_id["ACH-000012"]["status"] == "unmatched"
    assert by_id["ACH-000012"]["match_method"] == "unmatched_identifier"


def test_lettered_cellosaurus_ids_map_and_two_rrids_do_not():
    rows = build_model_xref(
        [
            cmp_row(model_id="SIDM00020", RRID="CVCL_A472", model_name="Lettered"),
            cmp_row(model_id="SIDM00021", RRID="CVCL_2717;CVCL_1888", BROAD_ID="ACH-000021", model_name="Split"),
        ],
        [depmap_row(ModelID="ACH-000020", RRID="CVCL_A472", CellLineName="Lettered")],
    )
    depmap = next(row for row in rows if row["source"] == "depmap" and row["alias_type"] == "primary")
    assert depmap["status"] == "mapped"
    assert depmap["cvcl"] == "CVCL_A472"
    split = next(row for row in rows if row["source_model_id"] == "SIDM00021" and row["alias_type"] == "primary")
    assert split["status"] == "ambiguous"
    assert split["cvcl"] is None


def test_cmp_sidm_maps_when_depmap_has_no_rrid():
    cmp = index_cmp([cmp_row(model_id="SIDM00003", RRID="CVCL_0003", BROAD_ID="ACH-000099")])
    decision = classify_depmap(depmap_row(ModelID="ACH-000003", SangerModelID="SIDM00003"), cmp)
    assert decision.status == "mapped"
    assert decision.cvcl == "CVCL_0003"
    assert decision.match_method == "cmp_sidm"
