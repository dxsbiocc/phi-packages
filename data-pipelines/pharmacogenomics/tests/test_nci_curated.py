from pgx.nci60 import curated_decision, collect_models, load_curated
from pgx.normalize.cell_lines import CellosaurusIndex
from pgx.process import connect


PANEL = """ID   COLO 205
AC   CVCL_0218
SY   COLO205
CC   Part of: NCI-60 cancer cell line panel.
//
ID   COLO205
AC   CVCL_F402
//
ID   HT-29
AC   CVCL_0320
SY   HT29
CC   Part of: NCI-60 cancer cell line panel.
//
ID   HT29
AC   CVCL_A8EZ
//
ID   M14
AC   CVCL_1395
CC   Part of: NCI-60 cancer cell line panel.
//
ID   M14K
AC   CVCL_0001
SY   M14
//
ID   MDA-MB-435
AC   CVCL_0417
CC   Problematic cell line: Contaminated. Misidentified. Derived from M14. ICLAC-00235.
CC   Part of: NCI-60 cancer cell line panel.
//
ID   PC-3
AC   CVCL_0035
CC   Part of: NCI-60 cancer cell line panel.
//
ID   PaCa-3
AC   CVCL_9999
SY   PC-3
//
ID   TK-10
AC   CVCL_1773
CC   Part of: NCI-60 cancer cell line panel.
//
ID   UCH-ICF-TK-10
AC   CVCL_E5GS
SY   TK-10
//
ID   A-204
AC   CVCL_1058
//
ID   Urd-B
AC   CVCL_1N55
SY   A-204
//
ID   NB4
AC   CVCL_0005
//
ID   SJNB-4
AC   CVCL_8821
SY   NB4
//
ID   ND-1
AC   CVCL_4781
CC   Problematic cell line: Contaminated. Misidentified. ICLAC-00244.
//
ID   UABMEL3
AC   CVCL_0002
//
"""

EXPECTED = {
    "A-204": "CVCL_1058",
    "COLO 205": "CVCL_0218",
    "HT29": "CVCL_0320",
    "M14": "CVCL_1395",
    "NB4": "CVCL_0005",
    "ND-1": "CVCL_4781",
    "PC-3": "CVCL_0035",
    "TK-10": "CVCL_1773",
}


def test_shipped_table_resolves_the_ambiguous_dtp_names():
    rows = load_curated()
    for name, cvcl in EXPECTED.items():
        assert rows[name]["cvcl"] == cvcl
    assert rows["MDA-MB-435"]["cvcl"] == "CVCL_0417"
    assert "ICLAC-00235" in rows["MDA-MB-435"]["notes"]
    assert "ICLAC-00244" in rows["ND-1"]["notes"]
    assert rows["M14"]["panel_role"] == "standard"
    assert rows["A-204"]["panel_role"] == "historical"


def test_curated_row_wins_and_keeps_the_cellosaurus_flag(tmp_path):
    path = tmp_path / "cellosaurus.txt"
    path.write_text(PANEL)
    index = CellosaurusIndex.load(path)
    assert index.decide_names(["COLO 205"]).status == "ambiguous"
    assert index.decide_names(["M14"]).status == "ambiguous"
    assert index.decide_names(["MDA-MB-435"]).cvcl == "CVCL_0417"
    assert index.flag("CVCL_0417") == "misidentified"
    assert index.flag("CVCL_1395") is None
    assert index.flag("CVCL_4781") == "misidentified"

    curated = load_curated()
    chosen = {name: curated[name] for name in [*EXPECTED, "MDA-MB-435"]}
    table = tmp_path / "cells.csv"
    names = [*EXPECTED, "MDA-MB-435", "UABMEL3"]
    table.write_text("CELL_NAME,PANEL_NAME\n" + "\n".join(f"{name},Panel" for name in names) + "\n")
    models = collect_models(connect(), [table], index, curated=chosen)
    primary = {row["source_model_id"]: row for row in models if row["alias_type"] == "primary"}
    for name, cvcl in EXPECTED.items():
        assert primary[name]["cvcl"] == cvcl
        assert primary[name]["match_method"] == "curated"
        assert primary[name]["status"] == "mapped"
    assert primary["M14"]["cellosaurus_flag"] is None
    assert primary["MDA-MB-435"]["cvcl"] == "CVCL_0417"
    assert primary["MDA-MB-435"]["match_method"] == "curated"
    assert primary["MDA-MB-435"]["cellosaurus_flag"] == "misidentified"
    assert primary["ND-1"]["cellosaurus_flag"] == "misidentified"
    assert primary["UABMEL3"]["cvcl"] == "CVCL_0002"
    assert primary["UABMEL3"]["match_method"] == "cellosaurus_name"


def test_curated_decision_keeps_a_missing_accession_mapped():
    decision = curated_decision({"cvcl": "CVCL_0000", "notes": "historical"}, cellosaurus=None)
    assert decision.cvcl == "CVCL_0000"
    assert decision.match_method == "curated"
    assert decision.status == "mapped"
    assert decision.cellosaurus_flag is None
