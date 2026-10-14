import zipfile

from pgx.ids import CompoundRecord
from pgx.nci60 import collect_compounds, collect_models, load_responses, nsc_id, open_table
from pgx.normalize.cell_lines import CellosaurusIndex
from pgx.normalize.compounds import iter_queries
from pgx.process import _write_rows, connect


CELLOSAURUS = """ID   NCI-H23
AC   CVCL_1547
SY   H23
//
ID   MDA-MB-435
AC   CVCL_0417
CC   Problematic cell line: Misidentified.
//
"""

GI50 = """RELEASE_DATE,EXPID,PREFIX,NSC,CONCENTRATION_UNIT,LOG_HI_CONCENTRATION,PANEL_NUMBER,CELL_NUMBER,PANEL_NAME,CELL_NAME,PANEL_CODE,COUNT,AVERAGE,STDDEV
20210223,0001MD02,S,740,M,-4,1,1,Non-Small Cell Lung Cancer,NCI-H23,LNS,1,-7.1391,0
20210223,0001MD02,S,740,M,-4,1,2,Melanoma,MDA-MB-435,MEL,1,-6.5,0
20210223,0001MD02,X,741,M,-4,1,1,Non-Small Cell Lung Cancer,NCI-H23,LNS,1,-1,0
"""

ONECONC = """RELEASE_DATE,EXPID,PREFIX,NSC,CONCENTRATION_UNIT,CONCENTRATION,PANEL_NUMBER,CELL_NUMBER,PANEL_NAME,CELL_NAME,PANEL_CODE,COUNT_GIPRCNT,AVERAGE_GIPRCNT,STDDEV_GIPRCNT
20210223,0001OM01,S,740,M,.000025,1,1,Non-Small Cell Lung Cancer,NCI-H23,LNS,1,4.3896,0
"""

DOSERESP = """RELEASE_DATE,EXPID,PREFIX,NSC,CONCENTRATION_UNIT,CONCENTRATION,CELL_NAME,AVERAGE_GIPRCNT,AVERAGE_PTC
20210223,0001DS01,S,740,M,-6,NCI-H23,80,90
20210223,0001DS01,S,740,M,-5,NCI-H23,20,40
"""


def _zip(path, name, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(name, text)


def test_nsc_lookup_precedes_the_stored_name():
    record = CompoundRecord(source="nci60", source_compound_id="NSC740", name="methotrexate")
    queries = iter_queries(record)
    assert [item[1] for item in queries[:2]] == ["NSC740", "NSC 740"]
    assert queries[0][2] == "pubchem_nsc"
    assert queries[0][3] == "high"
    assert "methotrexate" in [item[1] for item in queries[2:]]
    assert nsc_id("740") == "NSC740"
    assert nsc_id("NSC 740") == "NSC740"


def test_nci60_endpoints_oneconc_and_doseresp(tmp_path):
    raw = tmp_path / "raw"
    _zip(raw / "GI50.zip", "GI50.csv", GI50)
    _zip(raw / "ONECONC.zip", "ONECONC.csv", ONECONC)
    doseresp = raw / "DOSERESP.csv"
    doseresp.write_text(DOSERESP)
    cellosaurus_path = tmp_path / "cellosaurus.txt"
    cellosaurus_path.write_text(CELLOSAURUS)
    index = CellosaurusIndex.load(cellosaurus_path)
    con = connect()
    tables = [open_table(raw / "GI50.zip"), open_table(raw / "ONECONC.zip"), doseresp]
    models = collect_models(con, tables, index)
    primary = [row for row in models if row["alias_type"] == "primary"]
    by_name = {row["source_model_id"]: row for row in primary}
    assert by_name["NCI-H23"]["cvcl"] == "CVCL_1547"
    assert by_name["NCI-H23"]["status"] == "mapped"
    assert by_name["NCI-H23"]["match_method"] == "curated"
    assert by_name["MDA-MB-435"]["cvcl"] == "CVCL_0417"
    assert by_name["MDA-MB-435"]["match_method"] == "curated"
    assert by_name["MDA-MB-435"]["cellosaurus_flag"] == "misidentified"
    compounds = collect_compounds(con, tables)
    assert {record.source_compound_id for record in compounds} == {"NSC740"}

    xref = tmp_path / "model_xref"
    _write_rows(con, models, xref, list(models[0]))
    dest = tmp_path / "drug_response"
    count = load_responses(
        con,
        [("GI50.zip", "nci60_endpoint", tables[0]), ("ONECONC.zip", "nci60_oneconc", tables[1]), ("DOSERESP.csv", "nci60_doseresp", tables[2])],
        xref,
        dest,
    )
    rows = con.execute(
        """
        SELECT metric, source_compound_id, model_cvcl, value, unit, dose, response_class, combination_id
        FROM read_parquet(?)
        ORDER BY metric, dose NULLS FIRST
        """,
        [str(dest / "dataset=nci60" / "*.parquet")],
    ).fetchall()
    assert count == len(rows) == 7
    gi50 = [row for row in rows if row[0] == "GI50" and row[2] == "CVCL_1547"]
    assert len(gi50) == 1
    assert gi50[0][1] == "NSC740"
    assert gi50[0][3] == -7.1391
    assert gi50[0][4] == "log10_M"
    assert gi50[0][6] == "single"
    assert gi50[0][7] is None
    one = [row for row in rows if row[0] == "ONECONC_GIPRCNT"][0]
    assert abs(one[5] - 0.000025) < 1e-12
    doses = sorted(row[5] for row in rows if row[0] == "GIPRCNT")
    assert doses[0] == 1e-6
    assert {row[0] for row in rows if row[0] in {"GIPRCNT", "PTC"}} == {"GIPRCNT", "PTC"}
