import csv
import math
import shutil

import pytest

from pgx.normalize.cell_lines import CellosaurusIndex
from pgx.process import _write_rows, connect
from pgx.pset import collect_source, compound_key, extract_rds, load_source


CELLOSAURUS = """ID   MCF7
AC   CVCL_0031
//
ID   K562
AC   CVCL_0004
//
"""


def test_numeric_pset_ids_are_namespaced():
    assert compound_key("ctrpv2", "1788") == "ctrpv2:1788"
    assert compound_key("ctrpv2", "BRD-K46556387") == "BRD-K46556387"


def test_ctrp_uses_source_cellosaurus_and_one_pubchem_cid(tmp_path):
    raw = tmp_path / "raw"
    root = raw / "orcestra_ctrpv2" / "2.1"
    root.mkdir(parents=True)
    (root / "sample.tsv").write_text(
        "master_ccl_id\tCELL_NAME\tccle_primary_site\tCellosaurus.Accession\tIn.CTRPv2.1\n"
        "1\tMCF7\tbreast\tCVCL_0031\tTRUE\n"
        "2\tK562\thaematopoietic\tCVCL_0031\tTRUE\n"
        "9\tUNUSED\tlung\tCVCL_0004\tFALSE\n"
    )
    (root / "treatment.tsv").write_text(
        "master_cpd_id\tCompound.Name\tBroad.Compound.ID\tSMILES\tAnnotationDB.PubChem.CID\tAnnotationGx.PubChem.CID\tAnnotationGx.InChIKey\tAnnotationDB.ChEMBL.ID\n"
        "1788\tCIL55\tBRD-K46556387\tC\t\t6623618\tITESNYBWKVFAEV-UHFFFAOYSA-N\tCHEMBL1\n"
        "2000\tSPLIT\t\tC\t1\t2\t\t\n"
    )
    (root / "v2.1_auc_sensitivities.tsv").write_text(
        "experiment_id\tarea_under_curve\tmaster_cpd_id\tmaster_ccl_id\n"
        "7\t14.5\t1788\t1\n"
        "8\t3.0\t1788\t2\n"
        "9\t1.0\t2000\t1\n"
    )
    (root / "sensitivity.tsv").write_text(
        "response_id\tmaster_ccl_id\tdose_umol_1\tdose_umol_2\tviability_1\tviability_2\n"
        "CTRPv2.0_EXP_7_CPD_1788\t1\t0.01\t1\t90\t10\n"
    )
    (root / "recomputed_profiles.tsv").write_text(
        "response_id\tic50_recomputed\tauc_recomputed\n"
        "CTRPv2.0_EXP_7_CPD_1788\t0.2\t14.5\n"
    )
    cellosaurus = CellosaurusIndex.load(_write_cellosaurus(tmp_path))
    source = {
        "id": "orcestra_ctrpv2",
        "version": "2.1",
        "dataset": "ctrpv2",
        "enabled": True,
        "files": [
            {"name": "sample.tsv", "role": "pset_cells"},
            {"name": "treatment.tsv", "role": "pset_drugs"},
            {"name": "v2.1_auc_sensitivities.tsv", "role": "pset_profiles"},
            {"name": "sensitivity.tsv", "role": "pset_raw"},
            {"name": "recomputed_profiles.tsv", "role": "pset_recomputed"},
        ],
    }
    con = connect()
    models, compounds = collect_source(con, raw, source, cellosaurus)
    primary = {row["source_model_id"]: row for row in models if row["alias_type"] == "primary"}
    assert set(primary) == {"1", "2"}
    assert primary["1"]["cvcl"] == "CVCL_0031"
    assert primary["1"]["status"] == "mapped"
    assert primary["1"]["match_method"] == "pset_cellosaurus"
    assert primary["2"]["status"] == "ambiguous"
    assert primary["2"]["match_method"] == "conflict_identifier"
    by_id = {record.source_compound_id: record for record in compounds}
    assert by_id["ctrpv2:1788"].pubchem_cid == "6623618"
    assert by_id["ctrpv2:1788"].status == "mapped"
    assert by_id["ctrpv2:1788"].match_method == "pset_pubchem"
    assert by_id["ctrpv2:1788"].brd_id == "BRD-K46556387"
    assert by_id["ctrpv2:2000"].status == "ambiguous"
    assert by_id["ctrpv2:2000"].pubchem_cid is None

    xref = tmp_path / "model_xref"
    _write_rows(con, models, xref, list(models[0]))
    dest = tmp_path / "drug_response"
    count = load_source(con, raw, source, xref, dest)
    rows = con.execute(
        """
        SELECT metric, source_model_id, model_cvcl, source_compound_id, value, dose
        FROM read_parquet(?)
        ORDER BY metric, dose NULLS FIRST
        """,
        [str(dest / "dataset=ctrpv2" / "*.parquet")],
    ).fetchall()
    assert count == len(rows)
    auc = [row for row in rows if row[0] == "AUC" and row[4] == 14.5]
    assert auc[0][2] == "CVCL_0031"
    assert auc[0][3] == "ctrpv2:1788"
    viability = [row for row in rows if row[0] == "VIABILITY"]
    assert [row[5] for row in viability] == [0.01, 1.0]
    ic50 = [row for row in rows if row[0] == "IC50"]
    assert ic50[0][1] == "1"
    assert ic50[0][3] == "ctrpv2:1788"
    assert ic50[0][4] == 0.2
    conflict_auc = [row for row in rows if row[1] == "2"]
    assert conflict_auc[0][2] is None


def test_extract_csvs_load_without_r(tmp_path):
    raw = tmp_path / "raw"
    extract = raw / "orcestra_gcsi" / "2019" / "extract"
    extract.mkdir(parents=True)
    (extract / "cells.csv").write_text("source_model_id,cell_name,cellosaurus,tissue\nMCF7,MCF7,CVCL_0031,breast\n")
    (extract / "drugs.csv").write_text("source_compound_id,name,pubchem_cid,pubchem_cid_alt,smiles,inchikey,brd_id,chembl_id\n99,Paclitaxel,36314,,,,,\n")
    (extract / "profiles.csv").write_text("source_model_id,source_compound_id,metric,value,unit\nMCF7,99,auc,0.4,\n")
    (extract / "raw.csv").write_text("source_model_id,source_compound_id,dose,value,metric,unit\nMCF7,99,0.1,80,VIABILITY,percent\n")
    cellosaurus = CellosaurusIndex.load(_write_cellosaurus(tmp_path))
    source = {
        "id": "orcestra_gcsi",
        "version": "2019",
        "dataset": "gcsi",
        "enabled": True,
        "files": [],
    }
    con = connect()
    models, compounds = collect_source(con, raw, source, cellosaurus)
    assert models[0]["cvcl"] == "CVCL_0031"
    assert compounds[0].source_compound_id == "gcsi:99"
    assert compounds[0].pubchem_cid == "36314"
    xref = tmp_path / "model_xref"
    _write_rows(con, models, xref, list(models[0]))
    dest = tmp_path / "drug_response"
    load_source(con, raw, source, xref, dest)
    metrics = {row[0] for row in con.execute("SELECT metric FROM read_parquet(?)", [str(dest / "dataset=gcsi" / "*.parquet")]).fetchall()}
    assert metrics == {"AUC", "VIABILITY"}, "extractor column names are stored as the standard metric names"


def test_rds_list_roundtrip(tmp_path):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not installed")
    rds = tmp_path / "toy.rds"
    script = r"""
    cells <- data.frame(cellid=c('MCF7'), cellosaurus=c('CVCL_0031'), tissueid=c('breast'), stringsAsFactors=FALSE)
    drugs <- data.frame(drugid=c('99'), name=c('Paclitaxel'), pubchem_cid=c('36314'), stringsAsFactors=FALSE)
    profiles <- data.frame(source_model_id=c('MCF7'), source_compound_id=c('99'), metric=c('auc'), value=c(0.4), unit=c(NA), stringsAsFactors=FALSE)
    saveRDS(list(cells=cells, drugs=drugs, profiles=profiles), {path!r})
    """.format(path=str(rds))
    import subprocess
    subprocess.run(["Rscript", "-e", script], check=True)
    dest = tmp_path / "extract"
    assert extract_rds(rds, dest) == "extracted"
    text = (dest / "drugs.csv").read_text()
    assert "36314" in text
    assert "Paclitaxel" in text


def _write_cellosaurus(tmp_path):
    path = tmp_path / "cellosaurus.txt"
    path.write_text(CELLOSAURUS)
    return path

LEGACY_PSET = r"""
cell <- data.frame(`Cell line`=c('MCF7','T47D'), cellid=c('MCF7','T47D'), tissueid=c('breast','breast'), check.names=FALSE, stringsAsFactors=FALSE)
rownames(cell) <- cell$cellid
drug <- data.frame(drugid=c('Taxol','Gefitinib'), cid=c('36314','123631'), stringsAsFactors=FALSE)
rownames(drug) <- drug$drugid
ids <- c('Taxol_MCF7','Taxol_T47D','Gefitinib_MCF7')
# CTRPv2 keeps the replicate number in an experiment_id column; it is not the row key.
info <- data.frame(drugid=c('Taxol','Taxol','Gefitinib'), cellid=c('MCF7','T47D','MCF7'), experiment_id=c('1','1','1'), stringsAsFactors=FALSE, row.names=ids)
profiles <- data.frame(aac_recomputed=c(0.6, 0.2, NA), ic50_recomputed=c(1.5, Inf, 3), EC50=c(2, 3, 4), aac_published=c(0.5,0.1,0.3), row.names=ids)
raw <- array(NA_real_, dim=c(3, 2, 2), dimnames=list(ids, c('doses1','doses2'), c('Dose','Viability')))
raw[1, , 'Dose'] <- c(0.1, 1); raw[1, , 'Viability'] <- c(95, 40)
raw[2, 1, 'Dose'] <- 0.1; raw[2, 1, 'Viability'] <- 99
obj <- structure(list(), class='PharmacoSet')
attr(obj, 'cell') <- cell
attr(obj, 'drug') <- drug
attr(obj, 'sensitivity') <- list(info=info, profiles=profiles, raw=raw, n=NULL)
saveRDS(obj, {path!r})
"""


def _legacy_extract(tmp_path):
    import subprocess

    rds = tmp_path / "legacy.rds"
    subprocess.run(["Rscript", "-e", LEGACY_PSET.format(path=str(rds))], check=True)
    dest = tmp_path / "extract"
    assert extract_rds(rds, dest) == "extracted"
    return dest


def test_extractor_reads_a_pre_pharmacogx3_pharmacoset_without_pharmacogx(tmp_path):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not installed")
    dest = _legacy_extract(tmp_path)
    profiles = (dest / "profiles.csv").read_text()
    assert "aac_recomputed" in profiles and "ic50_recomputed" in profiles and "ec50" in profiles
    assert "aac_published" not in profiles, "published values are not stored beside the recomputed ones"
    with (dest / "raw.csv").open() as handle:
        raw = list(csv.DictReader(handle))
    assert sorted((row["source_model_id"], row["source_compound_id"], float(row["dose"]), float(row["value"])) for row in raw) == [
        ("MCF7", "Taxol", 0.1, 95.0),
        ("MCF7", "Taxol", 1.0, 40.0),
        ("T47D", "Taxol", 0.1, 99.0),
    ], "missing dose slots are not rows"
    assert {row["metric"] for row in raw} == {"VIABILITY"}
    assert "36314" in (dest / "drugs.csv").read_text()


def test_real_pharmacoset_columns_load_as_standard_metrics_and_drop_non_finite_values(tmp_path):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript is not installed")
    dest = _legacy_extract(tmp_path)
    raw_root = tmp_path / "raw"
    extract = raw_root / "orcestra_gcsi" / "2019" / "extract"
    extract.parent.mkdir(parents=True)
    dest.rename(extract)
    cellosaurus = CellosaurusIndex.load(_write_cellosaurus(tmp_path))
    source = {"id": "orcestra_gcsi", "version": "2019", "dataset": "gcsi", "enabled": True, "files": []}
    con = connect()
    models, _compounds = collect_source(con, raw_root, source, cellosaurus)
    xref = tmp_path / "model_xref"
    _write_rows(con, models, xref, list(models[0]))
    out = tmp_path / "drug_response"
    load_source(con, raw_root, source, xref, out)
    rows = con.execute(
        "SELECT metric, source_model_id, source_compound_id, value FROM read_parquet(?) ORDER BY 1, 2, 3",
        [str(out / "dataset=gcsi" / "profiles.parquet")],
    ).fetchall()
    metrics = {row[0] for row in rows}
    assert metrics == {"AAC", "IC50", "EC50"}
    assert sorted(row[3] for row in rows if row[0] == "IC50") == [1.5, 3.0], "an IC50 of Inf is not a value"
    assert sorted(row[3] for row in rows if row[0] == "AAC") == [0.2, 0.6], "an NA AAC is not a row either"
    assert all(math.isfinite(row[3]) for row in rows)
