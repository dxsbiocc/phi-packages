from pgx.download import load_sources
from pgx.process import process_all
from pgx.validate import validate


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def fixture_tree(tmp_path):
    raw = tmp_path / "raw"
    write(raw / "sanger_cmp/2026-10-14/model_list_latest.csv", """model_id,RRID,BROAD_ID,COSMIC_ID,model_name,synonyms,CCLE_ID,species,tissue,cancer_type
SIDM00001,CVCL_0001,ACH-000001,100,Alpha,,ALPHA,Homo sapiens,Lung,Carcinoma
SIDM00002,CVCL_0002,ACH-000002,200,Beta,,BETA,Homo sapiens,Lung,Carcinoma
SIDM00003,CVCL_0003,ACH-000003,300,Gamma,,GAMMA,Homo sapiens,Lung,Carcinoma
SIDM00004,CVCL_0004,ACH-000004,400,Delta,,DELTA,Homo sapiens,Lung,Carcinoma
SIDM00011,CVCL_0011,ACH-000011,111,Epsilon,,EPSILON,Homo sapiens,Lung,Carcinoma
SIDM00012,CVCL_0012,ACH-000012,112,Zeta,,ZETA,Homo sapiens,Lung,Carcinoma
SIDM00021,CVCL_0021,ACH-000021,121,Eta,,ETA,Homo sapiens,Lung,Carcinoma
SIDM00022,CVCL_0022,ACH-000022,122,Theta,,THETA,Homo sapiens,Lung,Carcinoma
""")
    write(raw / "depmap/24Q2/Model.csv", """ModelID,RRID,SangerModelID,COSMICID,CellLineName,StrippedCellLineName,CCLEName,OncotreeLineage,OncotreePrimaryDisease
ACH-000001,CVCL_0001,SIDM00001,100,Alpha,ALPHA,ALPHA,Lung,Carcinoma
ACH-000002,CVCL_0002,SIDM00002,200,Beta,BETA,BETA,Lung,Carcinoma
ACH-000003,CVCL_9999,SIDM00003,300,Gamma,GAMMA,GAMMA,Lung,Carcinoma
ACH-000099,,SIDM00999,,Mystery,MYSTERY,MYSTERY,Lung,Carcinoma
ACH-000011,CVCL_0011,SIDM00011,111,Epsilon,EPSILON,EPSILON,Lung,Carcinoma
ACH-000012,CVCL_0012,SIDM00012,112,Zeta,ZETA,ZETA,Lung,Carcinoma
ACH-000021,CVCL_0021,SIDM00021,121,Eta,ETA,ETA,Lung,Carcinoma
ACH-000022,CVCL_0022,SIDM00022,122,Theta,THETA,THETA,Lung,Carcinoma
""")
    write(raw / "depmap/24Q2/CRISPRGeneEffect.csv", """,TP53 (7157),MDM2 (4193)
ACH-000001,-1.2,0.1
ACH-000002,NA,-0.4
ACH-000099,-0.2,0.2
""")
    write(raw / "depmap/24Q2/OmicsSomaticMutations.csv", """Chrom,Pos,Ref,Alt,ProteinChange,HugoSymbol,MolecularConsequence,VepImpact,ModelID,EntrezGeneID
17,7676154,G,A,p.R175H,TP53,missense_variant,HIGH,ACH-000002,7157
17,100,A,T,p.Syn,TP53,synonymous_variant,LOW,ACH-000001,7157
12,1,A,T,p.X,MDM2,missense_variant,MODERATE,ACH-000011,4193
12,2,A,T,p.X,MDM2,missense_variant,MODERATE,ACH-000012,4193
17,3,G,A,p.R248Q,TP53,missense_variant,HIGH,ACH-000021,7157
17,4,G,A,p.R273H,TP53,missense_variant,HIGH,ACH-000022,7157
""")
    write(raw / "gdsc/8.4/screened_compounds_rel_8.4.csv", """DRUG_ID,SCREENING_SITE,DRUG_NAME,SYNONYMS,TARGET,TARGET_PATHWAY,PUBCHEM
1047,SANGER,Nutlin-3a (-),,MDM2,p53 pathway,11433190
1,MGH,Erlotinib,Tarceva,EGFR,EGFR signaling,
""")
    write(raw / "gdsc/8.4/GDSC1_fitted_dose_response_24Jul22.csv", """DATASET,COSMIC_ID,CELL_LINE_NAME,SANGER_MODEL_ID,DRUG_ID,DRUG_NAME,LN_IC50,AUC
GDSC1,100,Alpha,SIDM00001,1047,Nutlin-3a (-),0.2,0.10
GDSC1,200,Beta,SIDM00002,1047,Nutlin-3a (-),2.5,0.80
GDSC1,300,Gamma,SIDM00003,1047,Nutlin-3a (-),0.3,0.20
GDSC1,400,Delta,SIDM00004,1047,Nutlin-3a (-),3.0,0.90
GDSC1,111,Epsilon,SIDM00011,1047,Nutlin-3a (-),0.15,0.11
GDSC1,112,Zeta,SIDM00012,1047,Nutlin-3a (-),0.25,0.12
GDSC1,121,Eta,SIDM00021,1047,Nutlin-3a (-),2.2,0.70
GDSC1,122,Theta,SIDM00022,1047,Nutlin-3a (-),2.8,0.85
""")
    write(raw / "prism_19q4/19Q4/primary-screen-cell-line-info.csv", """depmap_id,ccle_name,primary_tissue
ACH-000001,ALPHA,lung
""")
    write(raw / "prism_19q4/19Q4/primary-screen-replicate-collapsed-treatment-info.csv", """column_name,broad_id,name,dose,screen_id,smiles
BRD-A12230535-001-06-7::2.5::HTS,BRD-A12230535-001-06-7,nutlin-3,2.5,HTS,C
""")
    write(raw / "prism_19q4/19Q4/primary-screen-replicate-collapsed-logfold-change.csv", """,BRD-A12230535-001-06-7::2.5::HTS
ACH-000001,-0.8
ACH-000002,-0.1
""")
    document = {
        "sources": [
            {"id": "sanger_cmp", "version": "2026-10-14", "enabled": True, "batch": 1, "license": "cmp", "citation": "cmp", "homepage": "https://cellmodelpassports.sanger.ac.uk/", "files": [
                {"name": "model_list_latest.csv", "role": "cmp_models"},
            ]},
            {"id": "depmap", "version": "24Q2", "enabled": True, "batch": 1, "license": "CC BY 4.0", "citation": "depmap", "homepage": "https://depmap.org/portal/", "files": [
                {"name": "Model.csv", "role": "depmap_models"},
                {"name": "CRISPRGeneEffect.csv", "role": "gene_matrix", "dataset": "depmap_crispr_24q2", "metric": "chronos_gene_effect"},
                {"name": "OmicsSomaticMutations.csv", "role": "depmap_mutations", "dataset": "depmap_mutations_24q2"},
            ]},
            {"id": "gdsc", "version": "8.4", "enabled": True, "batch": 1, "license": "gdsc", "citation": "gdsc", "homepage": "https://www.cancerrxgene.org/", "files": [
                {"name": "screened_compounds_rel_8.4.csv", "role": "gdsc_compounds"},
                {"name": "GDSC1_fitted_dose_response_24Jul22.csv", "role": "gdsc_dose_response"},
            ]},
            {"id": "prism_19q4", "version": "19Q4", "enabled": True, "batch": 1, "license": "CC BY 4.0", "citation": "prism", "homepage": "https://figshare.com/", "files": [
                {"name": "primary-screen-cell-line-info.csv", "role": "prism_cell_meta"},
                {"name": "primary-screen-replicate-collapsed-treatment-info.csv", "role": "prism_primary_compounds"},
                {"name": "primary-screen-replicate-collapsed-logfold-change.csv", "role": "prism_primary_lfc", "dataset": "prism_primary_19q4"},
            ]},
        ]
    }
    return document, raw


def test_process_harmonizes_and_nutlin_check_passes(tmp_path):
    document, raw = fixture_tree(tmp_path)
    counts = process_all(document, raw, tmp_path / "processed")
    assert counts["model_mapped"] >= 3
    assert counts["gene_effect"] == 5
    report = validate(tmp_path / "processed")
    depmap = report["coverage"]["models"]["depmap"]
    assert depmap["mapped"] == 6
    assert depmap["ambiguous"] == 1
    assert depmap["unmatched"] == 1
    assert report["coverage"]["compounds"]["gdsc"]["mapped"] == 1
    assert report["coverage"]["compounds"]["gdsc"]["unmatched"] == 1
    assert report["nutlin_tp53"]["status"] == "pass"
    assert report["coverage"]["tables"]["drug_response"]["with_cvcl"] == report["coverage"]["tables"]["drug_response"]["rows"]


def test_inverted_nutlin_fails(tmp_path):
    document, raw = fixture_tree(tmp_path)
    write(raw / "gdsc/8.4/GDSC1_fitted_dose_response_24Jul22.csv", """DATASET,COSMIC_ID,CELL_LINE_NAME,SANGER_MODEL_ID,DRUG_ID,DRUG_NAME,LN_IC50,AUC
GDSC1,100,Alpha,SIDM00001,1047,Nutlin-3a (-),4.0,0.90
GDSC1,200,Beta,SIDM00002,1047,Nutlin-3a (-),0.1,0.10
GDSC1,300,Gamma,SIDM00003,1047,Nutlin-3a (-),3.5,0.80
GDSC1,400,Delta,SIDM00004,1047,Nutlin-3a (-),0.2,0.20
GDSC1,111,Epsilon,SIDM00011,1047,Nutlin-3a (-),4.2,0.95
GDSC1,112,Zeta,SIDM00012,1047,Nutlin-3a (-),4.4,0.96
GDSC1,121,Eta,SIDM00021,1047,Nutlin-3a (-),0.05,0.05
GDSC1,122,Theta,SIDM00022,1047,Nutlin-3a (-),0.08,0.06
""")
    process_all(document, raw, tmp_path / "processed")
    from pgx.validate import ValidationError
    import pytest
    with pytest.raises(ValidationError):
        validate(tmp_path / "processed")


def test_packaged_sources_match_the_loader():
    document = load_sources(__import__("pathlib").Path(__file__).resolve().parents[1] / "sources.yaml")
    roles = {spec["role"] for source in document["sources"] for spec in source["files"]}
    assert {"depmap_models", "cmp_models", "gene_matrix", "gdsc_dose_response", "prism_primary_lfc", "prism_24q2_lfc"} <= roles
