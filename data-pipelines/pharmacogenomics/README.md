# Pharmacogenomics pipeline

Batch 1 turns DepMap, PRISM, GDSC, and Cell Model Passports releases into partitioned Parquet. The MCP server in `mcp-servers/pharmacogenomics` reads that directory. It does not download the releases itself.

```bash
pip install -e '.[test]'
pytest
phi-pgx download          # small metadata and crosswalks
phi-pgx download --live   # CRISPR, expression, mutations, copy number, dose-response matrices
phi-pgx process
phi-pgx validate
PHI_PGX_DATA_DIR=processed phi-omics-pharmacogenomics
```

`download` without `--live` skips files marked `live: true`. Those matrices are hundreds of megabytes to about 800 MB. `PHI_PGX_LIVE=1` is the same switch. Tests use fixtures and do not download them.

## Sources

| Source | Version | What batch 1 uses |
| --- | --- | --- |
| DepMap | 24Q2 | `Model.csv`, `CRISPRGeneEffect.csv`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, `OmicsCNGene.csv`, `OmicsSomaticMutations.csv` |
| PRISM | 19Q4 | Primary-screen LFC and secondary-screen dose-response curve parameters, plus cell and compound metadata. This is the Corsello 2020 release on Figshare, not an ORCESTRA PRISM PSet. |
| PRISM | Public 24Q2 | Collapsed LFC and compound metadata from the DepMap Figshare record |
| GDSC | release 8.4 | `screened_compounds_rel_8.4.csv`, GDSC1 and GDSC2 fitted dose response dated 24 Jul 2022. Those fitted files are single-agent matrices. URLs are pinned to `release-8.4`. |
| GDSC² combinations | Jaaks 2022 | Anchored combination screens, separate from GDSC1/GDSC2. Glossary, validation screen, and original screen from the Sanger Figshare records. The original screen is `live: true`. The license is internal research and education, not CC BY. |
| Cell Model Passports | model list | `model_list_latest.csv` for SIDM, BROAD/ACH, COSMIC, and RRID |
| HGNC | complete set | Approved symbols, previous symbols, aliases, Entrez ids, and Ensembl ids |
| Cellosaurus | release 56.0 | Accession, name, synonym, and problematic-cell-line comments |

DepMap's catalogue endpoints `https://depmap.org/portal/api/download/files` and `https://depmap.org/portal/api/no-captcha/download/files` returned an HTTP 200 Cloudflare HTML page from the network used to build this pipeline. Figshare article records answered, so `sources.yaml` stores those download URLs. Each DepMap and PRISM source also lists resolvers. `phi-pgx download --resolve-urls` tries the portal catalogue first and, if that body is HTML or has no URL, uses the Figshare article. A file's `url` is the reproducible default.

The downloader writes `raw/<source>/<version>/` and `sha256sum.json`. A finished file whose hash matches the manifest is skipped. A `.partial` file is resumed with HTTP `Range` when the server honors it. An HTML body is rejected, so a bot-check page is not stored as a CSV. Optional `sha256` values in `sources.yaml` fail the download on mismatch.

## Identifiers

Normalization in `pgx.normalize` runs before identifier mapping and is shared by DepMap, PRISM, GDSC, and Cell Model Passports. One module handles each entity. Unmatched and ambiguous rows are written to `review/` and are not guessed.

### Genes

DepMap headers such as `TP53 (7157)` are split into a source symbol and an Entrez id. A header can also carry no id, several ids (`7157;7158` or `7157, ENSG00000141510`), or a token that is not a gene symbol. The HGNC complete set (`sources.yaml`, `hgnc`) maps an Entrez id, then an approved symbol, then a previous symbol, then an alias. The approved symbol is stored in `symbol`. The header symbol stays in `source_symbol`, with `match_method` set to `hgnc_entrez`, `hgnc_symbol`, `hgnc_previous`, or `hgnc_alias`. An HGNC Ensembl id is kept when the header does not disagree with it. Several Entrez ids, or a symbol and an Entrez id that point at different genes, are `ambiguous`.

### Cell lines

Cell lines are keyed by a Cellosaurus accession (`CVCL_` plus four letters or digits). Names are case-folded and stripped of punctuation, spaces, hyphens, and underscores, so `NCI-H460` and `NCIH460` are the same key. `H460` is not: the normalizer does not drop an `NCI` prefix. `H460` matches NCI-H460 only when Cellosaurus lists it as a synonym (release 56 does, as `CVCL_0459`). A DepMap `CCLEName` tissue suffix such as `_LUNG` is removed for a second lookup key; the original name is kept. An explicit RRID wins over a name. A name or synonym is used only when the row has no RRID, SIDM, COSMIC, or ACH id that failed to resolve, and only when that key hits exactly one accession. Cellosaurus `Problematic cell line` comments set `cellosaurus_flag` to `misidentified`, `contaminated`, or `problematic` and the row stays mapped. Two accessions, or an RRID that disagrees with the name, stay `ambiguous`.

Without `cellosaurus.txt`, name matching falls back to Cell Model Passports names, with the same rule that a stable id blocks the name fallback.

### Compounds

Compounds are keyed by PubChem CID. Source names are whitespace- and Unicode-normalized for lookup. Stereo and enantiomer parentheticals (`Nutlin-3a (-)`) and trailing salt tokens are removed only from the lookup string; `name` keeps the source text. Synonym lists are split on `;`, `|`, and commas. GDSC `DRUG_ID` is stored as `gdsc_drug_id` and is not sent to PubChem as a name. An NCI NSC number is stored as `NSC740` (never the bare integer, which would collide with a GDSC drug id) and is queried as the PubChem synonyms `NSC740` and `NSC 740` after a structure and before the stored name. A PharmacoSet drug id that is only digits is stored as `{dataset}:{id}`, for example `ctrpv2:1788`. BRD ids are a last lookup when the name and any SMILES or InChIKey missed.

`phi-pgx resolve-compounds` prefers a source InChIKey, then SMILES, then an NSC synonym, then the exact name, then the stereo-stripped name, then the salt-stripped name, then a dose parenthetical stripped from the lookup only (`Bleomycin (50 uM)` is looked up as `Bleomycin`). A CID already present on the source row, including a PharmacoSet PubChem CID, is kept and is not sent to PubChem again. A comma-separated synonym list in the name (`Brivanib, BMS-540215`) is split for lookup. A comma between digits is a locant and is not split. The first query that returns one CID wins and records `match_method` plus `confidence` (`high` for structure, `medium` for an exact name or BRD, `low` for a stripped name). More than one CID is `ambiguous` and a looser query is not tried. Names that are still unmatched or ambiguous are looked up as an exact ChEMBL synonym, then an exact preferred name. One ChEMBL molecule with a standard InChIKey is resolved to a PubChem CID. More than one ChEMBL molecule stops. A parenthetical counterion such as `(TFA)` is not treated as a second drug.

### Drug combinations

A combination is not a compound. Two GDSC annotation names, `BAY-MPS-combo-1 (paclitaxel 5 uM)` and `BAY-MPS-combo 2 (paclitaxel 1 uM)`, and every GDSC² `COMBI_ID` (`anchor id:library id`) become a `drug_combination` record. The id is `gdsc.combination:<anchor>:<library>` or `gdsc.combination:annotation:<drug id>`. Components stay in anchor then library order. A `|` in a component id, such as `1032|1372` (`Afatinib|Trametinib`), is two drugs sharing that role, still in source order, and each is mapped on its own. The combination name is never sent to PubChem.

`drug_response` rows for those screens use `response_class=combination` and `combination_id`. `pubchem_cid` on those rows is null. The metrics are the screened deltas `SYNERGY_DELTA_AUC`, `SYNERGY_DELTA_EMAX`, `SYNERGY_DELTA_XMID`, `SYNERGY_DELTA_XMID_uM`, and `SYNERGY_DELTA_fAUC`, not the library single-agent AUC. The anchor and library concentrations for that row are stored beside the metric. Single-agent GDSC1/GDSC2 `AUC` and `LN_IC50` rows keep `combination_id` null. `compound_response` returns single-agent rows unless `combinations` is `include` or `only`, or the query is a combination id or an `Anchor + Library` name.

## Tables

Processed output:

- `model_xref/`, `compound_xref/`, and `combination_xref/`
- `gene_effect/`, `expression/`, `cnv/`, `mutation/`, `drug_response/` partitioned as `dataset=<name>/`
- `releases/` with version, license, and citation
- `gene_xref/` when the HGNC complete set and a gene matrix or mutation file are present
- `review/models/`, `review/compounds/`, and `review/genes/` for unresolved primary ids
- `validation/report.json` after `phi-pgx validate`

`drug_response.metric` is the published statistic: GDSC single-agent `AUC` and `LN_IC50` (natural log of micromolar; lower means a stronger response), GDSC² combination deltas (`SYNERGY_DELTA_AUC` and the other delta columns above), PRISM secondary `AUC`, `IC50`, and `EC50`, and PRISM primary / 24Q2 `LFC`. Combination rows are a separate response class from single-agent rows. Expression is `log2_tpm_plus_1` from the DepMap `TPMLogp1` matrix, not raw TPM. Copy number is the relative gene-level value in `OmicsCNGene.csv`. Mutations are the DepMap hg38 calls.

`validate` checks required columns and writes mapping coverage. The Nutlin-3a check expects a lower median GDSC `LN_IC50` and `AUC` in TP53 wild-type lines than in lines with a non-silent `ncbigene:7157` mutation. Groups smaller than 3 models are skipped. Lines absent from the mutation table are excluded. A failed check exits non-zero.

## Coverage on 2026-10-14

Gene headers were read from DepMap 24Q2 without downloading the matrices. Before normalization, the old `SYMBOL (digits)` regex kept every label. After the HGNC complete set (approved symbols only), Entrez ids that are not in that set stay unmatched, and a header whose symbol and Entrez id point at different genes stays ambiguous. Mapped rows whose header symbol is not the approved symbol keep the header text in `source_symbol`.

| Matrix | Labels | Before (regex kept) | After mapped | Unmatched | Ambiguous | Symbol updated |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CRISPRGeneEffect | 18443 | 18443 | 18437 | 5 | 1 | 206 |
| Expression TPMLogp1 | 19193 | 19193 | 19162 | 17 | 14 | 393 |

Examples of updated symbols: `AATK (9625)` to `LMTK1`, `BVES (11149)` to `POPDC1`. Unmatched CRISPR labels are LOC genes whose Entrez id is not an approved HGNC id, such as `LOC730098 (730098)`. `RSC1A1 (6248)` is ambiguous because that Entrez id is not the current id for the gene HGNC now calls `DDI2`.

Cell lines, primary rows, before Cellosaurus (Cell Model Passports crosswalk only) and after Cellosaurus 56.0. Problematic lines stay mapped and are counted in the flag columns.

| Source | Rows | Before mapped / ambiguous / unmatched | After mapped / ambiguous / unmatched | Misidentified | Contaminated | Problematic |
| --- | ---: | --- | --- | ---: | ---: | ---: |
| DepMap | 1959 | 1880 / 7 / 72 | 1884 / 8 / 67 | 61 | 8 | 62 |
| PRISM | 925 | 912 / 2 / 11 | 916 / 4 / 5 | 18 | 5 | 14 |
| Cell Model Passports | 2292 | 1926 / 1 / 365 | 1982 / 5 / 305 | 81 | 13 | 66 |

The extra DepMap mappings are Cellosaurus names or synonyms. Six PRISM rows matched by stripping a `_LUNG`-style CCLE suffix. Fifty-six CMP models with no RRID matched a Cellosaurus name. `H460` is not folded into `NCI-H460` by the normalizer; Cellosaurus lists `H460` as a synonym of `CVCL_0459`, so that synonym match is what resolves it.

Compounds had no CID column. Before resolution every row was unmatched: GDSC screened compounds 621, PRISM 19Q4 4739, PRISM 24Q2 6789. The GDSC² screen also names Galunisertib (`2265`), which is not in `screened_compounds`, so the compound table has 622 GDSC rows after the combination load. `phi-pgx resolve-compounds` ran against PubChem and, for GDSC names still open, exact ChEMBL synonyms.

| Source | Rows | Mapped | Ambiguous | Unmatched | Combination | With ChEMBL | With InChIKey |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GDSC | 622 | 501 | 5 | 114 | 2 | 472 | 501 |
| PRISM 19Q4 | 4739 | 4694 | 39 | 6 | 0 | 4334 | 4612 |
| PRISM 24Q2 | 6789 | 6655 | 119 | 15 | 0 | 6182 | 6598 |

Of the original 621 screened compounds, 500 are mapped, 5 ambiguous, 114 unmatched, and 2 are combinations. The extra mapped row is Galunisertib, PubChem CID 10090485 / `CHEMBL2364611`. GDSC `Nutlin-3a (-)` stayed stored under that name and mapped to PubChem CID 11433190 / `CHEMBL191334` by the stereo-stripped lookup `Nutlin-3a` (`pubchem_name_stereo`, confidence `low`). The second pass mapped `Brivanib, BMS-540215` via the synonym Brivanib (CID 11234052 / `CHEMBL377300`), `eEF2K Inhibitor, A-484954` via A-484954 (CID 14998470), `CAP-232, TT-232, TLN-232` via TT-232 (CID 6918265), and `FTY-720` via the ChEMBL synonym to fingolimod (CID 107970 / `CHEMBL314854`). `Bleomycin (50 uM)` and `Bleomycin (10 uM)` mapped from the dose-stripped name Bleomycin (CID 5360373). `ascorbate (vitamin C)` mapped from the alias ascorbate (CID 54679076).

The 114 unmatched screened compounds are left for review:

| Reason | Rows | What was left |
| --- | ---: | --- |
| internal_code | 100 | Sanger and lab series that were not queried: `KIN001-*`, `THZ-*`, `QL-*`, `TL-*`, `JW-*`, `JQ12`, `CDK9_5576` and other `TARGET_####` codes, `WIMM synthesis` (`THR-101`–`103`), pure numeric ids, `Genentech Cpd 10`, and similar closed codes |
| no_public_structure | 14 | Looked up and not mapped: `rTRAIL`, five `GSK*` salts, `BAY ACCi`, `BAY AKT1`, `BAY-MPS1`, `AZD7969`, `JNJ38877605`, `VNLG/124`, `CT7033-2`, `MetAP2 Inhibitor, A832234` |

Five names stayed ambiguous because more than one public record matched and ChEMBL did not select a single CID: Obatoclax Mesylate (two GDSC ids), CI-1033 (two ChEMBL molecules), ETP-45835 (two CIDs), and Sinularin (three CIDs).

GDSC² combinations are 2027 records: 2025 anchored-screen pairs from the original screen (131 of them also in the validation screen) and 2 annotation combinations. 2025 pairs are mapped, meaning every component has one PubChem CID. The two annotation rows are partial: paclitaxel is CID 36314 / `CHEMBL428647`, and `BAY-MPS` has no public structure. Fifty of the mapped pairs use a two-drug anchor or library, `Afatinib|Trametinib` (`1032|1372`), with both components mapped (Afatinib CID 10184653, Trametinib CID 11707110). Combination response rows: 1,483,535 from the original screen and 117,525 from the validation screen (1,601,060 total), all joined to a Cellosaurus id. They are stored as `response_class=combination` and are not mixed into single-agent AUC or LN_IC50.

PRISM 19Q4 used SMILES for 2341 of its mapped rows (`confidence` high). Two PRISM 24Q2 mebhydroline naphthalenedisulfonate rows mapped to the parent mebhydrolin CID 22530, the same CID as the 19Q4 SMILES row. A lookup that returned more than one CID was left ambiguous.

## Batch 2 drug datasets

NCI-60 and the frozen ORCESTRA PharmacoSets write the same `drug_response` columns as GDSC and PRISM. `dataset` is the source partition (`nci60`, `ctrpv2`, `gcsi`, `fimm`, `gray`, `uhnbreast`, `ccle`). Combinations stay on `drug_combination` and are not created for these screens. `compound_response` with `dataset` omitted returns every partition whose `source_compound_id` shares the queried compound.

`model_xref.source` for these screens is the dataset id. A PharmacoSet Cellosaurus accession is confirmed with `match_method=pset_cellosaurus`. A name that maps to a different accession is left `ambiguous` and is not overwritten. NCI-60 has cell-line names only. Those names go through `CellosaurusIndex.decide_names`. The canonical key strips punctuation (`NCI-H23` matches `NCIH23` when Cellosaurus lists that form). It does not fuzzy-match. Misidentified lines stay mapped and keep `cellosaurus_flag`.

### NCI-60 DTP, October 2026

The wiki page now says “FILES FOR DOWNLOAD October 2026 Release” (page 147193864, updated 7 Oct 2026). The April 2026 filenames are no longer published separately, so this pipeline pins the October versioned attachment URLs. License: US government work. Credit the National Cancer Institute. Not CC BY.

Downloaded on 2026-10-14, not behind `--live`:

| File | Bytes | sha256 |
| --- | ---: | --- |
| GI50.zip | 39212645 | `65212008012e0fc4648709f5f081c18e577986a40a8ba45ea271f8c3bb81f728` |
| TGI.zip | 34639387 | `f260a919f4f8c476f9398ed26aa6db6e8bd5a440896099b5da683f21541e0cc8` |
| LC50.zip | 30684621 | `0493ef82660e2f0f6d7a36ba56fbc212885ccb72193103ca50d21f3b4371beba` |
| IC50.zip | 37261424 | `9e5fd767f0738307b62b1bf83d3f9199d18757b3d23bcca213be65e419ad2ad2` |
| ONECONC.zip | 55646891 | `681291402e668f1b6edd95a0a6b4ef62ac50cfe997066c9382add32bf2d610ed` |
| public_nscs.csv | 2587451 | `34d74a29c94e546c8942e171ca5477bef43a48ca144437c609b19da6c536b125` |

`public_nscs.csv` is 334,952 NSC numbers and has no SMILES, names, or structures. No DTP structure file was pinned. Assay rows with `PREFIX` other than `S` are dropped. `S`-series rows in the downloaded zips: GI50 4,871,117, TGI 4,853,400, LC50 4,838,827, IC50 4,895,862, ONECONC 5,468,176. Endpoint `AVERAGE` is stored as `GI50`, `TGI`, `LC50`, or `IC50` with unit `log10_M` when `CONCENTRATION_UNIT` is `M`. `ONECONC` stores percent growth (`ONECONC_GIPRCNT`) and the prescreen concentration as `dose` (not log-transformed). `DOSERESP` stores `GIPRCNT` and `PTC` with `dose = 10^CONCENTRATION` when the unit is molar. The full matrices were counted from the CSVs. They were not written to `processed/` in this run; the loader is covered by fixture tests.

`DOSERESP.zip` (349,294,730 bytes, version 21) was not downloaded. Its zip method is Deflate64, which the Python zip reader cannot extract. `phi-pgx` calls `7z` when that binary is on `PATH`.

### ORCESTRA PharmacoSets

| Dataset | Pinned object | License | What was downloaded |
| --- | --- | --- | --- |
| ctrpv2 | [10.5281/zenodo.21794049](https://doi.org/10.5281/zenodo.21794049) CTRPv2.1, 2026-07-31 | CC BY 4.0 | `sample.tsv`, `treatment.tsv`, `v2.1_auc_sensitivities.tsv`, and `recomputed_profiles.tsv` (69,384,293 bytes, sha256 `23121c0d08299df0d090381e01eef213b47626da66180143c14952fae5fab4f1`). The RDS and `sensitivity.tsv` were not retrieved. |
| gcsi | [10.5281/zenodo.4737437](https://doi.org/10.5281/zenodo.4737437) `gCSI2.rds`, screen year 2019 | CC BY 4.0 | Not downloaded (about 1.7 GB, RDS only). |
| fimm | [10.5281/zenodo.3905448](https://doi.org/10.5281/zenodo.3905448) 2020-06-24 | CC0-1.0 | Not retrieved in this run. `phi-pgx download` fetches the 231,489-byte RDS because it is not marked `live`. An earlier CC BY 4.0 file is 10.5281/zenodo.3598728 and was not used. |
| gray | [10.5281/zenodo.3905458](https://doi.org/10.5281/zenodo.3905458) deposit 2020-06-24, screen year 2017 | CC0-1.0 | Not downloaded (155 MB RDS). |
| uhnbreast | [10.5281/zenodo.3601415](https://doi.org/10.5281/zenodo.3601415) `UHN_2019.rds` | CC BY 4.0 | Not downloaded (279 MB RDS). The later CC0 object 10.5281/zenodo.7826860 was not used. |
| ccle | [10.5281/zenodo.18759113](https://doi.org/10.5281/zenodo.18759113) CCLE PSet 2019 | CC BY 4.0 | Not downloaded (`CCLE_PharmacoSet.RDS`, about 3 GB). |

CTRPv2.1 already publishes the tables an R extract would emit. Those TSVs are the loaded path. A PharmacoSet CID is `match_method=pset_pubchem`. Two different CIDs on the same drug stay `ambiguous`. Numeric ids are prefixed with the dataset id.

RDS files are converted by `r/extract_pset.R` through `Rscript` (PharmacoGx for an S4 PharmacoSet, or a plain list for tests). This environment does not have R. The image that carries the dependency is `data-pipelines/pharmacogenomics/Dockerfile` (`rocker/r-ver:4.4.2` plus Bioconductor PharmacoGx):

```bash
docker build -t phi-pgx:batch2 data-pipelines/pharmacogenomics
docker run --rm -v "$PWD/raw:/data" --entrypoint Rscript phi-pgx:batch2 \
  /opt/phi-pgx/r/extract_pset.R /data/orcestra_gcsi/2019/gCSI2.rds /data/orcestra_gcsi/2019/extract
```

ORCESTRA lists NCI60_2026, but no public Zenodo file with that name was found. The 2021 deposit 10.5281/zenodo.5570629 is not pinned. A future PSet must use dataset `nci60_pset` so it does not overwrite DTP `nci60`. Do not point this pipeline at an ORCESTRA PRISM PSet.

### Coverage

Counts use Cellosaurus 56.0 and the files downloaded on 2026-10-14. NCI NSC numbers were not sent to PubChem. gCSI, GRAY, UHNBreast, CCLE, and FIMM have no sensitivity rows here because their RDS files were not extracted. CTRP recomputed profiles were downloaded and are loaded only together with `sensitivity.tsv`, which was not retrieved, so this run's CTRP response rows are the v2.1 AUC table.

| Dataset | Cell lines | Mapped | Ambiguous | Unmatched | Compounds | Mapped | Ambiguous | Unmatched |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| nci60 | 159 | 113 | 8 | 38 | 140,045 screened NSCs | not resolved | not resolved | 140,045 until `resolve-compounds` |
| ctrpv2 | 835 | 832 | 2 | 1 | 481 | 467 | 0 | 14 |

`public_nscs.csv` lists 334,952 NSC numbers. The 140,045 figure is the distinct `S`-series NSCs that appear in GI50, TGI, LC50, IC50, or ONECONC. Of the 113 mapped NCI lines, 9 are flagged misidentified and 2 problematic. The 8 ambiguous names are `A-204`, `COLO 205`, `HT29`, `M14`, `NB4`, `ND-1`, `PC-3`, and `TK-10`. The existing canonical key strips punctuation, and each of those keys hits more than one Cellosaurus accession (`COLO 205` is both `CVCL_0218` and `CVCL_F402`; `PC-3` also collides with synonyms that collapse to `pc3`). None was chosen. The 38 unmatched names have no Cellosaurus name or synonym, including engineered derivatives such as `HCT-116/P21/C` and `A-FOS 3`.

CTRPv2 `v2.1_auc_sensitivities.tsv` has 326,148 AUC rows. Its cells are the 835 `master_ccl_id` values in that table. 832 carry one confirmed Cellosaurus accession (18 of those flagged misidentified, 21 problematic, 3 contaminated). `HLF` is ambiguous: the PSet accession is `CVCL_2255` (HLF-a, problematic) and the name matches `CVCL_2947`. `HUH6` is ambiguous: the PSet accession is `CVCL_1296` (HuH-6 Clone 5) and the name matches `CVCL_4381`. `TTTHYROID` has no accession and no Cellosaurus name. 467 compounds have one PSet PubChem CID (`pset_pubchem`). The other 14 have a BRD id and no CID in the PSet, so they stay unmatched until `resolve-compounds`. Two different CIDs on one drug would stay ambiguous; this release had no such row among the screened compounds.

### Later

HMDB bulk XML remains planned and is not downloaded.

## Limitations

- NCI-60 April 2026 is not a separate download anymore. The pin is October 2026.
- No DTP SMILES or structure file is in the pinned set. NSC-to-CID mapping is a PubChem synonym lookup, and a missing or multi-CID synonym stays unresolved.
- `DOSERESP.zip` was not downloaded. It needs 7-Zip because the archive is Deflate64.
- gCSI, GRAY, UHNBreast, CCLE, and FIMM were not extracted. CTRP `sensitivity.tsv` was not retrieved, so recomputed IC50/AUC profiles (which need that file for the cell id) were not loaded. The AUC table does not need it.
- Recomputed CTRP profiles have a `response_id` and no cell id. They are loaded only when `sensitivity.tsv` is present to supply `master_ccl_id`.
- PharmacoGx raw sensitivity stored as an array, rather than a data.frame, is noted in `extract_notes.txt` and is not expanded.
- Cancelled NCI experiments are not in these zips. Rows whose `PREFIX` is not `S` are ignored anyway.
- These screens are research assays. They are not treatment recommendations.
