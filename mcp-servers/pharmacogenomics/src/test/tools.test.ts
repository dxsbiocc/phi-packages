import assert from "node:assert/strict";
import { mkdirSync } from "node:fs";
import { mkdtemp, rm } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import test from "node:test";
import { DuckDBInstance } from "@duckdb/node-api";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

const client: JsonClient = { async getJson() { return {}; }, async postJson() { return {}; } };

async function writeFixture(dir: string): Promise<void> {
  const connection = await (await DuckDBInstance.create(":memory:")).connect();
  const copy = async (sql: string, file: string) => {
    mkdirSync(path.dirname(file), { recursive: true });
    await connection.run(`COPY (${sql}) TO '${file}' (FORMAT PARQUET)`);
  };
  await copy(
    `SELECT * FROM (VALUES
      ('depmap','ACH-000001','ACH-000001','primary','CVCL_0001','depmap_rrid','mapped','depmap_rrid','Alpha','Lung','Carcinoma','Homo sapiens'),
      ('depmap','ACH-000001','CVCL_0001','rrid','CVCL_0001','depmap_rrid','mapped','depmap_rrid','Alpha','Lung','Carcinoma','Homo sapiens'),
      ('depmap','ACH-000001','Alpha','cell_line_name','CVCL_0001','depmap_rrid','mapped','depmap_rrid','Alpha','Lung','Carcinoma','Homo sapiens'),
      ('depmap','ACH-000002','ACH-000002','primary','CVCL_0002','depmap_rrid','mapped','depmap_rrid','Beta','Lung','Carcinoma','Homo sapiens')
    ) AS t(source, source_model_id, alias, alias_type, cvcl, match_method, status, evidence, display_name, tissue, disease, species)`,
    path.join(dir, "model_xref", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES
      ('gdsc','1047','Nutlin-3a (-)','','', '1047','11433190','CHEMBL191334','BDUHCSBCVGXTJM-WUFINQPMSA-N','', 'source_pubchem','mapped','PubChem CID column'),
      ('gdsc','1012','Vorinostat','','', '1012','5311','CHEMBL98','WAEXFXRVDQXREF-UHFFFAOYSA-N','', 'source_pubchem','mapped','PubChem CID column'),
      ('nci60','NSC740','NSC740','','', NULL,'11433190',NULL,NULL,NULL, 'pubchem_nsc','mapped','NSC synonym'),
      ('prism_19q4','BRD-A12230535-001-06-7','nutlin-3','','BRD-A12230535-001-06-7',NULL,NULL,NULL,NULL,NULL,'unmatched','unmatched','no PubChem CID')
    ) AS t(source, source_compound_id, name, synonyms, brd_id, gdsc_drug_id, pubchem_cid, chembl_id, inchikey, smiles, match_method, status, evidence)`,
    path.join(dir, "compound_xref", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES
      ('gdsc.combination:1012:1003','gdsc_combinations','gdsc_combination_validation','Vorinostat + Camptothecin','1012','1003','Vorinostat','Camptothecin','5311','24360','CHEMBL98','CHEMBL65','2','0.1','uM','mapped',NULL)
    ) AS t(combination_id, source, screen, name, anchor_drug_id, library_drug_id, anchor_name, library_name, anchor_pubchem_cid, library_pubchem_cid, anchor_chembl_id, library_chembl_id, anchor_concentrations, library_concentrations, concentration_unit, status, source_drug_id)`,
    path.join(dir, "combination_xref", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES
      ('depmap','24Q2','CC BY 4.0','DepMap 24Q2 Public','https://depmap.org/portal/'),
      ('gdsc','8.4','Sanger GDSC release 8.4','Yang et al. NAR 2013','https://www.cancerrxgene.org/')
    ) AS t(source, version, license, citation, homepage)`,
    path.join(dir, "releases", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES
      ('ACH-000001','CVCL_0001','7157','TP53',-1.2,'chronos_gene_effect'),
      ('ACH-000002','CVCL_0002','7157','TP53',-0.2,'chronos_gene_effect')
    ) AS t(source_model_id, model_cvcl, ncbi_gene_id, symbol, value, metric)`,
    path.join(dir, "gene_effect", "dataset=depmap_crispr_24q2", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES
      ('ACH-000001','CVCL_0001','7157','TP53',4.5,'log2_tpm_plus_1')
    ) AS t(source_model_id, model_cvcl, ncbi_gene_id, symbol, value, metric)`,
    path.join(dir, "expression", "dataset=depmap_expression_24q2", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES
      ('ACH-000001','CVCL_0001','7157','TP53',1.1,'relative_copy_number')
    ) AS t(source_model_id, model_cvcl, ncbi_gene_id, symbol, value, metric)`,
    path.join(dir, "cnv", "dataset=depmap_cn_24q2", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES
      ('ACH-000001','CVCL_0001','7157','TP53','chr17',7676154,'G','A','p.Syn','synonymous_variant','GRCh38','LOW')
    ) AS t(source_model_id, model_cvcl, ncbi_gene_id, symbol, chrom_norm, pos, ref, alt, protein_change, consequence, genome_build, vep_impact)`,
    path.join(dir, "mutation", "dataset=depmap_mutations_24q2", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES
      ('ACH-000001','CVCL_0001','1047','11433190','Nutlin-3a (-)', 'LN_IC50', 0.2, 'ln_uM', NULL::DOUBLE, NULL::VARCHAR, NULL::VARCHAR, NULL::VARCHAR, 'single', NULL::DOUBLE, NULL::DOUBLE),
      ('ACH-000002','CVCL_0002','1047','11433190','Nutlin-3a (-)', 'LN_IC50', 2.5, 'ln_uM', NULL::DOUBLE, NULL::VARCHAR, NULL::VARCHAR, NULL::VARCHAR, 'single', NULL::DOUBLE, NULL::DOUBLE),
      ('ACH-000001','CVCL_0001','1012','5311','Vorinostat', 'LN_IC50', 1.0, 'ln_uM', NULL::DOUBLE, NULL::VARCHAR, NULL::VARCHAR, NULL::VARCHAR, 'single', NULL::DOUBLE, NULL::DOUBLE),
      ('ACH-000001','CVCL_0001','1012:1003',NULL,'Vorinostat + Camptothecin', 'SYNERGY_DELTA_AUC', 0.05, NULL, NULL::DOUBLE, NULL::VARCHAR, NULL::VARCHAR, 'gdsc.combination:1012:1003', 'combination', 2.0, 0.1)
    ) AS t(source_model_id, model_cvcl, source_compound_id, pubchem_cid, compound_name, metric, value, unit, dose, screen_id, qc, combination_id, response_class, anchor_conc, library_conc)`,
    path.join(dir, "drug_response", "dataset=gdsc1", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES
      ('NCI-H23','CVCL_1547','NSC740','11433190','NSC740', 'GI50', -7.1, 'log10_M', NULL::DOUBLE, '0001MD02', NULL::VARCHAR, NULL::VARCHAR, 'single', NULL::DOUBLE, NULL::DOUBLE)
    ) AS t(source_model_id, model_cvcl, source_compound_id, pubchem_cid, compound_name, metric, value, unit, dose, screen_id, qc, combination_id, response_class, anchor_conc, library_conc)`,
    path.join(dir, "drug_response", "dataset=nci60", "part.parquet"),
  );
}

test("pharmacogenomics tools read fixture parquet", async () => {
  const dir = await mkdtemp(path.join(os.tmpdir(), "phi-pgx-"));
  process.env.PHI_PGX_DATA_DIR = dir;
  try {
    await writeFixture(dir);
    assert.equal(tools.length, 5);
    const calls: Record<string, Record<string, unknown>> = {
      usage_terms: {},
      resolve_cell_line: { query: "Alpha" },
      gene_dependency: { gene: "TP53", limit: 5 },
      compound_response: { compound: "pubchem:11433190", metric: "LN_IC50", limit: 5 },
      model_summary: { model: "ACH-000001", gene: "ncbigene:7157" },
    };
    for (const tool of tools) {
      const result = ToolResultSchema.parse(await tool.call(client, calls[tool.name] ?? {}));
      assert.equal(result.entity_type, "page", tool.name);
    }
    const resolved = ToolResultSchema.parse(await tools[1]!.call(client, { query: "Alpha" }));
    assert.equal(resolved.entity_type, "page");
    if (resolved.entity_type === "page") {
      assert.equal(resolved.items[0]?.id, "cellosaurus:CVCL_0001");
    }
    const dependency = ToolResultSchema.parse(await tools[2]!.call(client, { gene: "7157", limit: 1 }));
    assert.equal(dependency.entity_type, "page");
    if (dependency.entity_type === "page") {
      assert.equal(dependency.items[0]?.id, "cellosaurus:CVCL_0001");
      assert.equal(dependency.page.total, 2);
    }
    const response = ToolResultSchema.parse(await tools[3]!.call(client, { compound: "1047", limit: 5 }));
    assert.equal(response.entity_type, "page");
    if (response.entity_type === "page" && response.items[0]?.entity_type === "drug_response") {
      assert.equal(response.items[0].compound_id, "pubchem:11433190");
      assert.equal(response.items[0].combination_id, null);
      assert.equal(response.items[0].value, 0.2);
      assert.equal(response.items.every((item) => item.entity_type === "drug_response" && item.source_extra.response_class === "single"), true);
    }
    const singleAgent = ToolResultSchema.parse(await tools[3]!.call(client, { compound: "1012", limit: 5 }));
    assert.equal(singleAgent.entity_type, "page");
    if (singleAgent.entity_type === "page") {
      assert.equal(singleAgent.page.total, 1);
      assert.equal(singleAgent.items[0]?.entity_type === "drug_response" ? singleAgent.items[0].metric : null, "LN_IC50");
    }
    const included = ToolResultSchema.parse(await tools[3]!.call(client, { compound: "Vorinostat", combinations: "include", limit: 5 }));
    assert.equal(included.entity_type, "page");
    if (included.entity_type === "page") {
      assert.equal(included.page.total, 2);
      assert.equal(included.items[0]?.entity_type === "drug_response" ? included.items[0].combination_id : null, "gdsc.combination:1012:1003");
      assert.equal(included.items[0]?.entity_type === "drug_response" ? included.items[0].compound_id : "set", null);
    }
    const across = ToolResultSchema.parse(await tools[3]!.call(client, { compound: "pubchem:11433190", limit: 10 }));
    assert.equal(across.entity_type, "page");
    if (across.entity_type === "page") {
      const datasets = across.items.map((item) => item.entity_type === "drug_response" ? item.dataset_id : null);
      assert.equal(datasets.includes("pgx.dataset:gdsc1"), true);
      assert.equal(datasets.includes("pgx.dataset:nci60"), true);
      assert.equal(across.page.total, 3);
    }
    const nsc = ToolResultSchema.parse(await tools[3]!.call(client, { compound: "NSC740", limit: 5 }));
    assert.equal(nsc.entity_type, "page");
    if (nsc.entity_type === "page" && nsc.items[0]?.entity_type === "drug_response") {
      assert.equal(nsc.page.total, 1);
      assert.equal(nsc.items[0].dataset_id, "pgx.dataset:nci60");
      assert.equal(nsc.items[0].compound_id, "pubchem:11433190");
    }
    const only = ToolResultSchema.parse(await tools[3]!.call(client, { compound: "Vorinostat + Camptothecin", limit: 5 }));
    assert.equal(only.entity_type, "page");
    if (only.entity_type === "page" && only.items[0]?.entity_type === "drug_response") {
      assert.equal(only.page.total, 1);
      assert.equal(only.items[0].combination_id, "gdsc.combination:1012:1003");
      assert.equal(only.items[0].compound_id, null);
      assert.equal(only.items[0].metric, "SYNERGY_DELTA_AUC");
      assert.equal(only.items[0].source_extra.response_class, "combination");
    }
  } finally {
    delete process.env.PHI_PGX_DATA_DIR;
    await rm(dir, { recursive: true, force: true });
  }
});

test("missing pharmacogenomics data is an unavailable error", async () => {
  delete process.env.PHI_PGX_DATA_DIR;
  const result = ToolResultSchema.parse(await tools[1]!.call(client, { query: "Alpha" }));
  assert.equal(result.entity_type, "error");
  if (result.entity_type === "error") assert.equal(result.code, "unavailable");
});
