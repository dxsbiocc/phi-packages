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

async function writeFixture(dir: string, standardized: boolean): Promise<void> {
  const connection = await (await DuckDBInstance.create(":memory:")).connect();
  const copy = async (sql: string, file: string) => {
    mkdirSync(path.dirname(file), { recursive: true });
    await connection.run(`COPY (${sql}) TO '${file}' (FORMAT PARQUET)`);
  };
  await copy(
    `SELECT * FROM (VALUES
      ('gdsc','1047','Nutlin-3a (-)','','', '1047','11433190','CHEMBL191334','BDUHCSBCVGXTJM-WUFINQPMSA-N','', 'source_pubchem','mapped','PubChem CID column')
    ) AS t(source, source_compound_id, name, synonyms, brd_id, gdsc_drug_id, pubchem_cid, chembl_id, inchikey, smiles, match_method, status, evidence)`,
    path.join(dir, "compound_xref", "part.parquet"),
  );
  const extra = standardized ? ", 'lower_is_stronger', 1.0, 0.75, TRUE" : "";
  const extraNames = standardized ? ", response_direction, sensitivity_percentile, compound_percentile, beyond_tested_range" : "";
  await copy(
    `SELECT * FROM (VALUES
      ('ACH-000001','CVCL_0001','1047','11433190','Nutlin-3a (-)', 'LN_IC50', 0.2, 'ln_uM', NULL::DOUBLE, NULL::VARCHAR, NULL::VARCHAR, NULL::VARCHAR, 'single', NULL::DOUBLE, NULL::DOUBLE${extra})
    ) AS t(source_model_id, model_cvcl, source_compound_id, pubchem_cid, compound_name, metric, value, unit, dose, screen_id, qc, combination_id, response_class, anchor_conc, library_conc${extraNames})`,
    path.join(dir, "drug_response", "dataset=gdsc1", "part.parquet"),
  );
}

async function firstRow(standardized: boolean): Promise<Record<string, unknown>> {
  const dir = await mkdtemp(path.join(os.tmpdir(), "phi-pgx-std-"));
  process.env.PHI_PGX_DATA_DIR = dir;
  try {
    await writeFixture(dir, standardized);
    const tool = tools.find((item) => item.name === "compound_response")!;
    const result = ToolResultSchema.parse(await tool.call(client, { compound: "1047" }));
    assert.equal(result.entity_type, "page");
    if (result.entity_type !== "page") throw new Error("expected a page");
    return result.items[0]!.source_extra as Record<string, unknown>;
  } finally {
    delete process.env.PHI_PGX_DATA_DIR;
    await rm(dir, { recursive: true, force: true });
  }
}

test("compound_response returns the standardized percentiles when the pipeline wrote them", async () => {
  const extra = await firstRow(true);
  assert.equal(extra.sensitivity_percentile, 1);
  assert.equal(extra.compound_percentile, 0.75);
  assert.equal(extra.direction, "lower_is_stronger");
  assert.equal(extra.beyond_tested_range, true);
});

test("compound_response still works on parquet written before the percentile columns existed", async () => {
  const extra = await firstRow(false);
  assert.equal(extra.sensitivity_percentile, null);
  assert.equal(extra.compound_percentile, null);
  assert.equal(extra.direction, "lower_is_stronger");
});
