import assert from "node:assert/strict";
import { mkdirSync } from "node:fs";
import { mkdtemp, rm } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import test from "node:test";
import { DuckDBInstance } from "@duckdb/node-api";
import { ToolResultSchema, type Page, type ToolResult } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

const client: JsonClient = { async getJson() { return {}; }, async postJson() { return {}; } };

// Twelve models: 1-6 are Lung, 7-12 are Breast. TP53 expression rises 1..12. Models 1-4 carry a
// missense TP53 mutation, model 5 only a synonymous one, and models 5-12 are sequenced wild type.
const RESPONSE = [8, 7, 9, 6, 3, 4, 2, 1, 3.5, 0.5, 2.5, 1.5];

const cvcl = (index: number) => `CVCL_A${String(index).padStart(2, "0")}`;
const ach = (index: number) => `ACH-${String(index).padStart(6, "0")}`;
const range = (count: number) => Array.from({ length: count }, (_, index) => index + 1);

async function writeFixture(dir: string): Promise<void> {
  const connection = await (await DuckDBInstance.create(":memory:")).connect();
  const copy = async (sql: string, file: string) => {
    mkdirSync(path.dirname(file), { recursive: true });
    await connection.run(`COPY (${sql}) TO '${file}' (FORMAT PARQUET)`);
  };
  const values = (rowsList: string[]) => rowsList.join(",\n");
  await copy(
    `SELECT * FROM (VALUES ${values(range(12).map((i) => `('depmap','${ach(i)}','${ach(i)}','primary','${cvcl(i)}','depmap_rrid','mapped','depmap_rrid','Model ${i}','${i <= 6 ? "Lung" : "Breast"}','Carcinoma','Homo sapiens')`))})
     AS t(source, source_model_id, alias, alias_type, cvcl, match_method, status, evidence, display_name, tissue, disease, species)`,
    path.join(dir, "model_xref", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES ('gdsc','1047','Nutlin-3a (-)','','','1047','11433190','CHEMBL191334','BDUHCSBCVGXTJM-WUFINQPMSA-N','','source_pubchem','mapped','PubChem CID column'))
     AS t(source, source_compound_id, name, synonyms, brd_id, gdsc_drug_id, pubchem_cid, chembl_id, inchikey, smiles, match_method, status, evidence)`,
    path.join(dir, "compound_xref", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES ${values(range(12).map((i) => `('${ach(i)}','${cvcl(i)}','7157','TP53',${i},'log2_tpm_plus_1')`))})
     AS t(source_model_id, model_cvcl, ncbi_gene_id, symbol, value, metric)`,
    path.join(dir, "expression", "dataset=depmap_expression_24q2", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES ${values(range(12).map((i) => `('${ach(i)}','${cvcl(i)}','7157','TP53',${i % 3 === 0 ? 0.5 : 1.0},'relative_copy_number')`))})
     AS t(source_model_id, model_cvcl, ncbi_gene_id, symbol, value, metric)`,
    path.join(dir, "cnv", "dataset=depmap_cn_24q2", "part.parquet"),
  );
  const mutationRows = [
    ...range(4).map((i) => `('${ach(i)}','${cvcl(i)}','7157','TP53','chr17',${7676000 + i},'G','A','p.R175H','missense_variant','GRCh38','MODERATE')`),
    `('${ach(5)}','${cvcl(5)}','7157','TP53','chr17',7676200,'C','T','p.Syn','synonymous_variant','GRCh38','LOW')`,
    ...range(12).filter((i) => i >= 5).map((i) => `('${ach(i)}','${cvcl(i)}','9999','OTHER','chr1',${1000 + i},'A','C','p.K1E','missense_variant','GRCh38','MODERATE')`),
  ];
  await copy(
    `SELECT * FROM (VALUES ${values(mutationRows)})
     AS t(source_model_id, model_cvcl, ncbi_gene_id, symbol, chrom_norm, pos, ref, alt, protein_change, consequence, genome_build, vep_impact)`,
    path.join(dir, "mutation", "dataset=depmap_mutations_24q2", "part.parquet"),
  );
  await copy(
    `SELECT * FROM (VALUES ${values(range(12).map((i) => `('${ach(i)}','${cvcl(i)}','1047','11433190','Nutlin-3a (-)','LN_IC50',${RESPONSE[i - 1]},'ln_uM',NULL::DOUBLE,NULL::VARCHAR,NULL::VARCHAR,NULL::VARCHAR,'single',NULL::DOUBLE,NULL::DOUBLE)`))})
     AS t(source_model_id, model_cvcl, source_compound_id, pubchem_cid, compound_name, metric, value, unit, dose, screen_id, qc, combination_id, response_class, anchor_conc, library_conc)`,
    path.join(dir, "drug_response", "dataset=gdsc1", "part.parquet"),
  );
}

function tool(name: string) {
  const found = tools.find((item) => item.name === name);
  assert.ok(found, `tool ${name} is registered`);
  return found;
}

async function run(name: string, args: Record<string, unknown>): Promise<ToolResult> {
  return ToolResultSchema.parse(await tool(name).call(client, args));
}

function page(result: ToolResult): Page {
  assert.equal(result.entity_type, "page", JSON.stringify(result));
  if (result.entity_type !== "page") throw new Error("expected a page");
  return result;
}

async function withFixture(body: () => Promise<void>): Promise<void> {
  const dir = await mkdtemp(path.join(os.tmpdir(), "phi-pgx-omics-"));
  process.env.PHI_PGX_DATA_DIR = dir;
  try {
    await writeFixture(dir);
    await body();
  } finally {
    delete process.env.PHI_PGX_DATA_DIR;
    await rm(dir, { recursive: true, force: true });
  }
}

test("list_tissues counts mapped models per tissue", async () => {
  await withFixture(async () => {
    const result = page(await run("list_tissues", {}));
    assert.equal(result.page.total, 2);
    assert.deepEqual(result.items.map((item) => item.entity_type === "tissue" ? [item.name, item.source_extra.models] : null), [["Breast", 6], ["Lung", 6]]);
  });
});

test("gene_expression orders by value, filters by tissue, and summarizes the distribution", async () => {
  await withFixture(async () => {
    const all = page(await run("gene_expression", { gene: "TP53", limit: 3 }));
    assert.equal(all.page.total, 12);
    const first = all.items[0];
    assert.equal(first?.entity_type === "expression" ? first.value : null, 12);
    assert.equal(first?.entity_type === "expression" ? first.sample_id : null, "cellosaurus:CVCL_A12");
    assert.equal(first?.source_extra.tissue, "Breast");
    assert.equal((all.source_extra.distribution as Record<string, number>).median, 6.5);
    assert.equal((all.source_extra.distribution as Record<string, number>).n, 12);

    const lung = page(await run("gene_expression", { gene: "ncbigene:7157", tissue: "lung", order: "asc", limit: 2 }));
    assert.equal(lung.page.total, 6);
    assert.deepEqual(lung.items.map((item) => item.entity_type === "expression" ? item.value : null), [1, 2]);
    assert.equal((lung.source_extra.distribution as Record<string, number>).max, 6);
  });
});

test("a tissue with no models is a not_found error that points at list_tissues", async () => {
  await withFixture(async () => {
    const result = await run("gene_expression", { gene: "TP53", tissue: "Spleen" });
    assert.equal(result.entity_type, "error");
    if (result.entity_type === "error") {
      assert.equal(result.code, "not_found");
      assert.match(result.message, /list_tissues/);
    }
    const dependency = await run("gene_dependency", { gene: "TP53", tissue: "Spleen" });
    assert.equal(dependency.entity_type, "error");
  });
});

test("gene_copy_number returns sample rows with the value in source_extra", async () => {
  await withFixture(async () => {
    const result = page(await run("gene_copy_number", { gene: "TP53", tissue: "Breast", order: "asc", limit: 10 }));
    assert.equal(result.page.total, 6);
    assert.equal(result.items.every((item) => item.entity_type === "sample"), true);
    assert.equal(result.items[0]?.source_extra.copy_number, 0.5);
    assert.equal(result.items[0]?.source_extra.metric, "relative_copy_number");
  });
});

test("gene_mutations hides silent variants by default and reports frequency and hotspots", async () => {
  await withFixture(async () => {
    const result = page(await run("gene_mutations", { gene: "TP53" }));
    assert.equal(result.page.total, 4, "the synonymous call is left out");
    assert.equal(result.source_extra.mutated_models, 4);
    assert.equal(result.source_extra.sequenced_models, 12);
    assert.equal(result.source_extra.frequency, 4 / 12);
    assert.deepEqual(result.source_extra.top_protein_changes, [{ protein_change: "p.R175H", models: 4 }]);
    const first = result.items[0];
    assert.equal(first?.entity_type === "variant" ? first.consequence : null, "missense_variant");
    assert.equal(first?.entity_type === "variant" ? first.locus?.chrom : null, "chr17");

    const withSilent = page(await run("gene_mutations", { gene: "TP53", include_silent: true }));
    assert.equal(withSilent.page.total, 5);

    const byChange = page(await run("gene_mutations", { gene: "TP53", protein_change: "R175H", tissue: "Lung" }));
    assert.equal(byChange.page.total, 4);
    assert.equal(byChange.source_extra.sequenced_models, 6);
    assert.equal(byChange.source_extra.frequency, 4 / 6);

    const none = page(await run("gene_mutations", { gene: "TP53", consequence: "stop_gained" }));
    assert.equal(none.page.total, 0);

    const requestedSilent = page(await run("gene_mutations", { gene: "TP53", consequence: "synonymous_variant" }));
    assert.equal(requestedSilent.page.total, 1, "naming a silent consequence returns it without include_silent");
  });
});

test("biomarker_association tests mutation status with Mann-Whitney", async () => {
  await withFixture(async () => {
    const result = page(await run("biomarker_association", { compound: "1047", gene: "TP53", feature: "mutation" }));
    assert.equal(result.page.total, 1);
    const stats = result.items[0]?.source_extra as Record<string, number | string>;
    assert.equal(stats.drug_dataset, "gdsc1");
    assert.equal(stats.metric, "LN_IC50");
    assert.equal(stats.n_mutant, 4);
    assert.equal(stats.n_wild_type, 8, "the synonymous-only model counts as wild type");
    assert.equal(stats.median_mutant, 7.5);
    assert.equal(stats.median_wild_type, 2.25);
    assert.equal(stats.rank_biserial, 1);
    assert.ok(Math.abs(Number(stats.p_value) - 0.008474801892153832) < 1e-6, `p ${stats.p_value}`);
    assert.equal(stats.direction, "lower_is_stronger");
    assert.equal(result.source_extra.skipped_group_count, 0);
  });
});

test("biomarker_association correlates response with expression using Spearman", async () => {
  await withFixture(async () => {
    const result = page(await run("biomarker_association", { compound: "NUTLIN-3A (-)", gene: "ncbigene:7157", feature: "expression" }));
    assert.equal(result.page.total, 1);
    const stats = result.items[0]?.source_extra as Record<string, number>;
    assert.equal(stats.n_models, 12);
    assert.ok(Math.abs(stats.rho - -0.8251748251748254) < 1e-9, `rho ${stats.rho}`);
    assert.ok(Math.abs(stats.p_value - 0.0009513628555555001) < 1e-6, `p ${stats.p_value}`);
  });
});

test("biomarker_association skips groups that are too small and says why", async () => {
  await withFixture(async () => {
    const mutation = page(await run("biomarker_association", { compound: "1047", gene: "TP53", feature: "mutation", tissue: "Lung" }));
    assert.equal(mutation.page.total, 0);
    assert.equal(mutation.source_extra.skipped_group_count, 1);
    const skipped = mutation.source_extra.skipped_groups as { reason: string }[];
    assert.match(skipped[0]!.reason, /at least 3 models per group, found 4 mutant and 2 wild type/);

    const numeric = page(await run("biomarker_association", { compound: "1047", gene: "TP53", feature: "expression", tissue: "Lung" }));
    assert.equal(numeric.page.total, 0);
    assert.match((numeric.source_extra.skipped_groups as { reason: string }[])[0]!.reason, /at least 10 models/);
  });
});

test("biomarker_association rejects bad input and reports unknown compounds and missing data", async () => {
  await withFixture(async () => {
    await assert.rejects(tool("biomarker_association").call(client, { compound: "1047", gene: "TP53", feature: "methylation" }), /feature must be one of/);
    await assert.rejects(tool("biomarker_association").call(client, { compound: "1047", gene: "TP53", feature: "mutation", min_group: 1 }), /min_group/);
    const unknown = await run("biomarker_association", { compound: "no such drug", gene: "TP53", feature: "mutation" });
    assert.equal(unknown.entity_type, "error");
    if (unknown.entity_type === "error") assert.equal(unknown.code, "not_found");
  });
  delete process.env.PHI_PGX_DATA_DIR;
  const missing = await run("gene_mutations", { gene: "TP53" });
  assert.equal(missing.entity_type, "error");
  if (missing.entity_type === "error") assert.equal(missing.code, "unavailable");
});

test("compound_response accepts a tissue filter", async () => {
  await withFixture(async () => {
    assert.equal(page(await run("compound_response", { compound: "1047", limit: 10 })).page.total, 12);
    const response = page(await run("compound_response", { compound: "1047", tissue: "Breast", limit: 10 }));
    assert.equal(response.page.total, 6);
    const lung = page(await run("compound_response", { compound: "1047", tissue: "Lung", limit: 10 }));
    assert.equal(lung.page.total, 6);
    assert.equal(response.items.every((item) => item.entity_type === "drug_response"), true);
  });
});
