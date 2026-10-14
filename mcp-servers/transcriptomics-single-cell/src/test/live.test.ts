import assert from "node:assert/strict";
import { mkdtemp, rm } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import test from "node:test";
import { ToolResultSchema, type ToolResult } from "phi-omics-schemas";
import { UpstreamClient, UpstreamError, loadConfig } from "phi-omics-runtime";
import { tools } from "../server.js";

const dir = await mkdtemp(path.join(os.tmpdir(), "phi-omics-live-"));
process.env.PHI_OMICS_CACHE_DIR = dir;
process.env.PHI_OMICS_RATE_LIMIT_RPS = "2";
process.env.PHI_OMICS_TIMEOUT_MS = "30000";
const client = new UpstreamClient(loadConfig());

function tool(name: string) {
  const found = tools.find((item) => item.name === name);
  if (!found) throw new Error(name);
  return found;
}

async function live(t: test.TestContext, name: string, args: Record<string, unknown>, check: (result: ToolResult) => void) {
  if (process.env.PHI_MCP_SKIP_LIVE === "1") {
    t.skip("PHI_MCP_SKIP_LIVE=1");
    return;
  }
  try {
    const result = ToolResultSchema.parse(await tool(name).call(client, args));
    if (result.entity_type === "error" && (result.code === "unavailable" || result.retryable)) {
      t.skip(result.message);
      return;
    }
    check(result);
  } catch (error) {
    if (error instanceof UpstreamError && error.status == null) {
      t.skip(error.message);
      return;
    }
    throw error;
  }
}

test("live GTEx TP53 gene", async (t) => {
  await live(t, "gtex_search_genes", { gene: "TP53", limit: 1 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") {
      assert.equal(result.item_type, "gene");
      assert.equal(result.items[0]?.id, "ensembl:ENSG00000141510");
    }
  });
});

test("live CELLxGENE collection", async (t) => {
  await live(t, "cellxgene_get_collection", { collection_id: "af893e86-8e9f-41f1-a474-ef05359b1fb7" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "cellxgene.collection:af893e86-8e9f-41f1-a474-ef05359b1fb7");
  });
});

test("live HCA project", async (t) => {
  await live(t, "hca_get_project", { project_id: "74b6d569-3b11-42ef-b6b1-a0454522b4a0" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "hca.project:74b6d569-3b11-42ef-b6b1-a0454522b4a0");
  });
});

test("live HuBMAP dataset search", async (t) => {
  await live(t, "hubmap_search_datasets", { limit: 1 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.item_type, "dataset");
  });
});

test("live CellGuide markers", async (t) => {
  await live(t, "cellguide_cell_type", { cell_type_id: "CL:0000034" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "CL:0000034");
  });
});

test("live remaining transcriptomics tools", async (t) => {
  const calls: Array<[string, Record<string, unknown>]> = [
    ["gtex_list_tissues", { limit: 1 }],
    ["gtex_gene_expression", { gene: "TP53", limit: 1 }],
    ["gtex_list_egenes", { tissue: "Whole_Blood", limit: 1 }],
    ["gtex_list_eqtls", { gene: "ENSG00000227232.5", tissue: "Whole_Blood", limit: 1 }],
    ["cellxgene_list_collections", { limit: 1 }],
    ["cellxgene_list_datasets", { collection_id: "af893e86-8e9f-41f1-a474-ef05359b1fb7", limit: 1 }],
    ["cellxgene_get_dataset", { collection_id: "af893e86-8e9f-41f1-a474-ef05359b1fb7", dataset_id: "ed419b4e-db9b-40f1-8593-68fdf8dfb076" }],
    ["hca_search_projects", { limit: 1 }],
    ["hubmap_search_samples", { limit: 1 }],
    ["hubmap_search_donors", { limit: 1 }],
  ];
  for (const [name, args] of calls) {
    await live(t, name, args, (result) => {
      assert.equal(result.entity_type, "page", name);
      if (result.entity_type === "page") assert.ok(result.items.length > 0, name);
    });
  }
});

test.after(async () => {
  await rm(dir, { recursive: true, force: true });
});
