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

async function live(t: test.TestContext, name: string, args: Record<string, unknown>, check: (result: ToolResult) => void) {
  if (process.env.PHI_MCP_SKIP_LIVE === "1") {
    t.skip("PHI_MCP_SKIP_LIVE=1");
    return;
  }
  const found = tools.find((item) => item.name === name);
  if (!found) throw new Error(name);
  try {
    const result = ToolResultSchema.parse(await found.call(client, args));
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

test("live PRIDE project and file", async (t) => {
  await live(t, "pride_get_project", { accession: "PXD001357" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "pride:PXD001357");
  });
  await live(t, "pride_list_files", { accession: "PXD001357", limit: 1 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.item_type, "file");
  });
});

test("live Human Protein Atlas TP53", async (t) => {
  await live(t, "hpa_get_gene", { gene: "ENSG00000141510" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "ensembl:ENSG00000141510");
  });
  await live(t, "hpa_list_proteins", { gene: "ENSG00000141510", limit: 5 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "uniprot:P04637");
  });
});

test("live remaining proteomics tools", async (t) => {
  const calls: Array<[string, Record<string, unknown>, string]> = [
    ["pride_search_projects", { limit: 1 }, "pride:PXD001357"],
    ["hpa_search_genes", { query: "TP53", limit: 1 }, "ensembl:ENSG00000141510"],
  ];
  for (const [name, args, idPrefix] of calls) {
    await live(t, name, args, (result) => {
      assert.equal(result.entity_type, "page", name);
      if (result.entity_type === "page") assert.ok(result.items[0]?.id.startsWith(idPrefix), result.items[0]?.id);
    });
  }
});

test.after(async () => rm(dir, { recursive: true, force: true }));
