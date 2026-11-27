import assert from "node:assert/strict";
import { mkdtemp, rm } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import test from "node:test";
import { ToolResultSchema, type ToolResult } from "phi-omics-schemas";
import { UpstreamClient, UpstreamError, loadConfig } from "phi-omics-runtime";
import { tools } from "../server.js";

const dir = await mkdtemp(path.join(os.tmpdir(), "phi-structures-live-"));
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

test("live UniProtKB TP53", async (t) => {
  await live(t, "uniprot_search_proteins", { query: "TP53", reviewed_only: true, organism_id: 9606, limit: 2 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.ok(result.items.length > 0);
  });
  await live(t, "uniprot_get_protein", { accession: "P04637" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "uniprot:P04637");
  });
});

test("live AlphaFold DB P04637", async (t) => {
  await live(t, "alphafold_get_prediction", { uniprot_accession: "P04637" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.ok(result.items[0]?.id.startsWith("alphafold:AF-P04637"));
  });
  await live(t, "alphafold_list_files", { uniprot_accession: "P04637", limit: 5 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.item_type, "file");
  });
  await live(t, "alphafold_get_summary", { uniprot_accession: "P04637" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "uniprot:P04637");
  });
});

test.after(async () => rm(dir, { recursive: true, force: true }));
