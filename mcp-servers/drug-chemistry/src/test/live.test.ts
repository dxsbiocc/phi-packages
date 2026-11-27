import assert from "node:assert/strict";
import { mkdtemp, rm } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import test from "node:test";
import { ToolResultSchema, type ToolResult } from "phi-omics-schemas";
import { UpstreamClient, UpstreamError, loadConfig } from "phi-omics-runtime";
import { tools } from "../server.js";

const dir = await mkdtemp(path.join(os.tmpdir(), "phi-chemistry-live-"));
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

test("live ChEMBL CHEMBL941 imatinib tools", async (t) => {
  await live(t, "chembl_search_molecules", { query: "imatinib", limit: 2 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.ok(result.items.length > 0);
  });
  await live(t, "chembl_get_molecule", { chembl_id: "CHEMBL941" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "chembl:CHEMBL941");
  });
  await live(t, "chembl_get_activities", { molecule_chembl_id: "CHEMBL941", standard_type: "IC50", limit: 2 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.item_type, "drug_response");
  });
  await live(t, "chembl_get_mechanisms", { molecule_chembl_id: "CHEMBL941", limit: 2 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.ok(result.items.length > 0);
  });
});

test("live openFDA imatinib tools", async (t) => {
  await live(t, "openfda_search_labels", { drug: "imatinib", limit: 1 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.ok(result.items.length > 0);
  });
  await live(t, "openfda_adverse_events", { drug: "imatinib", limit: 5 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.ok(result.items.length > 0);
  });
});

test.after(async () => rm(dir, { recursive: true, force: true }));
