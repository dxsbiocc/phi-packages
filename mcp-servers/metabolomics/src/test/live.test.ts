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

test("live MetaboLights MTBLS1", async (t) => {
  await live(t, "metabolights_get_study", { study_id: "MTBLS1" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "mtbls:MTBLS1");
  });
  await live(t, "metabolights_list_assays", { study_id: "MTBLS1", limit: 5 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.item_type, "file");
  });
});

test("live HMDB 1-Methylhistidine", async (t) => {
  await live(t, "hmdb_get_metabolite", { hmdb_id: "HMDB0000001" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type !== "page" || result.items[0]?.entity_type !== "metabolite") return;
    assert.equal(result.items[0].id, "hmdb:HMDB0000001");
    assert.equal(result.items[0].formula, "C7H11N3O2");
    assert.equal(result.items[0].inchikey, "BRMWTNUJHUMWMS-LURJTMIESA-N");
    assert.ok(result.items[0].xrefs.includes("CHEBI:50599"));
  });
  await live(t, "hmdb_search_metabolites", { query: "1-Methylhistidine", limit: 5 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.ok(result.items.some((item) => item.id === "hmdb:HMDB0000001"), result.items.map((item) => item.id).join(","));
  });
});

test("live remaining metabolomics tools", async (t) => {
  const calls: Array<[string, Record<string, unknown>, string]> = [
    ["metabolights_list_studies", { limit: 1 }, "mtbls:MTBLS1"],
    ["metabolights_list_files", { study_id: "MTBLS1", limit: 1 }, "mtbls.file:MTBLS1:"],
  ];
  for (const [name, args, idPrefix] of calls) {
    await live(t, name, args, (result) => {
      assert.equal(result.entity_type, "page", name);
      if (result.entity_type === "page") assert.ok(result.items[0]?.id.startsWith(idPrefix), result.items[0]?.id);
    });
  }
});

test.after(async () => rm(dir, { recursive: true, force: true }));
