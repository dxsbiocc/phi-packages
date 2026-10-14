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

test("live GDC TP53", async (t) => {
  await live(t, "gdc_get_gene", { gene: "ENSG00000141510" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "ensembl:ENSG00000141510");
  });
});

test("live GDC controlled file metadata has no url", async (t) => {
  await live(t, "gdc_search_files", { access: "controlled", limit: 1 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page" && result.items[0]?.entity_type === "file") {
      assert.equal(result.items[0].access, "controlled");
      assert.equal(result.items[0].url, null);
    }
  });
});

test("live cBioPortal TP53 and brca_tcga", async (t) => {
  await live(t, "cbioportal_get_gene", { gene: "TP53" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "ncbigene:7157");
  });
  await live(t, "cbioportal_get_study", { study_id: "brca_tcga" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "cbioportal.study:brca_tcga");
  });
});

test("live PDC study", async (t) => {
  await live(t, "pdc_get_study", { study_id: "PDC000759" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "pdc.study:PDC000759");
  });
});

test("live remaining cancer tools", async (t) => {
  const calls: Array<[string, Record<string, unknown>]> = [
    ["gdc_search_cases", { limit: 1 }],
    ["gdc_get_case", { case_id: "TCGA-44-3918" }],
    ["gdc_search_mutations", { limit: 1 }],
    ["cbioportal_list_studies", { limit: 1 }],
    ["cbioportal_list_samples", { study_id: "brca_tcga", limit: 1 }],
    ["cbioportal_list_patients", { study_id: "brca_tcga", limit: 1 }],
    ["pdc_list_studies", { limit: 1 }],
  ];
  for (const [name, args] of calls) {
    await live(t, name, args, (result) => {
      assert.equal(result.entity_type, "page", name);
      if (result.entity_type === "page") assert.ok(result.items.length > 0, name);
    });
  }
});

test.after(async () => rm(dir, { recursive: true, force: true }));
