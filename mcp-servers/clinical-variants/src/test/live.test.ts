import assert from "node:assert/strict";
import { mkdtemp, rm } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import test from "node:test";
import { ToolResultSchema, type ToolResult } from "phi-omics-schemas";
import { UpstreamClient, UpstreamError, loadConfig } from "phi-omics-runtime";
import { tools } from "../server.js";

const dir = await mkdtemp(path.join(os.tmpdir(), "phi-variants-live-"));
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

test("live ClinVar tools", async (t) => {
  await live(t, "clinvar_search_variants", { gene: "TP53", clinical_significance: "pathogenic", limit: 2 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") {
      assert.equal(result.item_type, "variant");
      assert.ok(result.items.length > 0);
    }
  });
  await live(t, "clinvar_get_variant", { variation_id: "12351" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.items[0]?.id, "clinvar:12351");
  });
});

test("live gnomAD v4 tools", async (t) => {
  await live(t, "gnomad_get_gene_constraint", { gene: "TP53" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") {
      assert.equal(result.item_type, "gene");
      assert.equal(result.items[0]?.id, "ensembl:ENSG00000141510");
    }
  });
  await live(t, "gnomad_get_variant", { variant_id: "17-7674220-C-T" }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") assert.equal(result.item_type, "variant");
  });
});

test("live CIViC v2 tools", async (t) => {
  await live(t, "civic_search_variants", { gene: "BRAF", variant: "V600E", limit: 2 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") {
      assert.equal(result.item_type, "variant");
      assert.ok(result.items.length > 0);
    }
  });
  await live(t, "civic_search_evidence", { molecular_profile: "BRAF V600E", limit: 2 }, (result) => {
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") {
      assert.equal(result.item_type, "variant");
      assert.ok(result.items.length > 0);
    }
  });
});

test.after(async () => rm(dir, { recursive: true, force: true }));
