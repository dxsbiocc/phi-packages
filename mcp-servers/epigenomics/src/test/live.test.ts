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

async function live(t: test.TestContext, name: string, args: Record<string, unknown>, itemType: string) {
  if (process.env.PHI_MCP_SKIP_LIVE === "1") {
    t.skip("PHI_MCP_SKIP_LIVE=1");
    return;
  }
  const found = tools.find((item) => item.name === name);
  if (!found) throw new Error(name);
  try {
    const result: ToolResult = ToolResultSchema.parse(await found.call(client, args));
    if (result.entity_type === "error" && (result.code === "unavailable" || result.retryable)) {
      t.skip(result.message);
      return;
    }
    assert.equal(result.entity_type, "page");
    if (result.entity_type === "page") {
      assert.equal(result.item_type, itemType);
      assert.ok(result.items.length > 0);
    }
  } catch (error) {
    if (error instanceof UpstreamError && error.status == null) {
      t.skip(error.message);
      return;
    }
    throw error;
  }
}

test("live ENCODE experiment search", async (t) => live(t, "encode_search_experiments", { limit: 1 }, "dataset"));
test("live Roadmap filter", async (t) => live(t, "encode_search_roadmap", { limit: 1 }, "dataset"));
test("live ENCODE file metadata", async (t) => live(t, "encode_search_files", { limit: 1 }, "file"));
test("live ENCODE biosample", async (t) => live(t, "encode_get_biosample", { accession: "ENCBS997DQU" }, "sample"));
test("live ENCODE experiment and file records", async (t) => {
  await live(t, "encode_get_experiment", { accession: "ENCSR818VNL" }, "dataset");
  await live(t, "encode_get_file", { accession: "ENCFF372VWB" }, "file");
  await live(t, "encode_search_biosamples", { limit: 1 }, "sample");
});

test.after(async () => rm(dir, { recursive: true, force: true }));
