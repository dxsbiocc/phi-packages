import assert from "node:assert/strict";
import { mkdtemp, rm } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import { cacheDirFor } from "../config.js";
import { UpstreamClient } from "../http.js";
import { file, makePage, sourceOf } from "../page.js";
import { startHttp, toolResult, type ServerModule } from "../mcp.js";
import { z } from "zod";

test("cache directory follows the platform convention", () => {
  assert.equal(cacheDirFor("linux", {}, "/home/person"), "/home/person/.cache/phi-omics-mcp");
  assert.equal(cacheDirFor("darwin", { XDG_CACHE_HOME: "/tmp/cache" }, "/Users/person"), "/tmp/cache/phi-omics-mcp");
  assert.equal(
    cacheDirFor("win32", { LOCALAPPDATA: "C:\\Users\\person\\AppData\\Local" }, "C:\\Users\\person"),
    path.join("C:\\Users\\person\\AppData\\Local", "phi-omics-mcp", "cache"),
  );
  assert.equal(cacheDirFor("linux", { PHI_OMICS_CACHE_DIR: "/var/cache/omics" }, "/home/person"), "/var/cache/omics");
});

test("controlled file builder strips download URLs", () => {
  const result = file({
    id: "gdc.file:af0d3d92-f8b5-4612-805b-d6773595262b",
    access: "controlled",
    url: "https://api.gdc.cancer.gov/data/af0d3d92-f8b5-4612-805b-d6773595262b",
  });
  assert.equal(result.url, null);
});

test("tool results are validated before they are returned", () => {
  const page = makePage("dataset", [], 0, 1, 0, sourceOf("gtex", "GTEx", "https://gtexportal.org/"), { summary: "terms", url: null });
  const wrapped = toolResult(page);
  assert.equal(ToolResultSchema.parse(JSON.parse(wrapped.content[0].text)).entity_type, "page");
  assert.equal(wrapped.isError, false);
});

test("HTTP mode serves health and a Streamable HTTP initialize", async () => {
  const mod: ServerModule = {
    name: "runtime-test",
    version: "1.0.0",
    mcpPath: "/mcp",
    tools: [{
      name: "usage_terms",
      title: "Usage",
      description: "Return usage terms",
      inputSchema: z.object({}),
      call: async () => makePage("dataset", [], 0, 1, 0, sourceOf("gtex", "GTEx", "https://gtexportal.org/"), { summary: "terms", url: null }),
    }],
  };
  const dir = await mkdtemp(path.join(os.tmpdir(), "phi-omics-"));
  const client = new UpstreamClient({
    timeoutMs: 1000,
    retries: 0,
    rateLimitRps: 0,
    cacheTtlMs: 0,
    cacheDir: dir,
    userAgent: "test",
    maxBodyBytes: 1000,
    gtexDataset: "gtex_v10",
    gtexGencode: "v39",
    gtexGenome: "GRCh38/hg38",
    cellguideSnapshot: "1692975365",
  }, fetch);
  const listener = await startHttp([{ path: "/mcp", module: mod }], { host: "127.0.0.1", port: 0, basePath: "" }, client);
  try {
    const health = await fetch(`http://127.0.0.1:${listener.port}/health`);
    assert.equal(health.status, 200);
    const body = await health.json() as { ok: boolean };
    assert.equal(body.ok, true);
    const init = await fetch(`http://127.0.0.1:${listener.port}/mcp`, {
      method: "POST",
      headers: { "content-type": "application/json", accept: "application/json, text/event-stream" },
      body: JSON.stringify({
        jsonrpc: "2.0",
        id: 1,
        method: "initialize",
        params: { protocolVersion: "2025-03-26", capabilities: {}, clientInfo: { name: "test", version: "0" } },
      }),
    });
    const raw = await init.text();
    assert.equal(init.status, 200, raw);
    assert.match(raw, /runtime-test/);
  } finally {
    await listener.close();
    await rm(dir, { recursive: true, force: true });
  }
});
