import assert from "node:assert/strict";
import { mkdir, mkdtemp, readdir, rm, utimes, writeFile } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import test from "node:test";
import { z } from "zod";
import type { RuntimeConfig } from "../config.js";
import { withCallContext } from "../context.js";
import { InputError, sweepCache, UpstreamClient, UpstreamError } from "../http.js";
import { InboundLimiter, startHttp, withProvenance, type ServerModule } from "../mcp.js";
import { dataset, makePage, safeMap, sourceOf } from "../page.js";

function config(dir: string, overrides: Partial<RuntimeConfig> = {}): RuntimeConfig {
  return {
    timeoutMs: 2000,
    retries: 0,
    rateLimitRps: 0,
    cacheTtlMs: 0,
    cacheDir: dir,
    userAgent: "test",
    maxBodyBytes: 10_000,
    gtexDataset: "gtex_v10",
    gtexGencode: "v39",
    gtexGenome: "GRCh38/hg38",
    cellguideSnapshot: "1692975365",
    ...overrides,
  };
}

function json(body: unknown, init: ResponseInit = {}): Response {
  return new Response(JSON.stringify(body), { status: 200, headers: { "content-type": "application/json" }, ...init });
}

async function tempDir(): Promise<string> {
  return mkdtemp(path.join(os.tmpdir(), "phi-omics-hardening-"));
}

const source = sourceOf("test", "Test", "https://example.org/");
const emptyPage = () => makePage("dataset", [], 0, 20, 0, source, { summary: "terms", url: null });

test("retrieved_at is the time of the upstream response, not of the process or the call", async () => {
  const dir = await tempDir();
  try {
    const client = new UpstreamClient(config(dir, { cacheTtlMs: 60_000 }), async () => json({ ok: true }));
    await client.getJson("https://example.org/a");
    const firstFetch = Date.now();
    await new Promise((resolve) => setTimeout(resolve, 30));

    const { context } = await withCallContext(() => client.getJson("https://example.org/a"));
    const stamped = withProvenance(emptyPage(), context, Date.now());
    const retrieved = Date.parse(stamped.source.retrieved_at);
    assert.ok(retrieved <= firstFetch + 5, "a cached response must report when it was fetched");
    assert.ok(Date.now() - retrieved >= 25, "and that moment is before the second call");
  } finally {
    await rm(dir, { recursive: true, force: true });
  }
});

test("a call with no upstream fetch is stamped with the call time", () => {
  const now = Date.parse("2026-11-11T00:00:00Z");
  const stamped = withProvenance(emptyPage(), { fetchedAt: [], skipped: 0 }, now);
  assert.equal(stamped.source.retrieved_at, "2026-11-11T00:00:00.000Z");
});

test("the oldest of several upstream responses sets retrieved_at", () => {
  const stamped = withProvenance(emptyPage(), { fetchedAt: [Date.parse("2026-11-10T10:00:00Z"), Date.parse("2026-11-10T08:00:00Z")], skipped: 0 });
  assert.equal(stamped.source.retrieved_at, "2026-11-10T08:00:00.000Z");
});

test("rows that do not fit the schema are skipped and reported, and the cursor moves past them", async () => {
  const { value, context } = await withCallContext(async () => {
    const items = safeMap([{ id: "x:1" }, { id: "" }, { id: "x:3" }], (row) => dataset({ id: row.id }));
    return makePage("dataset", items, 0, 3, 10, source, { summary: "terms", url: null });
  });
  assert.equal(value.entity_type, "page");
  const stamped = withProvenance(value, context);
  if (stamped.entity_type !== "page") throw new Error("expected a page");
  assert.equal(stamped.items.length, 2);
  assert.equal(stamped.source_extra.skipped_rows, 1);
  assert.equal(stamped.page.returned, 2);
  assert.equal(stamped.page.next_offset, 3, "next page starts after the three rows that were read");
});

test("safeMap does not hide upstream failures", () => {
  assert.throws(() => safeMap([1], () => { throw new UpstreamError(500, "boom", true); }), UpstreamError);
});

test("identical concurrent requests share one upstream call", async () => {
  const dir = await tempDir();
  try {
    let calls = 0;
    const client = new UpstreamClient(config(dir), async () => {
      calls += 1;
      await new Promise((resolve) => setTimeout(resolve, 30));
      return json({ n: calls });
    });
    const results = await Promise.all([client.getJson("https://example.org/x"), client.getJson("https://example.org/x"), client.getJson("https://example.org/x")]);
    assert.equal(calls, 1);
    assert.deepEqual(results, [{ n: 1 }, { n: 1 }, { n: 1 }]);
  } finally {
    await rm(dir, { recursive: true, force: true });
  }
});

test("Retry-After lengthens the wait after a 429", async () => {
  const dir = await tempDir();
  try {
    const times: number[] = [];
    const client = new UpstreamClient(config(dir, { retries: 1 }), async () => {
      times.push(Date.now());
      return times.length === 1
        ? new Response("slow down", { status: 429, headers: { "retry-after": "0.4" } })
        : json({ ok: true });
    });
    assert.deepEqual(await client.getJson("https://example.org/limited"), { ok: true });
    assert.ok(times[1]! - times[0]! >= 380, `waited ${times[1]! - times[0]!}ms`);
  } finally {
    await rm(dir, { recursive: true, force: true });
  }
});

test("a body over the limit is refused while it streams, without retrying", async () => {
  const dir = await tempDir();
  try {
    let calls = 0;
    const client = new UpstreamClient(config(dir, { retries: 2, maxBodyBytes: 100 }), async () => {
      calls += 1;
      const stream = new ReadableStream<Uint8Array>({
        start(controller) {
          for (let i = 0; i < 10; i += 1) controller.enqueue(new TextEncoder().encode("x".repeat(50)));
          controller.close();
        },
      });
      return new Response(stream, { status: 200 });
    });
    await assert.rejects(client.getJson("https://example.org/big"), (error: unknown) => error instanceof UpstreamError && !error.retryable && /exceeded 100 bytes/.test(error.message));
    assert.equal(calls, 1);
  } finally {
    await rm(dir, { recursive: true, force: true });
  }
});

test("a declared Content-Length over the limit is refused before reading", async () => {
  const dir = await tempDir();
  try {
    const client = new UpstreamClient(config(dir, { maxBodyBytes: 100 }), async () => new Response("{}", { status: 200, headers: { "content-length": "5000" } }));
    await assert.rejects(client.getJson("https://example.org/declared"), /exceeded 100 bytes/);
  } finally {
    await rm(dir, { recursive: true, force: true });
  }
});

test("the cache sweep drops expired files, then the oldest ones above the size cap", async () => {
  const dir = await tempDir();
  try {
    const now = Date.now();
    const make = async (name: string, ageMs: number, size: number) => {
      const bucket = path.join(dir, name.slice(0, 2));
      await mkdir(bucket, { recursive: true });
      const file = path.join(bucket, `${name}.json`);
      await writeFile(file, "x".repeat(size));
      const when = new Date(now - ageMs);
      await utimes(file, when, when);
      return file;
    };
    await make("aa-expired", 2 * 3600_000, 10);
    await make("bb-oldest", 3000, 400);
    await make("cc-middle", 2000, 400);
    await make("dd-newest", 1000, 400);
    const result = await sweepCache(dir, 3600_000, 1000, now);
    assert.equal(result.removed, 2, "the expired file and the oldest of the rest");
    const left = (await Promise.all((await readdir(dir)).map((bucket) => readdir(path.join(dir, bucket))))).flat().sort();
    assert.deepEqual(left, ["cc-middle.json", "dd-newest.json"]);
  } finally {
    await rm(dir, { recursive: true, force: true });
  }
});

test("the inbound limiter blocks a client over its per-minute allowance and recovers", () => {
  const limiter = new InboundLimiter(2, 1000);
  assert.equal(limiter.check("a", 0), null);
  assert.equal(limiter.check("a", 10), null);
  assert.equal(limiter.check("a", 20), 1);
  assert.equal(limiter.check("b", 20), null, "another client is unaffected");
  assert.equal(limiter.check("a", 1100), null, "a new window starts");
  assert.equal(new InboundLimiter(0).check("a"), null, "zero disables the limit");
});

test("/ready reports a server that is missing its data without failing the probe", async () => {
  const dir = await tempDir();
  const mod: ServerModule = {
    name: "needs-data",
    version: "1.0.0",
    mcpPath: "/mcp",
    tools: [{ name: "usage_terms", title: "Usage", description: "terms", inputSchema: z.object({}), call: async () => emptyPage() }],
    ready: async () => "data directory is missing",
  };
  const previous = process.env.PHI_OMICS_ACCESS_LOG;
  process.env.PHI_OMICS_ACCESS_LOG = "0";
  const listener = await startHttp([{ path: "/mcp", module: mod }], { host: "127.0.0.1", port: 0, basePath: "" }, new UpstreamClient(config(dir)));
  try {
    const health = await fetch(`http://127.0.0.1:${listener.port}/health`);
    assert.equal(health.status, 200);
    const ready = await fetch(`http://127.0.0.1:${listener.port}/ready`);
    const body = await ready.json() as { ok: boolean; servers: { ready: boolean; reason: string | null }[] };
    assert.equal(ready.status, 200);
    assert.equal(body.servers[0]?.ready, false);
    assert.equal(body.servers[0]?.reason, "data directory is missing");
  } finally {
    await listener.close();
    if (previous === undefined) delete process.env.PHI_OMICS_ACCESS_LOG;
    else process.env.PHI_OMICS_ACCESS_LOG = previous;
    await rm(dir, { recursive: true, force: true });
  }
});

test("a tool call over HTTP returns provenance from the real fetch time and skipped-row counts", async () => {
  const dir = await tempDir();
  const mod: ServerModule = {
    name: "provenance",
    version: "1.0.0",
    mcpPath: "/mcp",
    tools: [{
      name: "list_things",
      title: "List",
      description: "list",
      inputSchema: z.object({}),
      call: async (client) => {
        const body = await client.getJson("https://example.org/things") as { ids: string[] };
        const items = safeMap(body.ids, (id) => dataset({ id }));
        return makePage("dataset", items, 0, 3, 3, source, { summary: "terms", url: null });
      },
    }],
  };
  const previous = process.env.PHI_OMICS_ACCESS_LOG;
  process.env.PHI_OMICS_ACCESS_LOG = "0";
  const client = new UpstreamClient(config(dir, { cacheTtlMs: 60_000 }), async () => json({ ids: ["x:1", "", "x:3"] }));
  await client.getJson("https://example.org/things");
  const fetchedAt = Date.now();
  await new Promise((resolve) => setTimeout(resolve, 40));
  const listener = await startHttp([{ path: "/mcp", module: mod }], { host: "127.0.0.1", port: 0, basePath: "" }, client);
  try {
    const response = await fetch(`http://127.0.0.1:${listener.port}/mcp`, {
      method: "POST",
      headers: { "content-type": "application/json", accept: "application/json, text/event-stream" },
      body: JSON.stringify({ jsonrpc: "2.0", id: 2, method: "tools/call", params: { name: "list_things", arguments: {} } }),
    });
    const raw = await response.text();
    assert.equal(response.status, 200, raw);
    const data = raw.split("\n").find((line) => line.startsWith("data:"))?.slice(5) ?? raw;
    const result = (JSON.parse(data) as { result: { structuredContent: { source: { retrieved_at: string }; source_extra: Record<string, unknown>; items: unknown[] } } }).result.structuredContent;
    assert.ok(Date.parse(result.source.retrieved_at) <= fetchedAt + 5, "reports the cached fetch time");
    assert.equal(result.items.length, 2);
    assert.equal(result.source_extra.skipped_rows, 1);
  } finally {
    await listener.close();
    if (previous === undefined) delete process.env.PHI_OMICS_ACCESS_LOG;
    else process.env.PHI_OMICS_ACCESS_LOG = previous;
    await rm(dir, { recursive: true, force: true });
  }
});

test("InputError keeps its name for the MCP wrapper", () => {
  assert.equal(new InputError("bad").name, "InputError");
});
