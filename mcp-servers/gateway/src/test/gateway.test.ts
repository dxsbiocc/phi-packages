import assert from "node:assert/strict";
import test from "node:test";
import { mkdtemp, rm } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { UpstreamClient, loadConfig, startHttp } from "phi-omics-runtime";
import { gatewayMounts } from "../mounts.js";

test("gateway health lists every server path", async () => {
  const dir = await mkdtemp(path.join(os.tmpdir(), "phi-omics-gw-"));
  const listener = await startHttp(gatewayMounts, { host: "127.0.0.1", port: 0, basePath: "" }, new UpstreamClient({ ...loadConfig(), cacheDir: dir, cacheTtlMs: 0, retries: 0 }));
  try {
    const response = await fetch(`http://127.0.0.1:${listener.port}/health`);
    const body = await response.json() as { ok: boolean; servers: Array<{ mcp_path: string }> };
    assert.equal(body.ok, true);
    assert.deepEqual(body.servers.map((server) => server.mcp_path), gatewayMounts.map((mount) => mount.path));
  } finally {
    await listener.close();
    await rm(dir, { recursive: true, force: true });
  }
});
