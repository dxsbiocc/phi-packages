#!/usr/bin/env node
import { loadConfig, readListenOptions, startHttp, UpstreamClient } from "phi-omics-runtime";
import { gatewayMounts as mounts } from "./mounts.js";

const selected = process.env.MCP_SERVER ?? "all";
const active = selected === "all"
  ? mounts
  : mounts.filter((mount) => mount.module.name === selected || mount.path.includes(`/${selected}/`));

if (active.length === 0) {
  process.stderr.write(`MCP_SERVER must be all or one of ${mounts.map((mount) => mount.module.name).join(", ")}\n`);
  process.exit(1);
}

const options = readListenOptions(process.argv.slice(2));
if (options.help) {
  process.stdout.write(`phi-omics-gateway serves Streamable HTTP MCP endpoints.

  /health
  /transcriptomics/mcp
  /cancer/mcp
  /epigenomics/mcp
  /proteomics/mcp
  /metabolomics/mcp

Set MCP_SERVER to one server name to serve only that server at /mcp.
Default mode is HTTP. --host defaults to 127.0.0.1; set HOST=0.0.0.0 when deploying.
`);
  process.exit(0);
}

const only = active[0];
if (!only) process.exit(1);
const listener = await startHttp(
  selected === "all" ? active : [{ path: "/mcp", module: only.module }],
  { host: options.host, port: options.port, basePath: options.basePath },
  new UpstreamClient(loadConfig()),
);
process.stderr.write(`phi-omics-gateway listening on http://${options.host}:${listener.port}${options.basePath}/health\n`);
