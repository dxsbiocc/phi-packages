import { access, constants, mkdir } from "node:fs/promises";
import { createServer, type IncomingMessage, type ServerResponse } from "node:http";
import { createMcpHandler, McpServer } from "@modelcontextprotocol/server";
import { toNodeHandler } from "@modelcontextprotocol/node";
import { serveStdio } from "@modelcontextprotocol/server/stdio";
import { ToolResultSchema, type ToolResult } from "phi-omics-schemas";
import type { ZodType } from "zod";
import { loadConfig } from "./config.js";
import { oldestFetch, withCallContext, type CallContext } from "./context.js";
import { InputError, UpstreamClient, UpstreamError, type JsonClient } from "./http.js";
import { failure, type ToolResult as Result } from "./page.js";

export interface OmicsTool {
  name: string;
  title: string;
  description: string;
  inputSchema: ZodType;
  call(client: JsonClient, args: Record<string, unknown>): Promise<Result>;
}

export interface ServerModule {
  name: string;
  version: string;
  mcpPath: string;
  tools: OmicsTool[];
  /** Optional readiness probe: null when usable, otherwise the reason this server cannot answer. */
  ready?: () => Promise<string | null>;
}

const serverSource = {
  id: "phi-omics-mcp",
  name: "Phi omics MCP",
  version: "1.0.0",
  homepage: "https://github.com/dxsbiocc/phi-packages",
  retrieved_at: new Date().toISOString(),
};

export function toolResult(result: ToolResult): {
  content: Array<{ type: "text"; text: string }>;
  structuredContent: ToolResult;
  isError: boolean;
} {
  const parsed = ToolResultSchema.safeParse(result);
  const value = parsed.success
    ? parsed.data
    : failure(serverSource, "unavailable", "Server produced a result that does not match the shared schema", false);
  return {
    content: [{ type: "text", text: JSON.stringify(value) }],
    structuredContent: value,
    isError: value.entity_type === "error",
  };
}

/**
 * Give a result honest provenance. `retrieved_at` becomes the time of the oldest upstream response
 * behind it (a cached response is older than the call; a call with no upstream fetch uses now), and
 * rows dropped for not fitting the shared schema are counted in `source_extra.skipped_rows`.
 */
export function withProvenance(result: Result, context: CallContext, now = Date.now()): Result {
  const retrievedAt = new Date(oldestFetch(context) ?? now).toISOString();
  const source = { ...result.source, retrieved_at: retrievedAt };
  const skipped = context.skipped;
  const sourceExtra = skipped > 0 ? { ...result.source_extra, skipped_rows: skipped } : result.source_extra;
  if (result.entity_type === "page") {
    // Rows were skipped, so the next page starts after them as well, not right after what was returned.
    const nextOffset = result.page.next_offset == null ? null : result.page.next_offset + skipped;
    return { ...result, source, source_extra: sourceExtra, page: { ...result.page, next_offset: nextOffset } };
  }
  return { ...result, source, source_extra: sourceExtra };
}

export function createServerFromTools(mod: ServerModule, client: JsonClient): McpServer {
  const server = new McpServer({ name: mod.name, version: mod.version });
  for (const tool of mod.tools) {
    server.registerTool(tool.name, {
      title: tool.title,
      description: tool.description,
      inputSchema: tool.inputSchema,
      annotations: { readOnlyHint: true, openWorldHint: true },
    }, async (args: unknown) => {
      try {
        const record = args && typeof args === "object" ? args as Record<string, unknown> : {};
        const started = Date.now();
        const { value, context } = await withCallContext(() => tool.call(client, record));
        const stamped = withProvenance(value, context);
        if (stamped.entity_type === "error") {
          log({ level: "warn", event: "tool_error", server: mod.name, tool: tool.name, code: stamped.code, retryable: stamped.retryable, ms: Date.now() - started });
        }
        return toolResult(stamped);
      } catch (error) {
        log({ level: "warn", event: "tool_exception", server: mod.name, tool: tool.name, message: error instanceof Error ? error.message : "tool failed" });
        if (error instanceof InputError) {
          return toolResult(failure(serverSource, "invalid_input", error.message, false));
        }
        if (error instanceof UpstreamError) {
          const code = error.status === 404 ? "not_found" : error.status == null ? "unavailable" : "upstream_http";
          return toolResult(failure(serverSource, code, error.message, error.retryable, { status: error.status }));
        }
        const message = error instanceof Error ? error.message : "tool failed";
        return toolResult(failure(serverSource, "unavailable", message, false));
      }
    });
  }
  return server;
}

export interface HttpMount {
  path: string;
  module: ServerModule;
}

export interface ListenOptions {
  host: string;
  port: number;
  basePath: string;
}

function normalizeBase(basePath: string): string {
  if (!basePath || basePath === "/") return "";
  const withSlash = basePath.startsWith("/") ? basePath : `/${basePath}`;
  return withSlash.replace(/\/$/, "");
}

function corsOrigin(): string {
  return process.env.PHI_OMICS_CORS_ORIGIN || "*";
}

function sendJson(res: ServerResponse, status: number, body: unknown, extra: Record<string, string> = {}): void {
  const payload = JSON.stringify(body);
  res.writeHead(status, {
    "content-type": "application/json",
    "content-length": Buffer.byteLength(payload),
    "access-control-allow-origin": corsOrigin(),
    "access-control-allow-headers": "content-type, accept, mcp-protocol-version, mcp-session-id",
    "access-control-allow-methods": "GET, POST, DELETE, OPTIONS",
    ...extra,
  });
  res.end(payload);
}

function log(entry: Record<string, unknown>): void {
  process.stderr.write(`${JSON.stringify({ level: "info", time: new Date().toISOString(), ...entry })}\n`);
}

/** Fixed-window request counter per client. A limit of 0 turns it off. */
export class InboundLimiter {
  private readonly windows = new Map<string, { start: number; count: number }>();
  constructor(private readonly perMinute: number, private readonly windowMs = 60_000) {}

  /** Returns null when the request may proceed, otherwise the seconds until the window resets. */
  check(client: string, now = Date.now()): number | null {
    if (this.perMinute <= 0) return null;
    if (this.windows.size > 10_000) {
      for (const [key, value] of this.windows) if (now - value.start >= this.windowMs) this.windows.delete(key);
    }
    const current = this.windows.get(client);
    if (!current || now - current.start >= this.windowMs) {
      this.windows.set(client, { start: now, count: 1 });
      return null;
    }
    current.count += 1;
    return current.count > this.perMinute ? Math.max(1, Math.ceil((current.start + this.windowMs - now) / 1000)) : null;
  }
}

function clientOf(req: IncomingMessage): string {
  if (process.env.PHI_OMICS_TRUST_PROXY === "1") {
    const forwarded = req.headers["x-forwarded-for"];
    const first = (Array.isArray(forwarded) ? forwarded[0] : forwarded)?.split(",")[0]?.trim();
    if (first) return first;
  }
  return req.socket.remoteAddress ?? "unknown";
}

async function cacheWritable(dir: string): Promise<string | null> {
  try {
    await mkdir(dir, { recursive: true });
    await access(dir, constants.W_OK);
    return null;
  } catch (error) {
    return error instanceof Error ? error.message : "cache directory is not writable";
  }
}

export function startHttp(mounts: HttpMount[], options: ListenOptions, client: JsonClient = new UpstreamClient(loadConfig())): Promise<{ port: number; close(): Promise<void> }> {
  const base = normalizeBase(options.basePath);
  const config = loadConfig();
  const limiter = new InboundLimiter(Number(process.env.PHI_OMICS_INBOUND_RPM || 0));
  const accessLog = process.env.PHI_OMICS_ACCESS_LOG !== "0";
  const handlers = new Map(mounts.map((mount) => {
    const handler = createMcpHandler(() => createServerFromTools(mount.module, client));
    return [`${base}${mount.path}`, toNodeHandler(handler)] as const;
  }));
  const httpServer = createServer((req: IncomingMessage, res: ServerResponse) => {
    const started = Date.now();
    const url = new URL(req.url ?? "/", `http://${req.headers.host ?? "localhost"}`);
    if (accessLog) {
      res.once("finish", () => log({ event: "request", method: req.method, path: url.pathname, status: res.statusCode, ms: Date.now() - started }));
    }
    if (req.method === "OPTIONS") {
      sendJson(res, 204, {});
      return;
    }
    if (req.method === "GET" && url.pathname === `${base}/health`) {
      sendJson(res, 200, {
        ok: true,
        servers: mounts.map((mount) => ({ name: mount.module.name, mcp_path: `${base}${mount.path}` })),
      });
      return;
    }
    if (req.method === "GET" && url.pathname === `${base}/ready`) {
      void (async () => {
        const cache = await cacheWritable(config.cacheDir);
        const servers = await Promise.all(mounts.map(async (mount) => {
          const reason = mount.module.ready ? await mount.module.ready().catch((error: unknown) => error instanceof Error ? error.message : "readiness check failed") : null;
          return { name: mount.module.name, mcp_path: `${base}${mount.path}`, ready: reason == null, reason };
        }));
        // A server that is missing its data is reported, but only an unusable cache fails the probe.
        sendJson(res, cache == null ? 200 : 503, { ok: cache == null, cache: { writable: cache == null, reason: cache }, servers });
      })();
      return;
    }
    const handler = handlers.get(url.pathname);
    if (!handler) {
      sendJson(res, 404, { ok: false, message: "not found" });
      return;
    }
    const wait = limiter.check(clientOf(req));
    if (wait != null) {
      sendJson(res, 429, { ok: false, message: "too many requests" }, { "retry-after": String(wait) });
      return;
    }
    void handler(req, res).catch((error: unknown) => {
      log({ level: "error", event: "mcp_error", path: url.pathname, message: error instanceof Error ? error.message : "mcp error" });
      if (!res.headersSent) sendJson(res, 500, { ok: false, message: error instanceof Error ? error.message : "mcp error" });
    });
  });
  return new Promise((resolve, reject) => {
    httpServer.once("error", reject);
    httpServer.listen(options.port, options.host, () => {
      const address = httpServer.address();
      const port = typeof address === "object" && address ? address.port : options.port;
      resolve({
        port,
        close: () => new Promise((done, fail) => httpServer.close((error) => error ? fail(error) : done())),
      });
    });
  });
}

export function startStdio(mod: ServerModule, client: JsonClient = new UpstreamClient(loadConfig())): void {
  serveStdio(() => createServerFromTools(mod, client));
}

export function readListenOptions(argv: string[], env: NodeJS.ProcessEnv = process.env): ListenOptions & { http: boolean; help: boolean } {
  let http = false;
  let help = false;
  let host = env.HOST || "127.0.0.1";
  let port = Number(env.PORT || 8080);
  let basePath = env.BASE_PATH || "";
  for (let i = 0; i < argv.length; i += 1) {
    const arg = argv[i];
    if (arg === "--http") http = true;
    else if (arg === "--stdio") http = false;
    else if (arg === "--help" || arg === "-h") help = true;
    else if (arg === "--host") host = argv[++i] ?? host;
    else if (arg === "--port") port = Number(argv[++i]);
    else if (arg === "--base-path") basePath = argv[++i] ?? "";
    else if (arg?.startsWith("-")) throw new Error(`unknown argument ${arg}`);
  }
  if (!Number.isInteger(port) || port < 0 || port > 65535) throw new Error("--port must be an integer from 0 to 65535");
  return { http, help, host, port, basePath };
}

export function helpText(bin: string, summary: string): string {
  return `${summary}

Usage:
  ${bin} [--stdio]              MCP over stdin/stdout (default; for desktop clients)
  ${bin} --http [--host H] [--port N] [--base-path P]
                                Streamable HTTP MCP on http://H:N/mcp plus /health

Environment:
  HOST, PORT, BASE_PATH
  PHI_OMICS_CACHE_DIR, PHI_OMICS_CACHE_TTL_MS, PHI_OMICS_RATE_LIMIT_RPS
  PHI_OMICS_TIMEOUT_MS, PHI_OMICS_RETRIES, PHI_OMICS_MAX_BODY_BYTES
  PHI_OMICS_USER_AGENT, GTEX_DATASET, GTEX_GENCODE, GTEX_GENOME, CELLGUIDE_SNAPSHOT
  PHI_PGX_DATA_DIR          pharmacogenomics server: processed Parquet directory
  NCBI_API_KEY, PHI_LITERATURE_EMAIL, S2_API_KEY
`;
}
