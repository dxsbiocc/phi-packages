import { createServer, type IncomingMessage, type ServerResponse } from "node:http";
import { createMcpHandler, McpServer } from "@modelcontextprotocol/server";
import { toNodeHandler } from "@modelcontextprotocol/node";
import { serveStdio } from "@modelcontextprotocol/server/stdio";
import { ToolResultSchema, type ToolResult } from "phi-omics-schemas";
import type { ZodType } from "zod";
import { loadConfig } from "./config.js";
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
        return toolResult(await tool.call(client, record));
      } catch (error) {
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

function sendJson(res: ServerResponse, status: number, body: unknown): void {
  const payload = JSON.stringify(body);
  res.writeHead(status, {
    "content-type": "application/json",
    "content-length": Buffer.byteLength(payload),
    "access-control-allow-origin": "*",
    "access-control-allow-headers": "content-type, accept, mcp-protocol-version, mcp-session-id",
    "access-control-allow-methods": "GET, POST, DELETE, OPTIONS",
  });
  res.end(payload);
}

export function startHttp(mounts: HttpMount[], options: ListenOptions, client: JsonClient = new UpstreamClient(loadConfig())): Promise<{ port: number; close(): Promise<void> }> {
  const base = normalizeBase(options.basePath);
  const handlers = new Map(mounts.map((mount) => {
    const handler = createMcpHandler(() => createServerFromTools(mount.module, client));
    return [`${base}${mount.path}`, toNodeHandler(handler)] as const;
  }));
  const httpServer = createServer((req: IncomingMessage, res: ServerResponse) => {
    const url = new URL(req.url ?? "/", `http://${req.headers.host ?? "localhost"}`);
    if (req.method === "OPTIONS") {
      sendJson(res, 204, {});
      return;
    }
    if (req.method === "GET" && (url.pathname === `${base}/health` || url.pathname === `${base}/ready`)) {
      sendJson(res, 200, {
        ok: true,
        servers: mounts.map((mount) => ({ name: mount.module.name, mcp_path: `${base}${mount.path}` })),
      });
      return;
    }
    const handler = handlers.get(url.pathname);
    if (!handler) {
      sendJson(res, 404, { ok: false, message: "not found" });
      return;
    }
    void handler(req, res).catch((error: unknown) => {
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
  NCBI_API_KEY, PHI_LITERATURE_EMAIL, S2_API_KEY
`;
}
