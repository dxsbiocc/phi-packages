import { createHash } from "node:crypto";
import { mkdir, readFile, writeFile } from "node:fs/promises";
import path from "node:path";
import type { RuntimeConfig } from "./config.js";

export class InputError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "InputError";
  }
}

export class UpstreamError extends Error {
  readonly status: number | null;
  readonly retryable: boolean;
  constructor(status: number | null, message: string, retryable: boolean) {
    super(message);
    this.name = "UpstreamError";
    this.status = status;
    this.retryable = retryable;
  }
}

export interface JsonClient {
  getJson(url: string): Promise<unknown>;
  postJson(url: string, body: unknown): Promise<unknown>;
}

const nextSlot = new Map<string, number>();

async function pace(host: string, rps: number): Promise<void> {
  if (rps <= 0) return;
  const gap = 1000 / rps;
  const now = Date.now();
  const slot = Math.max(now, nextSlot.get(host) ?? 0);
  nextSlot.set(host, slot + gap);
  const wait = slot - now;
  if (wait > 0) await new Promise((resolve) => setTimeout(resolve, wait));
}

function cacheFile(dir: string, key: string): string {
  const hash = createHash("sha256").update(key).digest("hex");
  return path.join(dir, hash.slice(0, 2), `${hash}.json`);
}

async function readCache(file: string): Promise<unknown | undefined> {
  try {
    const parsed = JSON.parse(await readFile(file, "utf8")) as { expires: number; body: unknown };
    if (parsed.expires < Date.now()) return undefined;
    return parsed.body;
  } catch {
    return undefined;
  }
}

async function writeCache(file: string, body: unknown, ttlMs: number): Promise<void> {
  if (ttlMs <= 0) return;
  await mkdir(path.dirname(file), { recursive: true });
  await writeFile(file, JSON.stringify({ expires: Date.now() + ttlMs, body }));
}

export class UpstreamClient implements JsonClient {
  constructor(
    private readonly config: RuntimeConfig,
    private readonly fetchImpl: typeof fetch = fetch,
  ) {}

  getJson(url: string): Promise<unknown> {
    return this.request(url, "GET");
  }

  postJson(url: string, body: unknown): Promise<unknown> {
    return this.request(url, "POST", body);
  }

  private async request(url: string, method: "GET" | "POST", body?: unknown): Promise<unknown> {
    const key = method === "GET" ? url : `${method} ${url} ${JSON.stringify(body)}`;
    const file = cacheFile(this.config.cacheDir, key);
    const cached = await readCache(file);
    if (cached !== undefined) return cached;

    const host = new URL(url).host;
    let last: UpstreamError | null = null;
    for (let attempt = 0; attempt <= this.config.retries; attempt += 1) {
      await pace(host, this.config.rateLimitRps);
      try {
        const response = await this.fetchImpl(url, {
          method,
          headers: {
            accept: "application/json",
            "user-agent": this.config.userAgent,
            ...(body === undefined ? {} : { "content-type": "application/json" }),
          },
          body: body === undefined ? undefined : JSON.stringify(body),
          signal: AbortSignal.timeout(this.config.timeoutMs),
        });
        const text = await response.text();
        if (text.length > this.config.maxBodyBytes) {
          throw new UpstreamError(response.status, `response exceeded ${this.config.maxBodyBytes} bytes`, false);
        }
        if (!response.ok) {
          const retryable = response.status === 429 || response.status >= 500;
          throw new UpstreamError(response.status, `${response.status} ${text.slice(0, 300)}`, retryable);
        }
        let parsed: unknown;
        try {
          parsed = text.length === 0 ? null : JSON.parse(text) as unknown;
        } catch {
          throw new UpstreamError(response.status, "upstream response was not JSON", false);
        }
        await writeCache(file, parsed, this.config.cacheTtlMs);
        return parsed;
      } catch (error) {
        const wrapped = error instanceof UpstreamError
          ? error
          : new UpstreamError(null, error instanceof Error ? error.message : "upstream request failed", true);
        last = wrapped;
        if (!wrapped.retryable || attempt === this.config.retries) break;
        await new Promise((resolve) => setTimeout(resolve, 200 * 2 ** attempt));
      }
    }
    throw last ?? new UpstreamError(null, "upstream request failed", true);
  }
}
