import { createHash, randomBytes } from "node:crypto";
import { mkdir, readdir, readFile, rename, rm, stat, writeFile } from "node:fs/promises";
import path from "node:path";
import type { RuntimeConfig } from "./config.js";
import { noteFetch } from "./context.js";

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

export interface UpstreamRequest {
  headers?: Record<string, string>;
  accept?: string;
}

export interface JsonClient {
  getJson(url: string, request?: UpstreamRequest): Promise<unknown>;
  postJson(url: string, body: unknown, request?: UpstreamRequest): Promise<unknown>;
  getText?(url: string, request?: UpstreamRequest): Promise<string>;
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

interface CacheHit {
  body: unknown;
  storedAt: number;
}

async function readCache(file: string, ttlMs: number): Promise<CacheHit | undefined> {
  try {
    const parsed = JSON.parse(await readFile(file, "utf8")) as { expires: number; stored?: number; body: unknown };
    if (parsed.expires < Date.now()) return undefined;
    // Entries written before `stored` existed: the best estimate is expiry minus the lifetime.
    return { body: parsed.body, storedAt: parsed.stored ?? parsed.expires - ttlMs };
  } catch {
    return undefined;
  }
}

async function writeCache(file: string, body: unknown, ttlMs: number): Promise<void> {
  if (ttlMs <= 0) return;
  await mkdir(path.dirname(file), { recursive: true });
  const now = Date.now();
  // Write then rename so a concurrent reader never sees half a file.
  const temp = `${file}.${randomBytes(6).toString("hex")}.tmp`;
  await writeFile(temp, JSON.stringify({ expires: now + ttlMs, stored: now, body }));
  await rename(temp, file);
}

const SWEEP_EVERY_MS = 10 * 60 * 1000;
let lastSweep = 0;

function warn(message: string, detail: unknown): void {
  const text = detail instanceof Error ? detail.message : String(detail);
  process.stderr.write(`${JSON.stringify({ level: "warn", message, detail: text })}\n`);
}

/**
 * Keep the on-disk cache bounded: drop files older than the TTL, then the oldest files until the
 * directory is under 80% of the cap. Runs at most every ten minutes and never blocks a request.
 */
export async function sweepCache(dir: string, ttlMs: number, maxBytes: number, now = Date.now()): Promise<{ removed: number; bytes: number }> {
  const files: { file: string; size: number; mtime: number }[] = [];
  for (const bucket of await readdir(dir, { withFileTypes: true }).catch(() => [])) {
    if (!bucket.isDirectory()) continue;
    const bucketDir = path.join(dir, bucket.name);
    for (const name of await readdir(bucketDir).catch(() => [] as string[])) {
      const file = path.join(bucketDir, name);
      const info = await stat(file).catch(() => null);
      if (info?.isFile()) files.push({ file, size: info.size, mtime: info.mtimeMs });
    }
  }
  let removed = 0;
  let kept: typeof files = [];
  for (const entry of files) {
    if (now - entry.mtime > ttlMs) {
      await rm(entry.file, { force: true });
      removed += 1;
    } else {
      kept.push(entry);
    }
  }
  let bytes = kept.reduce((sum, entry) => sum + entry.size, 0);
  if (bytes > maxBytes) {
    kept = kept.sort((left, right) => left.mtime - right.mtime);
    for (const entry of kept) {
      if (bytes <= maxBytes * 0.8) break;
      await rm(entry.file, { force: true });
      bytes -= entry.size;
      removed += 1;
    }
  }
  return { removed, bytes };
}

function maybeSweep(config: RuntimeConfig): void {
  const now = Date.now();
  if (config.cacheTtlMs <= 0 || now - lastSweep < SWEEP_EVERY_MS) return;
  lastSweep = now;
  void sweepCache(config.cacheDir, config.cacheTtlMs, config.cacheMaxBytes ?? DEFAULT_CACHE_MAX_BYTES).catch((error: unknown) => warn("cache sweep failed", error));
}

export const DEFAULT_CACHE_MAX_BYTES = 500_000_000;

function retryAfterMs(response: Response): number | null {
  const header = response.headers?.get?.("retry-after");
  if (!header) return null;
  const seconds = Number(header);
  const ms = Number.isFinite(seconds) ? seconds * 1000 : Date.parse(header) - Date.now();
  return Number.isFinite(ms) && ms > 0 ? Math.min(ms, 10_000) : null;
}

/** Read a response body without buffering more than `max` bytes. */
async function readBounded(response: Response, max: number): Promise<string> {
  const declared = Number(response.headers?.get?.("content-length"));
  if (Number.isFinite(declared) && declared > max) {
    await response.body?.cancel();
    throw new UpstreamError(response.status, `response exceeded ${max} bytes`, false);
  }
  if (!response.body) {
    const text = await response.text();
    if (Buffer.byteLength(text) > max) throw new UpstreamError(response.status, `response exceeded ${max} bytes`, false);
    return text;
  }
  const reader = response.body.getReader();
  const chunks: Uint8Array[] = [];
  let total = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    total += value.byteLength;
    if (total > max) {
      await reader.cancel();
      throw new UpstreamError(response.status, `response exceeded ${max} bytes`, false);
    }
    chunks.push(value);
  }
  return Buffer.concat(chunks).toString("utf8");
}

export class UpstreamClient implements JsonClient {
  constructor(
    private readonly config: RuntimeConfig,
    private readonly fetchImpl: typeof fetch = fetch,
  ) {}

  getJson(url: string, request?: UpstreamRequest): Promise<unknown> {
    return this.request(url, "GET", undefined, request, "json");
  }

  postJson(url: string, body: unknown, request?: UpstreamRequest): Promise<unknown> {
    return this.request(url, "POST", body, request, "json");
  }

  getText(url: string, request?: UpstreamRequest): Promise<string> {
    return this.request(url, "GET", undefined, request, "text").then((body) => String(body));
  }

  private readonly inFlight = new Map<string, Promise<CacheHit>>();

  private async request(url: string, method: "GET" | "POST", body?: unknown, request?: UpstreamRequest, as: "json" | "text" = "json"): Promise<unknown> {
    const headerKey = request?.headers && Object.keys(request.headers).length > 0 ? ` ${JSON.stringify(request.headers)}` : "";
    const acceptKey = request?.accept ? ` accept=${request.accept}` : "";
    const textKey = as === "text" ? "TEXT " : "";
    const key = method === "GET"
      ? `${textKey}${url}${headerKey}${acceptKey}`
      : `${textKey}${method} ${url} ${JSON.stringify(body)}${headerKey}${acceptKey}`;
    const file = cacheFile(this.config.cacheDir, key);
    const cached = await readCache(file, this.config.cacheTtlMs);
    if (cached !== undefined) {
      noteFetch(cached.storedAt);
      return cached.body;
    }

    // Identical concurrent requests share one upstream call.
    let pending = this.inFlight.get(key);
    if (!pending) {
      pending = this.fetchFresh(url, method, body, request, as, file).finally(() => this.inFlight.delete(key));
      this.inFlight.set(key, pending);
    }
    const hit = await pending;
    noteFetch(hit.storedAt);
    return hit.body;
  }

  private async fetchFresh(url: string, method: "GET" | "POST", body: unknown, request: UpstreamRequest | undefined, as: "json" | "text", file: string): Promise<CacheHit> {
    const host = new URL(url).host;
    let last: UpstreamError | null = null;
    let wait = 0;
    for (let attempt = 0; attempt <= this.config.retries; attempt += 1) {
      await pace(host, this.config.rateLimitRps);
      try {
        const response = await this.fetchImpl(url, {
          method,
          headers: {
            accept: request?.accept ?? (as === "json" ? "application/json" : "application/xml, text/xml, text/plain, */*"),
            "user-agent": this.config.userAgent,
            ...(body === undefined ? {} : { "content-type": "application/json" }),
            ...(request?.headers ?? {}),
          },
          body: body === undefined ? undefined : JSON.stringify(body),
          signal: AbortSignal.timeout(this.config.timeoutMs),
        });
        wait = retryAfterMs(response) ?? 0;
        const text = await readBounded(response, this.config.maxBodyBytes);
        if (!response.ok) {
          const retryable = response.status === 429 || response.status >= 500;
          throw new UpstreamError(response.status, `${response.status} ${text.slice(0, 300)}`, retryable);
        }
        let parsed: unknown = text;
        if (as === "json") {
          try {
            parsed = text.length === 0 ? null : JSON.parse(text) as unknown;
          } catch {
            throw new UpstreamError(response.status, "upstream response was not JSON", false);
          }
        }
        const storedAt = Date.now();
        try {
          await writeCache(file, parsed, this.config.cacheTtlMs);
          maybeSweep(this.config);
        } catch (error) {
          warn("cache write failed", error);
        }
        return { body: parsed, storedAt };
      } catch (error) {
        const wrapped = error instanceof UpstreamError
          ? error
          : new UpstreamError(null, error instanceof Error ? error.message : "upstream request failed", true);
        last = wrapped;
        if (!wrapped.retryable || attempt === this.config.retries) break;
        await new Promise((resolve) => setTimeout(resolve, Math.max(200 * 2 ** attempt, wait)));
      }
    }
    throw last ?? new UpstreamError(null, "upstream request failed", true);
  }
}
