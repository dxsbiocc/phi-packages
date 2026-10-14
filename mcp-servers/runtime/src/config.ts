import os from "node:os";
import path from "node:path";

export interface RuntimeConfig {
  timeoutMs: number;
  retries: number;
  rateLimitRps: number;
  cacheTtlMs: number;
  cacheDir: string;
  userAgent: string;
  maxBodyBytes: number;
  gtexDataset: string;
  gtexGencode: string;
  gtexGenome: string;
  cellguideSnapshot: string;
}

export function cacheDirFor(platform: NodeJS.Platform, env: NodeJS.ProcessEnv, home: string): string {
  if (env.PHI_OMICS_CACHE_DIR) return env.PHI_OMICS_CACHE_DIR;
  if (platform === "win32") {
    const base = env.LOCALAPPDATA || path.join(home, "AppData", "Local");
    return path.join(base, "phi-omics-mcp", "cache");
  }
  const base = env.XDG_CACHE_HOME || path.join(home, ".cache");
  return path.join(base, "phi-omics-mcp");
}

function integer(env: NodeJS.ProcessEnv, name: string, fallback: number): number {
  const raw = env[name];
  if (raw == null || raw === "") return fallback;
  const value = Number(raw);
  if (!Number.isFinite(value)) throw new Error(`${name} must be a number`);
  return value;
}

export function loadConfig(env: NodeJS.ProcessEnv = process.env): RuntimeConfig {
  return {
    timeoutMs: integer(env, "PHI_OMICS_TIMEOUT_MS", 30_000),
    retries: integer(env, "PHI_OMICS_RETRIES", 2),
    rateLimitRps: integer(env, "PHI_OMICS_RATE_LIMIT_RPS", 4),
    cacheTtlMs: integer(env, "PHI_OMICS_CACHE_TTL_MS", 6 * 60 * 60 * 1000),
    cacheDir: cacheDirFor(process.platform, env, os.homedir()),
    userAgent: env.PHI_OMICS_USER_AGENT || "phi-omics-mcp/1.0 (+https://github.com/dxsbiocc/phi-packages)",
    maxBodyBytes: integer(env, "PHI_OMICS_MAX_BODY_BYTES", 8_000_000),
    gtexDataset: env.GTEX_DATASET || "gtex_v10",
    gtexGencode: env.GTEX_GENCODE || "v39",
    gtexGenome: env.GTEX_GENOME || "GRCh38/hg38",
    cellguideSnapshot: env.CELLGUIDE_SNAPSHOT || "1692975365",
  };
}
