import { existsSync, readdirSync, readFileSync, realpathSync } from "node:fs";
import path from "node:path";
import { DuckDBInstance, type DuckDBConnection, type DuckDBValue } from "@duckdb/node-api";
import { CURIE_PATTERN, taxonCurie, type Access, type ToolResult } from "phi-omics-schemas";
import { failure, InputError, sourceOf } from "phi-omics-runtime";

export type { DuckDBValue };

export const dataSource = sourceOf(
  "omics-sniffer",
  "Omics Sniffer",
  "https://github.com/dxsbiocc/phi-packages",
);

const EXTRA_SPECIES: Record<string, string> = {
  "rattus norvegicus": "NCBITaxon:10116",
  rat: "NCBITaxon:10116",
  "danio rerio": "NCBITaxon:7955",
  zebrafish: "NCBITaxon:7955",
  "drosophila melanogaster": "NCBITaxon:7227",
  "caenorhabditis elegans": "NCBITaxon:6239",
  "saccharomyces cerevisiae": "NCBITaxon:4932",
  "arabidopsis thaliana": "NCBITaxon:3702",
  "sus scrofa": "NCBITaxon:9823",
  "macaca mulatta": "NCBITaxon:9544",
};

let connectionPromise: Promise<DuckDBConnection> | undefined;

async function getConnection(): Promise<DuckDBConnection> {
  if (!connectionPromise) {
    connectionPromise = (async () => {
      const instance = await DuckDBInstance.create(":memory:");
      return instance.connect();
    })();
  }
  return connectionPromise;
}

export interface ServingPaths {
  /** Base configured directory (e.g. /var/lib/phi-omics-sniffer/serving). */
  configuredRoot: string;
  /** Concrete resolved snapshot directory containing datasets.parquet (follows `current` symlink on every call). */
  snapshotDir: string;
}

/**
 * Resolve the active Parquet snapshot directory dynamically on every call so atomic
 * `serving/current -> releases/v_YYYYMMDD_HHMMSS` symlink swaps take effect live
 * without restarting the gateway process or container.
 */
export function resolveServingPaths(): ServingPaths | null {
  const raw = (process.env.PHI_OMICS_SNIFFER_DATA_DIR ?? process.env.OMICS_SNIFFER_SERVING_DIR ?? "").trim();
  if (!raw) return null;
  const configuredRoot = path.resolve(raw);
  if (!existsSync(configuredRoot)) return null;

  const currentCandidate = path.join(configuredRoot, "current");
  if (existsSync(path.join(currentCandidate, "datasets.parquet"))) {
    try {
      return { configuredRoot, snapshotDir: realpathSync(currentCandidate) };
    } catch {
      return null;
    }
  }
  const nestedServingCurrent = path.join(configuredRoot, "serving", "current");
  if (existsSync(path.join(nestedServingCurrent, "datasets.parquet"))) {
    try {
      return { configuredRoot: path.join(configuredRoot, "serving"), snapshotDir: realpathSync(nestedServingCurrent) };
    } catch {
      return null;
    }
  }
  if (existsSync(path.join(configuredRoot, "datasets.parquet"))) {
    try {
      return { configuredRoot, snapshotDir: realpathSync(configuredRoot) };
    } catch {
      return null;
    }
  }
  return null;
}

export async function queryRows(sql: string, params: DuckDBValue[] = []): Promise<Record<string, unknown>[]> {
  const conn = await getConnection();
  const reader = await conn.runAndReadAll(sql, params);
  return reader.getRowObjectsJson() as Record<string, unknown>[];
}

export function partitionFiles(snapshotDir: string, table: "samples" | "dataset_artifacts", sourceHint?: string | null): string | null {
  const tableDir = path.join(snapshotDir, table);
  if (!existsSync(tableDir)) return null;
  if (sourceHint && /^[A-Za-z0-9_.-]+$/.test(sourceHint)) {
    const specific = path.join(tableDir, `source=${sourceHint}`);
    if (existsSync(specific) && readdirSync(specific).some((f) => f.endsWith(".parquet"))) {
      return path.join(specific, "*.parquet");
    }
  }
  try {
    const subdirs = readdirSync(tableDir, { withFileTypes: true }).filter((d) => d.isDirectory() && d.name.startsWith("source="));
    const anyParquet = subdirs.some((d) => readdirSync(path.join(tableDir, d.name)).some((f) => f.endsWith(".parquet")));
    if (!anyParquet) return null;
  } catch {
    return null;
  }
  return path.join(tableDir, "source=*", "*.parquet");
}

export function normalizePrefix(rawPrefix: string): string {
  const cleaned = rawPrefix
    .trim()
    .replace(/[^A-Za-z0-9.]+/g, ".")
    .replace(/^\.+|\.+$/g, "")
    .replace(/\.+/g, ".");
  if (!cleaned) return "omics";
  return cleaned
    .split(".")
    .filter(Boolean)
    .map((seg) => (/^[A-Za-z]/.test(seg) ? seg : `x${seg}`))
    .join(".");
}

export function toCurie(rawId: string | null | undefined, defaultPrefix = "omics"): string {
  const trimmed = (rawId ?? "").trim();
  if (!trimmed) return `${defaultPrefix}:unknown`;
  const colon = trimmed.indexOf(":");
  if (colon > 0) {
    const pfx = normalizePrefix(trimmed.slice(0, colon));
    const local = trimmed.slice(colon + 1).trim().replace(/\s+/g, "_") || "unknown";
    const candidate = `${pfx}:${local}`;
    if (CURIE_PATTERN.test(candidate)) return candidate;
  }
  const safeLocal = trimmed.replace(/\s+/g, "_");
  const fallback = `${normalizePrefix(defaultPrefix)}:${safeLocal}`;
  return CURIE_PATTERN.test(fallback) ? fallback : "omics:unknown";
}

export function normalizeAccess(raw: unknown): Access {
  const val = typeof raw === "string" ? raw.trim().toLowerCase() : "";
  if (val === "public" || val === "open") return "open";
  if (val === "controlled" || val === "restricted" || val === "ega" || val === "dbgap") return "controlled";
  return "unknown";
}

export function mapTaxon(organism: string | null): string | null {
  if (!organism) return null;
  const direct = taxonCurie(organism);
  if (direct) return direct;
  return EXTRA_SPECIES[organism.trim().toLowerCase()] ?? null;
}

export function str(val: unknown): string | null {
  if (typeof val === "string") {
    const t = val.trim();
    return t.length > 0 ? t : null;
  }
  if (typeof val === "number" && Number.isFinite(val)) return String(val);
  if (typeof val === "bigint") return String(val);
  return null;
}

export function intOrNull(val: unknown): number | null {
  if (typeof val === "number" && Number.isFinite(val)) return Math.trunc(val);
  if (typeof val === "bigint") return Number(val);
  if (typeof val === "string" && val.trim() && Number.isFinite(Number(val))) return Math.trunc(Number(val));
  return null;
}

export function parseJsonSafe(raw: string | null, maxLen = 1500): unknown {
  if (!raw) return null;
  try {
    return JSON.parse(raw);
  } catch {
    return raw.length > maxLen ? `${raw.slice(0, maxLen)}...[truncated]` : raw;
  }
}

export function truncateText(raw: string | null, maxLen = 1000): string | null {
  if (!raw) return null;
  return raw.length <= maxLen ? raw : `${raw.slice(0, maxLen)}...[truncated]`;
}

export function bounds(args: { offset?: number; limit?: number }, maxLimit = 50, defaultLimit = 20): { offset: number; limit: number } {
  const offset = args.offset ?? 0;
  const limit = args.limit ?? defaultLimit;
  if (!Number.isInteger(offset) || offset < 0 || !Number.isInteger(limit) || limit < 1 || limit > maxLimit) {
    throw new InputError(`offset must be an integer >= 0 and limit an integer from 1 to ${maxLimit}`);
  }
  return { offset, limit };
}

export function readJsonIfExists(filePath: string): Record<string, unknown> | null {
  if (!existsSync(filePath)) return null;
  try {
    return JSON.parse(readFileSync(filePath, "utf8")) as Record<string, unknown>;
  } catch {
    return null;
  }
}

export function unavailable(message: string): ToolResult {
  return failure(dataSource, "unavailable", message, false);
}
