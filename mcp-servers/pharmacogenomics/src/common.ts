import { existsSync } from "node:fs";
import path from "node:path";
import { DuckDBInstance, type DuckDBConnection, type DuckDBValue } from "@duckdb/node-api";
import { curie, type ToolResult } from "phi-omics-schemas";
import { failure, InputError, sourceOf } from "phi-omics-runtime";

export type { DuckDBValue };

export const dataSource = sourceOf(
  "pharmacogenomics",
  "Pharmacogenomics",
  "https://github.com/dxsbiocc/phi-packages",
);

let cached: { dir: string; connection: DuckDBConnection } | undefined;

export function dataDir(): string | null {
  const configured = process.env.PHI_PGX_DATA_DIR;
  if (!configured) return null;
  const resolved = path.resolve(configured);
  return existsSync(resolved) ? resolved : null;
}

async function database(dir: string): Promise<DuckDBConnection> {
  if (cached?.dir === dir) return cached.connection;
  const instance = await DuckDBInstance.create(":memory:");
  const connection = await instance.connect();
  cached = { dir, connection };
  return connection;
}

export function tablePath(dir: string, table: string, partitioned: boolean): string | null {
  if (!/^[a-z_]+$/.test(table)) throw new InputError(`unknown table ${table}`);
  const root = path.resolve(dir, table);
  if (!root.startsWith(path.resolve(dir) + path.sep)) throw new InputError("invalid data path");
  const file = partitioned ? path.join(root, "**", "*.parquet") : path.join(root, "part.parquet");
  if (!existsSync(root)) return null;
  return file;
}

export async function rows(dir: string, sql: string, params: DuckDBValue[]): Promise<Record<string, unknown>[]> {
  const connection = await database(dir);
  const reader = await connection.runAndReadAll(sql, params);
  return reader.getRowObjectsJson() as Record<string, unknown>[];
}

export function text(value: unknown): string | null {
  if (typeof value === "string" && value.trim()) return value;
  if (typeof value === "number" && Number.isFinite(value)) return String(value);
  return null;
}

export function num(value: unknown): number | null {
  if (typeof value === "number" && Number.isFinite(value)) return value;
  if (typeof value === "bigint") return Number(value);
  if (typeof value === "string" && value.trim() && Number.isFinite(Number(value))) return Number(value);
  return null;
}

export function taxon(species: string | null): string | null {
  if (!species) return null;
  if (species.toLowerCase() === "homo sapiens") return "NCBITaxon:9606";
  return null;
}

export function modelCurie(cvcl: string | null, sourceModelId: string): string {
  if (cvcl) return curie("cellosaurus", cvcl);
  if (sourceModelId.startsWith("ACH-")) return curie("depmap", sourceModelId);
  if (sourceModelId.startsWith("SIDM")) return curie("sanger.model", sourceModelId);
  return curie("pgx.model", sourceModelId);
}

export function bounds(args: Record<string, unknown>): { offset: number; limit: number } {
  const offset = args.offset == null ? 0 : Number(args.offset);
  const limit = args.limit == null ? 20 : Number(args.limit);
  if (!Number.isInteger(offset) || offset < 0 || !Number.isInteger(limit) || limit < 1 || limit > 50) {
    throw new InputError("offset must be an integer >= 0 and limit an integer from 1 to 50");
  }
  return { offset, limit };
}

export function geneKey(value: string): { id: string | null; symbol: string | null } {
  const textValue = value.trim();
  const curieMatch = textValue.match(/^ncbigene:(\d+)$/i);
  if (curieMatch) return { id: curieMatch[1] ?? null, symbol: null };
  if (/^\d+$/.test(textValue)) return { id: textValue, symbol: null };
  return { id: null, symbol: textValue };
}

export function unavailable(message: string): ToolResult {
  return failure(dataSource, "unavailable", message, false);
}


/** Consequences the pipeline treats as silent. A missing consequence counts as non-silent. */
export const SILENT_CONSEQUENCES: readonly string[] = [
  "synonymous_variant",
  "intron_variant",
  "upstream_gene_variant",
  "downstream_gene_variant",
  "5_prime_utr_variant",
  "3_prime_utr_variant",
  "non_coding_transcript_exon_variant",
  "intergenic_variant",
];

/**
 * One row per mapped Cellosaurus model with its tissue and disease. Cell-line names differ across
 * sources, so tissue comes from whichever alias row has it. When the crosswalk is missing the
 * source is empty, so a tissue filter matches nothing rather than being ignored.
 */
export function modelsSource(dir: string): { sql: string; params: DuckDBValue[] } {
  const xref = tablePath(dir, "model_xref", false);
  if (!xref) {
    return { sql: "(SELECT NULL::VARCHAR AS cvcl, NULL::VARCHAR AS tissue, NULL::VARCHAR AS disease WHERE false)", params: [] };
  }
  return {
    sql: `(SELECT cvcl,
                  any_value(tissue) FILTER (WHERE tissue IS NOT NULL AND tissue <> '') AS tissue,
                  any_value(disease) FILTER (WHERE disease IS NOT NULL AND disease <> '') AS disease
           FROM read_parquet(?)
           WHERE status = 'mapped' AND cvcl IS NOT NULL
           GROUP BY cvcl)`,
    params: [xref],
  };
}

/** Number of mapped models whose tissue equals `tissue` (case-insensitive). */
export async function tissueModelCount(dir: string, tissue: string): Promise<number> {
  const models = modelsSource(dir);
  const found = await rows(dir, `SELECT count(*)::INTEGER AS n FROM ${models.sql} AS m WHERE lower(m.tissue) = lower(?)`, [...models.params, tissue]);
  return num(found[0]?.n) ?? 0;
}

export interface TissueFilter {
  sql: string;
  params: DuckDBValue[];
  /** True when a tissue was asked for but no mapped model has it. */
  unknown: boolean;
}

/** SQL fragment (starting with AND) limiting `model_cvcl` to one tissue, plus its parameters. */
export async function tissueFilter(dir: string, tissue: string | null): Promise<TissueFilter> {
  if (!tissue) return { sql: "", params: [], unknown: false };
  const models = modelsSource(dir);
  return {
    sql: `AND model_cvcl IN (SELECT m.cvcl FROM ${models.sql} AS m WHERE lower(m.tissue) = lower(?))`,
    params: [...models.params, tissue],
    unknown: (await tissueModelCount(dir, tissue)) === 0,
  };
}

export function unknownTissue(tissue: string): ToolResult {
  return failure(dataSource, "not_found", `No models with tissue ${tissue}. Call list_tissues for the tissue values in this data.`, false);
}

export function cvclCurie(cvcl: string | null): string | null {
  return cvcl ? curie("cellosaurus", cvcl) : null;
}

export function requireInt(value: unknown, name: string, min: number, max: number): number {
  const parsed = Number(value);
  if (!Number.isInteger(parsed) || parsed < min || parsed > max) throw new InputError(`${name} must be an integer from ${min} to ${max}`);
  return parsed;
}
