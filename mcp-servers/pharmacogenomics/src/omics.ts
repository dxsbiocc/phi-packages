import { curie, type ToolResult } from "phi-omics-schemas";
import { expression, InputError, makePage, safeMap, sample, sliceList, tissue as tissueEntity, variant } from "phi-omics-runtime";
import {
  bounds,
  cvclCurie,
  dataDir,
  dataSource,
  geneKey,
  modelsSource,
  num,
  rows,
  SILENT_CONSEQUENCES,
  tablePath,
  text,
  tissueFilter,
  unavailable,
  unknownTissue,
  type DuckDBValue,
} from "./common.js";

const NO_DATA = "Set PHI_PGX_DATA_DIR to the processed pharmacogenomics parquet directory.";
const SILENT_SQL = SILENT_CONSEQUENCES.map((name) => `'${name}'`).join(", ");

interface Scope {
  from: string;
  fromParams: DuckDBValue[];
  where: string;
  whereParams: DuckDBValue[];
}

type ScopeResult = { scope: Scope } | { error: ToolResult };

/** FROM and WHERE for one gene in one per-model table, optionally limited to a dataset and a tissue. */
async function geneScope(
  dir: string,
  table: string,
  args: Record<string, unknown>,
  extra: { sql: string; params: DuckDBValue[] } = { sql: "", params: [] },
): Promise<ScopeResult> {
  const file = tablePath(dir, table, true);
  if (!file) return { error: unavailable(`${table} parquet is not in PHI_PGX_DATA_DIR.`) };
  const geneArg = text(args.gene);
  if (!geneArg) throw new InputError("gene is required");
  const key = geneKey(geneArg);
  const tissue = text(args.tissue);
  const dataset = text(args.dataset);
  const tissueClause = await tissueFilter(dir, tissue);
  if (tissue && tissueClause.unknown) return { error: unknownTissue(tissue) };
  const models = modelsSource(dir);
  return {
    scope: {
      from: `read_parquet(?, hive_partitioning=true) AS f LEFT JOIN ${models.sql} AS mm ON mm.cvcl = f.model_cvcl`,
      fromParams: [file, ...models.params],
      where: `WHERE ((? IS NOT NULL AND f.ncbi_gene_id = ?) OR (? IS NOT NULL AND lower(f.symbol) = lower(?)))
                AND (? IS NULL OR f.dataset = ?) ${tissueClause.sql} ${extra.sql}`,
      whereParams: [key.id, key.id, key.symbol, key.symbol, dataset, dataset, ...tissueClause.params, ...extra.params],
    },
  };
}

async function countRows(dir: string, scope: Scope): Promise<number> {
  const found = await rows(dir, `SELECT count(*)::INTEGER AS n FROM ${scope.from} ${scope.where}`, [...scope.fromParams, ...scope.whereParams]);
  return num(found[0]?.n) ?? 0;
}

function direction(args: Record<string, unknown>): "ASC" | "DESC" {
  const order = text(args.order) ?? "desc";
  if (order !== "asc" && order !== "desc") throw new InputError("order must be asc or desc");
  return order === "asc" ? "ASC" : "DESC";
}

function modelId(row: Record<string, unknown>): string {
  return text(row.model_cvcl) ?? text(row.source_model_id) ?? "unknown";
}

export async function listTissues(args: Record<string, unknown>): Promise<ToolResult> {
  const dir = dataDir();
  if (!dir) return unavailable(NO_DATA);
  if (!tablePath(dir, "model_xref", false)) return unavailable("model_xref parquet is not in PHI_PGX_DATA_DIR.");
  const { offset, limit } = bounds(args);
  const models = modelsSource(dir);
  const found = await rows(
    dir,
    `SELECT m.tissue AS tissue, count(*)::INTEGER AS models FROM ${models.sql} AS m WHERE m.tissue IS NOT NULL GROUP BY 1 ORDER BY models DESC, tissue`,
    models.params,
  );
  const { rows: pageRows, total } = sliceList(found, offset, limit);
  const items = safeMap(pageRows, (row) => {
    const name = text(row.tissue)!;
    return tissueEntity({
      id: curie("pgx.tissue", name.replace(/\s+/g, "_")),
      name,
      source_extra: { models: num(row.models) },
    });
  });
  return makePage("tissue", items, offset, limit, total, dataSource, {
    summary: "Tissue values of mapped cell-line models, most models first. Pass a value as tissue to gene_dependency, gene_expression, gene_copy_number, gene_mutations, compound_response, or biomarker_association. Matching ignores case.",
    url: "https://depmap.org/portal/",
  });
}

async function distribution(dir: string, scope: Scope): Promise<Record<string, number | null>> {
  const found = await rows(
    dir,
    `SELECT count(f.value)::INTEGER AS n, median(f.value) AS median, quantile_cont(f.value, 0.25) AS q1,
            quantile_cont(f.value, 0.75) AS q3, min(f.value) AS min, max(f.value) AS max
     FROM ${scope.from} ${scope.where}`,
    [...scope.fromParams, ...scope.whereParams],
  );
  const row = found[0] ?? {};
  return { n: num(row.n), median: num(row.median), q1: num(row.q1), q3: num(row.q3), min: num(row.min), max: num(row.max) };
}

export async function geneExpression(args: Record<string, unknown>): Promise<ToolResult> {
  const dir = dataDir();
  if (!dir) return unavailable(NO_DATA);
  const built = await geneScope(dir, "expression", args);
  if ("error" in built) return built.error;
  const { scope } = built;
  const { offset, limit } = bounds(args);
  const order = direction(args);
  const total = await countRows(dir, scope);
  const found = await rows(
    dir,
    `SELECT f.dataset, f.source_model_id, f.model_cvcl, f.ncbi_gene_id, f.symbol, f.value, f.metric, mm.tissue, mm.disease
     FROM ${scope.from} ${scope.where}
     ORDER BY f.value ${order} NULLS LAST, f.dataset, f.source_model_id
     LIMIT ? OFFSET ?`,
    [...scope.fromParams, ...scope.whereParams, limit, offset],
  );
  const items = safeMap(found, (row) => {
    const geneId = text(row.ncbi_gene_id);
    return expression({
      id: curie("pgx.expression", `${text(row.dataset)}/${modelId(row)}/${geneId ?? "na"}`),
      gene_id: geneId ? curie("ncbigene", geneId) : null,
      sample_id: cvclCurie(text(row.model_cvcl)),
      value: num(row.value),
      unit: null,
      statistic: text(row.metric),
      source_extra: {
        dataset: text(row.dataset),
        symbol: text(row.symbol),
        source_model_id: text(row.source_model_id),
        tissue: text(row.tissue),
        disease: text(row.disease),
      },
    });
  });
  return makePage("expression", items, offset, limit, total, dataSource, {
    summary: `Gene expression across cell-line models, ${order === "DESC" ? "highest" : "lowest"} first. DepMap values are log2(TPM+1). Compare values within one dataset. The distribution covers every row that matches the filters, not only this page.`,
    url: "https://depmap.org/portal/",
  }, { distribution: await distribution(dir, scope) });
}

export async function geneCopyNumber(args: Record<string, unknown>): Promise<ToolResult> {
  const dir = dataDir();
  if (!dir) return unavailable(NO_DATA);
  const built = await geneScope(dir, "cnv", args);
  if ("error" in built) return built.error;
  const { scope } = built;
  const { offset, limit } = bounds(args);
  const order = direction(args);
  const total = await countRows(dir, scope);
  const found = await rows(
    dir,
    `SELECT f.dataset, f.source_model_id, f.model_cvcl, f.ncbi_gene_id, f.symbol, f.value, f.metric, mm.tissue, mm.disease
     FROM ${scope.from} ${scope.where}
     ORDER BY f.value ${order} NULLS LAST, f.dataset, f.source_model_id
     LIMIT ? OFFSET ?`,
    [...scope.fromParams, ...scope.whereParams, limit, offset],
  );
  const items = safeMap(found, (row) => {
    const geneId = text(row.ncbi_gene_id);
    const cvcl = text(row.model_cvcl);
    return sample({
      id: cvclCurie(cvcl) ?? curie("pgx.model", modelId(row)),
      name: cvcl,
      taxon_id: "NCBITaxon:9606",
      xrefs: geneId ? [curie("ncbigene", geneId)] : [],
      source_extra: {
        dataset: text(row.dataset),
        gene_id: geneId ? curie("ncbigene", geneId) : null,
        symbol: text(row.symbol),
        copy_number: num(row.value),
        metric: text(row.metric),
        source_model_id: text(row.source_model_id),
        tissue: text(row.tissue),
        disease: text(row.disease),
      },
    });
  });
  return makePage("sample", items, offset, limit, total, dataSource, {
    summary: `Gene copy number across cell-line models, ${order === "DESC" ? "highest" : "lowest"} first. The metric field names the transform. Compare values within one dataset. The distribution covers every row that matches the filters.`,
    url: "https://depmap.org/portal/",
  }, { distribution: await distribution(dir, scope) });
}

function proteinChange(value: string): string {
  return value.trim().replace(/^p\./i, "").toLowerCase();
}

export async function geneMutations(args: Record<string, unknown>): Promise<ToolResult> {
  const dir = dataDir();
  if (!dir) return unavailable(NO_DATA);
  const includeSilent = args.include_silent === true;
  const consequence = text(args.consequence);
  const change = text(args.protein_change);
  const silentSql = includeSilent ? "" : `AND lower(coalesce(f.consequence, '')) NOT IN (${SILENT_SQL})`;
  // Asking for a consequence by name overrides the silent default, so synonymous_variant can be requested.
  const filters: string[] = [consequence ? "" : silentSql];
  const filterParams: DuckDBValue[] = [];
  if (consequence) {
    // DepMap joins several annotations with "&", such as missense_variant&splice_region_variant.
    filters.push("AND list_contains(string_split(lower(coalesce(f.consequence, '')), '&'), lower(?))");
    filterParams.push(consequence);
  }
  if (change) {
    filters.push("AND lower(replace(coalesce(f.protein_change, ''), 'p.', '')) = ?");
    filterParams.push(proteinChange(change));
  }
  const filtered = await geneScope(dir, "mutation", args, { sql: filters.join(" "), params: filterParams });
  if ("error" in filtered) return filtered.error;
  const base = await geneScope(dir, "mutation", args, { sql: silentSql, params: [] });
  if ("error" in base) return base.error;
  const { offset, limit } = bounds(args);
  const total = await countRows(dir, filtered.scope);
  const found = await rows(
    dir,
    `SELECT f.dataset, f.source_model_id, f.model_cvcl, f.ncbi_gene_id, f.symbol, f.chrom_norm, f.pos, f.ref, f.alt,
            f.protein_change, f.consequence, f.genome_build, f.vep_impact, mm.tissue, mm.disease
     FROM ${filtered.scope.from} ${filtered.scope.where}
     ORDER BY f.model_cvcl NULLS LAST, f.source_model_id, f.pos, f.ref, f.alt, f.dataset
     LIMIT ? OFFSET ?`,
    [...filtered.scope.fromParams, ...filtered.scope.whereParams, limit, offset],
  );
  const items = safeMap(found, (row) => {
    const geneId = text(row.ncbi_gene_id);
    const chrom = text(row.chrom_norm);
    const position = num(row.pos);
    return variant({
      id: curie("pgx.variant", [text(row.dataset) ?? "mut", modelId(row), chrom ?? "na", String(position ?? 0), text(row.ref) ?? "N", text(row.alt) ?? "N"].join("/")),
      gene_id: geneId ? curie("ncbigene", geneId) : null,
      consequence: text(row.consequence),
      locus: chrom ? { genome_build: text(row.genome_build) ?? "GRCh38", chrom, start: position, end: position, strand: null } : null,
      source_extra: {
        dataset: text(row.dataset),
        symbol: text(row.symbol),
        protein_change: text(row.protein_change),
        vep_impact: text(row.vep_impact),
        sample_id: cvclCurie(text(row.model_cvcl)),
        source_model_id: text(row.source_model_id),
        tissue: text(row.tissue),
        disease: text(row.disease),
      },
    });
  });
  const summary = await mutationSummary(dir, base.scope, args);
  return makePage("variant", items, offset, limit, total, dataSource, {
    summary: `Somatic mutations of one gene across cell-line models. ${includeSilent ? "Silent variants are included." : "Silent variants (synonymous, intronic, UTR, up/downstream, intergenic) are left out unless include_silent is true."} Rows are grouped by model. frequency counts models with any non-silent call in the gene among models present in the mutation table, so a model with no call is read as wild type only if it was sequenced.`,
    url: "https://depmap.org/portal/",
  }, summary);
}

async function mutationSummary(dir: string, base: Scope, args: Record<string, unknown>): Promise<Record<string, unknown>> {
  const tissue = text(args.tissue);
  const dataset = text(args.dataset);
  const tissueClause = await tissueFilter(dir, tissue);
  const file = tablePath(dir, "mutation", true)!;
  const nonSilent = `AND lower(coalesce(f.consequence, '')) NOT IN (${SILENT_SQL})`;
  const mutated = await rows(
    dir,
    `SELECT count(DISTINCT f.model_cvcl)::INTEGER AS n FROM ${base.from} ${base.where} ${nonSilent}`,
    [...base.fromParams, ...base.whereParams],
  );
  const sequenced = await rows(
    dir,
    `SELECT count(DISTINCT model_cvcl)::INTEGER AS n FROM read_parquet(?, hive_partitioning=true)
     WHERE model_cvcl IS NOT NULL AND (? IS NULL OR dataset = ?) ${tissueClause.sql}`,
    [file, dataset, dataset, ...tissueClause.params],
  );
  const hotspots = await rows(
    dir,
    `SELECT f.protein_change AS protein_change, count(DISTINCT f.model_cvcl)::INTEGER AS models
     FROM ${base.from} ${base.where} ${nonSilent} AND f.protein_change IS NOT NULL AND f.protein_change <> ''
     GROUP BY 1 ORDER BY models DESC, protein_change LIMIT 10`,
    [...base.fromParams, ...base.whereParams],
  );
  const mutatedModels = num(mutated[0]?.n) ?? 0;
  const sequencedModels = num(sequenced[0]?.n) ?? 0;
  return {
    mutated_models: mutatedModels,
    sequenced_models: sequencedModels,
    frequency: sequencedModels > 0 ? mutatedModels / sequencedModels : null,
    top_protein_changes: hotspots.map((row) => ({ protein_change: text(row.protein_change), models: num(row.models) })),
  };
}
