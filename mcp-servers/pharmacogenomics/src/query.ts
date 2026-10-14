import { existsSync } from "node:fs";
import path from "node:path";
import { DuckDBInstance, type DuckDBConnection, type DuckDBValue } from "@duckdb/node-api";
import { curie, type ToolResult } from "phi-omics-schemas";
import {
  InputError,
  drugResponse,
  expression,
  failure,
  gene,
  makePage,
  sample,
  sourceOf,
  variant,
} from "phi-omics-runtime";

const dataSource = sourceOf(
  "pharmacogenomics",
  "Pharmacogenomics",
  "https://github.com/dxsbiocc/phi-packages",
);
const usage = {
  summary: "",
  url: "https://depmap.org/portal/about/",
};

export const USAGE_SUMMARY = [
  "Cell-line pharmacogenomics for research. These screens are not treatment recommendations.",
  "DepMap 24Q2 Public, PRISM Repurposing 19Q4 (primary and secondary screens, Corsello et al. 2020), and PRISM Repurposing Public 24Q2 are CC BY 4.0 on their Figshare records.",
  "GDSC1/GDSC2 release 8.4 is the 24 Jul 2022 Sanger single-agent dose-response release. GDSC² combination screens (Jaaks et al., Nature 2022) are a separate source and are not CC BY. Cell Model Passports supplies model identifiers.",
  "Cell lines are Cellosaurus CVCL accessions. Compounds are PubChem CIDs when a single CID is known. Unmatched and conflicting identifiers are left unresolved.",
  "NCI-60 growth-inhibition data are a US National Cancer Institute DTP release. Credit the National Cancer Institute. They are not CC BY. The pinned files are the October 2026 wiki attachments.",
  "ORCESTRA PharmacoSets are frozen Zenodo records. CTRPv2.1 (10.5281/zenodo.21794049), gCSI 2019 (10.5281/zenodo.4737437), UHNBreast 2019 (10.5281/zenodo.3601415), and CCLE 2019 (10.5281/zenodo.18759113) are CC BY 4.0. FIMM (10.5281/zenodo.3905448) and GRAY 2017 (10.5281/zenodo.3905458) are CC0.",
  "DepMap expression values are log2(TPM+1). CRISPR values are Chronos gene-effect scores. GDSC LN_IC50 is the natural log of micromolar IC50, and a lower AUC or LN_IC50 means a stronger response.",
].join(" ");

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

function tablePath(dir: string, table: string, partitioned: boolean): string | null {
  if (!/^[a-z_]+$/.test(table)) throw new InputError(`unknown table ${table}`);
  const root = path.resolve(dir, table);
  if (!root.startsWith(path.resolve(dir) + path.sep)) throw new InputError("invalid data path");
  const file = partitioned ? path.join(root, "**", "*.parquet") : path.join(root, "part.parquet");
  if (!existsSync(root)) return null;
  return file;
}

async function rows(dir: string, sql: string, params: DuckDBValue[]): Promise<Record<string, unknown>[]> {
  const connection = await database(dir);
  const reader = await connection.runAndReadAll(sql, params);
  return reader.getRowObjectsJson() as Record<string, unknown>[];
}

function text(value: unknown): string | null {
  if (typeof value === "string" && value.trim()) return value;
  if (typeof value === "number" && Number.isFinite(value)) return String(value);
  return null;
}

function num(value: unknown): number | null {
  if (typeof value === "number" && Number.isFinite(value)) return value;
  if (typeof value === "bigint") return Number(value);
  if (typeof value === "string" && value.trim() && Number.isFinite(Number(value))) return Number(value);
  return null;
}

function taxon(species: string | null): string | null {
  if (!species) return null;
  if (species.toLowerCase() === "homo sapiens") return "NCBITaxon:9606";
  return null;
}

function modelCurie(cvcl: string | null, sourceModelId: string): string {
  if (cvcl) return curie("cellosaurus", cvcl);
  if (sourceModelId.startsWith("ACH-")) return curie("depmap", sourceModelId);
  if (sourceModelId.startsWith("SIDM")) return curie("sanger.model", sourceModelId);
  return curie("pgx.model", sourceModelId);
}

function bounds(args: Record<string, unknown>): { offset: number; limit: number } {
  const offset = args.offset == null ? 0 : Number(args.offset);
  const limit = args.limit == null ? 20 : Number(args.limit);
  if (!Number.isInteger(offset) || offset < 0 || !Number.isInteger(limit) || limit < 1 || limit > 50) {
    throw new InputError("offset must be an integer >= 0 and limit an integer from 1 to 50");
  }
  return { offset, limit };
}

function geneKey(value: string): { id: string | null; symbol: string | null } {
  const textValue = value.trim();
  const curieMatch = textValue.match(/^ncbigene:(\d+)$/i);
  if (curieMatch) return { id: curieMatch[1] ?? null, symbol: null };
  if (/^\d+$/.test(textValue)) return { id: textValue, symbol: null };
  return { id: null, symbol: textValue };
}

function unavailable(message: string): ToolResult {
  return failure(dataSource, "unavailable", message, false);
}

export async function usageTerms(): Promise<ToolResult> {
  const dir = dataDir();
  let summary = USAGE_SUMMARY;
  if (!dir) {
    summary += " PHI_PGX_DATA_DIR is not set to a processed parquet directory, so the query tools are unavailable.";
  } else {
    const releases = tablePath(dir, "releases", false);
    if (releases) {
      const found = await rows(dir, "SELECT source, version, citation FROM read_parquet(?) ORDER BY source", [releases]);
      const lines = found.map((row) => `${text(row.source)} ${text(row.version)}: ${text(row.citation)}`);
      if (lines.length) summary += ` Loaded releases: ${lines.join(" | ")}`;
    }
  }
  return makePage("dataset", [], 0, 1, 0, dataSource, { summary, url: usage.url });
}

export async function resolveCellLine(args: Record<string, unknown>): Promise<ToolResult> {
  const dir = dataDir();
  if (!dir) return unavailable("Set PHI_PGX_DATA_DIR to the processed pharmacogenomics parquet directory.");
  const xref = tablePath(dir, "model_xref", false);
  if (!xref) return unavailable("model_xref parquet is not in PHI_PGX_DATA_DIR.");
  const query = text(args.query);
  if (!query) throw new InputError("query is required");
  const { offset, limit } = bounds(args);
  const bare = query.replace(/^cellosaurus:/i, "");
  const matched = await rows(
    dir,
    `SELECT source, source_model_id, alias, alias_type, cvcl, match_method, status, evidence,
            display_name, tissue, disease, species
     FROM read_parquet(?)
     WHERE lower(alias) = lower(?) OR cvcl = ? OR lower(alias) = lower(?)`,
    [xref, query, bare, bare],
  );
  if (!matched.length) return failure(dataSource, "not_found", `No cell line alias matched ${query}`, false);
  const grouped = new Map<string, typeof matched>();
  for (const row of matched) {
    const key = `${text(row.source)}:${text(row.source_model_id)}`;
    grouped.set(key, [...(grouped.get(key) ?? []), row]);
  }
  const byCvcl = new Map<string, Array<typeof matched>>();
  const separate: Array<typeof matched> = [];
  for (const group of grouped.values()) {
    const cvcl = text(group[0]?.cvcl);
    if (cvcl && text(group[0]?.status) === "mapped") {
      byCvcl.set(cvcl, [...(byCvcl.get(cvcl) ?? []), group]);
    } else {
      separate.push(group);
    }
  }
  const models = [...byCvcl.values(), ...separate.map((group) => [group])];
  const pageRows = models.slice(offset, offset + limit);
  const items = pageRows.map((groups) => {
    const group = groups.flat();
    const first = group[0]!;
    const cvcl = text(first.cvcl);
    const sourceModelId = text(first.source_model_id) ?? query;
    const xrefs = group.flatMap((row) => {
      const alias = text(row.alias);
      const kind = text(row.alias_type);
      if (!alias) return [];
      if (kind === "rrid") return [curie("cellosaurus", alias)];
      if (kind === "ach" || (kind === "primary" && alias.startsWith("ACH-"))) return [curie("depmap", alias)];
      if (kind === "sidm" || (kind === "primary" && alias.startsWith("SIDM"))) return [curie("sanger.model", alias)];
      if (kind === "cosmic") return [curie("cosmic", alias)];
      return [];
    });
    return sample({
      id: modelCurie(cvcl, sourceModelId),
      name: text(first.display_name),
      taxon_id: taxon(text(first.species)),
      xrefs: [...new Set(xrefs.filter((item) => item !== modelCurie(cvcl, sourceModelId)))],
      source_extra: {
        status: text(first.status),
        match_method: text(first.match_method),
        evidence: text(first.evidence),
        cvcl,
        source: text(first.source),
        source_model_id: sourceModelId,
        tissue: text(first.tissue),
        disease: text(first.disease),
      },
    });
  });
  const cvcls = new Set(models.map((groups) => text(groups[0]?.[0]?.cvcl)).filter((item): item is string => item != null));
  return makePage("sample", items, offset, limit, models.length, dataSource, {
    summary: cvcls.size > 1 ? "This alias matches more than one Cellosaurus accession. None was chosen automatically." : USAGE_SUMMARY,
    url: "https://www.cellosaurus.org/",
  }, { ambiguous: cvcls.size > 1 });
}

export async function geneDependency(args: Record<string, unknown>): Promise<ToolResult> {
  const dir = dataDir();
  if (!dir) return unavailable("Set PHI_PGX_DATA_DIR to the processed pharmacogenomics parquet directory.");
  const table = tablePath(dir, "gene_effect", true);
  if (!table) return unavailable("gene_effect parquet is not in PHI_PGX_DATA_DIR. Download CRISPRGeneEffect.csv with phi-pgx download --live and rerun process.");
  const geneArg = text(args.gene);
  if (!geneArg) throw new InputError("gene is required");
  const key = geneKey(geneArg);
  const dataset = text(args.dataset);
  const { offset, limit } = bounds(args);
  const where = `WHERE ((? IS NOT NULL AND ncbi_gene_id = ?) OR (? IS NOT NULL AND lower(symbol) = lower(?)))
                  AND (? IS NULL OR dataset = ?)`;
  const params = [key.id, key.id, key.symbol, key.symbol, dataset, dataset];
  const total = num((await rows(dir, `SELECT count(*)::INTEGER AS n FROM read_parquet(?, hive_partitioning=true) ${where}`, [table, ...params]))[0]?.n);
  const found = await rows(
    dir,
    `SELECT dataset, source_model_id, model_cvcl, ncbi_gene_id, symbol, value, metric
     FROM read_parquet(?, hive_partitioning=true) ${where}
     ORDER BY value ASC NULLS LAST
     LIMIT ? OFFSET ?`,
    [table, ...params, limit, offset],
  );
  const items = found.map((row) => {
    const cvcl = text(row.model_cvcl);
    const sourceModelId = text(row.source_model_id) ?? "unknown";
    const geneId = text(row.ncbi_gene_id);
    return sample({
      id: modelCurie(cvcl, sourceModelId),
      name: cvcl,
      taxon_id: "NCBITaxon:9606",
      xrefs: geneId ? [curie("ncbigene", geneId)] : [],
      source_extra: {
        dataset: text(row.dataset),
        gene_id: geneId ? curie("ncbigene", geneId) : null,
        symbol: text(row.symbol),
        gene_effect: num(row.value),
        metric: text(row.metric),
        source_model_id: sourceModelId,
        cvcl,
      },
    });
  });
  return makePage("sample", items, offset, limit, total, dataSource, {
    summary: "Chronos gene-effect scores. More negative values are stronger dependencies. Rows are ordered by value ascending.",
    url: "https://doi.org/10.25452/figshare.plus.25880521.v1",
  });
}

function combinationIdentity(query: string, row: Record<string, unknown>): boolean {
  const id = text(row.combination_id) ?? "";
  const local = id.startsWith("gdsc.combination:") ? id.slice("gdsc.combination:".length) : id;
  const name = (text(row.name) ?? "").toLowerCase();
  const folded = query.trim().toLowerCase();
  return folded === id.toLowerCase() || folded === local.toLowerCase() || folded === name;
}

function componentCurie(cid: string | null, chembl: string | null): string | null {
  if (cid) return curie("pubchem", cid);
  if (chembl) return curie("chembl", chembl);
  return null;
}

function splitRole(value: string | null): string[] {
  if (!value) return [];
  return value.split("|");
}

function roleComponents(
  role: "anchor" | "library",
  ids: string | null,
  names: string | null,
  cids: string | null,
  chembls: string | null,
  concentrations: number[],
  unit: string | null,
): Array<Record<string, unknown>> {
  const idParts = splitRole(ids);
  const nameParts = splitRole(names);
  const cidParts = splitRole(cids);
  const chemblParts = splitRole(chembls);
  const count = Math.max(idParts.length, nameParts.length, cidParts.length, 1);
  const items = [];
  for (let index = 0; index < count; index += 1) {
    const cid = cidParts[index]?.trim() || null;
    const chembl = chemblParts[index]?.trim() || null;
    items.push({
      role,
      name: nameParts[index]?.trim() || null,
      compound_id: componentCurie(cid, chembl),
      concentrations,
      unit,
    });
  }
  return items;
}

function concentrationsOf(value: unknown): number[] {
  const textValue = text(value);
  if (!textValue) return [];
  return textValue.split(";").map((part) => Number(part)).filter((part) => Number.isFinite(part));
}

export async function compoundResponse(args: Record<string, unknown>): Promise<ToolResult> {
  const dir = dataDir();
  if (!dir) return unavailable("Set PHI_PGX_DATA_DIR to the processed pharmacogenomics parquet directory.");
  const compounds = tablePath(dir, "compound_xref", false);
  const responses = tablePath(dir, "drug_response", true);
  if (!compounds || !responses) return unavailable("compound_xref or drug_response parquet is not in PHI_PGX_DATA_DIR.");
  const query = text(args.compound);
  if (!query) throw new InputError("compound is required");
  const mode = text(args.combinations) ?? "exclude";
  if (mode !== "exclude" && mode !== "include" && mode !== "only") {
    throw new InputError("combinations must be exclude, include, or only");
  }
  const dataset = text(args.dataset);
  const metric = text(args.metric);
  const { offset, limit } = bounds(args);
  const bare = query.replace(/^pubchem:/i, "");
  const matched = await rows(
    dir,
    `SELECT source, source_compound_id, name, pubchem_cid, chembl_id, inchikey, brd_id, gdsc_drug_id, status, match_method
     FROM read_parquet(?)
     WHERE status IS DISTINCT FROM 'combination'
       AND (
         lower(coalesce(name, '')) = lower(?)
         OR lower(source_compound_id) = lower(?)
         OR lower(coalesce(brd_id, '')) = lower(?)
         OR coalesce(gdsc_drug_id, '') = ?
         OR coalesce(pubchem_cid, '') = ?
         OR lower(coalesce(chembl_id, '')) = lower(?)
         OR lower(coalesce(inchikey, '')) = lower(?)
         OR list_contains(string_split(lower(coalesce(synonyms, '')), ';'), lower(?))
       )`,
    [compounds, query, query, query, query, bare, query, query, query],
  );
  const combinationFile = tablePath(dir, "combination_xref", false);
  const matchedCombinations = combinationFile
    ? await rows(
      dir,
      `SELECT combination_id, name, screen, status,
              anchor_drug_id, library_drug_id, anchor_name, library_name,
              anchor_pubchem_cid, library_pubchem_cid, anchor_chembl_id, library_chembl_id,
              anchor_concentrations, library_concentrations, concentration_unit
       FROM read_parquet(?)
       WHERE lower(coalesce(combination_id, '')) = lower(?)
          OR lower(replace(coalesce(combination_id, ''), 'gdsc.combination:', '')) = lower(?)
          OR lower(coalesce(name, '')) = lower(?)
          OR coalesce(anchor_drug_id, '') = ?
          OR coalesce(library_drug_id, '') = ?
          OR lower(coalesce(anchor_name, '')) = lower(?)
          OR lower(coalesce(library_name, '')) = lower(?)
          OR coalesce(anchor_pubchem_cid, '') = ?
          OR coalesce(library_pubchem_cid, '') = ?
          OR lower(coalesce(anchor_chembl_id, '')) = lower(?)
          OR lower(coalesce(library_chembl_id, '')) = lower(?)`,
      [combinationFile, query, query, query, query, query, query, query, bare, bare, query, query],
    )
    : [];
  if (!matched.length && !matchedCombinations.length) return failure(dataSource, "not_found", `No compound alias matched ${query}`, false);
  const forcedIds = matchedCombinations.filter((row) => combinationIdentity(query, row)).map((row) => text(row.combination_id)).filter((item): item is string => item != null);
  const combinationIds = mode === "exclude"
    ? [...new Set(forcedIds)]
    : [...new Set(matchedCombinations.map((row) => text(row.combination_id)).filter((item): item is string => item != null))];
  const singleIds = mode === "only"
    ? []
    : [...new Set(matched.map((row) => text(row.source_compound_id)).filter((item): item is string => item != null))];
  if (!singleIds.length && !combinationIds.length) return failure(dataSource, "not_found", `No compound alias matched ${query}`, false);
  const described = await rows(dir, "DESCRIBE SELECT * FROM read_parquet(?, hive_partitioning=true)", [responses]);
  const columns = new Set(described.map((row) => text(row.column_name)).filter((item): item is string => item != null));
  const hasClass = columns.has("response_class") && columns.has("combination_id");
  const classSelect = hasClass
    ? "combination_id, response_class, anchor_conc, library_conc"
    : "NULL::VARCHAR AS combination_id, 'single' AS response_class, NULL::DOUBLE AS anchor_conc, NULL::DOUBLE AS library_conc";
  const classWhere = hasClass
    ? `(
         (list_contains(string_split(?, ','), source_compound_id) AND coalesce(response_class, 'single') = 'single' AND combination_id IS NULL)
         OR (list_contains(string_split(?, ','), combination_id) AND coalesce(response_class, 'single') = 'combination')
       )`
    : "list_contains(string_split(?, ','), source_compound_id)";
  const cids = [...new Set(matched.map((row) => text(row.pubchem_cid)).filter((item): item is string => item != null))];
  const where = `WHERE ${classWhere}
                  AND (? IS NULL OR dataset = ?)
                  AND (? IS NULL OR metric = ?)
                  AND coalesce(qc, '') = ''`;
  const classParams: DuckDBValue[] = hasClass
    ? [singleIds.join(",") || "__none__", combinationIds.join(",") || "__none__"]
    : [singleIds.join(",") || "__none__"];
  const params: DuckDBValue[] = [...classParams, dataset, dataset, metric, metric];
  const total = num((await rows(dir, `SELECT count(*)::INTEGER AS n FROM read_parquet(?, hive_partitioning=true) ${where}`, [responses, ...params]))[0]?.n);
  const found = await rows(
    dir,
    `SELECT dataset, source_model_id, model_cvcl, source_compound_id, pubchem_cid, compound_name, metric, value, unit, dose, screen_id, ${classSelect}
     FROM read_parquet(?, hive_partitioning=true) ${where}
     ORDER BY value ASC NULLS LAST
     LIMIT ? OFFSET ?`,
    [responses, ...params, limit, offset],
  );
  const combinationById = new Map(matchedCombinations.map((row) => [text(row.combination_id), row]));
  const items = found.map((row) => {
    const responseClass = text(row.response_class) === "combination" || text(row.combination_id) ? "combination" : "single";
    const combinationId = text(row.combination_id);
    const combination = combinationId ? combinationById.get(combinationId) : undefined;
    const cid = responseClass === "combination" ? null : text(row.pubchem_cid);
    const cvcl = text(row.model_cvcl);
    const datasetName = text(row.dataset) ?? "unknown";
    const compoundId = text(row.source_compound_id) ?? query;
    const metricName = text(row.metric) ?? "value";
    const dose = num(row.dose);
    const anchorConc = num(row.anchor_conc);
    const libraryConc = num(row.library_conc);
    const idParts = [datasetName, cvcl ?? text(row.source_model_id) ?? "model", compoundId, metricName, dose == null ? "na" : String(dose), text(row.screen_id) ?? "na"];
    if (responseClass === "combination") idParts.push(anchorConc == null ? "na" : String(anchorConc), libraryConc == null ? "na" : String(libraryConc));
    return drugResponse({
      id: curie("pgx.response", idParts.join("/")),
      compound_id: cid ? curie("pubchem", cid) : null,
      combination_id: combinationId,
      sample_id: cvcl ? curie("cellosaurus", cvcl) : null,
      dataset_id: curie("pgx.dataset", datasetName),
      metric: metricName,
      value: num(row.value),
      unit: text(row.unit),
      xrefs: [
        compoundId.startsWith("BRD-") ? curie("brd", compoundId) : null,
        /^\d+$/.test(compoundId) ? curie("gdsc.drug", compoundId) : null,
        combinationId,
      ].filter((item): item is string => item != null),
      source_extra: {
        response_class: responseClass,
        source_compound_id: compoundId,
        compound_name: text(row.compound_name),
        dose,
        anchor_conc: anchorConc,
        library_conc: libraryConc,
        screen_id: text(row.screen_id),
        source_model_id: text(row.source_model_id),
        ambiguous_compound: responseClass === "single" && cids.length > 1,
        components: combination
          ? [
            ...roleComponents("anchor", text(combination.anchor_drug_id), text(combination.anchor_name), text(combination.anchor_pubchem_cid), text(combination.anchor_chembl_id), anchorConc == null ? concentrationsOf(combination.anchor_concentrations) : [anchorConc], text(combination.concentration_unit)),
            ...roleComponents("library", text(combination.library_drug_id), text(combination.library_name), text(combination.library_pubchem_cid), text(combination.library_chembl_id), libraryConc == null ? concentrationsOf(combination.library_concentrations) : [libraryConc], text(combination.concentration_unit)),
          ]
          : [],
      },
    });
  });
  const combinationSummary = matchedCombinations
    .filter((row) => combinationIds.includes(text(row.combination_id) ?? ""))
    .map((row) => ({
      combination_id: text(row.combination_id),
      name: text(row.name),
      status: text(row.status),
      components: [
        ...roleComponents("anchor", text(row.anchor_drug_id), text(row.anchor_name), text(row.anchor_pubchem_cid), text(row.anchor_chembl_id), concentrationsOf(row.anchor_concentrations), text(row.concentration_unit)),
        ...roleComponents("library", text(row.library_drug_id), text(row.library_name), text(row.library_pubchem_cid), text(row.library_chembl_id), concentrationsOf(row.library_concentrations), text(row.concentration_unit)),
      ],
    }));
  return makePage("drug_response", items, offset, limit, total, dataSource, {
    summary: combinationIds.length && !singleIds.length
      ? "Combination responses are separate from single-agent responses. Each component keeps its own PubChem or ChEMBL id. Delta metrics are the screened synergy values."
      : cids.length > 1
        ? "This alias matches more than one PubChem CID. Responses are returned for every matched source id and none was chosen automatically."
        : "Drug-response rows ordered by value ascending. Omitting dataset returns every loaded dataset for this compound. Lower GDSC AUC and LN_IC50 mean a stronger response. Combination rows are omitted unless the query is a combination or combinations is include or only.",
    url: "https://www.cancerrxgene.org/",
  }, { pubchem_cids: cids, match_count: matched.length, combinations: combinationSummary });
}

export async function modelSummary(args: Record<string, unknown>): Promise<ToolResult> {
  const resolved = await resolveCellLine({ query: args.model, limit: 5, offset: 0 });
  if (resolved.entity_type === "error") return resolved;
  const dir = dataDir();
  if (!dir) return unavailable("Set PHI_PGX_DATA_DIR to the processed pharmacogenomics parquet directory.");
  const geneArg = text(args.gene);
  const key = geneArg ? geneKey(geneArg) : null;
  const items = [];
  for (const item of resolved.items) {
    if (item.entity_type !== "sample") continue;
    const cvcl = text(item.source_extra.cvcl);
    if (!cvcl) {
      items.push(item);
      continue;
    }
    const counts: Record<string, number> = {};
    for (const table of ["gene_effect", "expression", "mutation", "cnv", "drug_response"] as const) {
      const file = tablePath(dir, table, true);
      counts[table] = file
        ? num((await rows(dir, "SELECT count(*)::INTEGER AS n FROM read_parquet(?, hive_partitioning=true) WHERE model_cvcl = ?", [file, cvcl]))[0]?.n) ?? 0
        : 0;
    }
    const extra: Record<string, unknown> = { ...item.source_extra, counts };
    if (key) {
      extra.gene = await geneSlice(dir, cvcl, key);
    }
    items.push(sample({ ...item, source_extra: extra }));
  }
  return makePage("sample", items, 0, Math.max(items.length, 1), items.length, dataSource, {
    summary: "Per-model counts across the local pharmacogenomics tables. A gene argument adds that gene's dependency, expression, copy number, and mutations.",
    url: usage.url,
  });
}

async function geneSlice(dir: string, cvcl: string, key: { id: string | null; symbol: string | null }): Promise<Record<string, unknown>> {
  const where = "WHERE model_cvcl = ? AND ((? IS NOT NULL AND ncbi_gene_id = ?) OR (? IS NOT NULL AND lower(symbol) = lower(?)))";
  const params = [cvcl, key.id, key.id, key.symbol, key.symbol];
  async function optional(table: string, sql: string): Promise<Record<string, unknown>[]> {
    const file = tablePath(dir, table, true);
    if (!file) return [];
    return rows(dir, sql, [file, ...params]);
  }
  const [effects, expressed, copies, mutations] = await Promise.all([
    optional("gene_effect", `SELECT dataset, value, metric, symbol, ncbi_gene_id FROM read_parquet(?, hive_partitioning=true) ${where} LIMIT 5`),
    optional("expression", `SELECT dataset, value, metric, symbol, ncbi_gene_id FROM read_parquet(?, hive_partitioning=true) ${where} LIMIT 5`),
    optional("cnv", `SELECT dataset, value, metric, symbol, ncbi_gene_id FROM read_parquet(?, hive_partitioning=true) ${where} LIMIT 5`),
    optional("mutation", `SELECT dataset, chrom_norm, pos, ref, alt, consequence, protein_change, symbol, ncbi_gene_id FROM read_parquet(?, hive_partitioning=true) ${where} LIMIT 10`),
  ]);
  const geneId = text(effects[0]?.ncbi_gene_id) ?? text(expressed[0]?.ncbi_gene_id) ?? key.id;
  return {
    gene: geneId ? gene({ id: curie("ncbigene", geneId), symbol: text(effects[0]?.symbol) ?? key.symbol, taxon_id: "NCBITaxon:9606" }) : null,
    expression: expressed.map((row) => expression({
      id: curie("pgx.expression", `${text(row.dataset)}/${cvcl}/${text(row.ncbi_gene_id) ?? "na"}`),
      gene_id: text(row.ncbi_gene_id) ? curie("ncbigene", text(row.ncbi_gene_id)!) : null,
      sample_id: curie("cellosaurus", cvcl),
      value: num(row.value),
      unit: null,
      statistic: text(row.metric),
    })),
    gene_effect: effects.map((row) => ({ dataset: text(row.dataset), value: num(row.value), metric: text(row.metric) })),
    copy_number: copies.map((row) => ({ dataset: text(row.dataset), value: num(row.value), metric: text(row.metric) })),
    variants: mutations.map((row) => variant({
      id: curie("pgx.variant", [text(row.dataset) ?? "mut", cvcl, text(row.chrom_norm) ?? "na", String(num(row.pos) ?? 0), text(row.ref) ?? "N", text(row.alt) ?? "N"].join("/")),
      gene_id: text(row.ncbi_gene_id) ? curie("ncbigene", text(row.ncbi_gene_id)!) : null,
      consequence: text(row.consequence),
      locus: text(row.chrom_norm) ? { genome_build: "GRCh38", chrom: text(row.chrom_norm), start: num(row.pos), end: num(row.pos), strand: null } : null,
      source_extra: { protein_change: text(row.protein_change) },
    })),
  };
}
