import { UpstreamError } from "./http.js";
import {
  CellTypeSchema,
  CompoundSchema,
  DatasetSchema,
  DrugCombinationSchema,
  DrugResponseSchema,
  DonorSchema,
  ErrorObjectSchema,
  ExpressionSchema,
  FileSchema,
  GeneSchema,
  MetaboliteSchema,
  PageSchema,
  ProteinSchema,
  PublicationSchema,
  SampleSchema,
  TissueSchema,
  VariantSchema,
  type Access,
  type CellType,
  type Compound,
  type Dataset,
  type DrugCombination,
  type DrugResponse,
  type Donor,
  type Entity,
  type EntityType,
  type ErrorObject,
  type Expression,
  type ExpressionUnit,
  type Gene,
  type Locus,
  type Metabolite,
  type OmicsFile,
  type Page,
  type Protein,
  type Publication,
  type Sample,
  type Source,
  type Tissue,
  type ToolResult,
  type Usage,
  type Variant,
} from "phi-omics-schemas";

export function sourceOf(id: string, name: string, homepage: string, version: string | null = null): Source {
  return { id, name, version, homepage, retrieved_at: new Date().toISOString() };
}

export function pageInfo(offset: number, limit: number, returned: number, total: number | null) {
  const more = total == null ? returned === limit : offset + returned < total;
  return {
    offset,
    limit,
    returned,
    total,
    next_offset: more ? offset + returned : null,
  };
}

export function makePage(itemType: EntityType, items: Entity[], offset: number, limit: number, total: number | null, source: Source, usage: Usage, sourceExtra: Record<string, unknown> = {}): Page {
  return PageSchema.parse({
    entity_type: "page",
    item_type: itemType,
    items,
    page: pageInfo(offset, limit, items.length, total),
    source,
    usage,
    source_extra: sourceExtra,
  });
}

export function failure(source: Source, code: ErrorObject["code"], message: string, retryable: boolean, sourceExtra: Record<string, unknown> = {}): ErrorObject {
  return ErrorObjectSchema.parse({
    entity_type: "error",
    code,
    message,
    source,
    retryable,
    source_extra: sourceExtra,
  });
}

export async function guard(source: Source, fn: () => Promise<ToolResult>): Promise<ToolResult> {
  try {
    return await fn();
  } catch (error) {
    if (error instanceof UpstreamError) {
      const code = error.status === 404 ? "not_found" : error.status == null ? "unavailable" : "upstream_http";
      return failure(source, code, error.message, error.retryable, { status: error.status });
    }
    throw error;
  }
}

export function usagePage(source: Source, summary: string, url: string | null = null): Page {
  return makePage("dataset", [], 0, 1, 0, source, { summary, url });
}

export function paging(args: { offset?: number; limit?: number }): { offset: number; limit: number } {
  return { offset: args.offset ?? 0, limit: args.limit ?? 20 };
}

export function sliceList<T>(items: T[], offset: number, limit: number): { rows: T[]; total: number } {
  return { rows: items.slice(offset, offset + limit), total: items.length };
}

export function locusOf(fields: {
  genomeBuild?: string | null;
  chrom?: string | null;
  start?: number | null;
  end?: number | null;
  strand?: Locus["strand"];
}): Locus | null {
  if (!fields.chrom && fields.start == null && fields.end == null) return null;
  return {
    genome_build: fields.genomeBuild ?? null,
    chrom: fields.chrom ?? null,
    start: fields.start ?? null,
    end: fields.end ?? null,
    strand: fields.strand ?? null,
  };
}

export function gene(fields: Partial<Gene> & Pick<Gene, "id">): Gene {
  return GeneSchema.parse({
    entity_type: "gene",
    xrefs: [],
    symbol: null,
    name: null,
    taxon_id: null,
    biotype: null,
    synonyms: [],
    locus: null,
    source_extra: {},
    ...fields,
  });
}

export function sample(fields: Partial<Sample> & Pick<Sample, "id">): Sample {
  return SampleSchema.parse({
    entity_type: "sample",
    xrefs: [],
    name: null,
    donor_id: null,
    tissue_id: null,
    cell_type_id: null,
    disease_ids: [],
    assay: null,
    taxon_id: null,
    source_extra: {},
    ...fields,
  });
}

export function donor(fields: Partial<Donor> & Pick<Donor, "id">): Donor {
  return DonorSchema.parse({
    entity_type: "donor",
    xrefs: [],
    sex: null,
    age_years: null,
    taxon_id: null,
    disease_ids: [],
    source_extra: {},
    ...fields,
  });
}

export function dataset(fields: Partial<Dataset> & Pick<Dataset, "id">): Dataset {
  const access: Access = fields.access ?? "unknown";
  return DatasetSchema.parse({
    entity_type: "dataset",
    xrefs: [],
    title: null,
    description: null,
    assay: null,
    taxon_id: null,
    tissue_ids: [],
    disease_ids: [],
    cell_type_ids: [],
    sample_count: null,
    source_extra: {},
    ...fields,
    access,
    url: access === "controlled" ? null : fields.url ?? null,
  });
}

export function file(fields: Partial<OmicsFile> & Pick<OmicsFile, "id">): OmicsFile {
  const access: Access = fields.access ?? "unknown";
  return FileSchema.parse({
    entity_type: "file",
    xrefs: [],
    name: null,
    format: null,
    data_type: null,
    size_bytes: null,
    dataset_id: null,
    source_extra: {},
    ...fields,
    access,
    url: access === "controlled" ? null : fields.url ?? null,
  });
}

export function cellType(fields: Partial<CellType> & Pick<CellType, "id">): CellType {
  return CellTypeSchema.parse({
    entity_type: "cell_type",
    xrefs: [],
    name: null,
    description: null,
    ontology_id: null,
    parent_ids: [],
    source_extra: {},
    ...fields,
  });
}

export function tissue(fields: Partial<Tissue> & Pick<Tissue, "id">): Tissue {
  return TissueSchema.parse({
    entity_type: "tissue",
    xrefs: [],
    name: null,
    ontology_id: null,
    source_extra: {},
    ...fields,
  });
}

export function protein(fields: Partial<Protein> & Pick<Protein, "id">): Protein {
  return ProteinSchema.parse({
    entity_type: "protein",
    xrefs: [],
    name: null,
    gene_id: null,
    gene_symbol: null,
    taxon_id: null,
    source_extra: {},
    ...fields,
  });
}

export function compound(fields: Partial<Compound> & Pick<Compound, "id">): Compound {
  return CompoundSchema.parse({
    entity_type: "compound",
    xrefs: [],
    name: null,
    pubchem_cid: null,
    chembl_id: null,
    inchikey: null,
    synonyms: [],
    source_extra: {},
    ...fields,
  });
}

export function drugCombination(fields: Partial<DrugCombination> & Pick<DrugCombination, "id">): DrugCombination {
  return DrugCombinationSchema.parse({
    entity_type: "drug_combination",
    xrefs: [],
    name: null,
    components: [],
    source_extra: {},
    ...fields,
  });
}

export function drugResponse(fields: Partial<DrugResponse> & Pick<DrugResponse, "id">): DrugResponse {
  return DrugResponseSchema.parse({
    entity_type: "drug_response",
    xrefs: [],
    compound_id: null,
    combination_id: null,
    sample_id: null,
    dataset_id: null,
    metric: null,
    value: null,
    unit: null,
    source_extra: {},
    ...fields,
  });
}

export function metabolite(fields: Partial<Metabolite> & Pick<Metabolite, "id">): Metabolite {
  return MetaboliteSchema.parse({
    entity_type: "metabolite",
    xrefs: [],
    name: null,
    formula: null,
    inchikey: null,
    source_extra: {},
    ...fields,
  });
}

export function expression(fields: Partial<Expression> & Pick<Expression, "id"> & { unit?: ExpressionUnit | null }): Expression {
  return ExpressionSchema.parse({
    entity_type: "expression",
    xrefs: [],
    gene_id: null,
    tissue_id: null,
    cell_type_id: null,
    sample_id: null,
    value: null,
    unit: null,
    statistic: null,
    source_extra: {},
    ...fields,
  });
}

export function publication(fields: Partial<Publication> & Pick<Publication, "id">): Publication {
  return PublicationSchema.parse({
    entity_type: "publication",
    xrefs: [],
    doi: null,
    title: null,
    authors: [],
    venue: null,
    abstract: null,
    type: null,
    is_oa: null,
    oa_urls: [],
    license: null,
    published_as: null,
    preprint_of: null,
    source_extra: {},
    ...fields,
    dates: { published: null, posted: null, updated: null, ...fields.dates },
  });
}

export function variant(fields: Partial<Variant> & Pick<Variant, "id">): Variant {
  return VariantSchema.parse({
    entity_type: "variant",
    xrefs: [],
    gene_id: null,
    locus: null,
    consequence: null,
    source_extra: {},
    ...fields,
  });
}

export function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === "object" && value !== null && !Array.isArray(value);
}

export function asRecord(value: unknown): Record<string, unknown> {
  return isRecord(value) ? value : {};
}

export function asArray(value: unknown): unknown[] {
  return Array.isArray(value) ? value : [];
}

export function text(value: unknown): string | null {
  if (typeof value === "string" && value.trim()) return value;
  if (typeof value === "number" && Number.isFinite(value)) return String(value);
  return null;
}

export function numberOrNull(value: unknown): number | null {
  if (typeof value === "number" && Number.isFinite(value)) return value;
  if (typeof value === "string" && value.trim() && Number.isFinite(Number(value))) return Number(value);
  return null;
}

export function labels(value: unknown): { ids: string[]; names: string[] } {
  const ids: string[] = [];
  const names: string[] = [];
  for (const item of asArray(value)) {
    const record = asRecord(item);
    const name = text(record.label) ?? text(record.term_name) ?? text(item);
    const id = text(record.ontology_term_id) ?? text(record.term_id);
    if (name) names.push(name);
    if (id) ids.push(id);
  }
  return { ids, names };
}

export type { ToolResult };
