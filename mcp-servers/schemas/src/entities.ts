import { z } from "zod";
import { CurieSchema } from "./ids.js";

export const AccessSchema = z.enum(["open", "controlled", "unknown"]);
export const StrandSchema = z.enum(["+", "-", "."]);
export const SexSchema = z.enum(["female", "male", "unknown"]);
export const ExpressionUnitSchema = z.enum([
  "TPM",
  "log2FC",
  "count",
  "UMI",
  "protein_level",
  "z_score",
]);

export const SourceExtraSchema = z.record(z.string(), z.unknown());

export const LocusSchema = z.object({
  genome_build: z.string().nullable(),
  chrom: z.string().nullable(),
  start: z.number().int().nullable(),
  end: z.number().int().nullable(),
  strand: StrandSchema.nullable(),
}).superRefine((locus, ctx) => {
  if (locus.start != null && locus.end != null && locus.end < locus.start) {
    ctx.addIssue({ code: "custom", message: "end is before start", path: ["end"] });
  }
  if (locus.chrom && !locus.chrom.startsWith("chr")) {
    ctx.addIssue({ code: "custom", message: "chrom must be chr-prefixed", path: ["chrom"] });
  }
});

const base = {
  id: CurieSchema,
  xrefs: z.array(CurieSchema),
  source_extra: SourceExtraSchema,
};

export const GeneSchema = z.object({
  entity_type: z.literal("gene"),
  ...base,
  symbol: z.string().nullable(),
  name: z.string().nullable(),
  taxon_id: CurieSchema.nullable(),
  biotype: z.string().nullable(),
  synonyms: z.array(z.string()),
  locus: LocusSchema.nullable(),
});

export const SampleSchema = z.object({
  entity_type: z.literal("sample"),
  ...base,
  name: z.string().nullable(),
  donor_id: CurieSchema.nullable(),
  tissue_id: CurieSchema.nullable(),
  cell_type_id: CurieSchema.nullable(),
  disease_ids: z.array(CurieSchema),
  assay: z.string().nullable(),
  taxon_id: CurieSchema.nullable(),
});

export const DonorSchema = z.object({
  entity_type: z.literal("donor"),
  ...base,
  sex: SexSchema.nullable(),
  age_years: z.number().nullable(),
  taxon_id: CurieSchema.nullable(),
  disease_ids: z.array(CurieSchema),
});

export const DatasetSchema = z.object({
  entity_type: z.literal("dataset"),
  ...base,
  title: z.string().nullable(),
  description: z.string().nullable(),
  assay: z.string().nullable(),
  taxon_id: CurieSchema.nullable(),
  tissue_ids: z.array(CurieSchema),
  disease_ids: z.array(CurieSchema),
  cell_type_ids: z.array(CurieSchema),
  sample_count: z.number().int().nullable(),
  access: AccessSchema,
  url: z.string().nullable(),
});

export const FileSchema = z.object({
  entity_type: z.literal("file"),
  ...base,
  name: z.string().nullable(),
  format: z.string().nullable(),
  data_type: z.string().nullable(),
  size_bytes: z.number().int().nullable(),
  access: AccessSchema,
  url: z.string().nullable(),
  dataset_id: CurieSchema.nullable(),
}).superRefine((file, ctx) => {
  if (file.access === "controlled" && file.url != null) {
    ctx.addIssue({
      code: "custom",
      message: "controlled-access files expose metadata only; url must be null",
      path: ["url"],
    });
  }
});

export const CellTypeSchema = z.object({
  entity_type: z.literal("cell_type"),
  ...base,
  name: z.string().nullable(),
  description: z.string().nullable(),
  ontology_id: CurieSchema.nullable(),
  parent_ids: z.array(CurieSchema),
});

export const TissueSchema = z.object({
  entity_type: z.literal("tissue"),
  ...base,
  name: z.string().nullable(),
  ontology_id: CurieSchema.nullable(),
});

export const ProteinSchema = z.object({
  entity_type: z.literal("protein"),
  ...base,
  name: z.string().nullable(),
  gene_id: CurieSchema.nullable(),
  gene_symbol: z.string().nullable(),
  taxon_id: CurieSchema.nullable(),
});

export const MetaboliteSchema = z.object({
  entity_type: z.literal("metabolite"),
  ...base,
  name: z.string().nullable(),
  formula: z.string().nullable(),
  inchikey: z.string().nullable(),
});

export const CompoundSchema = z.object({
  entity_type: z.literal("compound"),
  ...base,
  name: z.string().nullable(),
  pubchem_cid: z.string().nullable(),
  chembl_id: z.string().nullable(),
  inchikey: z.string().nullable(),
  synonyms: z.array(z.string()),
});

export const DrugResponseSchema = z.object({
  entity_type: z.literal("drug_response"),
  ...base,
  compound_id: CurieSchema.nullable(),
  sample_id: CurieSchema.nullable(),
  dataset_id: CurieSchema.nullable(),
  metric: z.string().nullable(),
  value: z.number().nullable(),
  unit: z.string().nullable(),
});

export const ExpressionSchema = z.object({
  entity_type: z.literal("expression"),
  ...base,
  gene_id: CurieSchema.nullable(),
  tissue_id: CurieSchema.nullable(),
  cell_type_id: CurieSchema.nullable(),
  sample_id: CurieSchema.nullable(),
  value: z.number().nullable(),
  unit: ExpressionUnitSchema.nullable(),
  statistic: z.string().nullable(),
});

export const VariantSchema = z.object({
  entity_type: z.literal("variant"),
  ...base,
  gene_id: CurieSchema.nullable(),
  locus: LocusSchema.nullable(),
  consequence: z.string().nullable(),
});

export const PUBLICATION_TYPES = [
  "journal-article",
  "preprint",
  "review",
  "book-chapter",
  "proceedings-article",
  "other",
] as const;

export const PublicationTypeSchema = z.enum(PUBLICATION_TYPES);

export const PublicationAuthorSchema = z.object({
  name: z.string().min(1),
  orcid: z.string().nullable(),
});

export const PublicationDatesSchema = z.object({
  published: z.string().nullable(),
  posted: z.string().nullable(),
  updated: z.string().nullable(),
});

export const OaUrlSchema = z.object({
  url: z.string().min(1),
  content_type: z.enum(["html", "pdf", "xml", "landing", "other"]),
  version: z.string().nullable(),
});

export const PublicationSchema = z.object({
  entity_type: z.literal("publication"),
  ...base,
  doi: z.string().nullable(),
  title: z.string().nullable(),
  authors: z.array(PublicationAuthorSchema),
  venue: z.string().nullable(),
  dates: PublicationDatesSchema,
  abstract: z.string().nullable(),
  type: PublicationTypeSchema.nullable(),
  is_oa: z.boolean().nullable(),
  oa_urls: z.array(OaUrlSchema),
  license: z.string().nullable(),
  published_as: CurieSchema.nullable(),
  preprint_of: CurieSchema.nullable(),
}).superRefine((publication, ctx) => {
  if (publication.doi) {
    if (/[A-Z]/.test(publication.doi) || /^https?:/i.test(publication.doi)) {
      ctx.addIssue({ code: "custom", message: "doi must be a lowercased bare DOI", path: ["doi"] });
    }
    if (publication.id !== `doi:${publication.doi}`) {
      ctx.addIssue({
        code: "custom",
        message: "id must be doi:<lowercased DOI> when a DOI is present",
        path: ["id"],
      });
    }
  } else if (!publication.id.startsWith("pmid:") && !publication.id.startsWith("arxiv:")) {
    ctx.addIssue({
      code: "custom",
      message: "without a DOI, id must be pmid: or arxiv:",
      path: ["id"],
    });
  }
});

export const EntitySchema = z.discriminatedUnion("entity_type", [
  GeneSchema,
  SampleSchema,
  DonorSchema,
  DatasetSchema,
  FileSchema,
  CellTypeSchema,
  TissueSchema,
  ProteinSchema,
  MetaboliteSchema,
  CompoundSchema,
  DrugResponseSchema,
  ExpressionSchema,
  VariantSchema,
  PublicationSchema,
]);

export const ENTITY_TYPES = [
  "gene",
  "sample",
  "donor",
  "dataset",
  "file",
  "cell_type",
  "tissue",
  "protein",
  "metabolite",
  "compound",
  "drug_response",
  "expression",
  "variant",
  "publication",
] as const;

export const EntityTypeSchema = z.enum(ENTITY_TYPES);
