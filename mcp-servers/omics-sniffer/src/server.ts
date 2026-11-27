import type { OmicsTool, ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { resolveServingPaths } from "./common.js";
import { getCatalogStats, getDataset, listArtifacts, listSamples, searchOmicsDatasets, usageTerms } from "./query.js";

const pageArgs = {
  offset: z.number().int().nonnegative().optional().describe("Pagination offset (>= 0, default 0)."),
  limit: z.number().int().min(1).max(50).optional().describe("Maximum rows to return (1..50, default 20)."),
};

export const tools: OmicsTool[] = [
  {
    name: "usage_terms",
    title: "Usage terms",
    description: "Return usage terms, repository coverage, and active Parquet snapshot status for Omics Sniffer.",
    inputSchema: z.object({}),
    call: async () => usageTerms(),
  },
  {
    name: "get_catalog_stats",
    title: "Catalog release and source coverage",
    description: "Return active Parquet snapshot release info, total indexed datasets, per-repository dataset counts, and incremental auto-sync freshness.",
    inputSchema: z.object({}),
    call: async () => getCatalogStats(),
  },
  {
    name: "search_omics_datasets",
    title: "Search public omics datasets",
    description: "Search curated public multi-omics datasets across GEO, SRA, ArrayExpress, PRIDE, MassIVE, MetaboLights, CELLxGENE, ENCODE, GTEx, GDC/TCGA, cBioPortal, CPTAC, DepMap, HuBMAP, HPA, ChIP-Atlas, GWAS Catalog, FinnGen, and 15+ other repositories.",
    inputSchema: z.object({
      query: z.string().optional().describe("Free-text keyword matched against dataset title or canonical ID (case-insensitive)."),
      organism: z.string().optional().describe("Organism substring filter, e.g. 'Homo sapiens' or 'Mus musculus'."),
      assay: z.string().optional().describe("Assay/modality substring filter, e.g. 'RNA-seq', 'scRNA-seq', 'ATAC', 'proteomics', 'Visium'."),
      disease: z.string().optional().describe("Disease keyword matched against disease_tags or title, e.g. 'HCC', 'MASLD', 'breast cancer'."),
      tissue: z.string().optional().describe("Tissue/organ substring filter, e.g. 'Liver', 'Brain', 'Lung'."),
      source: z.string().optional().describe("Repository source filter, e.g. 'geo', 'arrayexpress', 'pride', 'sra', 'cellxgene', 'encode', 'gdc_tcga'."),
      min_samples: z.number().int().nonnegative().optional().describe("Minimum sample count threshold."),
      access: z.enum(["open", "public", "controlled", "unknown"]).optional().describe("Access level filter ('open'/'public', 'controlled', or 'unknown')."),
      ...pageArgs,
    }),
    call: async (_client, args) => searchOmicsDatasets(args),
  },
  {
    name: "get_dataset",
    title: "Get dataset metadata, cross-links, and publications",
    description: "Fetch full metadata, cross-repository member accessions (GEO/SRA/BioProject/ArrayExpress/PRIDE), and linked PubMed/DOI publications for one dataset by canonical_id (e.g. 'geo:GSE164760') or bare accession (e.g. 'GSE164760', 'PXD000123').",
    inputSchema: z.object({
      canonical_id: z.string().min(1).describe("Dataset canonical_id (e.g. 'geo:GSE164760', 'ae:E-MTAB-1234') or bare repository accession ('GSE164760')."),
    }),
    call: async (_client, args) => getDataset(args),
  },
  {
    name: "list_samples",
    title: "List dataset sample metadata",
    description: "List paginated sample-level metadata (organism, tissue, disease, condition, sex, age, attributes, SRA run IDs) for a dataset.",
    inputSchema: z.object({
      dataset_id: z.string().min(1).describe("Target dataset canonical_id or bare accession."),
      tissue: z.string().optional().describe("Optional tissue substring filter."),
      disease: z.string().optional().describe("Optional disease substring filter."),
      condition: z.string().optional().describe("Optional experimental condition substring filter."),
      source: z.string().optional().describe("Optional repository source partition hint (e.g. 'geo', 'cbioportal')."),
      ...pageArgs,
    }),
    call: async (_client, args) => listSamples(args),
  },
  {
    name: "list_artifacts",
    title: "List dataset downloadable file artifacts",
    description: "List raw, processed, and supplementary file artifacts (counts, TPM, h5ad, FASTQ/FTP roots, clinical tables) for a dataset. Controlled-access artifacts omit download URLs.",
    inputSchema: z.object({
      dataset_id: z.string().min(1).describe("Target dataset canonical_id or bare accession."),
      kind: z.string().optional().describe("Artifact kind ('raw', 'processed', 'supplementary', 'paper_supplement')."),
      role: z.string().optional().describe("Artifact role substring (e.g. 'counts', 'tpm', 'fpkm', 'ftp_root', 'clinical_table')."),
      format: z.string().optional().describe("File format substring (e.g. 'tsv', 'h5ad', 'mtx', 'fastq', 'bam')."),
      source: z.string().optional().describe("Optional repository source partition hint."),
      ...pageArgs,
    }),
    call: async (_client, args) => listArtifacts(args),
  },
];

export const serverModule: ServerModule = {
  name: "omics-sniffer",
  version: "1.0.0",
  mcpPath: "/mcp",
  tools,
  ready: async () =>
    resolveServingPaths() == null
      ? "PHI_OMICS_SNIFFER_DATA_DIR is not set to an existing omics-sniffer Parquet serving directory"
      : null,
};
