import { curie, ensemblGeneCurie, humanChrom, ontologyCurie, strandOf } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  expression,
  failure,
  gene,
  guard,
  locusOf,
  makePage,
  numberOrNull,
  paging,
  sourceOf,
  text,
  tissue,
  variant,
  type JsonClient,
  type OmicsTool,
} from "phi-omics-runtime";
import { loadConfig } from "phi-omics-runtime";
import { z } from "zod";

const source = () => sourceOf("gtex", "GTEx Portal", "https://gtexportal.org/", loadConfig().gtexDataset);
const usage = {
  summary: "GTEx gene expression and eQTL summaries are open aggregate data. Cite GTEx and do not treat portal values as clinical measurements. Coordinates are 1-based and inclusive on the genome build named in each record.",
  url: "https://gtexportal.org/home/datasets",
};
const pageArgs = {
  offset: z.number().int().nonnegative().optional().describe("Items to skip. Default 0."),
  limit: z.number().int().min(1).max(50).optional().describe("Page size. Default 20, maximum 50."),
};

async function gtex(client: JsonClient, path: string, params: Record<string, string | number>): Promise<{ rows: Record<string, unknown>[]; total: number }> {
  const url = new URL(`https://gtexportal.org/api/v2/${path}`);
  for (const [key, value] of Object.entries(params)) url.searchParams.set(key, String(value));
  const body = asRecord(await client.getJson(url.toString()));
  const rows = asArray(body.data).map(asRecord);
  const total = numberOrNull(asRecord(body.paging_info).totalNumberOfItems) ?? rows.length;
  return { rows, total };
}

async function window(client: JsonClient, path: string, params: Record<string, string>, offset: number, limit: number) {
  const startPage = Math.floor(offset / limit);
  const skip = offset - startPage * limit;
  const first = await gtex(client, path, { ...params, page: startPage, itemsPerPage: limit });
  let rows = first.rows.slice(skip, skip + limit);
  if (skip > 0 && rows.length < limit && first.rows.length === limit) {
    const second = await gtex(client, path, { ...params, page: startPage + 1, itemsPerPage: limit });
    rows = rows.concat(second.rows.slice(0, limit - rows.length));
  }
  return { rows, total: first.total };
}

async function resolveGencode(client: JsonClient, geneQuery: string): Promise<string | null> {
  const config = loadConfig();
  const versioned = geneQuery.match(/ENSG\d+\.\d+/i);
  if (versioned) return versioned[0];
  const found = await gtex(client, "reference/gene", {
    geneId: geneQuery.replace(/^ensembl:/i, ""),
    gencodeVersion: config.gtexGencode,
    genomeBuild: config.gtexGenome,
    page: 0,
    itemsPerPage: 1,
  });
  return text(found.rows[0]?.gencodeId);
}

function mapGene(row: Record<string, unknown>) {
  const gencode = text(row.gencodeId);
  const ensembl = gencode ? ensemblGeneCurie(gencode) : null;
  const entrez = numberOrNull(row.entrezGeneId);
  const symbol = text(row.geneSymbol);
  return gene({
    id: ensembl?.id ?? curie("gtex.gene", symbol ?? "unknown"),
    symbol,
    name: text(row.description),
    taxon_id: "NCBITaxon:9606",
    biotype: text(row.geneType),
    xrefs: entrez == null ? [] : [`ncbigene:${entrez}`],
    locus: locusOf({
      genomeBuild: text(row.genomeBuild),
      chrom: humanChrom(text(row.chromosome)),
      start: numberOrNull(row.start),
      end: numberOrNull(row.end),
      strand: strandOf(row.strand),
    }),
    source_extra: { gencode_id: gencode, gencode_version: text(row.gencodeVersion), dataset_id: loadConfig().gtexDataset },
  });
}

export const gtexTools: OmicsTool[] = [
  {
    name: "gtex_search_genes",
    title: "Search GTEx genes",
    description: "Search GTEx reference genes by symbol or Ensembl id. Returns the shared gene entity.",
    inputSchema: z.object({ gene: z.string().describe("Gene symbol or Ensembl id, for example TP53."), ...pageArgs }),
    call: (client, args) => guard(source(), async () => {
      const { offset, limit } = paging(args);
      const config = loadConfig();
      const geneQuery = text(args.gene);
      if (!geneQuery) return failure(source(), "invalid_input", "gene is required", false);
      const { rows, total } = await window(client, "reference/gene", {
        geneId: geneQuery.replace(/^ensembl:/i, ""),
        gencodeVersion: config.gtexGencode,
        genomeBuild: config.gtexGenome,
      }, offset, limit);
      return makePage("gene", rows.map(mapGene), offset, limit, total, source(), usage);
    }),
  },
  {
    name: "gtex_list_tissues",
    title: "List GTEx tissues",
    description: "List GTEx tissue sites. Returns the shared tissue entity; ontology_id is UBERON when GTEx provides one.",
    inputSchema: z.object({ ...pageArgs }),
    call: (client, args) => guard(source(), async () => {
      const { offset, limit } = paging(args);
      const { rows } = await gtex(client, "dataset/tissueSiteDetail", { datasetId: loadConfig().gtexDataset });
      const page = rows.slice(offset, offset + limit).map((row) => tissue({
        id: curie("gtex.tissue", text(row.tissueSiteDetailId) ?? "unknown"),
        name: text(row.tissueSiteDetail) ?? text(row.tissueSite),
        ontology_id: ontologyCurie(text(row.ontologyId)),
        source_extra: { tissue_site: row.tissueSite, dataset_id: row.datasetId, rna_seq_sample_count: asRecord(row.rnaSeqSampleSummary).totalCount ?? null },
      }));
      return makePage("tissue", page, offset, limit, rows.length, source(), usage);
    }),
  },
  {
    name: "gtex_gene_expression",
    title: "GTEx median expression",
    description: "Median GTEx expression for one gene across tissues, in TPM.",
    inputSchema: z.object({
      gene: z.string().describe("Gene symbol, Ensembl id, or versioned GENCODE id."),
      tissue: z.string().optional().describe("Optional GTEx tissueSiteDetailId, for example Whole_Blood."),
      ...pageArgs,
    }),
    call: (client, args) => guard(source(), async () => {
      const { offset, limit } = paging(args);
      const geneQuery = text(args.gene);
      if (!geneQuery) return failure(source(), "invalid_input", "gene is required", false);
      const gencode = await resolveGencode(client, geneQuery);
      if (!gencode) return failure(source(), "not_found", `no GTEx gene matches ${geneQuery}`, false);
      const params: Record<string, string> = { gencodeId: gencode, datasetId: loadConfig().gtexDataset };
      const tissueId = text(args.tissue);
      if (tissueId) params.tissueSiteDetailId = tissueId;
      const { rows, total } = await window(client, "expression/medianGeneExpression", params, offset, limit);
      const ensembl = ensemblGeneCurie(gencode);
      const items = rows.map((row) => {
        const tissueName = text(row.tissueSiteDetailId) ?? "unknown";
        return expression({
          id: curie("gtex.expression", `${ensembl.id.slice("ensembl:".length)}:${tissueName}:median`),
          gene_id: ensembl.id,
          tissue_id: curie("gtex.tissue", tissueName),
          value: numberOrNull(row.median),
          unit: "TPM",
          statistic: "median",
          xrefs: ontologyCurie(text(row.ontologyId)) ? [ontologyCurie(text(row.ontologyId))!] : [],
          source_extra: { gencode_id: gencode, dataset_id: text(row.datasetId), unit_note: text(row.unit) },
        });
      });
      return makePage("expression", items, offset, limit, total, source(), usage);
    }),
  },
  {
    name: "gtex_list_egenes",
    title: "List GTEx eGenes",
    description: "List GTEx eGenes for a tissue. The value is the log2 allelic fold change.",
    inputSchema: z.object({
      tissue: z.string().describe("GTEx tissueSiteDetailId, for example Whole_Blood."),
      ...pageArgs,
    }),
    call: (client, args) => guard(source(), async () => {
      const { offset, limit } = paging(args);
      const tissueId = text(args.tissue);
      if (!tissueId) return failure(source(), "invalid_input", "tissue is required", false);
      const { rows, total } = await window(client, "association/egene", { tissueSiteDetailId: tissueId, datasetId: loadConfig().gtexDataset }, offset, limit);
      const items = rows.map((row) => {
        const gencode = text(row.gencodeId);
        const ensembl = gencode ? ensemblGeneCurie(gencode) : null;
        const tissueName = text(row.tissueSiteDetailId) ?? tissueId;
        return expression({
          id: curie("gtex.expression", `${ensembl?.id.slice("ensembl:".length) ?? "unknown"}:${tissueName}:egene`),
          gene_id: ensembl?.id ?? null,
          tissue_id: curie("gtex.tissue", tissueName),
          value: numberOrNull(row.log2AllelicFoldChange),
          unit: "log2FC",
          statistic: "log2_allelic_fold_change",
          source_extra: { gencode_id: gencode, p_value: row.pValue ?? null, q_value: row.qValue ?? null, dataset_id: text(row.datasetId), symbol: text(row.geneSymbol) },
        });
      });
      return makePage("expression", items, offset, limit, total, source(), usage);
    }),
  },
  {
    name: "gtex_list_eqtls",
    title: "List GTEx single-tissue eQTLs",
    description: "List single-tissue eQTL variants for one GENCODE gene and tissue. Returns the shared variant entity; effect size and p-value are in source_extra.",
    inputSchema: z.object({
      gene: z.string().describe("Versioned GENCODE id or a symbol/Ensembl id that GTEx can resolve."),
      tissue: z.string().describe("GTEx tissueSiteDetailId."),
      ...pageArgs,
    }),
    call: (client, args) => guard(source(), async () => {
      const { offset, limit } = paging(args);
      const geneQuery = text(args.gene);
      const tissueId = text(args.tissue);
      if (!geneQuery || !tissueId) return failure(source(), "invalid_input", "gene and tissue are required", false);
      const gencode = await resolveGencode(client, geneQuery);
      if (!gencode) return failure(source(), "not_found", `no GTEx gene matches ${geneQuery}`, false);
      const { rows, total } = await window(client, "association/singleTissueEqtl", {
        gencodeId: gencode,
        tissueSiteDetailId: tissueId,
        datasetId: loadConfig().gtexDataset,
      }, offset, limit);
      const ensembl = ensemblGeneCurie(gencode);
      const items = rows.map((row) => {
        const variantId = text(row.variantId) ?? "unknown";
        const snp = text(row.snpId);
        return variant({
          id: curie("gtex.variant", variantId),
          gene_id: ensembl.id,
          xrefs: snp ? [`dbsnp:${snp}`] : [],
          locus: locusOf({ genomeBuild: "GRCh38", chrom: humanChrom(text(row.chromosome)), start: numberOrNull(row.pos), end: numberOrNull(row.pos) }),
          source_extra: { nes: row.nes ?? null, p_value: row.pValue ?? null, tissue_id: tissueId, gencode_id: gencode, dataset_id: text(row.datasetId) },
        });
      });
      return makePage("variant", items, offset, limit, total, source(), usage);
    }),
  },
];
