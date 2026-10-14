import { curie, ensemblGeneCurie, humanChrom, strandOf } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  dataset,
  donor,
  failure,
  file,
  gene,
  guard,
  locusOf,
  makePage,
  numberOrNull,
  paging,
  sample,
  sourceOf,
  text,
  variant,
  type JsonClient,
  type OmicsTool,
} from "phi-omics-runtime";
import { z } from "zod";

const gdc = sourceOf("gdc", "NCI GDC", "https://gdc.cancer.gov/", "GRCh38");
const gdcUsage = {
  summary: "GDC case, gene, mutation, and file records are metadata. access=controlled files have a null url and are not downloaded, even if a token is available. Open-file URLs point at the GDC data endpoint but this server does not retrieve bytes.",
  url: "https://gdc.cancer.gov/access-data/data-access-policies",
};
const cbio = sourceOf("cbioportal", "cBioPortal", "https://www.cbioportal.org/");
const cbioUsage = {
  summary: "These tools use the public cBioPortal REST API for open study metadata. They do not replace the separate cBioPortal OAuth MCP connector and do not fetch controlled studies.",
  url: "https://www.cbioportal.org/api/swagger-ui/index.html",
};
const pdc = sourceOf("pdc", "NCI Proteomic Data Commons", "https://pdc.cancer.gov/");
const pdcUsage = {
  summary: "PDC study metadata comes from the public GraphQL API. Controlled proteomic files are not downloaded.",
  url: "https://pdc.cancer.gov/",
};
const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

function gdcUrl(path: string, params: Record<string, string | number>): string {
  const url = new URL(path, "https://api.gdc.cancer.gov");
  for (const [key, value] of Object.entries(params)) url.searchParams.set(key, String(value));
  return url.toString();
}

async function gdcHits(client: JsonClient, path: string, params: Record<string, string | number>) {
  const body = asRecord(await client.getJson(gdcUrl(path, params)));
  const data = asRecord(body.data);
  const pagination = asRecord(data.pagination);
  return {
    rows: asArray(data.hits).map(asRecord),
    total: numberOrNull(pagination.total),
    record: data,
  };
}

function accessOf(value: unknown): "open" | "controlled" | "unknown" {
  const level = text(value)?.toLowerCase();
  if (level === "open" || level === "public") return "open";
  if (level === "controlled" || level === "restricted" || level === "private") return "controlled";
  return "unknown";
}

function mapCase(row: Record<string, unknown>) {
  const id = text(row.case_id) ?? text(row.id) ?? "unknown";
  const demographic = asRecord(row.demographic);
  return donor({
    id: curie("gdc.case", id),
    age_years: numberOrNull(demographic.age_at_index),
    taxon_id: "NCBITaxon:9606",
    source_extra: { submitter_id: text(row.submitter_id), disease_type: text(row.disease_type), primary_site: text(row.primary_site) },
  });
}

function mapGdcFile(row: Record<string, unknown>) {
  const id = text(row.file_id) ?? text(row.id) ?? "unknown";
  const access = accessOf(row.access);
  const caseId = text(asRecord(asArray(row.cases)[0]).case_id);
  return file({
    id: curie("gdc.file", id),
    name: text(row.file_name),
    format: text(row.data_format),
    data_type: text(row.data_type),
    size_bytes: numberOrNull(row.file_size),
    access,
    url: access === "open" ? `https://api.gdc.cancer.gov/data/${id}` : null,
    source_extra: { case_id: caseId, submitter_id: text(asRecord(asArray(row.cases)[0]).submitter_id) },
  });
}

function mapGdcGene(row: Record<string, unknown>) {
  const geneId = text(row.gene_id);
  const ensembl = geneId ? ensemblGeneCurie(geneId) : null;
  return gene({
    id: ensembl?.id ?? curie("gdc.gene", text(row.symbol) ?? "unknown"),
    symbol: text(row.symbol),
    name: text(row.name),
    taxon_id: "NCBITaxon:9606",
    biotype: text(row.biotype),
    synonyms: asArray(row.synonyms).map(text).filter((item): item is string => item != null),
    locus: locusOf({
      genomeBuild: "GRCh38",
      chrom: humanChrom(text(row.gene_chromosome)),
      start: numberOrNull(row.gene_start),
      end: numberOrNull(row.gene_end),
      strand: strandOf(row.gene_strand),
    }),
    source_extra: {},
  });
}

async function cbioGet(client: JsonClient, path: string, params: Record<string, string | number> = {}) {
  const url = new URL(path.replace(/^\//, ""), "https://www.cbioportal.org/api/");
  for (const [key, value] of Object.entries(params)) url.searchParams.set(key, String(value));
  return client.getJson(url.toString());
}

async function pdcQuery(client: JsonClient, query: string, variables: Record<string, unknown>) {
  const body = asRecord(await client.postJson("https://proteomic.datacommons.cancer.gov/graphql", { query, variables }));
  if (asArray(body.errors).length > 0) {
    return failure(pdc, "upstream_http", text(asRecord(asArray(body.errors)[0]).message) ?? "PDC GraphQL error", false);
  }
  return body;
}

const studyFields = "pdc_study_id submitter_id_name disease_type analytical_fraction program_name";

function mapPdc(row: Record<string, unknown>) {
  const id = text(row.pdc_study_id) ?? "unknown";
  return dataset({
    id: curie("pdc.study", id),
    title: text(row.submitter_id_name),
    assay: text(row.analytical_fraction),
    taxon_id: "NCBITaxon:9606",
    access: "open",
    url: `https://pdc.cancer.gov/study/${id}`,
    source_extra: { disease_type: text(row.disease_type), program_name: text(row.program_name) },
  });
}

export const cancerTools: OmicsTool[] = [
  {
    name: "gdc_search_cases",
    title: "Search GDC cases",
    description: "Search GDC cases as donor entities. A query matches case_id or submitter_id exactly.",
    inputSchema: z.object({ query: z.string().optional(), ...pageArgs }),
    call: (client, args) => guard(gdc, async () => {
      const { offset, limit } = paging(args);
      const params: Record<string, string | number> = { from: offset, size: limit, fields: "case_id,submitter_id,disease_type,primary_site,demographic.age_at_index" };
      const query = text(args.query);
      if (query) {
        params.filters = JSON.stringify({ op: "or", content: ["case_id", "submitter_id"].map((field) => ({ op: "=", content: { field, value: query } })) });
      }
      const { rows, total } = await gdcHits(client, "/cases", params);
      return makePage("donor", rows.map(mapCase), offset, limit, total, gdc, gdcUsage);
    }),
  },
  {
    name: "gdc_get_case",
    title: "Get a GDC case",
    description: "Fetch one GDC case by case UUID or submitter id.",
    inputSchema: z.object({ case_id: z.string() }),
    call: (client, args) => guard(gdc, async () => {
      const id = text(args.case_id);
      if (!id) return failure(gdc, "invalid_input", "case_id is required", false);
      const field = id.includes("-") && id.length > 20 ? "case_id" : "submitter_id";
      const { rows } = await gdcHits(client, "/cases", {
        size: 1,
        fields: "case_id,submitter_id,disease_type,primary_site,demographic.age_at_index",
        filters: JSON.stringify({ op: "=", content: { field, value: id } }),
      });
      if (!rows[0]) return failure(gdc, "not_found", `case ${id} was not found`, false);
      return makePage("donor", [mapCase(rows[0])], 0, 1, 1, gdc, gdcUsage);
    }),
  },
  {
    name: "gdc_search_files",
    title: "Search GDC file metadata",
    description: "Search GDC file metadata. Controlled files are returned with url null.",
    inputSchema: z.object({ case_id: z.string().optional(), access: z.enum(["open", "controlled"]).optional(), ...pageArgs }),
    call: (client, args) => guard(gdc, async () => {
      const { offset, limit } = paging(args);
      const content = [];
      if (text(args.case_id)) content.push({ op: "=", content: { field: "cases.case_id", value: text(args.case_id) } });
      if (text(args.access)) content.push({ op: "=", content: { field: "access", value: text(args.access) } });
      const params: Record<string, string | number> = { from: offset, size: limit, fields: "file_id,file_name,data_type,data_format,access,file_size,cases.case_id,cases.submitter_id" };
      if (content.length === 1) params.filters = JSON.stringify(content[0]);
      if (content.length > 1) params.filters = JSON.stringify({ op: "and", content });
      const { rows, total } = await gdcHits(client, "/files", params);
      return makePage("file", rows.map(mapGdcFile), offset, limit, total, gdc, gdcUsage);
    }),
  },
  {
    name: "gdc_get_gene",
    title: "Get a GDC gene",
    description: "Fetch one GDC gene by Ensembl id or symbol.",
    inputSchema: z.object({ gene: z.string() }),
    call: (client, args) => guard(gdc, async () => {
      const query = text(args.gene);
      if (!query) return failure(gdc, "invalid_input", "gene is required", false);
      const bare = query.replace(/^ensembl:/i, "").replace(/\.\d+$/, "");
      if (/^ENSG\d+$/i.test(bare)) {
        const body = asRecord(await client.getJson(gdcUrl(`/genes/${bare}`, { fields: "gene_id,symbol,name,biotype,synonyms,gene_chromosome,gene_start,gene_end,gene_strand" })));
        const row = asRecord(body.data);
        if (!text(row.gene_id)) return failure(gdc, "not_found", `gene ${query} was not found`, false);
        return makePage("gene", [mapGdcGene(row)], 0, 1, 1, gdc, gdcUsage);
      }
      const { rows } = await gdcHits(client, "/genes", {
        size: 1,
        fields: "gene_id,symbol,name,biotype,synonyms,gene_chromosome,gene_start,gene_end,gene_strand",
        filters: JSON.stringify({ op: "=", content: { field: "symbol", value: query } }),
      });
      if (!rows[0]) return failure(gdc, "not_found", `gene ${query} was not found`, false);
      return makePage("gene", [mapGdcGene(rows[0])], 0, 1, 1, gdc, gdcUsage);
    }),
  },
  {
    name: "gdc_search_mutations",
    title: "Search GDC simple somatic mutations",
    description: "Search GDC SSM metadata as variant entities. An optional gene limits results to an Ensembl gene id.",
    inputSchema: z.object({ gene: z.string().optional(), ...pageArgs }),
    call: (client, args) => guard(gdc, async () => {
      const { offset, limit } = paging(args);
      const params: Record<string, string | number> = { from: offset, size: limit, fields: "ssm_id,chromosome,start_position,end_position,genomic_dna_change,mutation_type,mutation_subtype,ncbi_build" };
      const geneQuery = text(args.gene)?.replace(/^ensembl:/i, "").replace(/\.\d+$/, "");
      if (geneQuery) params.filters = JSON.stringify({ op: "=", content: { field: "consequence.transcript.gene.gene_id", value: geneQuery } });
      const { rows, total } = await gdcHits(client, "/ssms", params);
      const items = rows.map((row) => variant({
        id: curie("gdc.ssm", text(row.ssm_id) ?? text(row.id) ?? "unknown"),
        gene_id: geneQuery && /^ENSG\d+$/i.test(geneQuery) ? ensemblGeneCurie(geneQuery).id : null,
        locus: locusOf({ genomeBuild: text(row.ncbi_build) ?? "GRCh38", chrom: humanChrom(text(row.chromosome)), start: numberOrNull(row.start_position), end: numberOrNull(row.end_position) }),
        consequence: text(row.mutation_subtype) ?? text(row.mutation_type),
        source_extra: { genomic_dna_change: text(row.genomic_dna_change), mutation_type: text(row.mutation_type) },
      }));
      return makePage("variant", items, offset, limit, total, gdc, gdcUsage);
    }),
  },
  {
    name: "cbioportal_list_studies",
    title: "List cBioPortal studies",
    description: "List public cBioPortal studies as dataset entities.",
    inputSchema: z.object({ ...pageArgs }),
    call: (client, args) => guard(cbio, async () => {
      const { offset, limit } = paging(args);
      const start = Math.floor(offset / limit);
      const rows = asArray(await cbioGet(client, "/studies", { pageSize: limit, pageNumber: start })).map(asRecord);
      const items = rows.slice(offset - start * limit, offset - start * limit + limit).map((row) => {
        const access = row.publicStudy === false || row.readPermission === false ? "controlled" : "open";
        const id = text(row.studyId) ?? "unknown";
        return dataset({
          id: curie("cbioportal.study", id),
          title: text(row.name),
          description: text(row.description),
          sample_count: numberOrNull(row.allSampleCount),
          access,
          url: access === "open" ? `https://www.cbioportal.org/study/summary?id=${encodeURIComponent(id)}` : null,
          source_extra: { cancer_type_id: text(row.cancerTypeId), reference_genome: text(row.referenceGenome), pmid: text(row.pmid) },
        });
      });
      return makePage("dataset", items, offset, limit, rows.length < limit ? start * limit + rows.length : null, cbio, cbioUsage);
    }),
  },
  {
    name: "cbioportal_get_study",
    title: "Get a cBioPortal study",
    description: "Fetch one cBioPortal study by studyId.",
    inputSchema: z.object({ study_id: z.string() }),
    call: (client, args) => guard(cbio, async () => {
      const id = text(args.study_id);
      if (!id) return failure(cbio, "invalid_input", "study_id is required", false);
      const row = asRecord(await cbioGet(client, `/studies/${encodeURIComponent(id)}`));
      if (!text(row.studyId)) return failure(cbio, "not_found", `study ${id} was not found`, false);
      const access = row.publicStudy === false || row.readPermission === false ? "controlled" : "open";
      return makePage("dataset", [dataset({
        id: curie("cbioportal.study", id),
        title: text(row.name),
        description: text(row.description),
        sample_count: numberOrNull(row.allSampleCount),
        access,
        url: access === "open" ? `https://www.cbioportal.org/study/summary?id=${encodeURIComponent(id)}` : null,
        source_extra: { cancer_type_id: text(row.cancerTypeId), reference_genome: text(row.referenceGenome) },
      })], 0, 1, 1, cbio, cbioUsage);
    }),
  },
  {
    name: "cbioportal_list_samples",
    title: "List cBioPortal samples",
    description: "List samples for one cBioPortal study.",
    inputSchema: z.object({ study_id: z.string(), ...pageArgs }),
    call: (client, args) => guard(cbio, async () => {
      const { offset, limit } = paging(args);
      const studyId = text(args.study_id);
      if (!studyId) return failure(cbio, "invalid_input", "study_id is required", false);
      const rows = asArray(await cbioGet(client, `/studies/${encodeURIComponent(studyId)}/samples`, { pageSize: limit, pageNumber: Math.floor(offset / limit) })).map(asRecord);
      const items = rows.map((row) => sample({
        id: curie("cbioportal.sample", `${studyId}:${text(row.sampleId) ?? "unknown"}`),
        name: text(row.sampleId),
        donor_id: text(row.patientId) ? curie("cbioportal.patient", `${studyId}:${text(row.patientId)}`) : null,
        source_extra: { sample_type: text(row.sampleType), study_id: studyId },
      }));
      return makePage("sample", items, offset, limit, rows.length < limit ? offset + rows.length : null, cbio, cbioUsage);
    }),
  },
  {
    name: "cbioportal_list_patients",
    title: "List cBioPortal patients",
    description: "List patients for one cBioPortal study as donor entities.",
    inputSchema: z.object({ study_id: z.string(), ...pageArgs }),
    call: (client, args) => guard(cbio, async () => {
      const { offset, limit } = paging(args);
      const studyId = text(args.study_id);
      if (!studyId) return failure(cbio, "invalid_input", "study_id is required", false);
      const rows = asArray(await cbioGet(client, `/studies/${encodeURIComponent(studyId)}/patients`, { pageSize: limit, pageNumber: Math.floor(offset / limit) })).map(asRecord);
      const items = rows.map((row) => donor({
        id: curie("cbioportal.patient", `${studyId}:${text(row.patientId) ?? "unknown"}`),
        source_extra: { study_id: studyId },
      }));
      return makePage("donor", items, offset, limit, rows.length < limit ? offset + rows.length : null, cbio, cbioUsage);
    }),
  },
  {
    name: "cbioportal_get_gene",
    title: "Get a cBioPortal gene",
    description: "Fetch one cBioPortal gene by Hugo symbol. The id is an ncbigene CURIE.",
    inputSchema: z.object({ gene: z.string() }),
    call: (client, args) => guard(cbio, async () => {
      const symbol = text(args.gene);
      if (!symbol) return failure(cbio, "invalid_input", "gene is required", false);
      const row = asRecord(await cbioGet(client, `/genes/${encodeURIComponent(symbol)}`));
      const entrez = numberOrNull(row.entrezGeneId);
      if (entrez == null) return failure(cbio, "not_found", `gene ${symbol} was not found`, false);
      return makePage("gene", [gene({
        id: curie("ncbigene", String(entrez)),
        symbol: text(row.hugoGeneSymbol),
        biotype: text(row.type),
        taxon_id: "NCBITaxon:9606",
        source_extra: {},
      })], 0, 1, 1, cbio, cbioUsage);
    }),
  },
  {
    name: "pdc_list_studies",
    title: "List PDC studies",
    description: "List Proteomic Data Commons studies, including CPTAC and other programs present in PDC.",
    inputSchema: z.object({ ...pageArgs }),
    call: (client, args) => guard(pdc, async () => {
      const { offset, limit } = paging(args);
      const result = await pdcQuery(client, `query ($offset: Int, $limit: Int) { getPaginatedUIStudy(offset: $offset, limit: $limit) { total uiStudies { ${studyFields} } } }`, { offset, limit });
      if ("entity_type" in asRecord(result)) return result as ReturnType<typeof failure>;
      const payload = asRecord(asRecord(asRecord(result).data).getPaginatedUIStudy);
      const rows = asArray(payload.uiStudies).map(asRecord);
      return makePage("dataset", rows.map(mapPdc), offset, limit, numberOrNull(payload.total), pdc, pdcUsage);
    }),
  },
  {
    name: "pdc_get_study",
    title: "Get a PDC study",
    description: "Fetch one PDC study by PDC study id, for example PDC000759.",
    inputSchema: z.object({ study_id: z.string() }),
    call: (client, args) => guard(pdc, async () => {
      const id = text(args.study_id);
      if (!id || !/^PDC\d+$/.test(id)) return failure(pdc, "invalid_input", "study_id must look like PDC000759", false);
      const result = await pdcQuery(client, `query ($id: String) { getPaginatedUIStudy(offset: 0, limit: 1, pdc_study_id: $id) { total uiStudies { ${studyFields} } } }`, { id });
      if ("entity_type" in asRecord(result)) return result as ReturnType<typeof failure>;
      const row = asRecord(asArray(asRecord(asRecord(asRecord(result).data).getPaginatedUIStudy).uiStudies)[0]);
      if (text(row.pdc_study_id) !== id) return failure(pdc, "not_found", `study ${id} was not found`, false);
      return makePage("dataset", [mapPdc(row)], 0, 1, 1, pdc, pdcUsage);
    }),
  },
];
