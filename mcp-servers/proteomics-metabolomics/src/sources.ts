import { curie, ensemblGeneCurie, humanChrom } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  dataset,
  failure,
  file,
  gene,
  guard,
  locusOf,
  makePage,
  numberOrNull,
  paging,
  protein,
  sliceList,
  sourceOf,
  text,
  type JsonClient,
  type OmicsTool,
} from "phi-omics-runtime";
import { z } from "zod";

const pride = sourceOf("pride", "PRIDE Archive", "https://www.ebi.ac.uk/pride/");
const prideUsage = {
  summary: "PRIDE project and file metadata are public archive records. File URLs are HTTPS conversions of public FTP locations. This server does not download raw spectra.",
  url: "https://www.ebi.ac.uk/pride/markdownpage/pride-archive-license",
};
const hpa = sourceOf("hpa", "Human Protein Atlas", "https://www.proteinatlas.org/");
const hpaUsage = {
  summary: "Human Protein Atlas gene and protein summaries are open. Tissue specificity fields are copied into source_extra; this response is not a GTEx TPM matrix.",
  url: "https://www.proteinatlas.org/about/licence",
};
const mtbls = sourceOf("metabolights", "MetaboLights", "https://www.ebi.ac.uk/metabolights/");
const mtblsUsage = {
  summary: "MetaboLights study, assay, and file metadata are returned for public studies. Private studies omit file URLs. The studies API does not expose a parsed metabolite table, so this server does not invent metabolite rows.",
  url: "https://www.ebi.ac.uk/metabolights/about",
};
const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

function mapPrideProject(row: Record<string, unknown>) {
  const accession = text(row.accession) ?? "unknown";
  return dataset({
    id: curie("pride", accession),
    title: text(row.title),
    description: text(row.projectDescription),
    taxon_id: "NCBITaxon:9606",
    access: "open",
    url: `https://www.ebi.ac.uk/pride/archive/projects/${accession}`,
    source_extra: { submission_date: text(row.submissionDate) },
  });
}

function httpsFtp(value: string | null): string | null {
  if (!value) return null;
  return value.replace(/^ftp:\/\/ftp\.pride\.ebi\.ac\.uk\//, "https://ftp.pride.ebi.ac.uk/");
}

function mapPrideFile(row: Record<string, unknown>, project: string | null) {
  const locations = asArray(row.publicFileLocations).map(asRecord);
  const ftp = locations.map((item) => text(item.value)).find((item) => item?.startsWith("ftp://")) ?? null;
  const accession = text(row.accession) ?? text(row.fileName) ?? "unknown";
  const access = ftp ? "open" : "controlled";
  const datasetId = project ?? text(asArray(row.projectAccessions)[0]);
  return file({
    id: curie("pride.file", accession),
    name: text(row.fileName),
    format: text(row.fileName)?.split(".").pop() ?? null,
    data_type: text(asRecord(row.fileCategory).name),
    size_bytes: numberOrNull(row.fileSizeBytes),
    access,
    url: access === "open" ? httpsFtp(ftp) : null,
    dataset_id: datasetId ? curie("pride", datasetId) : null,
    source_extra: {},
  });
}

async function pridePage(client: JsonClient, path: string, offset: number, limit: number) {
  const url = new URL(path, "https://www.ebi.ac.uk");
  url.searchParams.set("pageSize", String(limit));
  url.searchParams.set("page", String(Math.floor(offset / limit)));
  const rows = asArray(await client.getJson(url.toString())).map(asRecord);
  return rows;
}

function mapHpaGene(row: Record<string, unknown>) {
  const ensembl = text(row.Ensembl);
  const parsed = ensembl ? ensemblGeneCurie(ensembl) : null;
  const uniprot = asArray(row.Uniprot).map(text).filter((item): item is string => item != null);
  const position = text(row.Position);
  const coords = position?.match(/^(\d+)-(\d+)$/);
  return gene({
    id: parsed?.id ?? curie("hgnc.symbol", text(row.Gene) ?? "unknown"),
    symbol: text(row.Gene),
    name: text(row["Gene description"]),
    taxon_id: "NCBITaxon:9606",
    synonyms: text(row["Gene synonym"]) ? [text(row["Gene synonym"])!] : [],
    xrefs: uniprot.map((id) => curie("uniprot", id)),
    locus: locusOf({
      genomeBuild: "GRCh38",
      chrom: humanChrom(text(row.Chromosome)),
      start: coords ? Number(coords[1]) : null,
      end: coords ? Number(coords[2]) : null,
    }),
    source_extra: {
      rna_tissue_specificity: text(row["RNA tissue specificity"]),
      rna_tissue_distribution: text(row["RNA tissue distribution"]),
    },
  });
}

function mapStudy(row: Record<string, unknown>, accession: string) {
  const isa = asRecord(row.isaInvestigation);
  const nested = asRecord(asArray(isa.studies)[0]);
  const status = text(row.studyStatus);
  const open = status == null || status.toLowerCase() === "public";
  return dataset({
    id: curie("mtbls", accession),
    title: text(nested.title) ?? text(isa.title) ?? text(row.title),
    description: text(nested.description) ?? text(isa.description) ?? text(row.description),
    access: open ? "open" : "controlled",
    url: open ? `https://www.ebi.ac.uk/metabolights/${accession}` : null,
    source_extra: { study_status: status, study_http_url: text(row.studyHttpUrl) },
  });
}

function mapMtblsFile(row: Record<string, unknown>, accession: string, isPrivate: boolean, dataType: string | null) {
  const name = text(row.file) ?? text(row.filename) ?? "unknown";
  const access = isPrivate ? "controlled" : "open";
  return file({
    id: curie("mtbls.file", `${accession}:${encodeURIComponent(name)}`),
    name,
    data_type: dataType ?? text(row.type),
    access,
    url: access === "open" ? `https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/${accession}/${encodeURIComponent(name)}` : null,
    dataset_id: curie("mtbls", accession),
    source_extra: { status: text(row.status), directory: row.directory ?? null },
  });
}

export const proteomicsTools: OmicsTool[] = [
  {
    name: "pride_search_projects",
    title: "Search PRIDE projects",
    description: "Search PRIDE Archive projects as dataset entities.",
    inputSchema: z.object({ keyword: z.string().optional(), ...pageArgs }),
    call: (client, args) => guard(pride, async () => {
      const { offset, limit } = paging(args);
      const url = new URL("https://www.ebi.ac.uk/pride/ws/archive/v2/projects");
      url.searchParams.set("pageSize", String(limit));
      url.searchParams.set("page", String(Math.floor(offset / limit)));
      if (text(args.keyword)) url.searchParams.set("keyword", text(args.keyword)!);
      const rows = asArray(await client.getJson(url.toString())).map(asRecord);
      return makePage("dataset", rows.map(mapPrideProject), offset, limit, rows.length < limit ? offset + rows.length : null, pride, prideUsage);
    }),
  },
  {
    name: "pride_get_project",
    title: "Get a PRIDE project",
    description: "Fetch one PRIDE project by accession, for example PXD001357.",
    inputSchema: z.object({ accession: z.string() }),
    call: (client, args) => guard(pride, async () => {
      const accession = text(args.accession);
      if (!accession || !/^PXD\d+$/i.test(accession)) return failure(pride, "invalid_input", "accession must look like PXD001357", false);
      const row = asRecord(await client.getJson(`https://www.ebi.ac.uk/pride/ws/archive/v2/projects/${accession}`));
      if (!text(row.accession)) return failure(pride, "not_found", `project ${accession} was not found`, false);
      return makePage("dataset", [mapPrideProject(row)], 0, 1, 1, pride, prideUsage);
    }),
  },
  {
    name: "pride_list_files",
    title: "List PRIDE project files",
    description: "List file metadata for one PRIDE project. URLs are present only for public files.",
    inputSchema: z.object({ accession: z.string(), ...pageArgs }),
    call: (client, args) => guard(pride, async () => {
      const { offset, limit } = paging(args);
      const accession = text(args.accession);
      if (!accession) return failure(pride, "invalid_input", "accession is required", false);
      const rows = await pridePage(client, `/pride/ws/archive/v2/projects/${encodeURIComponent(accession)}/files`, offset, limit);
      return makePage("file", rows.map((row) => mapPrideFile(row, accession)), offset, limit, rows.length < limit ? offset + rows.length : null, pride, prideUsage);
    }),
  },
  {
    name: "hpa_search_genes",
    title: "Search Human Protein Atlas genes",
    description: "Search Human Protein Atlas genes. Results use Ensembl CURIEs and UniProt cross-references.",
    inputSchema: z.object({ query: z.string(), ...pageArgs }),
    call: (client, args) => guard(hpa, async () => {
      const { offset, limit } = paging(args);
      const query = text(args.query);
      if (!query) return failure(hpa, "invalid_input", "query is required", false);
      const url = new URL("https://www.proteinatlas.org/api/search_download.php");
      url.searchParams.set("search", query);
      url.searchParams.set("format", "json");
      url.searchParams.set("columns", "g,gs,eg,up");
      url.searchParams.set("compress", "no");
      const { rows, total } = sliceList(asArray(await client.getJson(url.toString())).map(asRecord), offset, limit);
      return makePage("gene", rows.map(mapHpaGene), offset, limit, total, hpa, hpaUsage);
    }),
  },
  {
    name: "hpa_get_gene",
    title: "Get a Human Protein Atlas gene",
    description: "Fetch one Human Protein Atlas gene JSON by Ensembl id.",
    inputSchema: z.object({ gene: z.string().describe("Ensembl gene id, with or without the ensembl: prefix.") }),
    call: (client, args) => guard(hpa, async () => {
      const query = text(args.gene);
      if (!query) return failure(hpa, "invalid_input", "gene is required", false);
      let ensembl;
      try {
        ensembl = ensemblGeneCurie(query);
      } catch {
        return failure(hpa, "invalid_input", "gene must be an Ensembl id such as ENSG00000141510", false);
      }
      const row = asRecord(await client.getJson(`https://www.proteinatlas.org/${ensembl.id.slice("ensembl:".length)}.json`));
      if (!text(row.Ensembl) && !text(row.Gene)) return failure(hpa, "not_found", `gene ${query} was not found`, false);
      return makePage("gene", [mapHpaGene(row)], 0, 1, 1, hpa, hpaUsage);
    }),
  },
  {
    name: "hpa_list_proteins",
    title: "List Human Protein Atlas proteins",
    description: "List UniProt proteins linked from one Human Protein Atlas gene.",
    inputSchema: z.object({ gene: z.string(), ...pageArgs }),
    call: (client, args) => guard(hpa, async () => {
      const { offset, limit } = paging(args);
      const query = text(args.gene);
      if (!query) return failure(hpa, "invalid_input", "gene is required", false);
      let ensembl;
      try {
        ensembl = ensemblGeneCurie(query);
      } catch {
        return failure(hpa, "invalid_input", "gene must be an Ensembl id such as ENSG00000141510", false);
      }
      const row = asRecord(await client.getJson(`https://www.proteinatlas.org/${ensembl.id.slice("ensembl:".length)}.json`));
      const ids = asArray(row.Uniprot).map(text).filter((item): item is string => item != null);
      const { rows, total } = sliceList(ids, offset, limit);
      const items = rows.map((id) => protein({
        id: curie("uniprot", id),
        name: text(row["Gene description"]),
        gene_id: ensembl.id,
        gene_symbol: text(row.Gene),
        taxon_id: "NCBITaxon:9606",
        source_extra: { protein_class: row["Protein class"] ?? null },
      }));
      return makePage("protein", items, offset, limit, total, hpa, hpaUsage);
    }),
  },
  {
    name: "metabolights_list_studies",
    title: "List MetaboLights studies",
    description: "List public MetaboLights study accessions. Titles are filled by metabolights_get_study.",
    inputSchema: z.object({ ...pageArgs }),
    call: (client, args) => guard(mtbls, async () => {
      const { offset, limit } = paging(args);
      const body = asRecord(await client.getJson("https://www.ebi.ac.uk/metabolights/ws/studies"));
      const ids = asArray(body.content).map(text).filter((item): item is string => item != null);
      const { rows, total } = sliceList(ids, offset, limit);
      const items = rows.map((id) => dataset({ id: curie("mtbls", id), access: "open", url: `https://www.ebi.ac.uk/metabolights/${id}` }));
      return makePage("dataset", items, offset, limit, total, mtbls, mtblsUsage);
    }),
  },
  {
    name: "metabolights_get_study",
    title: "Get a MetaboLights study",
    description: "Fetch one MetaboLights study, for example MTBLS1.",
    inputSchema: z.object({ study_id: z.string() }),
    call: (client, args) => guard(mtbls, async () => {
      const id = text(args.study_id);
      if (!id || !/^MTBLS\d+$/i.test(id)) return failure(mtbls, "invalid_input", "study_id must look like MTBLS1", false);
      const body = asRecord(await client.getJson(`https://www.ebi.ac.uk/metabolights/ws/studies/${id}`));
      const study = asRecord(body.mtblsStudy);
      if (!text(study.studyStatus) && !text(study.studyHttpUrl)) return failure(mtbls, "not_found", `study ${id} was not found`, false);
      return makePage("dataset", [mapStudy(study, id)], 0, 1, 1, mtbls, mtblsUsage);
    }),
  },
  {
    name: "metabolights_list_files",
    title: "List MetaboLights files",
    description: "List files for one public MetaboLights study. Private studies omit URLs.",
    inputSchema: z.object({ study_id: z.string(), ...pageArgs }),
    call: (client, args) => guard(mtbls, async () => {
      const { offset, limit } = paging(args);
      const id = text(args.study_id);
      if (!id) return failure(mtbls, "invalid_input", "study_id is required", false);
      const body = asRecord(await client.getJson(`https://www.ebi.ac.uk/metabolights/ws/studies/${encodeURIComponent(id)}/files`));
      const { rows, total } = sliceList(asArray(body.study).map(asRecord).filter((row) => row.directory !== true), offset, limit);
      return makePage("file", rows.map((row) => mapMtblsFile(row, id, body.private === true, null)), offset, limit, total, mtbls, mtblsUsage);
    }),
  },
  {
    name: "metabolights_list_assays",
    title: "List MetaboLights assay files",
    description: "List assay file names for one MetaboLights study as file entities.",
    inputSchema: z.object({ study_id: z.string(), ...pageArgs }),
    call: (client, args) => guard(mtbls, async () => {
      const { offset, limit } = paging(args);
      const id = text(args.study_id);
      if (!id) return failure(mtbls, "invalid_input", "study_id is required", false);
      const body = asRecord(await client.getJson(`https://www.ebi.ac.uk/metabolights/ws/studies/${encodeURIComponent(id)}/assays`));
      const { rows, total } = sliceList(asArray(asRecord(body.data).assays).map(asRecord), offset, limit);
      return makePage("file", rows.map((row) => mapMtblsFile(row, id, false, "assay")), offset, limit, total, mtbls, mtblsUsage);
    }),
  },
];
