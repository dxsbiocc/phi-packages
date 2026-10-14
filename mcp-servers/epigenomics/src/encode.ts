import { curie, ontologyCurie, taxonCurie } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  dataset,
  failure,
  file,
  guard,
  makePage,
  numberOrNull,
  paging,
  sample,
  sourceOf,
  text,
  type JsonClient,
  type OmicsTool,
} from "phi-omics-runtime";
import { z } from "zod";

const source = sourceOf("encode", "ENCODE", "https://www.encodeproject.org/");
const usage = {
  summary: "ENCODE and Roadmap Epigenomics records are metadata from the ENCODE portal. Released, unrestricted file URLs are included; restricted or unreleased files keep url null. This server does not download file bytes. Roadmap experiments are the ENCODE award.rfa=Roadmap subset. ChIP-Atlas was not included: its experiment API returned AccessDenied.",
  url: "https://www.encodeproject.org/help/getting-started/",
};
const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

function accessionFrom(value: unknown): string | null {
  if (typeof value === "string") {
    const match = value.match(/ENC[A-Z]{2}\d+[A-Z]+/);
    return match?.[0] ?? null;
  }
  return text(asRecord(value).accession);
}

async function search(client: JsonClient, type: string, extra: Record<string, string>, offset: number, limit: number, fields: string[]) {
  // ENCODE's edge blocks the `from` query parameter, so request the window and slice it.
  const url = new URL("https://www.encodeproject.org/search/");
  url.searchParams.set("type", type);
  url.searchParams.set("status", "released");
  url.searchParams.set("format", "json");
  url.searchParams.set("limit", String(offset + limit));
  for (const [key, value] of Object.entries(extra)) url.searchParams.set(key, value);
  for (const field of fields) url.searchParams.append("field", field);
  const body = asRecord(await client.getJson(url.toString()));
  const rows = asArray(body["@graph"]).map(asRecord).slice(offset, offset + limit);
  return { rows, total: numberOrNull(body.total) };
}

function mapExperiment(row: Record<string, unknown>) {
  const accession = text(row.accession) ?? "unknown";
  const ontology = asRecord(row.biosample_ontology);
  const term = ontologyCurie(text(ontology.term_id));
  const classification = text(ontology.classification);
  const released = text(row.status) === "released" || text(row.status) == null;
  return dataset({
    id: curie("encode", accession),
    title: text(row.description) ?? text(row.biosample_summary),
    description: text(row.biosample_summary),
    assay: text(row.assay_title) ?? text(row.assay_term_name),
    taxon_id: "NCBITaxon:9606",
    tissue_ids: classification === "tissue" && term?.startsWith("UBERON:") ? [term] : [],
    cell_type_ids: classification !== "tissue" && term?.startsWith("CL:") ? [term] : [],
    access: released ? "open" : "unknown",
    url: released ? `https://www.encodeproject.org/experiments/${accession}/` : null,
    source_extra: { biosample_term: text(ontology.term_name), biosample_ontology: term, award_rfa: text(asRecord(row.award).rfa) },
  });
}

function mapFile(row: Record<string, unknown>) {
  const accession = text(row.accession) ?? "unknown";
  const restricted = row.restricted === true;
  const released = text(row.status) === "released";
  const access = restricted ? "controlled" : released ? "open" : "unknown";
  const href = text(row.href);
  const datasetAccession = accessionFrom(row.dataset);
  return file({
    id: curie("encode", accession),
    name: accession,
    format: text(row.file_format),
    data_type: text(row.output_type),
    size_bytes: numberOrNull(row.file_size),
    access,
    url: access === "open" && href ? new URL(href, "https://www.encodeproject.org").toString() : null,
    dataset_id: datasetAccession ? curie("encode", datasetAccession) : null,
    source_extra: { assay: text(row.assay_term_name), status: text(row.status) },
  });
}

function mapBiosample(row: Record<string, unknown>) {
  const accession = text(row.accession) ?? "unknown";
  const ontology = asRecord(row.biosample_ontology);
  const term = ontologyCurie(text(ontology.term_id));
  const organism = taxonCurie(text(asRecord(row.organism).scientific_name) ?? text(ontology.term_name));
  return sample({
    id: curie("encode", accession),
    name: text(row.summary) ?? text(ontology.term_name),
    donor_id: accessionFrom(row.donor) ? curie("encode", accessionFrom(row.donor)!) : null,
    tissue_id: term?.startsWith("UBERON:") ? term : null,
    cell_type_id: term && !term.startsWith("UBERON:") ? term : null,
    taxon_id: organism,
    source_extra: { classification: text(ontology.classification), term_name: text(ontology.term_name) },
  });
}

const experimentFields = ["accession", "assay_title", "assay_term_name", "description", "biosample_summary", "status", "biosample_ontology", "award"];
const fileFields = ["accession", "file_format", "file_size", "href", "status", "output_type", "restricted", "assay_term_name", "dataset"];
const biosampleFields = ["accession", "summary", "biosample_ontology", "donor", "organism"];

export const encodeTools: OmicsTool[] = [
  {
    name: "encode_search_experiments",
    title: "Search ENCODE experiments",
    description: "Search released ENCODE experiments as dataset entities. Optional search and assay_title narrow the query.",
    inputSchema: z.object({ search: z.string().optional(), assay_title: z.string().optional().describe("For example ChIP-seq or ATAC-seq."), ...pageArgs }),
    call: (client, args) => guard(source, async () => {
      const { offset, limit } = paging(args);
      if (offset + limit > 100) return failure(source, "invalid_input", "ENCODE search supports offset+limit <= 100", false);
      const extra: Record<string, string> = {};
      if (text(args.search)) extra.searchTerm = text(args.search)!;
      if (text(args.assay_title)) extra.assay_title = text(args.assay_title)!;
      const { rows, total } = await search(client, "Experiment", extra, offset, limit, experimentFields);
      return makePage("dataset", rows.map(mapExperiment), offset, limit, total, source, usage);
    }),
  },
  {
    name: "encode_get_experiment",
    title: "Get an ENCODE experiment",
    description: "Fetch one ENCODE experiment by accession, for example ENCSR818VNL.",
    inputSchema: z.object({ accession: z.string() }),
    call: (client, args) => guard(source, async () => {
      const accession = text(args.accession);
      if (!accession || !/^ENC/.test(accession)) return failure(source, "invalid_input", "accession must be an ENCODE accession", false);
      const row = asRecord(await client.getJson(`https://www.encodeproject.org/experiments/${encodeURIComponent(accession)}/?format=json`));
      if (!text(row.accession)) return failure(source, "not_found", `experiment ${accession} was not found`, false);
      return makePage("dataset", [mapExperiment(row)], 0, 1, 1, source, usage);
    }),
  },
  {
    name: "encode_search_roadmap",
    title: "Search Roadmap Epigenomics experiments",
    description: "Search ENCODE experiments funded under award.rfa=Roadmap.",
    inputSchema: z.object({ search: z.string().optional(), ...pageArgs }),
    call: (client, args) => guard(source, async () => {
      const { offset, limit } = paging(args);
      if (offset + limit > 100) return failure(source, "invalid_input", "ENCODE search supports offset+limit <= 100", false);
      const extra: Record<string, string> = { "award.rfa": "Roadmap" };
      if (text(args.search)) extra.searchTerm = text(args.search)!;
      const { rows, total } = await search(client, "Experiment", extra, offset, limit, experimentFields);
      return makePage("dataset", rows.map(mapExperiment), offset, limit, total, source, usage, { collection: "roadmap" });
    }),
  },
  {
    name: "encode_search_files",
    title: "Search ENCODE file metadata",
    description: "Search released ENCODE file metadata. Restricted files have url null.",
    inputSchema: z.object({ search: z.string().optional(), file_format: z.string().optional(), ...pageArgs }),
    call: (client, args) => guard(source, async () => {
      const { offset, limit } = paging(args);
      if (offset + limit > 100) return failure(source, "invalid_input", "ENCODE search supports offset+limit <= 100", false);
      const extra: Record<string, string> = {};
      if (text(args.search)) extra.searchTerm = text(args.search)!;
      if (text(args.file_format)) extra.file_format = text(args.file_format)!;
      const { rows, total } = await search(client, "File", extra, offset, limit, fileFields);
      return makePage("file", rows.map(mapFile), offset, limit, total, source, usage);
    }),
  },
  {
    name: "encode_get_file",
    title: "Get ENCODE file metadata",
    description: "Fetch metadata for one ENCODE file accession. Restricted files have url null.",
    inputSchema: z.object({ accession: z.string() }),
    call: (client, args) => guard(source, async () => {
      const accession = text(args.accession);
      if (!accession || !/^ENCFF/.test(accession)) return failure(source, "invalid_input", "accession must be an ENCFF file accession", false);
      const row = asRecord(await client.getJson(`https://www.encodeproject.org/files/${encodeURIComponent(accession)}/?format=json`));
      if (!text(row.accession)) return failure(source, "not_found", `file ${accession} was not found`, false);
      return makePage("file", [mapFile(row)], 0, 1, 1, source, usage);
    }),
  },
  {
    name: "encode_search_biosamples",
    title: "Search ENCODE biosamples",
    description: "Search released ENCODE biosamples as sample entities.",
    inputSchema: z.object({ search: z.string().optional(), ...pageArgs }),
    call: (client, args) => guard(source, async () => {
      const { offset, limit } = paging(args);
      if (offset + limit > 100) return failure(source, "invalid_input", "ENCODE search supports offset+limit <= 100", false);
      const extra: Record<string, string> = {};
      if (text(args.search)) extra.searchTerm = text(args.search)!;
      const { rows, total } = await search(client, "Biosample", extra, offset, limit, biosampleFields);
      return makePage("sample", rows.map(mapBiosample), offset, limit, total, source, usage);
    }),
  },
  {
    name: "encode_get_biosample",
    title: "Get an ENCODE biosample",
    description: "Fetch one released ENCODE biosample by ENCBS accession.",
    inputSchema: z.object({ accession: z.string() }),
    call: (client, args) => guard(source, async () => {
      const accession = text(args.accession);
      if (!accession || !/^ENCBS/.test(accession)) return failure(source, "invalid_input", "accession must be an ENCBS biosample accession", false);
      const row = asRecord(await client.getJson(`https://www.encodeproject.org/biosamples/${encodeURIComponent(accession)}/?format=json`));
      if (!text(row.accession)) return failure(source, "not_found", `biosample ${accession} was not found`, false);
      return makePage("sample", [mapBiosample(row)], 0, 1, 1, source, usage);
    }),
  },
];
