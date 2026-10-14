import { ontologyCurie, taxonCurie } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  dataset,
  failure,
  guard,
  labels,
  makePage,
  numberOrNull,
  paging,
  sliceList,
  sourceOf,
  text,
  type JsonClient,
  type OmicsTool,
} from "phi-omics-runtime";
import { z } from "zod";

const source = sourceOf("cellxgene", "CZ CELLxGENE Discover", "https://cellxgene.cziscience.com/");
const usage = {
  summary: "CELLxGENE collection and dataset metadata are open. H5AD assets can be large; this server returns metadata and explorer URLs, not matrix bytes. Controlled raw data, when linked, stays outside this response.",
  url: "https://cellxgene.cziscience.com/docs",
};
const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

function ontologyIds(value: unknown, prefix: string): string[] {
  return labels(value).ids.map((id) => ontologyCurie(id)).filter((id): id is string => id != null && id.startsWith(prefix));
}

function mapDataset(row: Record<string, unknown>, collectionId: string | null) {
  const id = text(row.dataset_id) ?? "unknown";
  const organisms = labels(row.organism);
  const assays = labels(row.assay);
  return dataset({
    id: `cellxgene.dataset:${id}`,
    title: text(row.title),
    assay: assays.names.join(", ") || null,
    taxon_id: taxonCurie(organisms.ids[0] ?? organisms.names[0]),
    tissue_ids: ontologyIds(row.tissue, "UBERON:"),
    disease_ids: ontologyIds(row.disease, "MONDO:").concat(ontologyIds(row.disease, "PATO:")),
    cell_type_ids: ontologyIds(row.cell_type, "CL:"),
    sample_count: numberOrNull(row.cell_count),
    access: "open",
    url: text(row.explorer_url),
    source_extra: {
      collection_id: collectionId,
      dataset_version_id: text(row.dataset_version_id),
      tissue_labels: labels(row.tissue).names,
      cell_type_labels: labels(row.cell_type).names,
      disease_labels: labels(row.disease).names,
      is_primary_data: row.is_primary_data ?? null,
    },
  });
}

function mapCollection(row: Record<string, unknown>) {
  const id = text(row.collection_id) ?? "unknown";
  return dataset({
    id: `cellxgene.collection:${id}`,
    title: text(row.name),
    description: text(row.description),
    access: text(row.visibility) === "PRIVATE" ? "controlled" : "open",
    url: text(row.collection_url),
    sample_count: asArray(row.datasets).length,
    source_extra: { doi: text(row.doi), collection_version_id: text(row.collection_version_id) },
  });
}

async function collections(client: JsonClient): Promise<Record<string, unknown>[]> {
  return asArray(await client.getJson("https://api.cellxgene.cziscience.com/curation/v1/collections")).map(asRecord);
}

export const cellxgeneTools: OmicsTool[] = [
  {
    name: "cellxgene_list_collections",
    title: "List CELLxGENE collections",
    description: "List CELLxGENE Discover collections as dataset entities. The full catalog is cached, then paged locally.",
    inputSchema: z.object({ ...pageArgs }),
    call: (client, args) => guard(source, async () => {
      const { offset, limit } = paging(args);
      const { rows, total } = sliceList(await collections(client), offset, limit);
      return makePage("dataset", rows.map(mapCollection), offset, limit, total, source, usage);
    }),
  },
  {
    name: "cellxgene_get_collection",
    title: "Get a CELLxGENE collection",
    description: "Fetch one CELLxGENE collection by id.",
    inputSchema: z.object({ collection_id: z.string().describe("Collection UUID.") }),
    call: (client, args) => guard(source, async () => {
      const id = text(args.collection_id);
      if (!id) return failure(source, "invalid_input", "collection_id is required", false);
      const row = asRecord(await client.getJson(`https://api.cellxgene.cziscience.com/curation/v1/collections/${encodeURIComponent(id)}`));
      if (!text(row.collection_id)) return failure(source, "not_found", `collection ${id} was not found`, false);
      return makePage("dataset", [mapCollection(row)], 0, 1, 1, source, usage);
    }),
  },
  {
    name: "cellxgene_list_datasets",
    title: "List datasets in a CELLxGENE collection",
    description: "List datasets embedded in one CELLxGENE collection.",
    inputSchema: z.object({ collection_id: z.string(), ...pageArgs }),
    call: (client, args) => guard(source, async () => {
      const { offset, limit } = paging(args);
      const id = text(args.collection_id);
      if (!id) return failure(source, "invalid_input", "collection_id is required", false);
      const row = asRecord(await client.getJson(`https://api.cellxgene.cziscience.com/curation/v1/collections/${encodeURIComponent(id)}`));
      const { rows, total } = sliceList(asArray(row.datasets).map(asRecord), offset, limit);
      return makePage("dataset", rows.map((item) => mapDataset(item, id)), offset, limit, total, source, usage);
    }),
  },
  {
    name: "cellxgene_get_dataset",
    title: "Get a CELLxGENE dataset",
    description: "Fetch one dataset from its parent collection. The standalone dataset URL is not public.",
    inputSchema: z.object({ collection_id: z.string(), dataset_id: z.string() }),
    call: (client, args) => guard(source, async () => {
      const collectionId = text(args.collection_id);
      const datasetId = text(args.dataset_id);
      if (!collectionId || !datasetId) return failure(source, "invalid_input", "collection_id and dataset_id are required", false);
      const row = asRecord(await client.getJson(`https://api.cellxgene.cziscience.com/curation/v1/collections/${encodeURIComponent(collectionId)}`));
      const found = asArray(row.datasets).map(asRecord).find((item) => item.dataset_id === datasetId);
      if (!found) return failure(source, "not_found", `dataset ${datasetId} was not found in collection ${collectionId}`, false);
      return makePage("dataset", [mapDataset(found, collectionId)], 0, 1, 1, source, usage);
    }),
  },
];
