import { curie, ontologyCurie, sexOf, taxonCurie } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  cellType,
  dataset,
  donor,
  failure,
  guard,
  loadConfig,
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

const hca = sourceOf("hca", "Human Cell Atlas", "https://data.humancellatlas.org/");
const hcaUsage = {
  summary: "HCA Azul returns project metadata. dataUseRestriction other than NRES is treated as controlled: the response keeps the description and omits a data URL. This server does not download sequence files.",
  url: "https://data.humancellatlas.org/",
};
const hubmap = sourceOf("hubmap", "HuBMAP", "https://portal.hubmapconsortium.org/");
const hubmapUsage = {
  summary: "HuBMAP public entity metadata is returned. data_access_level other than public is controlled, and no file URL is included. Donor contact fields from the search index are not copied.",
  url: "https://hubmapconsortium.org/",
};
const cellguide = () => sourceOf("cellguide", "CellGuide", "https://cellguide.cellxgene.cziscience.com/", loadConfig().cellguideSnapshot);
const cellguideUsage = {
  summary: "CellGuide canonical marker genes are read from the public snapshot named by CELLGUIDE_SNAPSHOT (default 1692975365, published 2023-08-25). The bucket has no stable latest pointer.",
  url: "https://cellguide.cellxgene.cziscience.com/",
};
const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

function hcaAccess(value: unknown): "open" | "controlled" | "unknown" {
  const textValue = text(value);
  if (!textValue) return "unknown";
  return textValue.toUpperCase() === "NRES" ? "open" : "controlled";
}

function mapProject(hit: Record<string, unknown>) {
  const project = asRecord(asArray(hit.projects)[0]);
  const organisms = asRecord(asArray(hit.donorOrganisms)[0]);
  const specimens = asRecord(asArray(hit.specimens)[0]);
  const id = text(project.projectId) ?? text(hit.entryId) ?? "unknown";
  const access = hcaAccess(project.dataUseRestriction);
  const species = asArray(organisms.genusSpecies).map(text).find((item): item is string => item != null) ?? null;
  return dataset({
    id: curie("hca.project", id),
    title: text(project.projectTitle),
    description: text(project.projectDescription),
    taxon_id: taxonCurie(species),
    tissue_ids: asArray(specimens.organ).map((organ) => text(organ)).filter((organ): organ is string => organ != null).map((organ) => ontologyCurie(organ)).filter((organ): organ is string => organ != null),
    sample_count: numberOrNull(project.estimatedCellCount),
    access,
    url: access === "open" ? `https://data.humancellatlas.org/explore/projects/${id}` : null,
    source_extra: {
      data_use_restriction: text(project.dataUseRestriction),
      genus_species: organisms.genusSpecies ?? null,
      organ: specimens.organ ?? null,
      organ_part: specimens.organPart ?? null,
    },
  });
}

async function hcaProjects(client: JsonClient, filters: Record<string, unknown> | null, size: number) {
  const url = new URL("https://service.azul.data.humancellatlas.org/index/projects");
  url.searchParams.set("size", String(size));
  if (filters) url.searchParams.set("filters", JSON.stringify(filters));
  return asRecord(await client.getJson(url.toString()));
}

function hubmapAccess(value: unknown): "open" | "controlled" | "unknown" {
  const level = text(value)?.toLowerCase();
  if (level === "public") return "open";
  if (!level) return "unknown";
  return "controlled";
}

async function hubmapSearch(client: JsonClient, entityType: string, query: string | null, offset: number, limit: number) {
  const must: unknown[] = [{ term: { "entity_type.keyword": entityType } }];
  if (query) must.push({ query_string: { query, lenient: true } });
  const body = await client.postJson("https://search.api.hubmapconsortium.org/v3/search", {
    from: offset,
    size: limit,
    _source: ["hubmap_id", "uuid", "entity_type", "data_access_level", "title", "description", "status", "dataset_type", "donor.hubmap_id", "origin_samples.hubmap_id"],
    query: { bool: { must } },
  });
  const hits = asRecord(asRecord(body).hits);
  const totalInfo = asRecord(hits.total);
  const total = text(totalInfo.relation) === "eq" ? numberOrNull(totalInfo.value) : null;
  return { rows: asArray(hits.hits).map((hit) => asRecord(asRecord(hit)._source)), total };
}

export const hcaTools: OmicsTool[] = [
  {
    name: "hca_search_projects",
    title: "Search HCA projects",
    description: "Search Human Cell Atlas projects. Optional organ uses an Azul organ filter. offset+limit must be at most 100.",
    inputSchema: z.object({ organ: z.string().optional().describe("Uberon label or name used as an Azul organ filter, for example brain."), project_id: z.string().optional(), ...pageArgs }),
    call: (client, args) => guard(hca, async () => {
      const { offset, limit } = paging(args);
      if (offset + limit > 100) return failure(hca, "invalid_input", "HCA search supports offset+limit <= 100", false);
      const filters: Record<string, unknown> = {};
      if (text(args.project_id)) filters.projectId = { is: [text(args.project_id)] };
      if (text(args.organ)) filters.organ = { is: [text(args.organ)] };
      const body = await hcaProjects(client, Object.keys(filters).length ? filters : null, offset + limit);
      const hits = asArray(body.hits).map(asRecord).slice(offset, offset + limit);
      const total = numberOrNull(asRecord(body.pagination).total);
      return makePage("dataset", hits.map(mapProject), offset, limit, total, hca, hcaUsage);
    }),
  },
  {
    name: "hca_get_project",
    title: "Get an HCA project",
    description: "Fetch one HCA project by project id.",
    inputSchema: z.object({ project_id: z.string() }),
    call: (client, args) => guard(hca, async () => {
      const id = text(args.project_id);
      if (!id) return failure(hca, "invalid_input", "project_id is required", false);
      const body = await hcaProjects(client, { projectId: { is: [id] } }, 1);
      const hit = asRecord(asArray(body.hits)[0]);
      if (!text(asRecord(asArray(hit.projects)[0]).projectId)) return failure(hca, "not_found", `project ${id} was not found`, false);
      return makePage("dataset", [mapProject(hit)], 0, 1, 1, hca, hcaUsage);
    }),
  },
];

export const hubmapTools: OmicsTool[] = [
  {
    name: "hubmap_search_datasets",
    title: "Search HuBMAP datasets",
    description: "Search public HuBMAP dataset metadata.",
    inputSchema: z.object({ query: z.string().optional().describe("Optional search string passed to the HuBMAP query_string query."), ...pageArgs }),
    call: (client, args) => guard(hubmap, async () => {
      const { offset, limit } = paging(args);
      const { rows, total } = await hubmapSearch(client, "Dataset", text(args.query), offset, limit);
      const items = rows.map((row) => {
        const access = hubmapAccess(row.data_access_level);
        const id = text(row.hubmap_id) ?? text(row.uuid) ?? "unknown";
        return dataset({
          id: curie("hubmap", id),
          title: text(row.title),
          description: text(row.description),
          assay: text(row.dataset_type),
          access,
          url: access === "open" && text(row.uuid) ? `https://portal.hubmapconsortium.org/browse/dataset/${text(row.uuid)}` : null,
          source_extra: { uuid: text(row.uuid), status: text(row.status), data_access_level: text(row.data_access_level) },
        });
      });
      return makePage("dataset", items, offset, limit, total, hubmap, hubmapUsage);
    }),
  },
  {
    name: "hubmap_search_samples",
    title: "Search HuBMAP samples",
    description: "Search HuBMAP sample metadata.",
    inputSchema: z.object({ query: z.string().optional(), ...pageArgs }),
    call: (client, args) => guard(hubmap, async () => {
      const { offset, limit } = paging(args);
      const { rows, total } = await hubmapSearch(client, "Sample", text(args.query), offset, limit);
      const items = rows.map((row) => sample({
        id: curie("hubmap", text(row.hubmap_id) ?? text(row.uuid) ?? "unknown"),
        donor_id: text(asRecord(row.donor).hubmap_id) ? curie("hubmap", text(asRecord(row.donor).hubmap_id)!) : null,
        source_extra: { uuid: text(row.uuid), data_access_level: text(row.data_access_level) },
      }));
      return makePage("sample", items, offset, limit, total, hubmap, hubmapUsage);
    }),
  },
  {
    name: "hubmap_search_donors",
    title: "Search HuBMAP donors",
    description: "Search HuBMAP donor metadata. Sex and age are included only when the public index provides them.",
    inputSchema: z.object({ query: z.string().optional(), ...pageArgs }),
    call: (client, args) => guard(hubmap, async () => {
      const { offset, limit } = paging(args);
      const { rows, total } = await hubmapSearch(client, "Donor", text(args.query), offset, limit);
      const items = rows.map((row) => donor({
        id: curie("hubmap", text(row.hubmap_id) ?? text(row.uuid) ?? "unknown"),
        sex: sexOf(row.sex),
        source_extra: { uuid: text(row.uuid), data_access_level: text(row.data_access_level) },
      }));
      return makePage("donor", items, offset, limit, total, hubmap, hubmapUsage);
    }),
  },
];

export const cellguideTools: OmicsTool[] = [
  {
    name: "cellguide_cell_type",
    title: "CellGuide canonical markers",
    description: "Look up one Cell Ontology id in the pinned CellGuide marker snapshot. Returns one cell_type entity; marker symbols are in source_extra.",
    inputSchema: z.object({ cell_type_id: z.string().describe("CL CURIE or CL_0000034 snapshot key.") }),
    call: (client, args) => guard(cellguide(), async () => {
      const raw = text(args.cell_type_id);
      if (!raw) return failure(cellguide(), "invalid_input", "cell_type_id is required", false);
      const ontology = ontologyCurie(raw.startsWith("CL") ? raw : `CL:${raw}`);
      if (!ontology?.startsWith("CL:")) return failure(cellguide(), "invalid_input", "cell_type_id must be a CL id", false);
      const key = ontology.replace(":", "_");
      const url = `https://cellguide.cellxgene.cziscience.com/${loadConfig().cellguideSnapshot}/canonical_marker_genes/${key}.json`;
      const markers = asArray(await client.getJson(url)).map(asRecord).slice(0, 50).map((row) => ({
        symbol: text(row.symbol),
        tissue: text(row.tissue),
        publication: text(row.publication),
      }));
      return makePage("cell_type", [cellType({
        id: ontology,
        ontology_id: ontology,
        source_extra: { snapshot: loadConfig().cellguideSnapshot, markers },
      })], 0, 1, 1, cellguide(), cellguideUsage);
    }),
  },
];
