import { curie } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  dataset,
  failure,
  file,
  guard,
  makePage,
  paging,
  sliceList,
  sourceOf,
  text,
  type OmicsTool,
} from "phi-omics-runtime";
import { z } from "zod";

const mtbls = sourceOf("metabolights", "MetaboLights", "https://www.ebi.ac.uk/metabolights/");
const mtblsUsage = {
  summary: "MetaboLights study, assay, and file metadata are returned for public studies. Private studies omit file URLs. The studies API does not expose a parsed metabolite table, so this server does not invent metabolite rows.",
  url: "https://www.ebi.ac.uk/metabolights/about",
};
const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

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

export const metabolomicsTools: OmicsTool[] = [
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
