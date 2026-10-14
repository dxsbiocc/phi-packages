import { sourceOf, usagePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { cellxgeneTools } from "./cellxgene.js";
import { gtexTools } from "./gtex.js";
import { cellguideTools, hcaTools, hubmapTools } from "./more.js";

const summary = [
  "GTEx: open aggregate expression and eQTL summaries. Cite GTEx. https://gtexportal.org/home/datasets",
  "CELLxGENE: collection and dataset metadata only; matrix files are not downloaded. https://cellxgene.cziscience.com/",
  "CellGuide: canonical markers from pinned public snapshot CELLGUIDE_SNAPSHOT. https://cellguide.cellxgene.cziscience.com/",
  "HCA: project metadata. Non-NRES data use restrictions are controlled and have no data URL. https://data.humancellatlas.org/",
  "HuBMAP: public search metadata. Non-public access levels omit portal URLs. https://portal.hubmapconsortium.org/",
].join("\n");

const usageTool: OmicsTool = {
  name: "usage_terms",
  title: "Usage terms",
  description: "Return usage terms for every source in this server. The page is empty; read usage.summary.",
  inputSchema: z.object({}),
  call: async () => usagePage(sourceOf("transcriptomics-single-cell", "Transcriptomics and single-cell", "https://github.com/dxsbiocc/phi-packages"), summary),
};

export const tools: OmicsTool[] = [usageTool, ...gtexTools, ...cellxgeneTools, ...cellguideTools, ...hcaTools, ...hubmapTools];

export const serverModule: ServerModule = {
  name: "transcriptomics-single-cell",
  version: "1.0.0",
  mcpPath: "/mcp",
  tools,
};
