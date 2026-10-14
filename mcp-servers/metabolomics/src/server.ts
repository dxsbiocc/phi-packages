import { sourceOf, usagePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { hmdbTools } from "./hmdb.js";
import { metabolomicsTools } from "./sources.js";

const summary = [
  "MetaboLights: public study, assay, and file metadata. Parsed metabolite tables are not exposed by the verified JSON API, so MetaboLights tools do not return metabolite rows. https://www.ebi.ac.uk/metabolights/",
  "HMDB: public per-metabolite XML at https://hmdb.ca/metabolites/HMDBxxxxxxx.xml, mapped to the shared metabolite schema. The contact-only API at https://hmdb.ca/simple/api is not used, and the multi-gigabyte bulk XML download is not fetched. Name search reads the first page of the public HTML search. HMDB is free for public use; commercial redistribution requires permission of the authors and acknowledgment of HMDB and its papers. https://hmdb.ca/citing",
].join("\n");

const usageTool: OmicsTool = {
  name: "usage_terms",
  title: "Usage terms",
  description: "Return usage terms for this metabolomics server.",
  inputSchema: z.object({}),
  call: async () => usagePage(sourceOf("metabolomics", "Metabolomics", "https://github.com/dxsbiocc/phi-packages"), summary),
};

export const tools: OmicsTool[] = [usageTool, ...metabolomicsTools, ...hmdbTools];
export const serverModule: ServerModule = { name: "metabolomics", version: "1.0.0", mcpPath: "/mcp", tools };
