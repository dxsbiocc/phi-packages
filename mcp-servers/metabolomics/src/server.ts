import { sourceOf, usagePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { metabolomicsTools } from "./sources.js";

const summary = [
  "MetaboLights: public study, assay, and file metadata. Parsed metabolite tables are not exposed by the verified JSON API, so no metabolite rows are returned. https://www.ebi.ac.uk/metabolights/",
].join("\n");

const usageTool: OmicsTool = {
  name: "usage_terms",
  title: "Usage terms",
  description: "Return usage terms for this metabolomics server.",
  inputSchema: z.object({}),
  call: async () => usagePage(sourceOf("metabolomics", "Metabolomics", "https://github.com/dxsbiocc/phi-packages"), summary),
};

export const tools: OmicsTool[] = [usageTool, ...metabolomicsTools];
export const serverModule: ServerModule = { name: "metabolomics", version: "1.0.0", mcpPath: "/mcp", tools };
