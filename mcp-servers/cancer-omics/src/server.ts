import { sourceOf, usagePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { cancerTools } from "./sources.js";

const summary = [
  "GDC: metadata for cases, genes, mutations, and files. Controlled files have url null and are not downloaded. https://gdc.cancer.gov/",
  "cBioPortal: public REST metadata for open studies. This is separate from the cBioPortal OAuth MCP connector. https://www.cbioportal.org/",
  "PDC: study metadata for CPTAC and other programs exposed by the Proteomic Data Commons. https://pdc.cancer.gov/",
  "DepMap was not included: the public models endpoint returned a bot-check page rather than JSON.",
].join("\n");

const usageTool: OmicsTool = {
  name: "usage_terms",
  title: "Usage terms",
  description: "Return usage terms for every source in this server.",
  inputSchema: z.object({}),
  call: async () => usagePage(sourceOf("cancer-omics", "Cancer omics", "https://github.com/dxsbiocc/phi-packages"), summary),
};

export const tools: OmicsTool[] = [usageTool, ...cancerTools];
export const serverModule: ServerModule = { name: "cancer-omics", version: "1.0.0", mcpPath: "/mcp", tools };
