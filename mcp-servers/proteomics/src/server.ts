import { sourceOf, usagePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { proteomicsTools } from "./sources.js";

const summary = [
  "PRIDE: public project and file metadata. https://www.ebi.ac.uk/pride/",
  "Human Protein Atlas: gene and UniProt protein summaries. https://www.proteinatlas.org/about/licence",
].join("\n");

const usageTool: OmicsTool = {
  name: "usage_terms",
  title: "Usage terms",
  description: "Return usage terms for this proteomics server.",
  inputSchema: z.object({}),
  call: async () => usagePage(sourceOf("proteomics", "Proteomics", "https://github.com/dxsbiocc/phi-packages"), summary),
};

export const tools: OmicsTool[] = [usageTool, ...proteomicsTools];
export const serverModule: ServerModule = { name: "proteomics", version: "1.0.0", mcpPath: "/mcp", tools };
