import { sourceOf, usagePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { pathwayTools } from "./sources.js";

const summary = [
  "Ensembl REST: gene lookup, canonical transcripts, comparative orthologs/paralogues, and Variant Effect Predictor (VEP) annotations. https://rest.ensembl.org/",
  "Reactome Content Service: biological pathway search, curated pathway hierarchy/literature details, and UniProt-to-pathway mappings (CC0). https://reactome.org/",
].join("\n");

const usageTool: OmicsTool = {
  name: "usage_terms",
  title: "Usage terms",
  description: "Return usage terms and attribution for this gene and pathway server.",
  inputSchema: z.object({}),
  call: async () => usagePage(sourceOf("gene-pathway", "Gene & Pathway", "https://github.com/dxsbiocc/phi-packages"), summary),
};

export const tools: OmicsTool[] = [usageTool, ...pathwayTools];
export const serverModule: ServerModule = { name: "gene-pathway", version: "1.0.0", mcpPath: "/mcp", tools };
