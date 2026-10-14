import { sourceOf, usagePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { encodeTools } from "./encode.js";

const summary = [
  "ENCODE: experiment, file, and biosample metadata. Restricted or unreleased files have url null. https://www.encodeproject.org/",
  "Roadmap Epigenomics: the same ENCODE experiment search filtered to award.rfa=Roadmap.",
  "ChIP-Atlas was not included: https://chip-atlas.dbcls.jp/api/v1/experiment returned AccessDenied, and no bounded JSON catalog was verified.",
].join("\n");

const usageTool: OmicsTool = {
  name: "usage_terms",
  title: "Usage terms",
  description: "Return usage terms for this epigenomics server.",
  inputSchema: z.object({}),
  call: async () => usagePage(sourceOf("epigenomics", "Epigenomics", "https://github.com/dxsbiocc/phi-packages"), summary),
};

export const tools: OmicsTool[] = [usageTool, ...encodeTools];
export const serverModule: ServerModule = { name: "epigenomics", version: "1.0.0", mcpPath: "/mcp", tools };
