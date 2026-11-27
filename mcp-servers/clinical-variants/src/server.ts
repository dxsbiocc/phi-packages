import { sourceOf, usagePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { variantTools } from "./sources.js";

const summary = [
  "ClinVar: NCBI variant clinical significance summaries via E-utilities (optional NCBI_API_KEY raises rate limit). https://www.ncbi.nlm.nih.gov/clinvar/",
  "gnomAD: population allele frequencies and gene constraint metrics via public GraphQL API (gnomAD v4 GRCh38, v3 GRCh38, v2.1 GRCh37). https://gnomad.broadinstitute.org/",
  "CIViC: Clinical Interpretation of Variants in Cancer variant records and clinical evidence items via GraphQL v2 API. https://civicdb.org/",
].join("\n");

const usageTool: OmicsTool = {
  name: "usage_terms",
  title: "Usage terms",
  description: "Return usage terms and data source attribution for this clinical variants server.",
  inputSchema: z.object({}),
  call: async () => usagePage(sourceOf("clinical-variants", "Clinical Variants", "https://github.com/dxsbiocc/phi-packages"), summary),
};

export const tools: OmicsTool[] = [usageTool, ...variantTools];
export const serverModule: ServerModule = { name: "clinical-variants", version: "1.0.0", mcpPath: "/mcp", tools };
