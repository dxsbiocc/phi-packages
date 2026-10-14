import type { OmicsTool, ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { compoundResponse, geneDependency, modelSummary, resolveCellLine, usageTerms } from "./query.js";

const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

export const tools: OmicsTool[] = [
  {
    name: "usage_terms",
    title: "Usage terms",
    description: "Return usage terms, licenses, and citations for the pharmacogenomics tables.",
    inputSchema: z.object({}),
    call: async () => usageTerms(),
  },
  {
    name: "resolve_cell_line",
    title: "Resolve a cell line",
    description: "Resolve a cell line by DepMap ACH id, Sanger SIDM, COSMIC id, Cellosaurus CVCL, or name. Ambiguous aliases are returned together and are not collapsed to one model.",
    inputSchema: z.object({ query: z.string().min(1), ...pageArgs }),
    call: async (_client, args) => resolveCellLine(args),
  },
  {
    name: "gene_dependency",
    title: "Gene dependency across models",
    description: "Chronos gene-effect scores for one gene across models, most dependent models first. Accepts an HGNC symbol or ncbigene id.",
    inputSchema: z.object({ gene: z.string().min(1), dataset: z.string().optional(), ...pageArgs }),
    call: async (_client, args) => geneDependency(args),
  },
  {
    name: "compound_response",
    title: "Compound response profile",
    description: "Drug-response rows for a compound across datasets. Accepts a PubChem CID, ChEMBL id, InChIKey, BRD id, GDSC drug id, or name. Rows are ordered by value ascending.",
    inputSchema: z.object({ compound: z.string().min(1), dataset: z.string().optional(), metric: z.string().optional(), ...pageArgs }),
    call: async (_client, args) => compoundResponse(args),
  },
  {
    name: "model_summary",
    title: "Model multi-omics summary",
    description: "Summarize one cell line: identifier crosswalk plus counts of dependency, expression, mutation, copy number, and drug-response rows. Optional gene adds that gene's measurements.",
    inputSchema: z.object({ model: z.string().min(1), gene: z.string().optional() }),
    call: async (_client, args) => modelSummary(args),
  },
];

export const serverModule: ServerModule = {
  name: "pharmacogenomics",
  version: "1.0.0",
  mcpPath: "/mcp",
  tools,
};
