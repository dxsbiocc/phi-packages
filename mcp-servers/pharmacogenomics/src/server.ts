import type { OmicsTool, ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { biomarkerAssociation } from "./biomarker.js";
import { geneCopyNumber, geneExpression, geneMutations, listTissues } from "./omics.js";
import { compoundResponse, dataDir, geneDependency, modelSummary, resolveCellLine, usageTerms } from "./query.js";

const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

const tissueArg = z.string().min(1).optional().describe("Limit to models of one tissue, such as Lung. Matching ignores case. Call list_tissues for the values.");
const orderArg = z.enum(["asc", "desc"]).optional().describe("Sort direction by value. Default desc.");

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
    description: "Chronos gene-effect scores for one gene across models, most dependent models first. Accepts an HGNC symbol or ncbigene id. Optional tissue limits the models.",
    inputSchema: z.object({ gene: z.string().min(1), dataset: z.string().optional(), tissue: tissueArg, ...pageArgs }),
    call: async (_client, args) => geneDependency(args),
  },
  {
    name: "compound_response",
    title: "Compound response profile",
    description: "Drug-response rows for a compound or a drug combination. Accepts a PubChem CID, ChEMBL id, InChIKey, BRD id, GDSC drug id, NSC id, combination id, or name. Omit dataset to compare the compound across every loaded dataset. Single-agent rows are the default. combinations=include adds combination rows, and combinations=only returns those rows. A combination id or an 'Anchor + Library' name returns combination rows either way. Optional tissue limits the models. Rows are ordered by value ascending.",
    inputSchema: z.object({
      compound: z.string().min(1),
      dataset: z.string().optional(),
      metric: z.string().optional(),
      tissue: tissueArg,
      combinations: z.enum(["exclude", "include", "only"]).optional(),
      ...pageArgs,
    }),
    call: async (_client, args) => compoundResponse(args),
  },
  {
    name: "list_tissues",
    title: "List model tissues",
    description: "Tissue values of the mapped cell-line models with the number of models in each. Use a value as the tissue filter of the other tools.",
    inputSchema: z.object({ ...pageArgs }),
    call: async (_client, args) => listTissues(args),
  },
  {
    name: "gene_expression",
    title: "Gene expression across models",
    description: "Expression of one gene across cell-line models (DepMap log2(TPM+1)), highest first by default, with the distribution over all matching models. Accepts an HGNC symbol or ncbigene id. Optional dataset and tissue limit the rows.",
    inputSchema: z.object({ gene: z.string().min(1), dataset: z.string().optional(), tissue: tissueArg, order: orderArg, ...pageArgs }),
    call: async (_client, args) => geneExpression(args),
  },
  {
    name: "gene_copy_number",
    title: "Gene copy number across models",
    description: "Copy number of one gene across cell-line models, highest first by default, with the distribution over all matching models. Accepts an HGNC symbol or ncbigene id. Optional dataset and tissue limit the rows.",
    inputSchema: z.object({ gene: z.string().min(1), dataset: z.string().optional(), tissue: tissueArg, order: orderArg, ...pageArgs }),
    call: async (_client, args) => geneCopyNumber(args),
  },
  {
    name: "gene_mutations",
    title: "Gene mutations across models",
    description: "Somatic mutations of one gene across cell-line models, with how many sequenced models carry a non-silent call and the most common protein changes. Silent variants are left out unless include_silent is true. Filter by consequence (a name such as synonymous_variant also matches inside combined annotations and overrides the silent default), protein_change (p.R175H or R175H), dataset, or tissue.",
    inputSchema: z.object({
      gene: z.string().min(1),
      dataset: z.string().optional(),
      tissue: tissueArg,
      consequence: z.string().min(1).optional(),
      protein_change: z.string().min(1).optional(),
      include_silent: z.boolean().optional(),
      ...pageArgs,
    }),
    call: async (_client, args) => geneMutations(args),
  },
  {
    name: "biomarker_association",
    title: "Compound response vs a gene feature",
    description: "Test whether a single-agent compound's response differs with one gene's mutation status (Mann-Whitney) or tracks its expression, copy number, or CRISPR dependency (Spearman), one row per response dataset and metric, most significant first. Statistics are descriptive and uncorrected for multiple testing. Set tissue to avoid tissue confounding.",
    inputSchema: z.object({
      compound: z.string().min(1),
      gene: z.string().min(1),
      feature: z.enum(["mutation", "expression", "copy_number", "dependency"]),
      dataset: z.string().optional().describe("Drug-response dataset, such as gdsc1."),
      metric: z.string().optional(),
      feature_dataset: z.string().optional().describe("Dataset of the gene feature, such as depmap_expression_24q2."),
      tissue: tissueArg,
      min_group: z.number().int().min(2).max(100).optional().describe("Smallest mutant or wild-type group to test. Default 3."),
      ...pageArgs,
    }),
    call: async (_client, args) => biomarkerAssociation(args),
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
  ready: async () => dataDir() == null ? "PHI_PGX_DATA_DIR is not set to an existing processed parquet directory" : null,
};
