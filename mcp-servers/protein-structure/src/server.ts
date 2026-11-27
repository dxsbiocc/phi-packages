import { sourceOf, usagePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { structureTools } from "./sources.js";

const summary = [
  "UniProtKB: protein sequences, domains, functional sites, disease annotations, and PDB cross-references via REST API (CC-BY 4.0). https://www.uniprot.org/",
  "AlphaFold DB: predicted protein 3D structure models (mmCIF, PDB, BinaryCIF, PAE) and per-residue pLDDT confidence scores (CC-BY 4.0). https://alphafold.ebi.ac.uk/",
].join("\n");

const usageTool: OmicsTool = {
  name: "usage_terms",
  title: "Usage terms",
  description: "Return usage terms and license attribution for this protein structure server.",
  inputSchema: z.object({}),
  call: async () => usagePage(sourceOf("protein-structure", "Protein Structure", "https://github.com/dxsbiocc/phi-packages"), summary),
};

export const tools: OmicsTool[] = [usageTool, ...structureTools];
export const serverModule: ServerModule = { name: "protein-structure", version: "1.0.0", mcpPath: "/mcp", tools };
