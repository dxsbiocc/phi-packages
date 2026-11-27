import { sourceOf, usagePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { chemistryTools } from "./sources.js";

const summary = [
  "ChEMBL: bioactive drug-like small molecules, physicochemical properties, IC50/Ki/Kd bioactivities, and mechanisms of action via EMBL-EBI REST API (CC BY-SA 3.0). https://www.ebi.ac.uk/chembl/",
  "OpenFDA: FDA-approved prescription drug labels and FAERS spontaneous adverse event reaction counts via api.fda.gov (FAERS counts reflect reporting volume, not incidence rates). https://open.fda.gov/",
].join("\n");

const usageTool: OmicsTool = {
  name: "usage_terms",
  title: "Usage terms",
  description: "Return usage terms and caveats for this drug chemistry server.",
  inputSchema: z.object({}),
  call: async () => usagePage(sourceOf("drug-chemistry", "Drug Chemistry", "https://github.com/dxsbiocc/phi-packages"), summary),
};

export const tools: OmicsTool[] = [usageTool, ...chemistryTools];
export const serverModule: ServerModule = { name: "drug-chemistry", version: "1.0.0", mcpPath: "/mcp", tools };
