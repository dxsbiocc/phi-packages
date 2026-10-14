import { makePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import { z } from "zod";
import { literatureTools } from "./tools.js";
import { literatureSource } from "./util.js";

const summary = [
  "Phi literature searches public bibliographic APIs and returns one publication record.",
  "Default search is PubMed (MeSH and field tags, newest first). Preprint topics go to Europe PMC. bioRxiv and medRxiv are browsed by date. DOI lookup uses OpenAlex, then Crossref. A 10.1101 DOI uses the bioRxiv/medRxiv API so published_as can carry the journal DOI.",
  "Search abstracts are truncated to 400 characters unless include_abstract is true. Open full text is inlined only for Creative Commons or public-domain XML from Europe PMC; otherwise the response contains links.",
  "Optional NCBI_API_KEY raises the NCBI rate limit. PHI_LITERATURE_EMAIL is the contact address for the Crossref, OpenAlex, and Unpaywall polite pools. S2_API_KEY is optional and only used when sources includes semantic_scholar. No paid API key is required.",
].join(" ");

const usageTool: OmicsTool = {
  name: "usage_terms",
  title: "Usage terms",
  description: "Return usage terms, source routing, and optional environment variables for this literature server.",
  inputSchema: z.object({}),
  call: async () => makePage("publication", [], 0, 1, 0, literatureSource, {
    summary,
    url: "https://github.com/dxsbiocc/phi-packages",
  }),
};

export const tools: OmicsTool[] = [usageTool, ...literatureTools];
export const serverModule: ServerModule = { name: "literature", version: "1.0.0", mcpPath: "/mcp", tools };
