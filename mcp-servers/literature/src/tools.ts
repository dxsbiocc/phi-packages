import type { ToolResult } from "phi-omics-schemas";
import { failure, type OmicsTool } from "phi-omics-runtime";
import { z } from "zod";
import { getFullText, getPublication, relatedPublications, resolveReferences, searchPublications } from "./service.js";
import { literatureSource, LITERATURE_SOURCES, toFailure } from "./util.js";

const pageArgs = {
  offset: z.number().int().nonnegative().optional().describe("Zero-based offset. Default 0."),
  limit: z.number().int().min(1).max(50).optional().describe("Page size. Default 20, maximum 50."),
};
const sourceArg = z.enum(LITERATURE_SOURCES);

function issue(error: z.ZodError): string {
  return error.issues[0]?.message ?? "invalid input";
}

async function run(fn: () => Promise<ToolResult>): Promise<ToolResult> {
  try {
    return await fn();
  } catch (error) {
    return toFailure(error);
  }
}

export const literatureTools: OmicsTool[] = [
  {
    name: "search_publications",
    title: "Search publications",
    description: "Find papers on a topic, gene, or disease, newest first. Default source is PubMed, so MeSH terms and field tags such as [ti], [au], and [pdat] work. Abstracts are truncated to 400 characters unless include_abstract is true. Set type to preprint to search Europe PMC preprints, including a bioRxiv or medRxiv journal filter via server. Omit query and pass server plus date_from/date_to to browse bioRxiv or medRxiv by date. sources overrides that routing and accepts at most 3 databases; the server does not query every API.",
    inputSchema: z.object({
      query: z.string().optional().describe("PubMed query, or free text for the selected source. Required except for bioRxiv/medRxiv date browse."),
      sources: z.array(sourceArg).max(3).optional().describe("Override routing. pubmed, europepmc, biorxiv, medrxiv, crossref, openalex, semantic_scholar, arxiv."),
      date_from: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional().describe("Inclusive start date, YYYY, YYYY-MM, or YYYY-MM-DD."),
      date_to: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional().describe("Inclusive end date."),
      type: z.enum(["journal-article", "preprint", "review", "any"]).optional().describe("Publication type filter. preprint switches the default source to Europe PMC, or bioRxiv when query is empty."),
      server: z.enum(["biorxiv", "medrxiv"]).optional().describe("Preprint server for date browse, or a journal filter on Europe PMC preprint search."),
      category: z.string().optional().describe("bioRxiv/medRxiv subject category, matched within the fetched date page. On Europe PMC it is added as a query term."),
      sort: z.enum(["date", "relevance"]).optional().describe("Default date, newest first."),
      include_abstract: z.boolean().optional().describe("Default false. Search hits include a 400-character abstract unless this is true."),
      ...pageArgs,
    }),
    call: (client, args) => run(async () => {
      const parsed = z.object({
        query: z.string().optional(),
        sources: z.array(sourceArg).max(3).optional(),
        date_from: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional(),
        date_to: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional(),
        type: z.enum(["journal-article", "preprint", "review", "any"]).optional(),
        server: z.enum(["biorxiv", "medrxiv"]).optional(),
        category: z.string().optional(),
        sort: z.enum(["date", "relevance"]).optional(),
        include_abstract: z.boolean().optional(),
        offset: z.number().int().nonnegative().optional(),
        limit: z.number().int().min(1).max(50).optional(),
      }).safeParse(args);
      if (!parsed.success) return failure(literatureSource, "invalid_input", issue(parsed.error), false);
      return searchPublications(client, parsed.data);
    }),
  },
  {
    name: "get_publication",
    title: "Get one publication",
    description: "Look up one paper by DOI, PMID, PMCID, arXiv id, OpenAlex id, Semantic Scholar id, or a title fragment or citation string. A DOI is the primary id even when there is no PMID. A 10.1101 DOI returns the preprint and published_as when bioRxiv or medRxiv reports a journal DOI. The abstract is included unless include_abstract is false. Identifier crosswalk and license live on this record.",
    inputSchema: z.object({
      id: z.string().optional().describe("DOI, PMID, PMCID, arXiv, openalex:W…, or s2: id. A bare DOI or PMID is accepted."),
      citation: z.string().optional().describe("Title fragment or citation string when you do not have an id. Returns the closest matches."),
      include_abstract: z.boolean().optional().describe("Default true."),
    }),
    call: (client, args) => run(async () => {
      const parsed = z.object({
        id: z.string().optional(),
        citation: z.string().optional(),
        include_abstract: z.boolean().optional(),
      }).safeParse(args);
      if (!parsed.success) return failure(literatureSource, "invalid_input", issue(parsed.error), false);
      return getPublication(client, parsed.data);
    }),
  },
  {
    name: "get_fulltext",
    title: "Get open full text",
    description: "Return the best open text for one paper, for reading or quoting methods. A DOI or PMID with no PMCID is checked in Europe PMC, then JATS XML is inlined only when the license is Creative Commons or public domain, clipped to max_chars (default 4000). Otherwise the tool returns open HTML, PDF, and XML links and does not paste the body. section selects abstract, methods, results, discussion, or all.",
    inputSchema: z.object({
      id: z.string().describe("DOI, PMID, PMCID, or other publication id."),
      section: z.enum(["abstract", "methods", "results", "discussion", "all"]).optional().describe("Default all. abstract can be returned from metadata without a full-text license."),
      max_chars: z.number().int().min(200).max(20000).optional().describe("Maximum characters of inlined text. Default 4000."),
    }),
    call: (client, args) => run(async () => {
      const parsed = z.object({
        id: z.string(),
        section: z.enum(["abstract", "methods", "results", "discussion", "all"]).optional(),
        max_chars: z.number().int().min(200).max(20000).optional(),
      }).safeParse(args);
      if (!parsed.success) return failure(literatureSource, "invalid_input", issue(parsed.error), false);
      return getFullText(client, parsed.data);
    }),
  },
  {
    name: "related_publications",
    title: "Related publications",
    description: "Expand from one paper. relation=similar uses PubMed related citations when a PMID is known, otherwise OpenAlex related works. citations and references use OpenAlex, with OpenCitations as a fallback. One provider is queried, not every citation API.",
    inputSchema: z.object({
      id: z.string().describe("DOI, PMID, or other publication id."),
      relation: z.enum(["similar", "citations", "references"]).optional().describe("Default similar."),
      ...pageArgs,
    }),
    call: (client, args) => run(async () => {
      const parsed = z.object({
        id: z.string(),
        relation: z.enum(["similar", "citations", "references"]).optional(),
        offset: z.number().int().nonnegative().optional(),
        limit: z.number().int().min(1).max(50).optional(),
      }).safeParse(args);
      if (!parsed.success) return failure(literatureSource, "invalid_input", issue(parsed.error), false);
      return relatedPublications(client, parsed.data);
    }),
  },
  {
    name: "resolve_references",
    title: "Resolve a reference list",
    description: "Batch-resolve up to 25 references. Each item may be a DOI, PMID, PMCID, arXiv id, or a citation string. DOIs and PMIDs are looked up in batches. Unresolved lines are returned in source_extra.unresolved instead of failing the whole list. Abstracts are omitted unless include_abstract is true.",
    inputSchema: z.object({
      references: z.array(z.string().min(1)).min(1).max(25).describe("Reference strings, in order."),
      include_abstract: z.boolean().optional().describe("Default false."),
    }),
    call: (client, args) => run(async () => {
      const parsed = z.object({
        references: z.array(z.string().min(1)).min(1).max(25),
        include_abstract: z.boolean().optional(),
      }).safeParse(args);
      if (!parsed.success) return failure(literatureSource, "invalid_input", issue(parsed.error), false);
      return resolveReferences(client, parsed.data);
    }),
  },
];
