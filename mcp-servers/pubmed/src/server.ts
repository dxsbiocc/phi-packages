import type { Source, ToolResult } from "phi-omics-schemas";
import { failure, makePage, type OmicsTool, type ServerModule } from "phi-omics-runtime";
import {
  getEntityAnnotations,
  getFullText,
  getPublication,
  relatedPublications,
  resolveReferences,
  searchPublications,
} from "phi-literature/service";
import { toFailure } from "phi-literature/util";
import { z } from "zod";

export const pubmedSource: Source = {
  id: "pubmed",
  name: "NCBI PubMed, PMC & PubTator 3.0",
  version: null,
  homepage: "https://pubmed.ncbi.nlm.nih.gov/",
  retrieved_at: "1970-01-01T00:00:00.000Z",
};

const USAGE_SUMMARY = [
  "Self-hosted PubMed, PubMed Central (PMC), and PubTator 3.0 MCP server operated by Phi Packages over official NCBI E-utilities, Europe PMC open-access JATS XML, and PubTator 3.0 REST APIs.",
  "Supports full PubMed query syntax (MeSH headings and field tags such as [ti], [au], [ta], [pdat]) as well as PubTator 3.0 entity queries (@GENE_, @DISEASE_, @CHEMICAL_, @VARIANT_, relations:).",
  "Search abstracts are clipped to 400 characters unless include_abstract is true. Full-text sections are inlined only when a Creative Commons or public-domain license permits redistribution; otherwise open links are returned.",
  "Follow NCBI usage guidelines and cite original articles and NCBI/NLM resources.",
].join(" ");

const pageArgs = {
  offset: z.number().int().nonnegative().optional().describe("Zero-based result offset (default 0)."),
  limit: z.number().int().min(1).max(50).optional().describe("Maximum records to return (1..50, default 20)."),
};

function rebrand(result: ToolResult): ToolResult {
  return {
    ...result,
    source: {
      ...pubmedSource,
      retrieved_at: result.source.retrieved_at,
    },
  };
}

async function safeCall(fn: () => Promise<ToolResult>): Promise<ToolResult> {
  try {
    return rebrand(await fn());
  } catch (err) {
    return rebrand(toFailure(err));
  }
}

export const tools: OmicsTool[] = [
  {
    name: "usage_terms",
    title: "Usage terms",
    description: "Return usage terms, NCBI E-utilities / PMC / PubTator 3.0 coverage, and rate-limit guidance.",
    inputSchema: z.object({}),
    call: async () =>
      makePage("publication", [], 0, 1, 0, pubmedSource, {
        summary: USAGE_SUMMARY,
        url: "https://www.ncbi.nlm.nih.gov/home/about/policies/",
      }),
  },
  {
    name: "pubmed_search_articles",
    title: "Search PubMed biomedical articles",
    description:
      "Search PubMed/MEDLINE biomedical literature using PubMed query syntax (MeSH terms, [ti], [au], [ta], [dp]) or PubTator 3.0 entity queries (@GENE_TP53, @DISEASE_, @CHEMICAL_, @VARIANT_, relations:). Supports date windows, review filtering, and relevance or date sorting.",
    inputSchema: z.object({
      query: z.string().min(1).describe("PubMed search query (supports MeSH and field tags) or PubTator 3.0 entity query."),
      date_from: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional().describe("Inclusive start publication date (YYYY, YYYY-MM, or YYYY-MM-DD)."),
      date_to: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional().describe("Inclusive end publication date (YYYY, YYYY-MM, or YYYY-MM-DD)."),
      article_type: z.enum(["any", "journal-article", "review"]).optional().describe("Optional publication type filter (default 'any')."),
      sort: z.enum(["date", "relevance"]).optional().describe("Sort order: 'date' (newest first, default) or 'relevance'."),
      include_abstract: z.boolean().optional().describe("Include full abstracts in search hits instead of 400-character previews (default false)."),
      ...pageArgs,
    }),
    call: (client, args) =>
      safeCall(async () => {
        const parsed = z
          .object({
            query: z.string().min(1),
            date_from: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional(),
            date_to: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional(),
            article_type: z.enum(["any", "journal-article", "review"]).optional(),
            sort: z.enum(["date", "relevance"]).optional(),
            include_abstract: z.boolean().optional(),
            offset: z.number().int().nonnegative().optional(),
            limit: z.number().int().min(1).max(50).optional(),
          })
          .safeParse(args);
        if (!parsed.success) {
          return failure(pubmedSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
        }
        const isPubtator = /@(?:GENE|DISEASE|CHEMICAL|VARIANT|SPECIES|CELLLINE)_|\brelations:/i.test(parsed.data.query);
        return searchPublications(client, {
          query: parsed.data.query,
          sources: [isPubtator ? "pubtator3" : "pubmed"],
          date_from: parsed.data.date_from,
          date_to: parsed.data.date_to,
          type: parsed.data.article_type,
          sort: parsed.data.sort,
          include_abstract: parsed.data.include_abstract,
          offset: parsed.data.offset,
          limit: parsed.data.limit,
        });
      }),
  },
  {
    name: "pubmed_get_article",
    title: "Get PubMed article metadata and abstract",
    description:
      "Retrieve complete bibliographic metadata, structured abstract, ORCID-linked authors, MeSH/publication details, PMID/PMCID/DOI identifier crosswalks, and open-access URLs for one article by PMID, PMCID, DOI, or citation string.",
    inputSchema: z.object({
      id: z.string().optional().describe("PMID (e.g. '38000001' or 'pmid:38000001'), PMCID ('PMC...'), or DOI ('10.1038/...')."),
      citation: z.string().optional().describe("Title fragment or formatted citation string when no identifier is known."),
      include_abstract: z.boolean().optional().describe("Include full abstract text (default true)."),
    }),
    call: (client, args) =>
      safeCall(async () => {
        const parsed = z
          .object({
            id: z.string().optional(),
            citation: z.string().optional(),
            include_abstract: z.boolean().optional(),
          })
          .safeParse(args);
        if (!parsed.success) {
          return failure(pubmedSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
        }
        return getPublication(client, parsed.data);
      }),
  },
  {
    name: "pubmed_get_fulltext",
    title: "Get open-access PMC full text",
    description:
      "Fetch open-access full text (or a targeted section: abstract, methods, results, discussion) from PubMed Central / Europe PMC JATS XML when Creative Commons or public-domain licensed, or open HTML/PDF links when redistribution is restricted.",
    inputSchema: z.object({
      id: z.string().min(1).describe("PMID, PMCID, or DOI of the target article."),
      section: z.enum(["all", "abstract", "methods", "results", "discussion"]).optional().describe("Target section to extract (default 'all')."),
      max_chars: z.number().int().min(200).max(20000).optional().describe("Maximum characters of inlined full text (default 4000)."),
    }),
    call: (client, args) =>
      safeCall(async () => {
        const parsed = z
          .object({
            id: z.string().min(1),
            section: z.enum(["all", "abstract", "methods", "results", "discussion"]).optional(),
            max_chars: z.number().int().min(200).max(20000).optional(),
          })
          .safeParse(args);
        if (!parsed.success) {
          return failure(pubmedSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
        }
        return getFullText(client, parsed.data);
      }),
  },
  {
    name: "pubmed_related_articles",
    title: "Find related, citing, or referenced articles",
    description:
      "Discover neighbor articles from a seed PMID or DOI: 'similar' uses NCBI PubMed neighbor_score related citations; 'citations' and 'references' traverse open citation graphs.",
    inputSchema: z.object({
      id: z.string().min(1).describe("Seed PMID, PMCID, or DOI."),
      relation: z.enum(["similar", "citations", "references"]).optional().describe("Relationship type (default 'similar')."),
      ...pageArgs,
    }),
    call: (client, args) =>
      safeCall(async () => {
        const parsed = z
          .object({
            id: z.string().min(1),
            relation: z.enum(["similar", "citations", "references"]).optional(),
            offset: z.number().int().nonnegative().optional(),
            limit: z.number().int().min(1).max(50).optional(),
          })
          .safeParse(args);
        if (!parsed.success) {
          return failure(pubmedSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
        }
        return relatedPublications(client, parsed.data);
      }),
  },
  {
    name: "pubmed_lookup_citations",
    title: "Batch-resolve citations and convert PMID/PMCID/DOI",
    description:
      "Batch-resolve up to 25 PMIDs, PMCIDs, DOIs, or raw reference strings into normalized publication records with complete PMID/PMCID/DOI crosswalks.",
    inputSchema: z.object({
      references: z.array(z.string().min(1)).min(1).max(25).describe("List of PMIDs, PMCIDs, DOIs, or citation strings (max 25)."),
      include_abstract: z.boolean().optional().describe("Include abstracts on resolved records (default false)."),
    }),
    call: (client, args) =>
      safeCall(async () => {
        const parsed = z
          .object({
            references: z.array(z.string().min(1)).min(1).max(25),
            include_abstract: z.boolean().optional(),
          })
          .safeParse(args);
        if (!parsed.success) {
          return failure(pubmedSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
        }
        return resolveReferences(client, parsed.data);
      }),
  },
  {
    name: "pubmed_entity_annotations",
    title: "Extract PubTator 3.0 biomedical entities",
    description:
      "Retrieve NCBI PubTator 3.0 BioCJSON entity annotations (genes, diseases, chemicals, variants, species, and cell lines) for an article by PMID, PMCID, or DOI.",
    inputSchema: z.object({
      id: z.string().min(1).describe("PMID, PMCID, or DOI."),
    }),
    call: (client, args) =>
      safeCall(async () => {
        const parsed = z.object({ id: z.string().min(1) }).safeParse(args);
        if (!parsed.success) {
          return failure(pubmedSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
        }
        return getEntityAnnotations(client, parsed.data);
      }),
  },
];

export const serverModule: ServerModule = {
  name: "pubmed",
  version: "1.0.0",
  mcpPath: "/mcp",
  tools,
};
