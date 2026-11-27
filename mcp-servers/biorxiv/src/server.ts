import { normalizeDoi, publicationId, type Publication, type Source, type ToolResult } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  failure,
  makePage,
  numberOrNull,
  publication,
  text,
  type JsonClient,
  type OmicsTool,
  type ServerModule,
} from "phi-omics-runtime";
import { searchPublications } from "phi-literature/service";
import { toFailure } from "phi-literature/util";
import { z } from "zod";

export const biorxivSource: Source = {
  id: "biorxiv",
  name: "bioRxiv & medRxiv Preprint Server",
  version: null,
  homepage: "https://www.biorxiv.org/",
  retrieved_at: "1970-01-01T00:00:00.000Z",
};

const USAGE_SUMMARY = [
  "Self-hosted bioRxiv and medRxiv preprint MCP server operated by Phi Packages over the official Cold Spring Harbor Laboratory (CSHL) bioRxiv/medRxiv REST API (api.biorxiv.org) and Europe PMC preprint index.",
  "Supports both full-text/keyword search across bioRxiv and medRxiv (biorxiv_search_preprints) and official CSHL date-window / subject-category browsing, version history inspection, and peer-reviewed journal DOI tracking.",
  "Preprints are not peer-reviewed; check published_as / biorxiv_get_published to verify whether a peer-reviewed journal version has appeared, and respect each preprint's Creative Commons or author license.",
].join(" ");

const pageArgs = {
  offset: z.number().int().nonnegative().optional().describe("Zero-based offset (default 0)."),
  limit: z.number().int().min(1).max(50).optional().describe("Page size (1..50, default 20)."),
};

function rebrand(result: ToolResult): ToolResult {
  return {
    ...result,
    source: {
      ...biorxivSource,
      retrieved_at: result.source.retrieved_at,
    },
  };
}

function cleanPreprintDoi(raw: string): string | null {
  const trimmed = raw.trim().replace(/^doi:/i, "").replace(/^https?:\/\/(?:dx\.)?doi\.org\//i, "");
  return normalizeDoi(trimmed);
}

function utcOffsetDate(days: number): string {
  const d = new Date();
  d.setUTCDate(d.getUTCDate() + days);
  return d.toISOString().slice(0, 10);
}

async function getPreprintDetails(
  client: JsonClient,
  rawDoi: string,
  preferredServer?: "biorxiv" | "medrxiv",
): Promise<ToolResult> {
  const doi = cleanPreprintDoi(rawDoi);
  if (!doi) return failure(biorxivSource, "invalid_input", `Invalid preprint DOI: '${rawDoi}'`, false);

  const servers: Array<"biorxiv" | "medrxiv"> = preferredServer
    ? [preferredServer]
    : ["biorxiv", "medrxiv"];

  for (const srv of servers) {
    const body = asRecord(await client.getJson(`https://api.biorxiv.org/details/${srv}/${doi}/na/json`));
    const collection = asArray(body.collection).map(asRecord);
    if (collection.length === 0) continue;

    const latest = collection[collection.length - 1]!;
    const latestVersion = text(latest.version) ?? "1";
    const publishedDoi = normalizeDoi(text(latest.published));
    const versions = collection.map((v) => ({
      version: text(v.version),
      date: text(v.date),
      category: text(v.category),
      type: text(v.type),
      jatsxml: text(v.jatsxml),
    }));

    const item: Publication = publication({
      id: publicationId({ doi }),
      doi,
      xrefs: publishedDoi ? [`doi:${publishedDoi}`] : [],
      title: text(latest.title) ?? doi,
      authors: (text(latest.authors) ?? "")
        .split(/;\s*/)
        .filter(Boolean)
        .map((name) => ({ name, orcid: null })),
      venue: srv,
      dates: {
        published: null,
        posted: text(latest.date),
        updated: collection.length > 1 ? text(latest.date) : null,
      },
      abstract: text(latest.abstract),
      type: "preprint",
      is_oa: true,
      license: text(latest.license),
      oa_urls: [
        { url: `https://www.${srv}.org/content/${doi}v${latestVersion}`, content_type: "landing", version: latestVersion },
        { url: `https://www.${srv}.org/content/${doi}v${latestVersion}.full.pdf`, content_type: "pdf", version: latestVersion },
        ...(text(latest.jatsxml)
          ? [{ url: text(latest.jatsxml)!, content_type: "xml" as const, version: latestVersion }]
          : []),
      ],
      published_as: publishedDoi ? `doi:${publishedDoi}` : null,
      preprint_of: null,
      source_extra: {
        provider: srv,
        category: text(latest.category),
        version: latestVersion,
        version_count: collection.length,
        versions,
        author_corresponding: text(latest.author_corresponding),
        author_corresponding_institution: text(latest.author_corresponding_institution),
        funder: latest.funder ?? null,
      },
    });

    return makePage("publication", [item], 0, 1, 1, biorxivSource, {
      summary: USAGE_SUMMARY,
      url: `https://www.${srv}.org/content/${doi}v${latestVersion}`,
    });
  }

  return failure(biorxivSource, "not_found", `No bioRxiv or medRxiv preprint found for DOI '${doi}'.`, false);
}

async function getPublishedLinks(
  client: JsonClient,
  args: {
    doi?: string;
    server?: "biorxiv" | "medrxiv";
    date_from?: string;
    date_to?: string;
    offset?: number;
    limit?: number;
  },
): Promise<ToolResult> {
  const offset = args.offset ?? 0;
  const limit = args.limit ?? 20;

  if (args.doi?.trim()) {
    const doi = cleanPreprintDoi(args.doi);
    if (!doi) return failure(biorxivSource, "invalid_input", `Invalid DOI '${args.doi}'`, false);
    const servers: Array<"biorxiv" | "medrxiv"> = args.server ? [args.server] : ["biorxiv", "medrxiv"];
    for (const srv of servers) {
      const body = asRecord(await client.getJson(`https://api.biorxiv.org/pubs/${srv}/${doi}/na/json`));
      const rows = asArray(body.collection).map(asRecord);
      if (rows.length === 0) continue;
      const hit = rows[0]!;
      const pubDoi = normalizeDoi(text(hit.published_doi));
      const preDoi = normalizeDoi(text(hit.biorxiv_doi) ?? text(hit.preprint_doi) ?? doi);
      if (!pubDoi) continue;
      const item = publication({
        id: publicationId({ doi: pubDoi }),
        doi: pubDoi,
        xrefs: preDoi ? [`doi:${preDoi}`] : [],
        title: text(hit.preprint_title) ?? pubDoi,
        authors: (text(hit.preprint_authors) ?? "")
          .split(/;\s*/)
          .filter(Boolean)
          .map((name) => ({ name, orcid: null })),
        venue: text(hit.published_journal),
        dates: {
          published: text(hit.published_date),
          posted: text(hit.preprint_date),
          updated: null,
        },
        abstract: text(hit.preprint_abstract),
        type: "journal-article",
        is_oa: null,
        license: null,
        oa_urls: [{ url: `https://doi.org/${pubDoi}`, content_type: "landing", version: null }],
        published_as: null,
        preprint_of: preDoi ? `doi:${preDoi}` : null,
        source_extra: {
          provider: srv,
          preprint_doi: preDoi,
          published_doi: pubDoi,
          published_journal: text(hit.published_journal),
          published_citation_count: numberOrNull(hit.published_citation_count),
          preprint_category: text(hit.preprint_category),
        },
      });
      return makePage("publication", [item], 0, 1, 1, biorxivSource, {
        summary: USAGE_SUMMARY,
        url: `https://doi.org/${pubDoi}`,
      });
    }
    return failure(
      biorxivSource,
      "not_found",
      `No peer-reviewed journal publication recorded yet on bioRxiv/medRxiv for '${doi}'.`,
      false,
    );
  }

  const srv = args.server ?? "biorxiv";
  const dateTo = args.date_to ?? utcOffsetDate(0);
  const dateFrom = args.date_from ?? utcOffsetDate(-14);
  const body = asRecord(await client.getJson(`https://api.biorxiv.org/pubs/${srv}/${dateFrom}/${dateTo}/${offset}/json`));
  const msg = asRecord(asArray(body.messages)[0]);
  const total = numberOrNull(msg.total) ?? numberOrNull(msg.count) ?? 0;
  const rows = asArray(body.collection).map(asRecord).slice(0, limit);

  const items: Publication[] = [];
  for (const hit of rows) {
    const pubDoi = normalizeDoi(text(hit.published_doi));
    const preDoi = normalizeDoi(text(hit.biorxiv_doi) ?? text(hit.preprint_doi));
    if (!pubDoi) continue;
    items.push(
      publication({
        id: publicationId({ doi: pubDoi }),
        doi: pubDoi,
        xrefs: preDoi ? [`doi:${preDoi}`] : [],
        title: text(hit.preprint_title) ?? pubDoi,
        authors: (text(hit.preprint_authors) ?? "")
          .split(/;\s*/)
          .filter(Boolean)
          .map((name) => ({ name, orcid: null })),
        venue: text(hit.published_journal),
        dates: {
          published: text(hit.published_date),
          posted: text(hit.preprint_date),
          updated: null,
        },
        abstract: null,
        type: "journal-article",
        is_oa: null,
        license: null,
        oa_urls: [{ url: `https://doi.org/${pubDoi}`, content_type: "landing", version: null }],
        published_as: null,
        preprint_of: preDoi ? `doi:${preDoi}` : null,
        source_extra: {
          provider: srv,
          preprint_doi: preDoi,
          published_doi: pubDoi,
          published_journal: text(hit.published_journal),
          preprint_category: text(hit.preprint_category),
        },
      }),
    );
  }

  return makePage("publication", items, offset, limit, total || items.length, biorxivSource, {
    summary: USAGE_SUMMARY,
    url: `https://api.biorxiv.org/pubs/${srv}/${dateFrom}/${dateTo}/${offset}`,
  });
}

export const tools: OmicsTool[] = [
  {
    name: "usage_terms",
    title: "Usage terms",
    description: "Return usage terms, preprint licensing notes, and routing details for bioRxiv and medRxiv.",
    inputSchema: z.object({}),
    call: async () =>
      makePage("publication", [], 0, 1, 0, biorxivSource, {
        summary: USAGE_SUMMARY,
        url: "https://api.biorxiv.org/",
      }),
  },
  {
    name: "biorxiv_search_preprints",
    title: "Search bioRxiv & medRxiv preprints by keyword",
    description:
      "Search bioRxiv and medRxiv preprints by topic, gene, disease, method, or author via Europe PMC's indexed preprint corpus, with optional server ('biorxiv' or 'medrxiv'), subject category, and posting date filters.",
    inputSchema: z.object({
      query: z.string().min(1).describe("Keyword, gene, disease, author, or boolean query string."),
      server: z.enum(["biorxiv", "medrxiv"]).optional().describe("Restrict search to 'biorxiv' or 'medrxiv' (omit to search both)."),
      category: z.string().optional().describe("Optional subject category keyword, e.g. 'bioinformatics', 'genomics', 'cancer biology'."),
      date_from: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional().describe("Inclusive start posting date (YYYY, YYYY-MM, or YYYY-MM-DD)."),
      date_to: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional().describe("Inclusive end posting date."),
      sort: z.enum(["date", "relevance"]).optional().describe("Sort by 'date' (newest first, default) or 'relevance'."),
      include_abstract: z.boolean().optional().describe("Include full abstract instead of 400-character preview (default false)."),
      ...pageArgs,
    }),
    call: async (client, args) => {
      try {
        const parsed = z
          .object({
            query: z.string().min(1),
            server: z.enum(["biorxiv", "medrxiv"]).optional(),
            category: z.string().optional(),
            date_from: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional(),
            date_to: z.string().regex(/^\d{4}(-\d{2}(-\d{2})?)?$/).optional(),
            sort: z.enum(["date", "relevance"]).optional(),
            include_abstract: z.boolean().optional(),
            offset: z.number().int().nonnegative().optional(),
            limit: z.number().int().min(1).max(50).optional(),
          })
          .safeParse(args);
        if (!parsed.success) {
          return failure(biorxivSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
        }
        const effectiveQuery = parsed.data.server
          ? parsed.data.query
          : `(${parsed.data.query}) AND (JOURNAL:"bioRxiv" OR JOURNAL:"medRxiv")`;
        return rebrand(
          await searchPublications(client, {
            query: effectiveQuery,
            sources: ["europepmc"],
            type: "preprint",
            server: parsed.data.server,
            category: parsed.data.category,
            date_from: parsed.data.date_from,
            date_to: parsed.data.date_to,
            sort: parsed.data.sort,
            include_abstract: parsed.data.include_abstract,
            offset: parsed.data.offset,
            limit: parsed.data.limit,
          }),
        );
      } catch (err) {
        return rebrand(toFailure(err));
      }
    },
  },
  {
    name: "biorxiv_browse_recent",
    title: "Browse recent bioRxiv or medRxiv preprints by date",
    description:
      "Browse newly posted or revised preprints on bioRxiv or medRxiv over a date window via the official CSHL api.biorxiv.org/details endpoint (defaults to the last 7 days), optionally filtered by subject category.",
    inputSchema: z.object({
      server: z.enum(["biorxiv", "medrxiv"]).optional().describe("Preprint server to browse ('biorxiv' default, or 'medrxiv')."),
      date_from: z.string().regex(/^\d{4}-\d{2}-\d{2}$/).optional().describe("Start date YYYY-MM-DD (defaults to 7 days ago)."),
      date_to: z.string().regex(/^\d{4}-\d{2}-\d{2}$/).optional().describe("End date YYYY-MM-DD (defaults to today UTC)."),
      category: z.string().optional().describe("Exact CSHL subject category filter, e.g. 'bioinformatics', 'genomics', 'immunology'."),
      include_abstract: z.boolean().optional().describe("Include full abstracts (default false)."),
      ...pageArgs,
    }),
    call: async (client, args) => {
      try {
        const parsed = z
          .object({
            server: z.enum(["biorxiv", "medrxiv"]).optional(),
            date_from: z.string().regex(/^\d{4}-\d{2}-\d{2}$/).optional(),
            date_to: z.string().regex(/^\d{4}-\d{2}-\d{2}$/).optional(),
            category: z.string().optional(),
            include_abstract: z.boolean().optional(),
            offset: z.number().int().nonnegative().optional(),
            limit: z.number().int().min(1).max(50).optional(),
          })
          .safeParse(args);
        if (!parsed.success) {
          return failure(biorxivSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
        }
        const srv = parsed.data.server ?? "biorxiv";
        return rebrand(
          await searchPublications(client, {
            sources: [srv],
            type: "preprint",
            server: srv,
            date_from: parsed.data.date_from,
            date_to: parsed.data.date_to,
            category: parsed.data.category,
            include_abstract: parsed.data.include_abstract,
            offset: parsed.data.offset,
            limit: parsed.data.limit,
          }),
        );
      } catch (err) {
        return rebrand(toFailure(err));
      }
    },
  },
  {
    name: "biorxiv_get_preprint",
    title: "Get bioRxiv / medRxiv preprint details and versions",
    description:
      "Retrieve complete preprint metadata, full abstract, version history (v1..vN), corresponding author/institution, license, PDF/JATS XML links, and linked peer-reviewed journal DOI (published_as) for one bioRxiv or medRxiv DOI.",
    inputSchema: z.object({
      doi: z.string().min(1).describe("Preprint DOI, e.g. '10.1101/2024.05.01.592000' or 'doi:10.1101/...'."),
      server: z.enum(["biorxiv", "medrxiv"]).optional().describe("Optional server hint ('biorxiv' or 'medrxiv')."),
    }),
    call: async (client, args) => {
      try {
        const parsed = z
          .object({
            doi: z.string().min(1),
            server: z.enum(["biorxiv", "medrxiv"]).optional(),
          })
          .safeParse(args);
        if (!parsed.success) {
          return failure(biorxivSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
        }
        return await getPreprintDetails(client, parsed.data.doi, parsed.data.server);
      } catch (err) {
        return rebrand(toFailure(err));
      }
    },
  },
  {
    name: "biorxiv_get_published",
    title: "Resolve peer-reviewed journal versions of preprints",
    description:
      "Check whether a bioRxiv or medRxiv preprint has been published in a peer-reviewed journal (by preprint DOI), or browse recent preprint-to-journal publications over a date window via api.biorxiv.org/pubs.",
    inputSchema: z.object({
      doi: z.string().optional().describe("Optional preprint DOI (e.g. '10.1101/...') to check journal publication status."),
      server: z.enum(["biorxiv", "medrxiv"]).optional().describe("Preprint server ('biorxiv' default, or 'medrxiv')."),
      date_from: z.string().regex(/^\d{4}-\d{2}-\d{2}$/).optional().describe("Start date YYYY-MM-DD when browsing published preprints."),
      date_to: z.string().regex(/^\d{4}-\d{2}-\d{2}$/).optional().describe("End date YYYY-MM-DD when browsing published preprints."),
      ...pageArgs,
    }),
    call: async (client, args) => {
      try {
        const parsed = z
          .object({
            doi: z.string().optional(),
            server: z.enum(["biorxiv", "medrxiv"]).optional(),
            date_from: z.string().regex(/^\d{4}-\d{2}-\d{2}$/).optional(),
            date_to: z.string().regex(/^\d{4}-\d{2}-\d{2}$/).optional(),
            offset: z.number().int().nonnegative().optional(),
            limit: z.number().int().min(1).max(50).optional(),
          })
          .safeParse(args);
        if (!parsed.success) {
          return failure(biorxivSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
        }
        return await getPublishedLinks(client, parsed.data);
      } catch (err) {
        return rebrand(toFailure(err));
      }
    },
  },
];

export const serverModule: ServerModule = {
  name: "biorxiv",
  version: "1.0.0",
  mcpPath: "/mcp",
  tools,
};
