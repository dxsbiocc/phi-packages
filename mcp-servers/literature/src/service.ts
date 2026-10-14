import type { Publication, ToolResult } from "phi-omics-schemas";
import {
  InputError,
  failure,
  makePage,
  type JsonClient,
} from "phi-omics-runtime";
import { clipSections, parseJats, redistributableLicense } from "./jats.js";
import {
  europePmcXml,
  fetchPubmed,
  getArxiv,
  getBiorxiv,
  getCrossrefWork,
  getEuropePmcByQuery,
  getOpenAlexWork,
  getSemanticScholar,
  hydrateOpenAlex,
  openAlexByFilter,
  openCitationDois,
  pubmedNeighbors,
  searchArxiv,
  searchBiorxiv,
  searchCrossref,
  searchEuropePmc,
  searchOpenAlex,
  searchPubmed,
  searchSemanticScholar,
  unpaywall,
  type HitPage,
  type Query,
} from "./providers.js";
import { combinePublications, mergePublications, present, tryPublication } from "./record.js";
import {
  classifyReference,
  contactEmail,
  dateRank,
  literatureSource,
  providerMeta,
  toFailure,
  type LiteratureSourceName,
  type ParsedId,
} from "./util.js";

export interface SearchInput {
  query?: string;
  sources?: LiteratureSourceName[];
  date_from?: string;
  date_to?: string;
  type?: "journal-article" | "preprint" | "review" | "any";
  server?: "biorxiv" | "medrxiv";
  category?: string;
  sort?: "date" | "relevance";
  include_abstract?: boolean;
  offset?: number;
  limit?: number;
}

interface Loaded {
  publication: Publication;
  provider: LiteratureSourceName;
  openAlexId: string | null;
  referencedWorks: string[];
  relatedWorks: string[];
}

function queryOf(input: SearchInput): Query {
  return {
    query: input.query,
    dateFrom: input.date_from,
    dateTo: input.date_to,
    type: input.type,
    server: input.server,
    category: input.category,
    sort: input.sort ?? "date",
    offset: input.offset ?? 0,
    limit: input.limit ?? 20,
  };
}

function route(input: SearchInput): LiteratureSourceName[] {
  if (input.sources && input.sources.length > 0) {
    if (input.sources.length > 3) throw new InputError("Pass at most 3 sources. The default queries one database chosen for the task.");
    return input.sources;
  }
  if (input.type === "preprint" || input.server) {
    if (!input.query?.trim()) return [input.server ?? "biorxiv"];
    return ["europepmc"];
  }
  return ["pubmed"];
}

async function searchOne(client: JsonClient, source: LiteratureSourceName, query: Query): Promise<HitPage> {
  if (source === "pubmed") return searchPubmed(client, query);
  if (source === "europepmc") return searchEuropePmc(client, query);
  if (source === "openalex") return searchOpenAlex(client, query);
  if (source === "crossref") return searchCrossref(client, query);
  if (source === "semantic_scholar") return searchSemanticScholar(client, query);
  if (source === "arxiv") return searchArxiv(client, query);
  return searchBiorxiv(client, source, query);
}

function show(items: Publication[], fullAbstract: boolean): Publication[] {
  return items.map((item) => present(item, fullAbstract ? 12_000 : 400, fullAbstract ? 30 : 8));
}

export async function searchPublications(client: JsonClient, input: SearchInput): Promise<ToolResult> {
  const sources = route(input);
  const query = queryOf(input);
  const fullAbstract = input.include_abstract === true;
  if (sources.length === 1) {
    const only = sources[0]!;
    const hit = await searchOne(client, only, query);
    const items = show(hit.items, fullAbstract);
    return makePage("publication", items, query.offset, query.limit, hit.total, hit.source, hit.usage, hit.extra);
  }
  const hits: HitPage[] = [];
  const warnings: string[] = [];
  let last: unknown;
  for (const source of sources) {
    try {
      hits.push(await searchOne(client, source, query));
    } catch (error) {
      last = error;
      warnings.push(`${source}: ${error instanceof Error ? error.message : "failed"}`);
    }
  }
  if (hits.length === 0) return toFailure(last);
  const items = show(mergePublications(hits.map((hit) => hit.items)).sort((left, right) => dateRank(right).localeCompare(dateRank(left))).slice(0, query.limit), fullAbstract);
  return makePage("publication", items, query.offset, query.limit, null, literatureSource, {
    summary: "Merged publication records from the requested sources. Duplicates match on DOI, then PMID, then PMCID.",
    url: "https://github.com/dxsbiocc/phi-packages",
  }, { providers: hits.map((hit) => hit.extra.provider), warnings });
}

async function loadParsed(client: JsonClient, parsed: ParsedId): Promise<Loaded | null> {
  if (parsed.kind === "doi" && parsed.value.startsWith("10.1101/")) {
    const publication = await getBiorxiv(client, parsed.value);
    if (publication) return { publication, provider: publication.venue === "medrxiv" ? "medrxiv" : "biorxiv", openAlexId: null, referencedWorks: [], relatedWorks: [] };
  }
  if (parsed.kind === "doi") {
    const openAlex = await getOpenAlexWork(client, `works/https://doi.org/${parsed.value}`);
    if (openAlex) return { publication: openAlex.publication, provider: "openalex", openAlexId: openAlex.openAlexId, referencedWorks: openAlex.referencedWorks, relatedWorks: openAlex.relatedWorks };
    const crossref = await getCrossrefWork(client, parsed.value);
    if (crossref) return { publication: crossref, provider: "crossref", openAlexId: null, referencedWorks: [], relatedWorks: [] };
    return null;
  }
  if (parsed.kind === "pmid") {
    const publication = await getEuropePmcByQuery(client, `EXT_ID:${parsed.value} AND SRC:MED`) ?? (await fetchPubmed(client, [parsed.value]))[0] ?? null;
    return publication ? { publication, provider: publication.source_extra.provider === "pubmed" ? "pubmed" : "europepmc", openAlexId: null, referencedWorks: [], relatedWorks: [] } : null;
  }
  if (parsed.kind === "pmcid") {
    const publication = await getEuropePmcByQuery(client, `PMCID:${parsed.value}`);
    return publication ? { publication, provider: "europepmc", openAlexId: null, referencedWorks: [], relatedWorks: [] } : null;
  }
  if (parsed.kind === "arxiv") {
    const publication = await getArxiv(client, parsed.value);
    return publication ? { publication, provider: "arxiv", openAlexId: null, referencedWorks: [], relatedWorks: [] } : null;
  }
  if (parsed.kind === "openalex") {
    const openAlex = await getOpenAlexWork(client, `works/${parsed.value}`);
    return openAlex ? { publication: openAlex.publication, provider: "openalex", openAlexId: openAlex.openAlexId, referencedWorks: openAlex.referencedWorks, relatedWorks: openAlex.relatedWorks } : null;
  }
  const publication = await getSemanticScholar(client, parsed.value);
  return publication ? { publication, provider: "semantic_scholar", openAlexId: null, referencedWorks: [], relatedWorks: [] } : null;
}

export async function getPublication(client: JsonClient, input: { id?: string; citation?: string; include_abstract?: boolean }): Promise<ToolResult> {
  const raw = input.id?.trim() || input.citation?.trim();
  if (!raw) return failure(literatureSource, "invalid_input", "Pass id or citation", false);
  const parsed = classifyReference(raw);
  const fullAbstract = input.include_abstract !== false;
  if (parsed.kind === "citation") {
    const pubmed = await searchPubmed(client, { query: raw, sort: "relevance", offset: 0, limit: 5 });
    const hit = pubmed.items.length > 0 ? pubmed : await searchOpenAlex(client, { query: raw, sort: "relevance", offset: 0, limit: 5 });
    const items = show(hit.items, fullAbstract);
    return makePage("publication", items, 0, Math.max(items.length, 1), items.length, hit.source, hit.usage, {
      ...hit.extra,
      match: items.length === 1 ? "top" : "ambiguous",
      note: "Citation or title lookup returns the closest records. Pass a DOI or PMID to get_publication when you already have one.",
    });
  }
  const loaded = await loadParsed(client, parsed);
  if (!loaded) return failure(providerMeta[parsed.kind === "doi" && parsed.value.startsWith("10.1101/") ? "biorxiv" : parsed.kind === "pmid" ? "europepmc" : "openalex"].source, "not_found", `No publication found for ${raw}`, false);
  const meta = providerMeta[loaded.provider];
  return makePage("publication", show([loaded.publication], fullAbstract), 0, 1, 1, meta.source, { summary: meta.summary, url: meta.url }, { provider: loaded.provider });
}

export async function getFullText(client: JsonClient, input: { id: string; section?: "abstract" | "methods" | "results" | "discussion" | "all"; max_chars?: number }): Promise<ToolResult> {
  const loaded = await loadParsed(client, classifyReference(input.id));
  if (!loaded) return failure(literatureSource, "not_found", `No publication found for ${input.id}`, false);
  const section = input.section ?? "all";
  const maxChars = input.max_chars ?? 4000;
  let publication = await withEuropePmcFullText(client, loaded.publication);
  const pmcid = publication.xrefs.find((id) => id.startsWith("pmcid:"))?.slice(6) ?? null;
  let xmlLicense: string | null = null;
  let parsedSections: { heading: string; text: string }[] = [];
  if (section === "abstract" && publication.abstract) {
    parsedSections = [{ heading: "Abstract", text: publication.abstract }];
  } else if (pmcid) {
    const xml = await europePmcXml(client, pmcid);
    if (xml) {
      const jats = parseJats(xml);
      xmlLicense = jats.license;
      parsedSections = jats.sections;
    }
  }
  const license = redistributableLicense(publication.license, xmlLicense) ?? publication.license ?? xmlLicense;
  const canInline = section === "abstract" || (parsedSections.length > 0 && (license?.startsWith("cc-by") || license === "cc0"));
  const clipped = canInline ? clipSections(parsedSections, maxChars, section === "abstract" ? "abstract" : section) : { sections: [], truncated: false };
  let links = publication.oa_urls;
  if (clipped.sections.length === 0 && publication.doi) {
    const open = await unpaywall(client, publication.doi);
    if (open) {
      links = [...links, ...open.urls.filter((url) => !links.some((existing) => existing.url === url.url))];
      publication = tryPublication({
        doi: publication.doi,
        pmid: publication.xrefs.find((id) => id.startsWith("pmid:"))?.slice(5),
        pmcid,
        arxiv: publication.xrefs.find((id) => id.startsWith("arxiv:"))?.slice(6),
        openalex: publication.xrefs.find((id) => id.startsWith("openalex:"))?.slice(9),
        s2: publication.xrefs.find((id) => id.startsWith("s2:"))?.slice(3),
        title: publication.title,
        authors: publication.authors,
        venue: publication.venue,
        published: publication.dates.published,
        posted: publication.dates.posted,
        updated: publication.dates.updated,
        abstract: publication.abstract,
        type: publication.type,
        is_oa: publication.is_oa ?? open.isOa,
        oa_urls: links,
        license: publication.license ?? open.license,
        published_as: publication.published_as,
        preprint_of: publication.preprint_of,
        source_extra: publication.source_extra,
      }) ?? publication;
    }
  }
  if (license && publication.license !== license) {
    publication = tryPublication({
      doi: publication.doi,
      pmid: publication.xrefs.find((id) => id.startsWith("pmid:"))?.slice(5),
      pmcid,
      arxiv: publication.xrefs.find((id) => id.startsWith("arxiv:"))?.slice(6),
      openalex: publication.xrefs.find((id) => id.startsWith("openalex:"))?.slice(9),
      s2: publication.xrefs.find((id) => id.startsWith("s2:"))?.slice(3),
      title: publication.title,
      authors: publication.authors,
      venue: publication.venue,
      published: publication.dates.published,
      posted: publication.dates.posted,
      updated: publication.dates.updated,
      abstract: publication.abstract,
      type: publication.type,
      is_oa: publication.is_oa ?? (license.startsWith("cc-by") || license === "cc0"),
      oa_urls: publication.oa_urls,
      license,
      published_as: publication.published_as,
      preprint_of: publication.preprint_of,
      source_extra: publication.source_extra,
    }) ?? publication;
  }
  const meta = providerMeta[loaded.provider];
  const mode = clipped.sections.length > 0 ? "sections" : "links";
  return makePage("publication", show([publication], false), 0, 1, 1, meta.source, { summary: meta.summary, url: meta.url }, {
    provider: loaded.provider,
    fulltext: mode,
    sections: clipped.sections,
    truncated: clipped.truncated,
    license: publication.license ?? license,
    links: (publication.oa_urls.length > 0 ? publication.oa_urls : links).map((link) => ({ url: link.url, content_type: link.content_type })),
    note: mode === "sections"
      ? "Open full text is clipped. Closed or unknown licenses are not copied."
      : contactEmail()
        ? "No redistributable full text was inlined. Use the open links."
        : "No redistributable full text was inlined. Set PHI_LITERATURE_EMAIL to enable Unpaywall link lookup.",
  });
}

async function withEuropePmcFullText(client: JsonClient, publication: Publication): Promise<Publication> {
  if (publication.xrefs.some((id) => id.startsWith("pmcid:"))) return publication;
  const pmid = publication.id.startsWith("pmid:") ? publication.id.slice(5) : publication.xrefs.find((id) => id.startsWith("pmid:"))?.slice(5);
  const query = publication.doi ? `DOI:${publication.doi}` : pmid ? `EXT_ID:${pmid} AND SRC:MED` : null;
  if (!query) return publication;
  const europe = await getEuropePmcByQuery(client, query);
  if (!europe?.xrefs.some((id) => id.startsWith("pmcid:"))) return publication;
  return combinePublications(publication, europe);
}

async function ensureOpenAlex(client: JsonClient, loaded: Loaded): Promise<Loaded> {
  if (loaded.openAlexId || !loaded.publication.doi) return loaded;
  const openAlex = await getOpenAlexWork(client, `works/https://doi.org/${loaded.publication.doi}`);
  if (!openAlex) return loaded;
  return { ...loaded, openAlexId: openAlex.openAlexId, referencedWorks: openAlex.referencedWorks, relatedWorks: openAlex.relatedWorks };
}

export async function relatedPublications(client: JsonClient, input: { id: string; relation?: "similar" | "citations" | "references"; offset?: number; limit?: number }): Promise<ToolResult> {
  const parsed = classifyReference(input.id);
  if (parsed.kind === "citation") return failure(literatureSource, "invalid_input", "related_publications needs a DOI, PMID, PMCID, arXiv, OpenAlex, or Semantic Scholar id", false);
  const loaded = await loadParsed(client, parsed);
  if (!loaded) return failure(literatureSource, "not_found", `No publication found for ${input.id}`, false);
  const offset = input.offset ?? 0;
  const limit = input.limit ?? 20;
  const relation = input.relation ?? "similar";
  const pmid = loaded.publication.id.startsWith("pmid:") ? loaded.publication.id.slice(5) : loaded.publication.xrefs.find((id) => id.startsWith("pmid:"))?.slice(5);
  if (relation === "similar" && pmid) {
    const ids = await pubmedNeighbors(client, pmid);
    const items = show(await fetchPubmed(client, ids.slice(offset, offset + limit)), false);
    const meta = providerMeta.pubmed;
    return makePage("publication", items, offset, limit, ids.length, meta.source, { summary: meta.summary, url: meta.url }, { provider: "pubmed", relation });
  }
  const graph = await ensureOpenAlex(client, loaded);
  if (relation === "citations" && graph.openAlexId) {
    const hit = await openAlexByFilter(client, `cites:${graph.openAlexId}`, offset, limit);
    return makePage("publication", show(hit.items, false), offset, limit, hit.total, hit.source, hit.usage, { ...hit.extra, relation });
  }
  if (relation === "references" && graph.referencedWorks.length > 0) {
    const items = show(await hydrateOpenAlex(client, graph.referencedWorks, offset, limit), false);
    const meta = providerMeta.openalex;
    return makePage("publication", items, offset, limit, graph.referencedWorks.length, meta.source, { summary: meta.summary, url: meta.url }, { provider: "openalex", relation });
  }
  if (relation === "similar" && graph.relatedWorks.length > 0) {
    const items = show(await hydrateOpenAlex(client, graph.relatedWorks, offset, limit), false);
    const meta = providerMeta.openalex;
    return makePage("publication", items, offset, limit, graph.relatedWorks.length, meta.source, { summary: meta.summary, url: meta.url }, { provider: "openalex", relation });
  }
  const doi = loaded.publication.doi;
  if (doi && (relation === "citations" || relation === "references")) {
    try {
      const dois = await openCitationDois(client, doi, relation);
      const slice = dois.slice(offset, offset + limit);
      const hit = slice.length === 0 ? null : await openAlexByFilter(client, `doi:${slice.join("|")}`, 0, slice.length);
      const items = show(hit?.items ?? [], false);
      const meta = hit ? providerMeta.openalex : providerMeta.openalex;
      return makePage("publication", items, offset, limit, dois.length, meta.source, { summary: "OpenCitations DOI list, hydrated from OpenAlex when the records exist.", url: "https://opencitations.net/" }, { provider: hit ? "openalex" : "opencitations", relation, unresolved_dois: slice.filter((item) => !items.some((publication) => publication.doi === item)) });
    } catch (error) {
      return toFailure(error, providerMeta.openalex.source);
    }
  }
  return failure(providerMeta.openalex.source, "not_found", `No ${relation} were found for ${input.id}`, false);
}

export async function resolveReferences(client: JsonClient, input: { references: string[]; include_abstract?: boolean }): Promise<ToolResult> {
  const fullAbstract = input.include_abstract === true;
  const parsed = input.references.map((reference, index) => ({ index, reference, parsed: classifyReference(reference) }));
  const dois = parsed.filter((item) => item.parsed.kind === "doi").map((item) => item.parsed.kind === "doi" ? item.parsed.value : "");
  const pmids = parsed.filter((item) => item.parsed.kind === "pmid").map((item) => item.parsed.kind === "pmid" ? item.parsed.value : "");
  const byDoi = new Map<string, Publication>();
  const byPmid = new Map<string, Publication>();
  if (dois.length > 0) {
    const hit = await openAlexByFilter(client, `doi:${dois.join("|")}`, 0, dois.length);
    for (const item of hit.items) if (item.doi) byDoi.set(item.doi, item);
    for (const doi of dois) {
      if (byDoi.has(doi)) continue;
      if (doi.startsWith("10.1101/")) {
        const preprint = await getBiorxiv(client, doi);
        if (preprint?.doi) byDoi.set(preprint.doi, preprint);
      }
    }
  }
  if (pmids.length > 0) {
    for (const item of await fetchPubmed(client, pmids)) {
      const pmid = item.id.startsWith("pmid:") ? item.id.slice(5) : item.xrefs.find((id) => id.startsWith("pmid:"))?.slice(5);
      if (pmid) byPmid.set(pmid, item);
    }
  }
  const items: Publication[] = [];
  const unresolved: Array<{ index: number; reference: string; reason: string }> = [];
  for (const item of parsed) {
    if (item.parsed.kind === "doi") {
      const found = byDoi.get(item.parsed.value);
      if (found) items.push(found);
      else unresolved.push({ index: item.index, reference: item.reference, reason: "not_found" });
      continue;
    }
    if (item.parsed.kind === "pmid") {
      const found = byPmid.get(item.parsed.value);
      if (found) items.push(found);
      else unresolved.push({ index: item.index, reference: item.reference, reason: "not_found" });
      continue;
    }
    if (item.parsed.kind === "pmcid" || item.parsed.kind === "arxiv" || item.parsed.kind === "openalex" || item.parsed.kind === "s2") {
      const loaded = await loadParsed(client, item.parsed);
      if (loaded) items.push(loaded.publication);
      else unresolved.push({ index: item.index, reference: item.reference, reason: "not_found" });
      continue;
    }
    const hit = await searchPubmed(client, { query: item.reference, sort: "relevance", offset: 0, limit: 1 });
    const top = hit.items[0];
    if (top) items.push(top);
    else unresolved.push({ index: item.index, reference: item.reference, reason: "not_found" });
  }
  return makePage("publication", show(items, fullAbstract), 0, Math.max(input.references.length, 1), items.length, literatureSource, {
    summary: "Batch resolution of a reference list. DOI and PMID lookups are batched. Unresolved lines are listed in source_extra.unresolved.",
    url: "https://github.com/dxsbiocc/phi-packages",
  }, { requested: input.references.length, unresolved });
}
