import type { Source, Usage } from "phi-omics-schemas";
import {
  asArray,
  asRecord as runtimeRecord,
  InputError,
  numberOrNull,
  text,
  type JsonClient,
} from "phi-omics-runtime";
import type { Publication } from "phi-omics-schemas";
import { parseXml } from "./jats.js";
import { tryPublication, type AuthorDraft, type OaDraft } from "./record.js";
import {
  asRecord,
  contactEmail,
  dateFromParts,
  isoDate,
  list,
  ncbiUrl,
  normalizeLicense,
  normalizePmcid,
  openAlexWindow,
  orcidOf,
  plain,
  providerMeta,
  readText,
  scholarRequest,
  stripTags,
  unique,
  utcDay,
  withMailto,
  type LiteratureSourceName,
} from "./util.js";

export interface HitPage {
  items: Publication[];
  total: number | null;
  source: Source;
  usage: Usage;
  extra: Record<string, unknown>;
}

export interface Query {
  query?: string;
  dateFrom?: string;
  dateTo?: string;
  type?: "journal-article" | "preprint" | "review" | "any";
  server?: "biorxiv" | "medrxiv";
  category?: string;
  sort: "date" | "relevance";
  offset: number;
  limit: number;
}

export interface OpenAlexHit {
  publication: Publication;
  openAlexId: string | null;
  referencedWorks: string[];
  relatedWorks: string[];
}

const SEARCH_SELECT = "id,doi,display_name,publication_date,publication_year,type,ids,authorships,primary_location,best_oa_location,locations,open_access,abstract_inverted_index,cited_by_count";
const DETAIL_SELECT = `${SEARCH_SELECT},referenced_works,related_works`;

function usageOf(name: LiteratureSourceName): { source: Source; usage: Usage } {
  const meta = providerMeta[name];
  return { source: meta.source, usage: { summary: meta.summary, url: meta.url } };
}

function page(name: LiteratureSourceName, items: Publication[], total: number | null, extra: Record<string, unknown> = {}): HitPage {
  const meta = usageOf(name);
  return { items, total, source: meta.source, usage: meta.usage, extra: { provider: name, ...extra } };
}

function mapType(value: string | null | undefined): Publication["type"] {
  if (!value) return null;
  const token = value.toLowerCase().replace(/_/g, "-");
  if (token.includes("preprint") || token === "posted-content") return "preprint";
  if (token.includes("review")) return "review";
  if (token.includes("book-chapter") || token === "chapter") return "book-chapter";
  if (token.includes("proceedings")) return "proceedings-article";
  if (token === "article" || token.includes("journal")) return "journal-article";
  return "other";
}

function invertedAbstract(value: unknown): string | null {
  const pairs: Array<[number, string]> = [];
  for (const [word, positions] of Object.entries(asRecord(value))) {
    if (!Array.isArray(positions)) continue;
    for (const position of positions) {
      if (typeof position === "number") pairs.push([position, word]);
    }
  }
  if (pairs.length === 0) return null;
  pairs.sort((left, right) => left[0] - right[0]);
  return pairs.map((pair) => pair[1]).join(" ");
}

function oaLink(url: string | null, contentType: OaDraft["content_type"], version: string | null = null): OaDraft | null {
  if (!url || !/^https?:\/\//i.test(url)) return null;
  return { url, content_type: contentType, version };
}

export async function fetchPubmed(client: JsonClient, ids: string[]): Promise<Publication[]> {
  if (ids.length === 0) return [];
  const xml = await readText(client, ncbiUrl("efetch.fcgi", { db: "pubmed", id: ids.join(","), retmode: "xml" }));
  const byPmid = new Map<string, Publication>();
  for (const article of pubmedArticles(xml)) {
    const item = fromPubmedArticle(article);
    const pmid = item?.id.startsWith("pmid:") ? item.id.slice(5) : item?.xrefs.find((id) => id.startsWith("pmid:"))?.slice(5);
    if (item && pmid) byPmid.set(pmid, item);
  }
  return ids.map((id) => byPmid.get(id)).filter((item): item is Publication => item != null);
}

function pubmedArticles(xml: string): unknown[] {
  const root = asRecord(parseXml(xml));
  const set = asRecord(root.PubmedArticleSet ?? root);
  return list(set.PubmedArticle);
}

function fromPubmedArticle(node: unknown): Publication | null {
  const article = asRecord(node);
  const medline = asRecord(article.MedlineCitation);
  const art = asRecord(medline.Article);
  const pmid = plain(medline.PMID).match(/\d+/)?.[0] ?? null;
  const journal = asRecord(art.Journal);
  const ids = list(asRecord(asRecord(article.PubmedData).ArticleIdList).ArticleId).map((id) => ({
    type: text(asRecord(id)["@_IdType"]),
    value: plain(id).trim(),
  }));
  const doi = ids.find((id) => id.type === "doi")?.value ?? null;
  const pmcid = ids.find((id) => id.type === "pmc")?.value ?? null;
  const pubTypes = list(asRecord(art.PublicationTypeList).PublicationType).map((item) => plain(item)).filter(Boolean);
  const articleDates = list(art.ArticleDate);
  const electronic = articleDates
    .map(asRecord)
    .find((row) => text(row["@_DateType"]) === "Electronic");
  const type = pubTypes.some((item) => /review/i.test(item))
    ? "review"
    : pubTypes.some((item) => /preprint/i.test(item))
      ? "preprint"
      : "journal-article";
  return tryPublication({
    doi,
    pmid,
    pmcid,
    title: plain(art.ArticleTitle),
    authors: pubmedAuthors(art.AuthorList),
    venue: plain(journal.Title) || null,
    published: dateFromParts(asRecord(asRecord(journal.JournalIssue).PubDate)),
    posted: electronic ? dateFromParts(electronic) : null,
    abstract: pubmedAbstract(art.Abstract),
    type,
    oa_urls: doi ? [oaLink(`https://doi.org/${doi.toLowerCase()}`, "landing")] : [],
    source_extra: { provider: "pubmed", pub_types: pubTypes.slice(0, 8) },
  });
}

function pubmedAuthors(value: unknown): AuthorDraft[] {
  return list(asRecord(value).Author).flatMap((node) => {
    const author = asRecord(node);
    const name = plain(author.CollectiveName) || [plain(author.ForeName), plain(author.LastName)].filter(Boolean).join(" ");
    if (!name.trim()) return [];
    const orcid = list(author.Identifier)
      .map((id) => text(asRecord(id)["@_Source"])?.toUpperCase() === "ORCID" ? orcidOf(plain(id)) : null)
      .find((id) => id != null) ?? null;
    return [{ name: name.trim(), orcid }];
  });
}

function pubmedAbstract(value: unknown): string | null {
  const parts = list(asRecord(value).AbstractText);
  if (parts.length === 0) return plain(value).trim() || null;
  const joined = parts.map((part) => {
    const label = text(asRecord(part)["@_Label"]);
    const body = plain(part).replace(/\s+/g, " ").trim();
    return label && body ? `${label}: ${body}` : body;
  }).filter(Boolean).join("\n");
  return joined || null;
}

export async function searchPubmed(client: JsonClient, query: Query): Promise<HitPage> {
  if (query.type === "preprint") {
    throw new InputError("PubMed does not index most preprints. Search with type=preprint and no sources to use Europe PMC.");
  }
  const term = pubmedTerm(query);
  const body = runtimeRecord(await client.getJson(ncbiUrl("esearch.fcgi", {
    db: "pubmed",
    term,
    retstart: String(query.offset),
    retmax: String(query.limit),
    sort: query.sort === "date" ? "pub_date" : "relevance",
    retmode: "json",
  })));
  const result = runtimeRecord(body.esearchresult);
  const ids = asArray(result.idlist).map(text).filter((id): id is string => id != null);
  const items = await fetchPubmed(client, ids);
  return page("pubmed", items, numberOrNull(result.count), { query_translation: text(result.querytranslation) });
}

function pubmedTerm(query: Query): string {
  let term = query.query?.trim() ?? "";
  if (!term) throw new InputError("query is required for PubMed search");
  if (query.type === "review" && !/\[pt\]/i.test(term)) term = `(${term}) AND review[pt]`;
  if (query.type === "journal-article" && !/\[pt\]/i.test(term)) term = `(${term}) AND journal article[pt]`;
  if ((query.dateFrom || query.dateTo) && !/\[(?:pdat|dp|edat)\]/i.test(term)) {
    const start = (query.dateFrom ?? "1900").replace(/-/g, "/");
    const end = (query.dateTo ?? "3000").replace(/-/g, "/");
    term = `(${term}) AND ${start}:${end}[pdat]`;
  }
  return term;
}

export async function pubmedNeighbors(client: JsonClient, pmid: string): Promise<string[]> {
  const body = runtimeRecord(await client.getJson(ncbiUrl("elink.fcgi", {
    dbfrom: "pubmed",
    db: "pubmed",
    id: pmid,
    cmd: "neighbor_score",
    retmode: "json",
  })));
  const linkset = runtimeRecord(asArray(body.linksets)[0]);
  const databases = asArray(linkset.linksetdbs).map(runtimeRecord);
  const links = asArray(databases.find((row) => text(row.linkname)?.includes("pubmed"))?.links ?? databases[0]?.links);
  return links.map((link) => text(runtimeRecord(link).id) ?? text(link)).filter((id): id is string => id != null && id !== pmid);
}

export async function searchEuropePmc(client: JsonClient, query: Query): Promise<HitPage> {
  let textQuery = query.query?.trim() ?? "";
  if (!textQuery) throw new InputError("query is required for Europe PMC search");
  if (query.type === "preprint" && !/SRC:PPR/i.test(textQuery)) textQuery = `(${textQuery}) AND SRC:PPR`;
  if (query.server === "biorxiv") textQuery = `(${textQuery}) AND JOURNAL:"bioRxiv"`;
  if (query.server === "medrxiv") textQuery = `(${textQuery}) AND JOURNAL:"medRxiv"`;
  if (query.category && query.type === "preprint") textQuery = `(${textQuery}) AND (${query.category})`;
  if (query.type === "review" && !/PUB_TYPE:/i.test(textQuery)) textQuery = `(${textQuery}) AND PUB_TYPE:review`;
  if (query.dateFrom || query.dateTo) {
    textQuery = `(${textQuery}) AND FIRST_PDATE:[${query.dateFrom ?? "1900-01-01"} TO ${query.dateTo ?? "3000-01-01"}]`;
  }
  const fetchSize = query.offset + query.limit;
  if (fetchSize > 1000) throw new InputError("offset + limit must be at most 1000 for Europe PMC");
  const url = new URL("https://www.ebi.ac.uk/europepmc/webservices/rest/search");
  url.searchParams.set("query", textQuery);
  url.searchParams.set("format", "json");
  url.searchParams.set("resultType", "core");
  url.searchParams.set("pageSize", String(fetchSize));
  url.searchParams.set("cursorMark", "*");
  if (query.sort === "date") url.searchParams.set("sort", "P_PDATE_D desc");
  const body = runtimeRecord(await client.getJson(url.toString()));
  const results = asArray(runtimeRecord(body.resultList).result).slice(query.offset, query.offset + query.limit);
  return page("europepmc", results.map(fromEuropePmc).filter((item): item is Publication => item != null), numberOrNull(body.hitCount));
}

export async function getEuropePmcByQuery(client: JsonClient, textQuery: string): Promise<Publication | null> {
  const url = new URL("https://www.ebi.ac.uk/europepmc/webservices/rest/search");
  url.searchParams.set("query", textQuery);
  url.searchParams.set("format", "json");
  url.searchParams.set("resultType", "core");
  url.searchParams.set("pageSize", "1");
  const body = runtimeRecord(await client.getJson(url.toString()));
  return fromEuropePmc(asArray(runtimeRecord(body.resultList).result)[0]);
}

function fromEuropePmc(node: unknown): Publication | null {
  const row = runtimeRecord(node);
  if (!text(row.id) && !text(row.doi)) return null;
  const authors = asArray(runtimeRecord(row.authorList).author).flatMap((authorNode) => {
    const author = runtimeRecord(authorNode);
    const name = text(author.fullName) ?? [text(author.firstName), text(author.lastName)].filter(Boolean).join(" ");
    if (!name) return [];
    const ids = list(author.authorId);
    const orcid = ids.map((id) => {
      const record = runtimeRecord(id);
      return text(record.type)?.toUpperCase() === "ORCID" ? orcidOf(text(record.value)) : null;
    }).find((id) => id != null) ?? null;
    return [{ name, orcid }];
  });
  const preprint = text(row.source) === "PPR" || asArray(runtimeRecord(row.pubTypeList).pubType).some((item) => /preprint/i.test(String(item)));
  const posted = text(row.firstPublicationDate);
  const urls = asArray(runtimeRecord(row.fullTextUrlList).fullTextUrl).flatMap((linkNode) => {
    const link = runtimeRecord(linkNode);
    const href = text(link.url);
    const free = text(link.availabilityCode) === "F" || text(link.availabilityCode) === "OA";
    if (!href || !free) return [];
    const style = text(link.documentStyle)?.toLowerCase();
    const contentType: OaDraft["content_type"] = style === "pdf" ? "pdf" : style === "html" ? "html" : style === "doi" ? "landing" : "other";
    return [oaLink(href, contentType)];
  });
  return tryPublication({
    doi: text(row.doi),
    pmid: text(row.pmid) ?? (text(row.source) === "MED" ? text(row.id) : null),
    pmcid: text(row.pmcid),
    title: text(row.title),
    authors,
    venue: text(runtimeRecord(runtimeRecord(row.journalInfo).journal).title),
    published: preprint ? null : posted,
    posted: preprint ? posted : text(row.electronicPublicationDate),
    updated: isoDate(text(row.dateOfRevision)),
    abstract: text(row.abstractText),
    type: preprint ? "preprint" : mapType(text(asArray(runtimeRecord(row.pubTypeList).pubType).find((item) => /review/i.test(String(item))) ?? "journal-article")),
    is_oa: text(row.isOpenAccess) === "Y" ? true : text(row.isOpenAccess) === "N" ? false : null,
    license: text(row.license),
    oa_urls: urls,
    source_extra: {
      provider: "europepmc",
      citation_count: numberOrNull(row.citedByCount),
      epmc_source: text(row.source),
    },
  });
}

function openAlexUrl(path: string, params: Record<string, string> = {}): string {
  const url = withMailto(new URL(`https://api.openalex.org/${path}`));
  for (const [key, value] of Object.entries(params)) url.searchParams.set(key, value);
  return url.toString();
}

export function fromOpenAlex(node: unknown): OpenAlexHit | null {
  const work = runtimeRecord(node);
  if (text(work.error) || text(work.message)?.includes("not found")) return null;
  const ids = runtimeRecord(work.ids);
  const doi = text(work.doi) ?? text(ids.doi);
  const openAlexId = text(ids.openalex)?.match(/W\d+/i)?.[0]?.toUpperCase() ?? text(work.id)?.match(/W\d+/i)?.[0]?.toUpperCase() ?? null;
  const locationDois = unique([work.primary_location, work.best_oa_location, ...asArray(work.locations)].map((location) => {
    const record = runtimeRecord(location);
    return text(record.landing_page_url);
  }));
  const normalizedDois = unique(locationDois.map((value) => {
    const match = value?.match(/10\.\d{4,9}\/\S+/i);
    return match?.[0] ?? null;
  }));
  const ownDoi = (doi?.match(/10\.\d{4,9}\/\S+/i)?.[0] ?? null);
  const other = normalizedDois.find((value) => value && ownDoi && value.toLowerCase() !== ownDoi.toLowerCase()) ?? null;
  const preprint = ownDoi?.toLowerCase().startsWith("10.1101/") || text(work.type) === "preprint";
  const best = runtimeRecord(work.best_oa_location);
  const primary = runtimeRecord(work.primary_location);
  const publication = tryPublication({
    doi: ownDoi,
    pmid: text(ids.pmid)?.match(/(\d{1,9})$/)?.[1] ?? null,
    pmcid: text(ids.pmcid)?.match(/(PMC\d+)/i)?.[1] ?? null,
    openalex: openAlexId,
    title: text(work.display_name),
    authors: asArray(work.authorships).flatMap((row) => {
      const author = runtimeRecord(runtimeRecord(row).author);
      const name = text(author.display_name);
      return name ? [{ name, orcid: orcidOf(text(author.orcid)) }] : [];
    }),
    venue: text(runtimeRecord(primary.source).display_name),
    published: preprint ? null : text(work.publication_date) ?? (work.publication_year ? String(work.publication_year) : null),
    posted: preprint ? text(work.publication_date) : null,
    abstract: invertedAbstract(work.abstract_inverted_index),
    type: mapType(text(work.type)),
    is_oa: runtimeRecord(work.open_access).is_oa === true ? true : runtimeRecord(work.open_access).is_oa === false ? false : null,
    license: text(best.license),
    oa_urls: [
      oaLink(text(best.url_for_pdf), "pdf", text(best.version)),
      oaLink(text(best.url) ?? text(runtimeRecord(work.open_access).oa_url), "html", text(best.version)),
      oaLink(text(primary.landing_page_url), "landing", null),
    ],
    published_as: preprint && other && !other.toLowerCase().startsWith("10.1101/") ? other : null,
    preprint_of: !preprint && other?.toLowerCase().startsWith("10.1101/") ? other : null,
    source_extra: { provider: "openalex", citation_count: numberOrNull(work.cited_by_count), openalex_id: openAlexId },
  });
  if (!publication) return null;
  return {
    publication,
    openAlexId,
    referencedWorks: asArray(work.referenced_works).map(text).filter((id): id is string => id != null),
    relatedWorks: asArray(work.related_works).map(text).filter((id): id is string => id != null),
  };
}

export async function getOpenAlexWork(client: JsonClient, path: string): Promise<OpenAlexHit | null> {
  try {
    return fromOpenAlex(await client.getJson(openAlexUrl(path, { select: DETAIL_SELECT })));
  } catch (error) {
    if (error instanceof Error && "status" in error && error.status === 404) return null;
    throw error;
  }
}

export async function searchOpenAlex(client: JsonClient, query: Query): Promise<HitPage> {
  if (!query.query?.trim()) throw new InputError("query is required for OpenAlex search");
  const filters: string[] = [];
  if (query.dateFrom) filters.push(`from_publication_date:${query.dateFrom.length === 4 ? `${query.dateFrom}-01-01` : query.dateFrom}`);
  if (query.dateTo) filters.push(`to_publication_date:${query.dateTo.length === 4 ? `${query.dateTo}-12-31` : query.dateTo}`);
  if (query.type === "preprint") filters.push("type:preprint");
  if (query.type === "review") filters.push("type:review");
  if (query.type === "journal-article") filters.push("type:article");
  const window = openAlexWindow(query.offset, query.limit);
  const params: Record<string, string> = {
    search: query.query,
    "per-page": String(window.perPage),
    page: String(window.page),
    select: SEARCH_SELECT,
    sort: query.sort === "date" ? "publication_date:desc" : "relevance_score:desc",
  };
  if (filters.length > 0) params.filter = filters.join(",");
  const body = runtimeRecord(await client.getJson(openAlexUrl("works", params)));
  const items = asArray(body.results).slice(window.skip, window.skip + query.limit).map((row) => fromOpenAlex(row)?.publication).filter((item): item is Publication => item != null);
  return page("openalex", items, numberOrNull(runtimeRecord(body.meta).count));
}

export async function hydrateOpenAlex(client: JsonClient, openAlexIds: string[], offset: number, limit: number): Promise<Publication[]> {
  const slice = openAlexIds.slice(offset, offset + limit);
  const shorts = slice.map((id) => id.match(/W\d+/i)?.[0]).filter((id): id is string => id != null);
  if (shorts.length === 0) return [];
  const body = runtimeRecord(await client.getJson(openAlexUrl("works", {
    filter: `openalex:${shorts.join("|")}`,
    "per-page": String(shorts.length),
    select: SEARCH_SELECT,
  })));
  const byId = new Map<string, Publication>();
  for (const row of asArray(body.results)) {
    const hit = fromOpenAlex(row);
    if (hit?.openAlexId && hit.publication) byId.set(hit.openAlexId, hit.publication);
  }
  return shorts.map((id) => byId.get(id.toUpperCase())).filter((item): item is Publication => item != null);
}

export async function openAlexByFilter(client: JsonClient, filter: string, offset: number, limit: number): Promise<HitPage> {
  const window = openAlexWindow(offset, limit);
  const body = runtimeRecord(await client.getJson(openAlexUrl("works", {
    filter,
    "per-page": String(window.perPage),
    page: String(window.page),
    select: SEARCH_SELECT,
    sort: "publication_date:desc",
  })));
  const items = asArray(body.results).slice(window.skip, window.skip + limit).map((row) => fromOpenAlex(row)?.publication).filter((item): item is Publication => item != null);
  return page("openalex", items, numberOrNull(runtimeRecord(body.meta).count));
}

export async function getCrossrefWork(client: JsonClient, doi: string): Promise<Publication | null> {
  try {
    const body = runtimeRecord(await client.getJson(crossrefUrl(`works/${encodeURIComponent(doi)}`)));
    return fromCrossref(body.message);
  } catch (error) {
    if (error instanceof Error && "status" in error && error.status === 404) return null;
    throw error;
  }
}

export async function searchCrossref(client: JsonClient, query: Query): Promise<HitPage> {
  if (!query.query?.trim()) throw new InputError("query is required for Crossref search");
  const filters: string[] = [];
  if (query.dateFrom) filters.push(`from-pub-date:${query.dateFrom}`);
  if (query.dateTo) filters.push(`until-pub-date:${query.dateTo}`);
  if (query.type === "preprint") filters.push("type:posted-content");
  if (query.type === "journal-article" || query.type === "review") filters.push("type:journal-article");
  const params: Record<string, string> = {
    "query.bibliographic": query.query,
    rows: String(query.limit),
    offset: String(query.offset),
    sort: query.sort === "date" ? "published" : "relevance",
    order: "desc",
  };
  if (filters.length > 0) params.filter = filters.join(",");
  const body = runtimeRecord(await client.getJson(crossrefUrl("works", params)));
  const message = runtimeRecord(body.message);
  const items = asArray(message.items).map(fromCrossref).filter((item): item is Publication => item != null);
  return page("crossref", items, numberOrNull(message["total-results"]));
}

function crossrefUrl(path: string, params: Record<string, string> = {}): string {
  const url = withMailto(new URL(`https://api.crossref.org/${path}`));
  for (const [key, value] of Object.entries(params)) url.searchParams.set(key, value);
  return url.toString();
}

function fromCrossref(node: unknown): Publication | null {
  const work = runtimeRecord(node);
  const doi = text(work.DOI);
  if (!doi) return null;
  const relation = runtimeRecord(work.relation);
  const publishedAs = relationDoi(relation["is-preprint-of"]);
  const preprintOf = relationDoi(relation["has-preprint"]);
  const issued = asArray(runtimeRecord(work.issued)["date-parts"])[0];
  const posted = asArray(runtimeRecord(work.posted)["date-parts"])[0];
  const preprint = text(work.type) === "posted-content";
  const licenseUrl = text(runtimeRecord(asArray(work.license)[0]).URL);
  return tryPublication({
    doi,
    title: stripTags(text(asArray(work.title)[0])),
    authors: asArray(work.author).flatMap((row) => {
      const author = runtimeRecord(row);
      const name = text(author.name) ?? [text(author.given), text(author.family)].filter(Boolean).join(" ");
      return name ? [{ name, orcid: orcidOf(text(author.ORCID)) }] : [];
    }),
    venue: text(asArray(work["container-title"])[0]),
    published: preprint ? null : partsToIso(issued),
    posted: preprint ? partsToIso(posted) ?? partsToIso(issued) : null,
    abstract: stripTags(text(work.abstract)),
    type: mapType(text(work.type)),
    is_oa: licenseUrl && /creativecommons|cc0|publicdomain/i.test(licenseUrl) ? true : null,
    license: licenseUrl,
    oa_urls: [
      oaLink(doi ? `https://doi.org/${doi.toLowerCase()}` : null, "landing"),
      ...asArray(work.link).map((link) => {
        const record = runtimeRecord(link);
        const content = text(record["content-type"]) ?? "";
        const contentType: OaDraft["content_type"] = content.includes("pdf") ? "pdf" : content.includes("html") ? "html" : "other";
        return oaLink(text(record.URL), contentType);
      }),
    ],
    published_as: preprint ? publishedAs : null,
    preprint_of: preprint ? null : preprintOf,
    source_extra: { provider: "crossref" },
  });
}

function relationDoi(value: unknown): string | null {
  return text(runtimeRecord(asArray(value)[0]).id);
}

function partsToIso(value: unknown): string | null {
  if (!Array.isArray(value) || value.length === 0) return null;
  const [year, month, day] = value;
  if (typeof year !== "number") return null;
  if (typeof month !== "number") return String(year);
  const mm = String(month).padStart(2, "0");
  if (typeof day !== "number") return `${year}-${mm}`;
  return `${year}-${mm}-${String(day).padStart(2, "0")}`;
}

export async function getBiorxiv(client: JsonClient, doi: string): Promise<Publication | null> {
  for (const server of ["biorxiv", "medrxiv"] as const) {
    const body = runtimeRecord(await client.getJson(`https://api.biorxiv.org/details/${server}/${doi}/na/json`));
    const rows = asArray(body.collection).map(runtimeRecord);
    const latest = rows[rows.length - 1];
    if (!latest) continue;
    return fromBiorxiv(latest, server);
  }
  return null;
}

export async function searchBiorxiv(client: JsonClient, server: "biorxiv" | "medrxiv", query: Query): Promise<HitPage> {
  if (query.query?.trim()) {
    throw new InputError(`${server} is browsed by date. For preprints on a topic, call search_publications with type=preprint and no sources.`);
  }
  const dateTo = query.dateTo ?? utcDay(0);
  const dateFrom = query.dateFrom ?? utcDay(-7);
  const first = runtimeRecord(await client.getJson(`https://api.biorxiv.org/details/${server}/${dateFrom}/${dateTo}/0/json`));
  const message = runtimeRecord(asArray(first.messages)[0]);
  const total = numberOrNull(message.total) ?? numberOrNull(message.count);
  const cursor = Math.max(0, (total ?? 0) - query.offset - query.limit);
  const firstRows = asArray(first.collection).map(runtimeRecord);
  const body = cursor > 0 && cursor + query.limit > firstRows.length
    ? runtimeRecord(await client.getJson(`https://api.biorxiv.org/details/${server}/${dateFrom}/${dateTo}/${cursor}/json`))
    : first;
  let rows = asArray(body.collection).map(runtimeRecord);
  if (cursor > 0 && body === first) rows = rows.slice(cursor, cursor + query.limit);
  else rows = rows.slice(0, query.limit);
  if (query.category) {
    const category = query.category.toLowerCase();
    rows = rows.filter((row) => text(row.category)?.toLowerCase() === category);
  }
  rows.reverse();
  const items = rows.map((row) => fromBiorxiv(row, server)).filter((item): item is Publication => item != null);
  return page(server, items, total, { interval: `${dateFrom}/${dateTo}`, category: query.category ?? null });
}

function fromBiorxiv(row: Record<string, unknown>, server: "biorxiv" | "medrxiv"): Publication | null {
  const doi = text(row.doi);
  const version = text(row.version);
  const published = text(row.published);
  return tryPublication({
    doi,
    title: text(row.title),
    authors: (text(row.authors) ?? "").split(/;\s*/).filter(Boolean).map((name) => ({ name, orcid: null })),
    venue: server,
    posted: text(row.date),
    abstract: text(row.abstract),
    type: "preprint",
    is_oa: true,
    license: text(row.license),
    oa_urls: [
      oaLink(doi ? `https://www.${server}.org/content/${doi}` : null, "landing", version),
      oaLink(doi ? `https://www.${server}.org/content/${doi}${version ? `v${version}` : ""}.full.pdf` : null, "pdf", version),
      oaLink(text(row.jatsxml), "xml", version),
    ],
    published_as: published,
    source_extra: { provider: server, category: text(row.category), version },
  });
}

const S2_FIELDS = "paperId,title,abstract,authors,externalIds,url,openAccessPdf,publicationDate,venue,journal,publicationTypes";

export async function searchSemanticScholar(client: JsonClient, query: Query): Promise<HitPage> {
  if (!query.query?.trim()) throw new InputError("query is required for Semantic Scholar search");
  const url = new URL("https://api.semanticscholar.org/graph/v1/paper/search");
  url.searchParams.set("query", query.query);
  url.searchParams.set("offset", String(query.offset));
  url.searchParams.set("limit", String(query.limit));
  url.searchParams.set("fields", S2_FIELDS);
  if (query.dateFrom || query.dateTo) url.searchParams.set("year", `${query.dateFrom?.slice(0, 4) ?? ""}-${query.dateTo?.slice(0, 4) ?? ""}`);
  const body = runtimeRecord(await client.getJson(url.toString(), scholarRequest()));
  const items = asArray(body.data).map(fromSemanticScholar).filter((item): item is Publication => item != null);
  return page("semantic_scholar", items, numberOrNull(body.total));
}

export async function getSemanticScholar(client: JsonClient, paperId: string): Promise<Publication | null> {
  const url = `https://api.semanticscholar.org/graph/v1/paper/${encodeURIComponent(paperId)}?fields=${encodeURIComponent(S2_FIELDS)}`;
  try {
    return fromSemanticScholar(await client.getJson(url, scholarRequest()));
  } catch (error) {
    if (error instanceof Error && "status" in error && error.status === 404) return null;
    throw error;
  }
}

function fromSemanticScholar(node: unknown): Publication | null {
  const paper = runtimeRecord(node);
  const external = runtimeRecord(paper.externalIds);
  const types = asArray(paper.publicationTypes).map(text).filter((item): item is string => item != null);
  const journal = runtimeRecord(paper.journal);
  return tryPublication({
    doi: text(external.DOI),
    pmid: text(external.PubMed),
    pmcid: text(external.PubMedCentral),
    arxiv: text(external.ArXiv),
    s2: text(paper.paperId),
    title: text(paper.title),
    authors: asArray(paper.authors).flatMap((row) => {
      const name = text(runtimeRecord(row).name);
      return name ? [{ name, orcid: null }] : [];
    }),
    venue: text(journal.name) ?? text(paper.venue),
    published: text(paper.publicationDate),
    abstract: text(paper.abstract),
    type: types.some((item) => /review/i.test(item)) ? "review" : types.some((item) => /preprint/i.test(item)) ? "preprint" : mapType(types[0] ?? "article"),
    is_oa: text(runtimeRecord(paper.openAccessPdf).url) ? true : null,
    oa_urls: [
      oaLink(text(runtimeRecord(paper.openAccessPdf).url), "pdf"),
      oaLink(text(paper.url), "landing"),
    ],
    source_extra: { provider: "semantic_scholar" },
  });
}

export async function searchArxiv(client: JsonClient, query: Query): Promise<HitPage> {
  const search = query.query?.trim()
    ? (/\bcat:/i.test(query.query) ? query.query : `cat:q-bio* AND all:${query.query}`)
    : "cat:q-bio*";
  const url = new URL("https://export.arxiv.org/api/query");
  url.searchParams.set("search_query", search);
  url.searchParams.set("start", String(query.offset));
  url.searchParams.set("max_results", String(query.limit));
  url.searchParams.set("sortBy", query.sort === "date" ? "submittedDate" : "relevance");
  url.searchParams.set("sortOrder", "descending");
  return arxivFeed(client, url);
}

export async function getArxiv(client: JsonClient, id: string): Promise<Publication | null> {
  const url = new URL("https://export.arxiv.org/api/query");
  url.searchParams.set("id_list", id);
  const hit = await arxivFeed(client, url);
  return hit.items[0] ?? null;
}

async function arxivFeed(client: JsonClient, url: URL): Promise<HitPage> {
  const root = asRecord(parseXml(await readText(client, url.toString())));
  const feed = asRecord(root.feed);
  const total = numberOrNull(plain(feed.totalResults));
  const items = list(feed.entry).map(fromArxiv).filter((item): item is Publication => item != null);
  return page("arxiv", items, total);
}

function fromArxiv(node: unknown): Publication | null {
  const entry = asRecord(node);
  const id = plain(entry.id);
  const arxiv = id.match(/arxiv\.org\/abs\/([^v\s]+)(?:v\d+)?/i)?.[1] ?? null;
  const rawDoi = plain(entry.doi);
  const ownDoi = rawDoi.toLowerCase().startsWith("10.48550/") ? rawDoi : null;
  const journalDoi = rawDoi && !ownDoi ? rawDoi : null;
  const links = list(entry.link);
  const pdf = links.map(asRecord).find((link) => text(link["@_title"]) === "pdf" || text(link["@_type"])?.includes("pdf"));
  return tryPublication({
    doi: ownDoi,
    arxiv,
    title: plain(entry.title),
    authors: list(entry.author).flatMap((author) => {
      const name = plain(asRecord(author).name);
      return name ? [{ name, orcid: null }] : [];
    }),
    venue: "arXiv",
    posted: isoDate(plain(entry.published)),
    updated: isoDate(plain(entry.updated)),
    abstract: plain(entry.summary),
    type: "preprint",
    is_oa: true,
    license: plain(entry.license) || null,
    oa_urls: [
      oaLink(id.startsWith("http") ? id : null, "landing"),
      oaLink(text(pdf?.["@_href"]), "pdf"),
    ],
    published_as: journalDoi,
    source_extra: { provider: "arxiv", category: text(asRecord(entry.primary_category)["@_term"]) },
  });
}

export async function unpaywall(client: JsonClient, doi: string): Promise<{ license: string | null; isOa: boolean | null; urls: OaDraft[] } | null> {
  const email = contactEmail();
  if (!email) return null;
  try {
    const body = runtimeRecord(await client.getJson(`https://api.unpaywall.org/v2/${encodeURIComponent(doi)}?email=${encodeURIComponent(email)}`));
    const locations = [body.best_oa_location, ...asArray(body.oa_locations)].map(runtimeRecord);
    const urls = locations.flatMap((location) => [
      oaLink(text(location.url_for_pdf), "pdf", text(location.version)),
      oaLink(text(location.url), "html", text(location.version)),
    ]).filter((item): item is OaDraft => item != null);
    return {
      license: normalizeLicense(text(runtimeRecord(body.best_oa_location).license)),
      isOa: body.is_oa === true ? true : body.is_oa === false ? false : null,
      urls,
    };
  } catch (error) {
    if (error instanceof Error && "status" in error && (error.status === 404 || error.status === 422)) return null;
    throw error;
  }
}

export async function openCitationDois(client: JsonClient, doi: string, kind: "citations" | "references"): Promise<string[]> {
  const path = kind === "citations" ? "citations" : "references";
  const body = await client.getJson(`https://opencitations.net/index/api/v1/${path}/doi:${encodeURIComponent(doi)}`);
  const field = kind === "citations" ? "citing" : "cited";
  return unique(asArray(body).map((row) => text(runtimeRecord(row)[field])));
}

export async function europePmcXml(client: JsonClient, pmcid: string): Promise<string | null> {
  const id = normalizePmcid(pmcid) ?? (pmcid.toUpperCase().startsWith("PMC") ? pmcid.toUpperCase() : `PMC${pmcid}`);
  try {
    const xml = await readText(client, `https://www.ebi.ac.uk/europepmc/webservices/rest/${id}/fullTextXML`);
    return xml.includes("<article") ? xml : null;
  } catch (error) {
    if (error instanceof Error && "status" in error && error.status === 404) return null;
    throw error;
  }
}
