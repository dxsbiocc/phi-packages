import {
  normalizeArxiv,
  normalizeDoi,
  type Publication,
  type Source,
} from "phi-omics-schemas";
import {
  InputError,
  failure,
  sourceOf,
  type JsonClient,
  type UpstreamRequest,
} from "phi-omics-runtime";

export const LITERATURE_SOURCES = [
  "pubmed",
  "europepmc",
  "biorxiv",
  "medrxiv",
  "crossref",
  "openalex",
  "semantic_scholar",
  "arxiv",
] as const;

export type LiteratureSourceName = (typeof LITERATURE_SOURCES)[number];

export const literatureSource = sourceOf("literature", "Phi literature", "https://github.com/dxsbiocc/phi-packages");

export const providerMeta: Record<LiteratureSourceName, { source: Source; summary: string; url: string }> = {
  pubmed: {
    source: sourceOf("pubmed", "PubMed", "https://pubmed.ncbi.nlm.nih.gov/"),
    summary: "PubMed citations via NCBI E-utilities. Use MeSH and field tags. Abstracts in search hits are truncated.",
    url: "https://www.ncbi.nlm.nih.gov/home/about/policies/",
  },
  europepmc: {
    source: sourceOf("europepmc", "Europe PMC", "https://europepmc.org/"),
    summary: "Europe PMC abstracts, preprint records, and open full text XML.",
    url: "https://europepmc.org/About",
  },
  biorxiv: {
    source: sourceOf("biorxiv", "bioRxiv", "https://www.biorxiv.org/"),
    summary: "bioRxiv preprints by date interval. A journal DOI is returned in published_as when bioRxiv reports one.",
    url: "https://www.biorxiv.org/about/FAQ",
  },
  medrxiv: {
    source: sourceOf("medrxiv", "medRxiv", "https://www.medrxiv.org/"),
    summary: "medRxiv preprints by date interval. A journal DOI is returned in published_as when medRxiv reports one.",
    url: "https://www.medrxiv.org/about/FAQ",
  },
  crossref: {
    source: sourceOf("crossref", "Crossref", "https://www.crossref.org/"),
    summary: "Crossref works, used when a DOI is missing from OpenAlex.",
    url: "https://www.crossref.org/documentation/",
  },
  openalex: {
    source: sourceOf("openalex", "OpenAlex", "https://openalex.org/"),
    summary: "OpenAlex works, including DOI-only online-first articles, citations, and references.",
    url: "https://openalex.org/about",
  },
  semantic_scholar: {
    source: sourceOf("semantic-scholar", "Semantic Scholar", "https://www.semanticscholar.org/"),
    summary: "Semantic Scholar graph API. Opt-in. Set S2_API_KEY for the higher rate limit.",
    url: "https://www.semanticscholar.org/product/api/terms",
  },
  arxiv: {
    source: sourceOf("arxiv", "arXiv", "https://arxiv.org/"),
    summary: "arXiv q-bio metadata. Opt-in. Journal DOIs are linked with published_as and are not merged into the preprint.",
    url: "https://info.arxiv.org/help/api/index.html",
  },
};

const MONTHS: Record<string, string> = {
  jan: "01", feb: "02", mar: "03", apr: "04", may: "05", jun: "06",
  jul: "07", aug: "08", sep: "09", oct: "10", nov: "11", dec: "12",
};

export function contactEmail(): string | null {
  const email = process.env.PHI_LITERATURE_EMAIL || process.env.NCBI_EMAIL || process.env.UNPAYWALL_EMAIL || "";
  return email.trim() || null;
}

export function ncbiUrl(script: string, params: Record<string, string>): string {
  const url = new URL(`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/${script}`);
  url.searchParams.set("tool", "phi-literature");
  const email = contactEmail();
  if (email) url.searchParams.set("email", email);
  const key = process.env.NCBI_API_KEY;
  if (key) url.searchParams.set("api_key", key);
  for (const [keyName, value] of Object.entries(params)) url.searchParams.set(keyName, value);
  return url.toString();
}

export function withMailto(url: URL): URL {
  const email = contactEmail();
  if (email && !url.searchParams.has("mailto")) url.searchParams.set("mailto", email);
  return url;
}

export function scholarRequest(): UpstreamRequest {
  const key = process.env.S2_API_KEY || process.env.SEMANTIC_SCHOLAR_API_KEY;
  return key ? { headers: { "x-api-key": key } } : {};
}

export async function readText(client: JsonClient, url: string, request?: UpstreamRequest): Promise<string> {
  if (!client.getText) throw new InputError("This HTTP client cannot fetch XML or plain text");
  return client.getText(url, request);
}

export function clipText(value: string | null, max: number): { text: string | null; truncated: boolean } {
  if (!value) return { text: null, truncated: false };
  const clean = value.replace(/\s+/g, " ").trim();
  if (clean.length <= max) return { text: clean, truncated: false };
  return { text: `${clean.slice(0, Math.max(0, max - 1)).trimEnd()}…`, truncated: true };
}

export function isoDate(value: unknown): string | null {
  if (typeof value === "number" && value >= 1000 && value <= 3000) return String(value);
  if (typeof value !== "string") return null;
  const text = value.trim();
  const full = text.match(/^(\d{4})-(\d{2})-(\d{2})/);
  if (full) return `${full[1]}-${full[2]}-${full[3]}`;
  const month = text.match(/^(\d{4})-(\d{2})$/);
  if (month) return `${month[1]}-${month[2]}`;
  if (/^\d{4}$/.test(text)) return text;
  const named = text.match(/(\d{4})\s+([A-Za-z]{3,9})(?:\s+(\d{1,2}))?/);
  if (named) {
    const mm = MONTHS[named[2]!.toLowerCase().slice(0, 3)];
    if (!mm) return named[1]!;
    if (!named[3]) return `${named[1]}-${mm}`;
    return `${named[1]}-${mm}-${named[3].padStart(2, "0")}`;
  }
  const year = text.match(/(\d{4})/);
  return year?.[1] ?? null;
}

export function dateFromParts(node: unknown): string | null {
  const record = asRecord(node);
  const year = plain(record.Year);
  if (!/^\d{4}$/.test(year)) return isoDate(plain(node)) ?? isoDate(node);
  const monthRaw = plain(record.Month).toLowerCase();
  const month = MONTHS[monthRaw.slice(0, 3)] ?? (/^\d{1,2}$/.test(monthRaw) ? monthRaw.padStart(2, "0") : "");
  if (!month) return year;
  const day = plain(record.Day);
  if (!/^\d{1,2}$/.test(day)) return `${year}-${month}`;
  return `${year}-${month}-${day.padStart(2, "0")}`;
}

export function utcDay(offsetDays = 0): string {
  return new Date(Date.now() + offsetDays * 86_400_000).toISOString().slice(0, 10);
}

export function normalizeLicense(value: string | null | undefined): string | null {
  if (!value) return null;
  const raw = value.trim().toLowerCase();
  if (!raw || raw === "na" || raw === "none") return null;
  if (/publicdomain\/zero|\bcc0\b|public domain|cc-zero/.test(raw)) return "cc0";
  const path = raw.match(/creativecommons\.org\/licenses\/([a-z0-9-]+)/);
  const token = (path?.[1] ?? raw).replace(/_/g, "-").replace(/\s+/g, "-");
  const body = token.replace(/^cc-?/, "");
  if (/^by(?:-(?:nc|nd|sa))*$/.test(body)) return `cc-${body}`;
  const loose = raw.replace(/_/g, " ").match(/\bcc(?:\s|-)?by(?:(?:\s|-)(?:nc|nd|sa))*/);
  if (loose) return loose[0].replace(/\s+/g, "-").replace("ccby", "cc-by");
  return raw.length > 120 ? raw.slice(0, 120) : raw;
}

export function allowsFullText(license: string | null | undefined): boolean {
  if (!license) return false;
  return license === "cc0" || license.startsWith("cc-by");
}

export function orcidOf(value: unknown): string | null {
  if (typeof value !== "string") return null;
  const match = value.match(/(\d{4}-\d{4}-\d{4}-\d{3}[\dX])/i);
  return match?.[1] ? `https://orcid.org/${match[1]}` : null;
}

export function stripTags(value: string | null): string | null {
  if (!value) return null;
  const text = value.replace(/<[^>]+>/g, " ").replace(/&nbsp;/g, " ").replace(/\s+/g, " ").trim();
  return text || null;
}

export function unique(values: Array<string | null | undefined>): string[] {
  const seen = new Set<string>();
  const out: string[] = [];
  for (const value of values) {
    if (!value || seen.has(value)) continue;
    seen.add(value);
    out.push(value);
  }
  return out;
}

export function list(value: unknown): unknown[] {
  if (value == null) return [];
  return Array.isArray(value) ? value : [value];
}

export function asRecord(value: unknown): Record<string, unknown> {
  return typeof value === "object" && value !== null && !Array.isArray(value) ? value as Record<string, unknown> : {};
}

export function plain(node: unknown): string {
  if (node == null) return "";
  if (typeof node === "string" || typeof node === "number") return String(node);
  if (Array.isArray(node)) return node.map(plain).filter(Boolean).join(" ");
  if (typeof node === "object") {
    const record = node as Record<string, unknown>;
    if (typeof record["#text"] === "string" || typeof record["#text"] === "number") {
      const extra = Object.entries(record)
        .filter(([key]) => key !== "#text" && !key.startsWith("@"))
        .map(([, value]) => plain(value))
        .filter(Boolean);
      return [String(record["#text"]), ...extra].filter(Boolean).join(" ");
    }
    return Object.entries(record)
      .filter(([key]) => !key.startsWith("@"))
      .map(([, value]) => plain(value))
      .filter(Boolean)
      .join(" ");
  }
  return "";
}

export function extractDoi(value: string): string | null {
  const match = value.match(/10\.\d{4,9}\/[^\s<>"']+/i);
  return match ? normalizeDoi(match[0]) : null;
}

export type ParsedId =
  | { kind: "doi"; value: string }
  | { kind: "pmid"; value: string }
  | { kind: "pmcid"; value: string }
  | { kind: "arxiv"; value: string }
  | { kind: "openalex"; value: string }
  | { kind: "s2"; value: string }
  | { kind: "citation"; value: string };

export function classifyReference(input: string): ParsedId {
  const text = input.trim();
  const doi = extractDoi(text);
  if (doi && (text.length < doi.length + 24 || /doi\.org|^doi:/i.test(text) || text.toLowerCase().startsWith("10."))) {
    return { kind: "doi", value: doi };
  }
  const pmid = text.match(/^(?:pmid|pubmed)\s*:\s*(\d{1,9})$/i)?.[1] ?? (text.match(/^\d{1,9}$/)?.[0] ?? null);
  if (pmid) return { kind: "pmid", value: pmid };
  const pmc = text.match(/^(?:pmcid\s*:\s*)?(PMC\d+)$/i)?.[1];
  if (pmc) return { kind: "pmcid", value: pmc.toUpperCase() };
  const arxivRaw = text.match(/^arxiv:\s*(.+)$/i)?.[1] ?? text;
  const arxiv = /^arxiv:/i.test(text) || /^\d{4}\.\d{4,5}(v\d+)?$/i.test(text) ? normalizeArxiv(arxivRaw) : null;
  if (arxiv && (/^arxiv:/i.test(text) || /^\d{4}\.\d{4,5}/.test(text))) return { kind: "arxiv", value: arxiv };
  const openalex = text.match(/^(?:openalex:)?(W\d+)$/i)?.[1];
  if (openalex && (/^openalex:/i.test(text) || /^W\d+$/i.test(text))) return { kind: "openalex", value: openalex.toUpperCase() };
  const s2 = text.match(/^(?:s2|corpusid)\s*:\s*(\S+)$/i)?.[1];
  if (s2) return { kind: "s2", value: s2 };
  if (doi) return { kind: "doi", value: doi };
  const embeddedPmid = text.match(/\b(?:pmid|pubmed)\s*[:\s]\s*(\d{1,9})\b/i)?.[1];
  if (embeddedPmid && !doi) return { kind: "pmid", value: embeddedPmid };
  return { kind: "citation", value: text };
}

export function normalizePmcid(value: string | null | undefined): string | null {
  if (!value) return null;
  const match = value.toUpperCase().match(/(PMC\d+)/);
  return match?.[1] ?? null;
}

export function openAlexWindow(offset: number, limit: number): { perPage: number; page: number; skip: number } {
  if (offset + limit <= 200) return { perPage: offset + limit, page: 1, skip: offset };
  if (offset % limit === 0) return { perPage: limit, page: offset / limit + 1, skip: 0 };
  throw new InputError("offset + limit must be at most 200, or offset must be a multiple of limit, for this source");
}

export function toFailure(error: unknown, source: Source = literatureSource) {
  if (error instanceof InputError) return failure(source, "invalid_input", error.message, false);
  if (error && typeof error === "object" && "status" in error && "retryable" in error) {
    const status = typeof error.status === "number" ? error.status : null;
    const retryable = error.retryable === true;
    const message = error instanceof Error ? error.message : "upstream request failed";
    const code = status === 404 ? "not_found" : status == null ? "unavailable" : "upstream_http";
    return failure(source, code, message, retryable, { status });
  }
  return failure(source, "unavailable", error instanceof Error ? error.message : "tool failed", false);
}

export function dateRank(item: Publication): string {
  return item.dates.published ?? item.dates.posted ?? item.dates.updated ?? "";
}
