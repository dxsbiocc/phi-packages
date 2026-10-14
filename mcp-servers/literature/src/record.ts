import {
  normalizeArxiv,
  normalizeDoi,
  normalizePmid,
  publicationId,
  type Publication,
} from "phi-omics-schemas";
import { publication } from "phi-omics-runtime";
import { clipText, normalizeLicense, normalizePmcid, orcidOf, unique } from "./util.js";

export interface AuthorDraft {
  name: string;
  orcid: string | null;
}

export interface OaDraft {
  url: string;
  content_type: "html" | "pdf" | "xml" | "landing" | "other";
  version: string | null;
}

export interface PublicationDraft {
  doi?: string | null;
  pmid?: string | null;
  pmcid?: string | null;
  arxiv?: string | null;
  openalex?: string | null;
  s2?: string | null;
  title?: string | null;
  authors?: AuthorDraft[];
  venue?: string | null;
  published?: string | null;
  posted?: string | null;
  updated?: string | null;
  abstract?: string | null;
  type?: Publication["type"];
  is_oa?: boolean | null;
  oa_urls?: Array<OaDraft | null>;
  license?: string | null;
  published_as?: string | null;
  preprint_of?: string | null;
  source_extra?: Record<string, unknown>;
}

function linkId(value: string | null | undefined): string | null {
  if (!value) return null;
  const doi = normalizeDoi(value);
  if (doi) return `doi:${doi}`;
  const pmid = normalizePmid(value.replace(/^pmid:/i, ""));
  if (/^pmid:/i.test(value) && pmid) return `pmid:${pmid}`;
  const arxiv = normalizeArxiv(value);
  if (arxiv && /^arxiv:/i.test(value)) return `arxiv:${arxiv}`;
  return null;
}

export function tryPublication(draft: PublicationDraft): Publication | null {
  const doi = normalizeDoi(draft.doi);
  const pmid = normalizePmid(draft.pmid);
  const arxiv = normalizeArxiv(draft.arxiv);
  let id: string;
  try {
    id = publicationId({ doi, pmid, arxiv });
  } catch {
    return null;
  }
  const pmcid = normalizePmcid(draft.pmcid);
  const openalex = draft.openalex?.match(/W\d+/i)?.[0]?.toUpperCase() ?? null;
  const publishedAs = linkId(draft.published_as);
  const preprintOf = linkId(draft.preprint_of);
  const xrefs = unique([
    pmid ? `pmid:${pmid}` : null,
    pmcid ? `pmcid:${pmcid}` : null,
    arxiv ? `arxiv:${arxiv}` : null,
    openalex ? `openalex:${openalex}` : null,
    draft.s2 ? `s2:${draft.s2}` : null,
  ]).filter((xref) => xref !== id);
  const authors = (draft.authors ?? [])
    .map((author) => ({ name: author.name.trim(), orcid: orcidOf(author.orcid) ?? (author.orcid?.startsWith("https://orcid.org/") ? author.orcid : null) }))
    .filter((author) => author.name.length > 0);
  try {
    return publication({
      id,
      doi,
      xrefs,
      title: draft.title?.replace(/\s+/g, " ").trim() || null,
      authors,
      venue: draft.venue?.trim() || null,
      dates: {
        published: draft.published ?? null,
        posted: draft.posted ?? null,
        updated: draft.updated ?? null,
      },
      abstract: draft.abstract?.trim() || null,
      type: draft.type ?? null,
      is_oa: draft.is_oa ?? null,
      oa_urls: (draft.oa_urls ?? []).filter((item): item is OaDraft => item != null && /^https?:\/\//i.test(item.url)),
      license: normalizeLicense(draft.license),
      published_as: publishedAs && publishedAs !== id ? publishedAs : null,
      preprint_of: preprintOf && preprintOf !== id ? preprintOf : null,
      source_extra: draft.source_extra ?? {},
    });
  } catch {
    return null;
  }
}

export function present(item: Publication, maxAbstract: number, maxAuthors: number): Publication {
  const clipped = clipText(item.abstract, maxAbstract);
  const authors = item.authors.slice(0, maxAuthors);
  const sourceExtra = { ...item.source_extra };
  if (clipped.truncated) sourceExtra.abstract_truncated = true;
  if (item.authors.length > authors.length) sourceExtra.author_count = item.authors.length;
  return publication({ ...item, abstract: clipped.text, authors, source_extra: sourceExtra });
}

function score(item: Publication): number {
  return (item.abstract ? 4 : 0) + (item.license ? 2 : 0) + (item.is_oa != null ? 1 : 0) + Math.min(item.authors.length, 5) + item.xrefs.length;
}

export function combinePublications(left: Publication, right: Publication): Publication {
  const [primary, secondary] = score(left) >= score(right) ? [left, right] : [right, left];
  const doi = primary.doi ?? secondary.doi;
  const pmid = [...primary.xrefs, ...secondary.xrefs, primary.id, secondary.id].find((id) => id.startsWith("pmid:"))?.slice(5) ?? null;
  const arxiv = [...primary.xrefs, ...secondary.xrefs, primary.id, secondary.id].find((id) => id.startsWith("arxiv:"))?.slice(6) ?? null;
  const pmcid = [...primary.xrefs, ...secondary.xrefs].find((id) => id.startsWith("pmcid:"))?.slice(6) ?? null;
  const openalex = [...primary.xrefs, ...secondary.xrefs].find((id) => id.startsWith("openalex:"))?.slice(9) ?? null;
  const s2 = [...primary.xrefs, ...secondary.xrefs].find((id) => id.startsWith("s2:"))?.slice(3) ?? null;
  const providers = unique([
    typeof primary.source_extra.provider === "string" ? primary.source_extra.provider : null,
    typeof secondary.source_extra.provider === "string" ? secondary.source_extra.provider : null,
  ]);
  return tryPublication({
    doi,
    pmid,
    pmcid,
    arxiv,
    openalex,
    s2,
    title: primary.title ?? secondary.title,
    authors: primary.authors.length >= secondary.authors.length ? primary.authors : secondary.authors,
    venue: primary.venue ?? secondary.venue,
    published: primary.dates.published ?? secondary.dates.published,
    posted: primary.dates.posted ?? secondary.dates.posted,
    updated: primary.dates.updated ?? secondary.dates.updated,
    abstract: (primary.abstract?.length ?? 0) >= (secondary.abstract?.length ?? 0) ? primary.abstract : secondary.abstract,
    type: primary.type ?? secondary.type,
    is_oa: primary.is_oa ?? secondary.is_oa,
    oa_urls: [...primary.oa_urls, ...secondary.oa_urls],
    license: primary.license ?? secondary.license,
    published_as: primary.published_as ?? secondary.published_as,
    preprint_of: primary.preprint_of ?? secondary.preprint_of,
    source_extra: { ...secondary.source_extra, ...primary.source_extra, provider: providers.join("+") },
  }) ?? primary;
}

export function mergePublications(groups: Publication[][]): Publication[] {
  const order: string[] = [];
  const byKey = new Map<string, Publication>();
  for (const group of groups) {
    for (const item of group) {
      const key = item.doi ? `doi:${item.doi}` : item.xrefs.find((id) => id.startsWith("pmid:") || id.startsWith("pmcid:")) ?? item.id;
      const existing = byKey.get(key);
      if (!existing) {
        order.push(key);
        byKey.set(key, item);
      } else {
        byKey.set(key, combinePublications(existing, item));
      }
    }
  }
  return order.map((key) => byKey.get(key)!);
}
