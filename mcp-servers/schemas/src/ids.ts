import { z } from "zod";

/** prefix:local, with dotted prefixes such as gdc.case or gtex.tissue. */
export const CURIE_PATTERN = /^[A-Za-z][A-Za-z0-9]*(\.[A-Za-z][A-Za-z0-9]*)*:[^\s]+$/;

export const CurieSchema = z.string().regex(CURIE_PATTERN, "expected a CURIE (prefix:local)");

const ONTOLOGY_PREFIX: Record<string, string> = {
  UBERON: "UBERON",
  CL: "CL",
  EFO: "EFO",
  MONDO: "MONDO",
  PATO: "PATO",
  NCBITAXON: "NCBITaxon",
  CHEBI: "CHEBI",
  HGNC: "HGNC",
  GO: "GO",
  DOID: "DOID",
  HP: "HP",
};

const SPECIES: Record<string, string> = {
  "homo sapiens": "NCBITaxon:9606",
  human: "NCBITaxon:9606",
  "mus musculus": "NCBITaxon:10090",
  mouse: "NCBITaxon:10090",
};

export function curie(prefix: string, local: string): string {
  const id = `${prefix}:${local}`;
  if (!CURIE_PATTERN.test(id)) {
    throw new Error(`invalid CURIE ${id}`);
  }
  return id;
}

export function ontologyCurie(term: string | null | undefined): string | null {
  if (!term) return null;
  const match = term.trim().replace(/_/g, ":").match(/^([A-Za-z]+):(\d+)$/);
  if (!match) return null;
  const prefix = ONTOLOGY_PREFIX[match[1].toUpperCase()];
  if (!prefix) return null;
  return `${prefix}:${match[2]}`;
}

export function taxonCurie(label: string | null | undefined): string | null {
  if (!label) return null;
  const ontology = ontologyCurie(label);
  if (ontology?.startsWith("NCBITaxon:")) return ontology;
  return SPECIES[label.trim().toLowerCase()] ?? null;
}

export function ensemblGeneCurie(id: string): { id: string; version: string | null } {
  const match = id.trim().match(/^(?:ensembl:)?(ENSG\d+)(?:\.(\d+))?$/i);
  if (!match) throw new Error(`not an Ensembl gene id: ${id}`);
  return { id: `ensembl:${match[1]}`, version: match[2] ?? null };
}

/** Human chromosome with a chr prefix. Coordinates elsewhere are 1-based and inclusive. */
export function humanChrom(value: string | number | null | undefined): string | null {
  if (value == null || value === "") return null;
  let text = String(value).trim();
  if (/^chr/i.test(text)) text = text.slice(3);
  if (!text) return null;
  const body = text.toUpperCase() === "MT" ? "M" : text.toUpperCase();
  return `chr${body}`;
}

export function strandOf(value: unknown): "+" | "-" | "." | null {
  if (value === "+" || value === "-" || value === ".") return value;
  if (value === 1 || value === "1") return "+";
  if (value === -1 || value === "-1") return "-";
  return null;
}

/** Bare lowercased DOI, without a resolver URL or `doi:` prefix. */
export function normalizeDoi(value: string | null | undefined): string | null {
  if (!value) return null;
  const text = value.trim()
    .replace(/^https?:\/\/(?:dx\.)?doi\.org\//i, "")
    .replace(/^doi:/i, "")
    .replace(/[.,;]+$/, "");
  if (!/^10\.\d{4,9}\/\S+$/i.test(text)) return null;
  return text.toLowerCase();
}

/** PMID digits, without a `pmid:` prefix. */
export function normalizePmid(value: string | null | undefined): string | null {
  if (!value) return null;
  const match = value.trim().match(/^(?:pmid:)?(\d{1,9})$/i);
  return match?.[1] ?? null;
}

/** arXiv id without a version suffix or URL. */
export function normalizeArxiv(value: string | null | undefined): string | null {
  if (!value) return null;
  const text = value.trim()
    .replace(/^https?:\/\/arxiv\.org\/abs\//i, "")
    .replace(/^arxiv:/i, "")
    .replace(/v\d+$/i, "");
  if (!/^(\d{4}\.\d{4,5}|[a-z-]+(?:\.[a-z]{2})?\/\d{7})$/i.test(text)) return null;
  return text.toLowerCase();
}

/**
 * Publication primary key. DOI wins, including online-first articles that have
 * no PMID. Otherwise PMID, then arXiv. Preprint and published version stay
 * separate records.
 */
export function publicationId(ids: { doi?: string | null; pmid?: string | null; arxiv?: string | null }): string {
  const doi = normalizeDoi(ids.doi);
  if (doi) return curie("doi", doi);
  const pmid = normalizePmid(ids.pmid);
  if (pmid) return curie("pmid", pmid);
  const arxiv = normalizeArxiv(ids.arxiv);
  if (arxiv) return curie("arxiv", arxiv);
  throw new Error("publication id requires a DOI, PMID, or arXiv id");
}

export function sexOf(value: unknown): "female" | "male" | "unknown" | null {
  if (typeof value !== "string" || !value.trim()) return null;
  const text = value.trim().toLowerCase();
  if (text === "female" || text === "f") return "female";
  if (text === "male" || text === "m") return "male";
  if (text === "unknown" || text === "unspecified" || text === "not reported") return "unknown";
  return null;
}
