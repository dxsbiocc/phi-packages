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

export function sexOf(value: unknown): "female" | "male" | "unknown" | null {
  if (typeof value !== "string" || !value.trim()) return null;
  const text = value.trim().toLowerCase();
  if (text === "female" || text === "f") return "female";
  if (text === "male" || text === "m") return "male";
  if (text === "unknown" || text === "unspecified" || text === "not reported") return "unknown";
  return null;
}
