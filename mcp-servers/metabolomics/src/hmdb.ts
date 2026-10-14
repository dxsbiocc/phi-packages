import { curie, ontologyCurie } from "phi-omics-schemas";
import {
  failure,
  guard,
  makePage,
  metabolite,
  paging,
  sliceList,
  sourceOf,
  text,
  UpstreamError,
  type JsonClient,
  type OmicsTool,
} from "phi-omics-runtime";
import { z } from "zod";

const hmdb = sourceOf("hmdb", "Human Metabolome Database", "https://hmdb.ca/");
const hmdbUsage = {
  summary: "HMDB records are read from the public per-metabolite XML (https://hmdb.ca/metabolites/HMDBxxxxxxx.xml). There is no public JSON API; https://hmdb.ca/simple/api is contact-only and is not used. Name search reads the public HTML metabolite search and does not download the multi-gigabyte bulk XML. HMDB is free to use as a public resource. Commercial redistribution requires permission of the authors and acknowledgment of HMDB and its papers. Cite HMDB when a substantial portion of the database is used. https://hmdb.ca/citing",
  url: "https://hmdb.ca/citing",
};
const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

export function hmdbAccession(value: string): string | null {
  const match = value.trim().match(/^(?:hmdb:)?HMDB(\d+)$/i);
  if (!match?.[1]) return null;
  const digits = match[1];
  return digits.length > 7 ? `HMDB${digits}` : `HMDB${digits.padStart(7, "0")}`;
}

function decodeXml(value: string): string {
  return value
    .replace(/<!\[CDATA\[([\s\S]*?)\]\]>/g, "$1")
    .replace(/&#(\d+);/g, (_, digits: string) => String.fromCodePoint(Number(digits)))
    .replace(/&#x([0-9a-fA-F]+);/g, (_, hex: string) => String.fromCodePoint(parseInt(hex, 16)))
    .replace(/&lt;/g, "<")
    .replace(/&gt;/g, ">")
    .replace(/&quot;/g, "\"")
    .replace(/&apos;/g, "'")
    .replace(/&amp;/g, "&")
    .trim();
}

function firstScalar(xml: string, tag: string): string | null {
  const match = xml.match(new RegExp(`<(?:[\\w.-]+:)?${tag}(?:\\s[^>]*)?>([\\s\\S]*?)</(?:[\\w.-]+:)?${tag}>`, "i"));
  if (!match?.[1]) return null;
  const inner = match[1].trim();
  if (!inner || inner.startsWith("<")) return null;
  const decoded = decodeXml(inner);
  return decoded.length ? decoded : null;
}

function secondaryAccessions(xml: string): string[] {
  const block = xml.match(/<(?:[\w.-]+:)?secondary_accessions(?:\s[^>]*)?>([\s\S]*?)<\/(?:[\w.-]+:)?secondary_accessions>/i)?.[1] ?? "";
  return [...block.matchAll(/<(?:[\w.-]+:)?accession(?:\s[^>]*)?>([^<]+)<\/(?:[\w.-]+:)?accession>/gi)]
    .map((match) => decodeXml(match[1] ?? ""))
    .filter((item) => item.length > 0)
    .slice(0, 20);
}

function plainDescription(xml: string): string | null {
  const match = xml.match(/<(?:[\w.-]+:)?description(?:\s[^>]*)?>([\s\S]*?)<\/(?:[\w.-]+:)?description>/i);
  if (!match?.[1]) return null;
  const plain = decodeXml(match[1]).replace(/<[^>]+>/g, " ").replace(/\s+/g, " ").trim();
  if (!plain) return null;
  return plain.length > 500 ? `${plain.slice(0, 500)}…` : plain;
}

function xref(prefix: string, value: string | null): string | null {
  if (!value || value.toUpperCase() === "NA") return null;
  if (prefix === "CHEBI") {
    const local = value.replace(/^CHEBI:/i, "");
    return /^\d+$/.test(local) ? ontologyCurie(`CHEBI:${local}`) : null;
  }
  if (/\s/.test(value)) return null;
  return curie(prefix, value);
}

export function metaboliteFromHmdbXml(xml: string) {
  const accession = hmdbAccession(firstScalar(xml, "accession") ?? "");
  if (!accession) return null;
  const xrefs = [
    xref("CHEBI", firstScalar(xml, "chebi_id")),
    xref("kegg", firstScalar(xml, "kegg_id")),
    xref("pubchem", firstScalar(xml, "pubchem_compound_id")),
    xref("chemspider", firstScalar(xml, "chemspider_id")),
    xref("drugbank", firstScalar(xml, "drugbank_id")),
    xref("cas", firstScalar(xml, "cas_registry_number")),
  ].filter((item): item is string => item != null);
  return metabolite({
    id: curie("hmdb", accession),
    name: firstScalar(xml, "name"),
    formula: firstScalar(xml, "chemical_formula"),
    inchikey: firstScalar(xml, "inchikey"),
    xrefs,
    source_extra: {
      status: firstScalar(xml, "status"),
      iupac_name: firstScalar(xml, "iupac_name"),
      smiles: firstScalar(xml, "smiles"),
      inchi: firstScalar(xml, "inchi"),
      average_molecular_weight: firstScalar(xml, "average_molecular_weight"),
      monoisotopic_molecular_weight: firstScalar(xml, "monisotopic_molecular_weight"),
      secondary_accessions: secondaryAccessions(xml),
      description: plainDescription(xml),
    },
  });
}

function challenged(body: string): boolean {
  return /just a moment|cf-mitigated|challenges\.cloudflare\.com/i.test(body.slice(0, 2500));
}

async function readText(client: JsonClient, url: string): Promise<string> {
  if (!client.getText) throw new UpstreamError(null, "this client cannot fetch text", false);
  return client.getText(url);
}

function challengeOrThrow(error: unknown): never {
  if (error instanceof UpstreamError && (error.status === 403 || /just a moment|cloudflare/i.test(error.message))) {
    throw new UpstreamError(403, "HMDB returned a Cloudflare challenge instead of the public record", true);
  }
  throw error;
}

export const hmdbTools: OmicsTool[] = [
  {
    name: "hmdb_get_metabolite",
    title: "Get an HMDB metabolite",
    description: "Fetch one HMDB metabolite from its public XML record, for example HMDB0000001. Results use the shared metabolite schema.",
    inputSchema: z.object({
      hmdb_id: z.string().describe("HMDB accession, with or without the hmdb: prefix. Short legacy ids such as HMDB1 are padded to 7 digits."),
    }),
    call: (client, args) => guard(hmdb, async () => {
      const raw = text(args.hmdb_id);
      const accession = raw ? hmdbAccession(raw) : null;
      if (!accession) return failure(hmdb, "invalid_input", "hmdb_id must look like HMDB0000001", false);
      let xml: string;
      try {
        xml = await readText(client, `https://hmdb.ca/metabolites/${accession}.xml`);
      } catch (error) {
        challengeOrThrow(error);
      }
      if (challenged(xml)) return failure(hmdb, "unavailable", "HMDB returned a Cloudflare challenge instead of XML", true);
      const item = metaboliteFromHmdbXml(xml);
      if (!item) return failure(hmdb, "not_found", `metabolite ${accession} was not found`, false);
      return makePage("metabolite", [item], 0, 1, 1, hmdb, hmdbUsage);
    }),
  },
  {
    name: "hmdb_search_metabolites",
    title: "Search HMDB metabolites",
    description: "Search HMDB metabolites by name using the public HTML search. Each hit is a metabolite id; call hmdb_get_metabolite for formula and InChIKey. Only the first HTML result page is read.",
    inputSchema: z.object({ query: z.string(), ...pageArgs }),
    call: (client, args) => guard(hmdb, async () => {
      const { offset, limit } = paging(args);
      const query = text(args.query);
      if (!query) return failure(hmdb, "invalid_input", "query is required", false);
      const url = new URL("https://hmdb.ca/unearth/q");
      url.searchParams.set("query", query);
      url.searchParams.set("searcher", "metabolites");
      let html: string;
      try {
        html = await readText(client, url.toString());
      } catch (error) {
        challengeOrThrow(error);
      }
      if (challenged(html)) return failure(hmdb, "unavailable", "HMDB returned a Cloudflare challenge instead of search results", true);
      const seen = new Set<string>();
      const found: Array<{ accession: string; name: string | null }> = [];
      for (const match of html.matchAll(/href="\/metabolites\/(HMDB\d+)"[^>]*>([^<]*)</gi)) {
        const accession = hmdbAccession(match[1] ?? "");
        if (!accession || seen.has(accession)) continue;
        seen.add(accession);
        const name = decodeXml(match[2] ?? "").trim();
        found.push({ accession, name: name.length ? name : null });
      }
      const { rows, total } = sliceList(found, offset, limit);
      const items = rows.map((row) => metabolite({
        id: curie("hmdb", row.accession),
        name: row.name,
        source_extra: { search_query: query },
      }));
      return makePage("metabolite", items, offset, limit, total, hmdb, hmdbUsage, {
        search: "html",
        note: "Formula and InChIKey are filled by hmdb_get_metabolite. HMDB does not publish a JSON search API.",
      });
    }),
  },
];
