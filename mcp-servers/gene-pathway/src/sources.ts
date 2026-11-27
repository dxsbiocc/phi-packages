import { curie, ensemblGeneCurie, humanChrom, strandOf, taxonCurie } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  dataset,
  failure,
  gene,
  guard,
  type JsonClient,
  locusOf,
  makePage,
  numberOrNull,
  type OmicsTool,
  paging,
  safeMap,
  sliceList,
  sourceOf,
  text,
  variant,
} from "phi-omics-runtime";
import { z } from "zod";

const ensembl = sourceOf("ensembl", "Ensembl REST", "https://rest.ensembl.org/", "GRCh38");
const ensemblUsage = {
  summary: "Ensembl REST gene lookups, comparative homology (orthologues/paralogues), and Variant Effect Predictor (VEP) annotations.",
  url: "https://rest.ensembl.org/documentation/info/terms_and_conditions",
};
const reactome = sourceOf("reactome", "Reactome", "https://reactome.org/", "ContentService");
const reactomeUsage = {
  summary: "Reactome curated human and model-organism biological pathways, events, literature links, and UniProt mappings (CC0).",
  url: "https://reactome.org/license",
};

const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

function stripHighlights(value: string | null): string | null {
  if (!value) return null;
  return value.replace(/<[^>]+>/g, "").replace(/\s+/g, " ").trim() || null;
}

function safeEnsemblId(id: string): string {
  try {
    return ensemblGeneCurie(id).id;
  } catch {
    return curie("ensembl", id.trim().replace(/\s+/g, "_"));
  }
}

function speciesTaxon(species: string | null): string | null {
  if (!species) return null;
  const clean = species.replace(/_/g, " ");
  return taxonCurie(clean);
}

function ensemblJson(client: JsonClient, url: string) {
  return client.getJson(url, { accept: "application/json" });
}

export const pathwayTools: OmicsTool[] = [
  {
    name: "ensembl_get_gene",
    title: "Get one Ensembl gene",
    description: "Look up an Ensembl gene by Ensembl ID (e.g. ENSG00000141510) or gene symbol + species (default homo_sapiens), including transcripts and genomic coordinates.",
    inputSchema: z.object({
      id: z.string().optional().describe("Ensembl ID, e.g. ENSG00000141510."),
      symbol: z.string().optional().describe("Gene symbol, e.g. TP53."),
      species: z.string().optional().describe("Ensembl species name for symbol lookup, default homo_sapiens."),
    }),
    call: (client, args) => guard(ensembl, async () => {
      const id = text(args.id)?.replace(/^ensembl:/i, "");
      const symbol = text(args.symbol);
      const species = (text(args.species) ?? "homo_sapiens").toLowerCase().replace(/\s+/g, "_");
      if (!id && !symbol) return failure(ensembl, "invalid_input", "Pass id or symbol", false);
      const url = id
        ? `https://rest.ensembl.org/lookup/id/${encodeURIComponent(id)}?expand=1`
        : `https://rest.ensembl.org/lookup/symbol/${encodeURIComponent(species)}/${encodeURIComponent(symbol!)}?expand=1`;
      const row = asRecord(await ensemblJson(client, url));
      const gid = text(row.id);
      if (!gid) return failure(ensembl, "not_found", `Gene ${id ?? symbol} was not found`, false);
      const start = numberOrNull(row.start);
      const end = numberOrNull(row.end);
      const transcripts = asArray(row.Transcript).map(asRecord).slice(0, 10).map((t) => ({
        id: text(t.id),
        biotype: text(t.biotype),
        is_canonical: t.is_canonical === 1 || t.is_canonical === true,
        start: numberOrNull(t.start),
        end: numberOrNull(t.end),
      }));
      const item = gene({
        id: safeEnsemblId(gid),
        symbol: text(row.display_name) ?? symbol ?? null,
        name: text(row.description),
        taxon_id: speciesTaxon(text(row.species) ?? species),
        biotype: text(row.biotype),
        locus: start != null ? locusOf({
          genomeBuild: text(row.assembly_name),
          chrom: humanChrom(text(row.seq_region_name)),
          start,
          end: end != null && end >= start ? end : start,
          strand: strandOf(row.strand),
        }) : null,
        source_extra: {
          version: numberOrNull(row.version),
          canonical_transcript: text(row.canonical_transcript),
          species: text(row.species),
          transcripts,
        },
      });
      return makePage("gene", [item], 0, 1, 1, ensembl, ensemblUsage);
    }),
  },
  {
    name: "ensembl_get_homologs",
    title: "Get Ensembl homologs",
    description: "Fetch orthologues and/or paralogues for one Ensembl gene ID across species (Ensembl Compara).",
    inputSchema: z.object({
      gene_id: z.string().describe("Ensembl gene ID, e.g. ENSG00000141510."),
      species: z.string().optional().describe("Source species, default homo_sapiens."),
      target_species: z.string().optional().describe("Optional target species, e.g. mus_musculus."),
      type: z.enum(["orthologues", "paralogues", "all"]).optional().describe("Default orthologues."),
      ...pageArgs,
    }),
    call: (client, args) => guard(ensembl, async () => {
      const { offset, limit } = paging(args);
      const geneId = text(args.gene_id)?.replace(/^ensembl:/i, "");
      if (!geneId) return failure(ensembl, "invalid_input", "gene_id is required", false);
      const species = (text(args.species) ?? "homo_sapiens").toLowerCase().replace(/\s+/g, "_");
      const type = text(args.type) ?? "orthologues";
      const url = new URL(`https://rest.ensembl.org/homology/id/${encodeURIComponent(species)}/${encodeURIComponent(geneId)}`);
      url.searchParams.set("type", type);
      url.searchParams.set("format", "condensed");
      const targetSpecies = text(args.target_species);
      if (targetSpecies) url.searchParams.set("target_species", targetSpecies.toLowerCase().replace(/\s+/g, "_"));
      const body = asRecord(await ensemblJson(client, url.toString()));
      const entry = asRecord(asArray(body.data)[0]);
      const homologies = asArray(entry.homologies).map(asRecord);
      const { rows, total } = sliceList(homologies, offset, limit);
      const items = safeMap(rows, (row) => {
        const target = asRecord(row.target);
        const tid = text(target.id) ?? "unknown";
        const sp = text(target.species);
        return gene({
          id: safeEnsemblId(tid),
          taxon_id: speciesTaxon(sp),
          xrefs: text(target.protein_id) ? [curie("ensembl.protein", text(target.protein_id)!)] : [],
          source_extra: {
            homology_type: text(row.type),
            taxonomy_level: text(row.taxonomy_level),
            species: sp,
            protein_id: text(target.protein_id),
            perc_id: numberOrNull(target.perc_id),
            source_gene: safeEnsemblId(geneId),
          },
        });
      });
      return makePage("gene", items, offset, limit, total, ensembl, ensemblUsage);
    }),
  },
  {
    name: "ensembl_vep_hgvs",
    title: "Predict variant consequences via Ensembl VEP",
    description: "Run Ensembl Variant Effect Predictor (VEP) for an HGVS notation (e.g. ENST00000269305.9:c.818G>A) or dbSNP rsID (e.g. rs121912651).",
    inputSchema: z.object({
      variant: z.string().describe("HGVS expression or dbSNP rsID."),
      species: z.string().optional().describe("Default homo_sapiens."),
    }),
    call: (client, args) => guard(ensembl, async () => {
      const notation = text(args.variant);
      if (!notation) return failure(ensembl, "invalid_input", "variant is required", false);
      const species = (text(args.species) ?? "homo_sapiens").toLowerCase().replace(/\s+/g, "_");
      const isRsid = /^rs\d+$/i.test(notation);
      const url = isRsid
        ? `https://rest.ensembl.org/vep/${encodeURIComponent(species)}/id/${encodeURIComponent(notation)}?canonical=1&SIFT=b&PolyPhen=b`
        : `https://rest.ensembl.org/vep/${encodeURIComponent(species)}/hgvs/${encodeURIComponent(notation)}?canonical=1&SIFT=b&PolyPhen=b`;
      const rows = asArray(await ensemblJson(client, url)).map(asRecord);
      if (rows.length === 0) return failure(ensembl, "not_found", `No VEP predictions found for ${notation}`, false);
      const items = safeMap(rows, (row) => {
        const consequences = asArray(row.transcript_consequences).map(asRecord);
        const canonical = consequences.find((c) => c.canonical === 1) ?? consequences[0];
        const colocated = asArray(row.colocated_variants).map(asRecord);
        const rsids = colocated.map((c) => text(c.id)).filter((s): s is string => s != null && /^rs\d+$/i.test(s));
        const start = numberOrNull(row.start);
        const end = numberOrNull(row.end);
        const cleanKey = encodeURIComponent((text(row.id) ?? notation).replace(/\s+/g, "_"));
        return variant({
          id: curie("ensembl.vep", cleanKey),
          gene_id: text(canonical?.gene_id) ? safeEnsemblId(text(canonical?.gene_id)!) : null,
          locus: start != null ? locusOf({
            genomeBuild: text(row.assembly_name),
            chrom: humanChrom(text(row.seq_region_name)),
            start,
            end: end != null && end >= start ? end : start,
            strand: strandOf(row.strand),
          }) : null,
          consequence: text(row.most_severe_consequence),
          xrefs: rsids.map((r) => curie("dbsnp", r)),
          source_extra: {
            input: text(row.input) ?? notation,
            allele_string: text(row.allele_string),
            variant_class: text(row.variant_class),
            transcript_consequences: consequences.slice(0, 12).map((c) => ({
              gene_id: text(c.gene_id),
              gene_symbol: text(c.gene_symbol),
              transcript_id: text(c.transcript_id),
              biotype: text(c.biotype),
              impact: text(c.impact),
              canonical: c.canonical === 1,
              consequence_terms: asArray(c.consequence_terms).map(text).filter(Boolean),
              sift_prediction: text(c.sift_prediction),
              polyphen_prediction: text(c.polyphen_prediction),
              amino_acids: text(c.amino_acids),
              protein_start: numberOrNull(c.protein_start),
            })),
          },
        });
      });
      return makePage("variant", items, 0, Math.max(items.length, 1), items.length, ensembl, ensemblUsage);
    }),
  },
  {
    name: "reactome_search_pathways",
    title: "Search Reactome pathways",
    description: "Search Reactome pathways by keyword, gene, or process (default species: Homo sapiens).",
    inputSchema: z.object({
      query: z.string().describe("Search keyword, e.g. TP53 or apoptosis."),
      species: z.string().optional().describe("Default Homo sapiens."),
      ...pageArgs,
    }),
    call: (client, args) => guard(reactome, async () => {
      const { offset, limit } = paging(args);
      const q = text(args.query);
      if (!q) return failure(reactome, "invalid_input", "query is required", false);
      const sp = text(args.species) ?? "Homo sapiens";
      const url = new URL("https://reactome.org/ContentService/search/query");
      url.searchParams.set("query", q);
      url.searchParams.set("species", sp);
      url.searchParams.set("types", "Pathway");
      url.searchParams.set("cluster", "true");
      const body = asRecord(await client.getJson(url.toString()));
      const entries = asArray(body.results).flatMap((group) => asArray(asRecord(group).entries).map(asRecord));
      const { rows, total } = sliceList(entries, offset, limit);
      const items = safeMap(rows, (row) => {
        const stId = text(row.stId) ?? text(row.id) ?? "unknown";
        const speciesNames = asArray(row.species).map(text).filter((s): s is string => s != null);
        return dataset({
          id: curie("reactome", stId),
          title: stripHighlights(text(row.name)),
          description: stripHighlights(text(row.summation)),
          taxon_id: speciesTaxon(speciesNames[0] ?? sp),
          access: "open",
          url: `https://reactome.org/content/detail/${encodeURIComponent(stId)}`,
          source_extra: {
            species: speciesNames,
            compartments: asArray(row.compartmentNames).map(text).filter(Boolean),
            in_disease: row.isDisease ?? null,
          },
        });
      });
      return makePage("dataset", items, offset, limit, total, reactome, reactomeUsage);
    }),
  },
  {
    name: "reactome_get_pathway",
    title: "Get one Reactome pathway",
    description: "Fetch detailed Reactome pathway metadata, sub-events, literature references (PMIDs), and GO biological process by stable ID (e.g. R-HSA-3700989).",
    inputSchema: z.object({
      pathway_id: z.string().describe("Reactome stable identifier, e.g. R-HSA-3700989."),
    }),
    call: (client, args) => guard(reactome, async () => {
      const id = text(args.pathway_id)?.replace(/^reactome:/i, "");
      if (!id || !/^R-[A-Z]{3}-\d+$/i.test(id)) {
        return failure(reactome, "invalid_input", "pathway_id must look like R-HSA-3700989", false);
      }
      const row = asRecord(await client.getJson(`https://reactome.org/ContentService/data/query/enhanced/${encodeURIComponent(id)}`));
      const stId = text(row.stId) ?? id;
      if (!text(row.displayName) && !text(row.stId)) {
        return failure(reactome, "not_found", `Pathway ${id} was not found`, false);
      }
      const summations = asArray(row.summation).map((s) => text(asRecord(s).text)).filter((s): s is string => s != null);
      const pmids = asArray(row.literatureReference)
        .map((r) => numberOrNull(asRecord(r).pubMedIdentifier))
        .filter((n): n is number => n != null);
      const events = asArray(row.hasEvent).map(asRecord).slice(0, 20).map((e) => ({
        stId: text(e.stId),
        displayName: text(e.displayName),
        schemaClass: text(e.schemaClass),
      }));
      const speciesName = text(asRecord(asArray(row.species)[0]).displayName);
      return makePage("dataset", [dataset({
        id: curie("reactome", stId),
        title: text(row.displayName),
        description: summations[0] ?? null,
        taxon_id: speciesTaxon(speciesName),
        access: "open",
        url: `https://reactome.org/content/detail/${encodeURIComponent(stId)}`,
        xrefs: pmids.slice(0, 15).map((p) => curie("pmid", String(p))),
        source_extra: {
          schema_class: text(row.schemaClass),
          species: speciesName,
          has_diagram: row.hasDiagram ?? null,
          is_in_disease: row.isInDisease ?? null,
          go_biological_process: text(asRecord(row.goBiologicalProcess).displayName),
          sub_events: events,
        },
      })], 0, 1, 1, reactome, reactomeUsage);
    }),
  },
  {
    name: "reactome_protein_pathways",
    title: "List Reactome pathways for a UniProt protein",
    description: "Map a UniProt accession (e.g. P04637) to its curated Reactome pathways.",
    inputSchema: z.object({
      uniprot_accession: z.string().describe("UniProt accession, e.g. P04637."),
      species: z.string().optional().describe("NCBI Taxonomy ID or species name, default 9606."),
      ...pageArgs,
    }),
    call: (client, args) => guard(reactome, async () => {
      const { offset, limit } = paging(args);
      const acc = text(args.uniprot_accession)?.replace(/^uniprot:/i, "").toUpperCase();
      if (!acc) return failure(reactome, "invalid_input", "uniprot_accession is required", false);
      const species = text(args.species) ?? "9606";
      const url = `https://reactome.org/ContentService/data/mapping/UniProt/${encodeURIComponent(acc)}/pathways?species=${encodeURIComponent(species)}`;
      const rowsRaw = asArray(await client.getJson(url)).map(asRecord);
      const { rows, total } = sliceList(rowsRaw, offset, limit);
      const items = safeMap(rows, (row) => {
        const stId = text(row.stId) ?? "unknown";
        return dataset({
          id: curie("reactome", stId),
          title: text(row.displayName),
          taxon_id: taxonCurie(text(row.speciesName)) ?? (/^\d+$/.test(species) ? curie("NCBITaxon", species) : null),
          access: "open",
          url: `https://reactome.org/content/detail/${encodeURIComponent(stId)}`,
          xrefs: [curie("uniprot", acc)],
          source_extra: {
            has_diagram: row.hasDiagram ?? null,
            is_in_disease: row.isInDisease ?? null,
            species: text(row.speciesName),
          },
        });
      });
      return makePage("dataset", items, offset, limit, total, reactome, reactomeUsage);
    }),
  },
];
