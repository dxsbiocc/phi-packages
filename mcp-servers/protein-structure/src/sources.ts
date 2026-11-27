import { curie } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  failure,
  file,
  guard,
  makePage,
  numberOrNull,
  type OmicsTool,
  paging,
  protein,
  safeMap,
  sourceOf,
  text,
} from "phi-omics-runtime";
import { z } from "zod";

const uniprot = sourceOf("uniprot", "UniProtKB", "https://www.uniprot.org/", "REST");
const uniprotUsage = {
  summary: "UniProtKB protein entries, function, domains, active sites, disease links, and PDB cross-references (CC-BY 4.0).",
  url: "https://www.uniprot.org/help/license",
};
const alphafold = sourceOf("alphafold", "AlphaFold DB", "https://alphafold.ebi.ac.uk/", "v4");
const alphafoldUsage = {
  summary: "AlphaFold Protein Structure Database predicted models and per-residue confidence metrics (CC-BY 4.0).",
  url: "https://alphafold.ebi.ac.uk/terms-of-use",
};

const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

function cleanAccession(raw: string): string {
  return raw.trim().replace(/^uniprot:/i, "").toUpperCase();
}

function mapUniprotEntry(row: Record<string, unknown>, includeSequence = false) {
  const acc = text(row.primaryAccession) ?? "unknown";
  const protDesc = asRecord(row.proteinDescription);
  const recName = asRecord(protDesc.recommendedName);
  const subName = asRecord(asArray(protDesc.submissionNames)[0]);
  const fullName = text(asRecord(recName.fullName).value) ?? text(asRecord(subName.fullName).value) ?? text(row.uniProtkbId);
  const genes = asArray(row.genes).map(asRecord);
  const geneSymbol = text(asRecord(genes[0]?.geneName).value);
  const organism = asRecord(row.organism);
  const taxonId = numberOrNull(organism.taxonId);

  const xrefsRaw = asArray(row.uniProtKBCrossReferences).map(asRecord);
  const pdbIds = xrefsRaw
    .filter((x) => text(x.database) === "PDB" && text(x.id))
    .map((x) => text(x.id)!);
  const ensemblGene = xrefsRaw.find((x) => text(x.database) === "Ensembl")
    ? (() => {
        const ref = xrefsRaw.find((x) => text(x.database) === "Ensembl")!;
        const props = asArray(ref.properties).map(asRecord);
        const gid = props.find((p) => text(p.key) === "GeneId")?.value;
        const clean = text(gid)?.split(".")[0];
        return clean && /^ENSG\d+/i.test(clean) ? `ensembl:${clean}` : null;
      })()
    : null;

  const comments = asArray(row.comments).map(asRecord);
  const functionComments = comments
    .filter((c) => text(c.commentType) === "FUNCTION")
    .flatMap((c) => asArray(c.texts).map((t) => text(asRecord(t).value)))
    .filter((s): s is string => s != null);
  const diseaseComments = comments
    .filter((c) => text(c.commentType) === "DISEASE")
    .map((c) => {
      const dis = asRecord(c.disease);
      return {
        disease_id: text(dis.diseaseId),
        acronym: text(dis.acronym),
        description: text(dis.description),
      };
    });

  const features = asArray(row.features).map(asRecord);
  const domainsAndRegions = features
    .filter((f) => ["Domain", "Region", "Motif", "Repeat"].includes(text(f.type) ?? ""))
    .slice(0, 20)
    .map((f) => ({
      type: text(f.type),
      description: text(f.description),
      start: numberOrNull(asRecord(asRecord(f.location).start).value),
      end: numberOrNull(asRecord(asRecord(f.location).end).value),
    }));
  const activeAndBindingSites = features
    .filter((f) => ["Active site", "Binding site", "Site"].includes(text(f.type) ?? ""))
    .slice(0, 20)
    .map((f) => ({
      type: text(f.type),
      description: text(f.description),
      start: numberOrNull(asRecord(asRecord(f.location).start).value),
      end: numberOrNull(asRecord(asRecord(f.location).end).value),
    }));

  const seqObj = asRecord(row.sequence);
  const xrefs: string[] = [
    ...pdbIds.slice(0, 15).map((id) => curie("pdb", id)),
  ];
  if (geneSymbol) xrefs.unshift(curie("hgnc.symbol", geneSymbol));

  return protein({
    id: curie("uniprot", acc),
    name: fullName,
    gene_id: ensemblGene ?? (geneSymbol ? curie("hgnc.symbol", geneSymbol) : null),
    gene_symbol: geneSymbol,
    taxon_id: taxonId != null ? curie("NCBITaxon", String(taxonId)) : null,
    xrefs,
    source_extra: {
      uniprot_id: text(row.uniProtkbId),
      entry_type: text(row.entryType),
      organism: text(organism.scientificName),
      sequence_length: numberOrNull(seqObj.length),
      function_comments: functionComments.slice(0, 5),
      disease_associations: diseaseComments.slice(0, 10),
      domains_and_regions: domainsAndRegions,
      active_and_binding_sites: activeAndBindingSites,
      pdb_ids: pdbIds.slice(0, 20),
      ...(includeSequence && text(seqObj.value) ? { sequence: text(seqObj.value) } : {}),
    },
  });
}

export const structureTools: OmicsTool[] = [
  {
    name: "uniprot_search_proteins",
    title: "Search UniProtKB proteins",
    description: "Search UniProtKB proteins by gene, protein name, or query syntax. Optional reviewed_only and organism_id (e.g. 9606 for human).",
    inputSchema: z.object({
      query: z.string().describe("UniProt query, e.g. TP53 or kinase."),
      reviewed_only: z.boolean().optional().describe("Default false. Set true for Swiss-Prot reviewed entries only."),
      organism_id: z.number().int().positive().optional().describe("NCBI Taxonomy ID filter, e.g. 9606."),
      ...pageArgs,
    }),
    call: (client, args) => guard(uniprot, async () => {
      const { offset, limit } = paging(args);
      const rawQuery = text(args.query);
      if (!rawQuery) return failure(uniprot, "invalid_input", "query is required", false);
      if (offset + limit > 100) {
        return failure(uniprot, "invalid_input", "offset + limit must be at most 100 for UniProt search", false);
      }
      const clauses = [`(${rawQuery})`];
      if (args.reviewed_only === true) clauses.push("reviewed:true");
      if (typeof args.organism_id === "number") clauses.push(`organism_id:${args.organism_id}`);
      const url = new URL("https://rest.uniprot.org/uniprotkb/search");
      url.searchParams.set("query", clauses.join(" AND "));
      url.searchParams.set("format", "json");
      url.searchParams.set("size", String(offset + limit));
      const body = asRecord(await client.getJson(url.toString()));
      const rows = asArray(body.results).map(asRecord).slice(offset, offset + limit);
      const items = safeMap(rows, (row) => mapUniprotEntry(row, false));
      const total = rows.length < limit ? offset + rows.length : null;
      return makePage("protein", items, offset, limit, total, uniprot, uniprotUsage, {
        query: clauses.join(" AND "),
      });
    }),
  },
  {
    name: "uniprot_get_protein",
    title: "Get one UniProtKB protein",
    description: "Fetch one UniProtKB entry by accession (e.g. P04637), returning function, disease links, domains, active sites, and PDB cross-references.",
    inputSchema: z.object({
      accession: z.string().describe("UniProt accession, e.g. P04637 or uniprot:P04637."),
      include_sequence: z.boolean().optional().describe("Default false. Include amino acid sequence when true."),
    }),
    call: (client, args) => guard(uniprot, async () => {
      const raw = text(args.accession);
      if (!raw) return failure(uniprot, "invalid_input", "accession is required", false);
      const acc = cleanAccession(raw);
      if (!/^[A-Z0-9]{6,10}$/.test(acc)) {
        return failure(uniprot, "invalid_input", "accession must be a valid UniProt accession like P04637", false);
      }
      const body = asRecord(await client.getJson(`https://rest.uniprot.org/uniprotkb/${encodeURIComponent(acc)}.json`));
      if (!text(body.primaryAccession)) {
        return failure(uniprot, "not_found", `UniProt entry ${acc} was not found`, false);
      }
      return makePage("protein", [mapUniprotEntry(body, args.include_sequence === true)], 0, 1, 1, uniprot, uniprotUsage);
    }),
  },
  {
    name: "alphafold_get_prediction",
    title: "Get AlphaFold prediction metadata",
    description: "Fetch AlphaFold DB structure prediction metadata, global/per-category pLDDT confidence scores, and coordinate URLs for a UniProt accession.",
    inputSchema: z.object({
      uniprot_accession: z.string().describe("UniProt accession, e.g. P04637."),
    }),
    call: (client, args) => guard(alphafold, async () => {
      const raw = text(args.uniprot_accession);
      if (!raw) return failure(alphafold, "invalid_input", "uniprot_accession is required", false);
      const acc = cleanAccession(raw);
      const rows = asArray(await client.getJson(`https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(acc)}`)).map(asRecord);
      if (rows.length === 0) return failure(alphafold, "not_found", `No AlphaFold predictions found for ${acc}`, false);
      const items = safeMap(rows, (row) => {
        const entryId = text(row.entryId) ?? `AF-${acc}-F1`;
        const geneSym = text(row.gene);
        const taxId = numberOrNull(row.taxId);
        return protein({
          id: curie("alphafold", entryId),
          name: text(row.uniprotDescription),
          gene_id: geneSym ? curie("hgnc.symbol", geneSym) : null,
          gene_symbol: geneSym,
          taxon_id: taxId != null ? curie("NCBITaxon", String(taxId)) : null,
          xrefs: [curie("uniprot", text(row.uniprotAccession) ?? acc)],
          source_extra: {
            model_entity_id: text(row.modelEntityId),
            global_metric_value: numberOrNull(row.globalMetricValue),
            fraction_plddt_very_high: numberOrNull(row.fractionPlddtVeryHigh),
            fraction_plddt_confident: numberOrNull(row.fractionPlddtConfident),
            fraction_plddt_low: numberOrNull(row.fractionPlddtLow),
            fraction_plddt_very_low: numberOrNull(row.fractionPlddtVeryLow),
            sequence_start: numberOrNull(row.uniprotStart),
            sequence_end: numberOrNull(row.uniprotEnd),
            latest_version: numberOrNull(row.latestVersion),
            cif_url: text(row.cifUrl),
            pdb_url: text(row.pdbUrl),
            bcif_url: text(row.bcifUrl),
            pae_doc_url: text(row.paeDocUrl),
            pae_image_url: text(row.paeImageUrl),
          },
        });
      });
      return makePage("protein", items, 0, Math.max(items.length, 1), items.length, alphafold, alphafoldUsage);
    }),
  },
  {
    name: "alphafold_list_files",
    title: "List AlphaFold structure files",
    description: "List downloadable AlphaFold DB structure and confidence files (mmCIF, PDB, BinaryCIF, PAE JSON/PNG) as open file entities.",
    inputSchema: z.object({
      uniprot_accession: z.string().describe("UniProt accession, e.g. P04637."),
      ...pageArgs,
    }),
    call: (client, args) => guard(alphafold, async () => {
      const { offset, limit } = paging(args);
      const raw = text(args.uniprot_accession);
      if (!raw) return failure(alphafold, "invalid_input", "uniprot_accession is required", false);
      const acc = cleanAccession(raw);
      const rows = asArray(await client.getJson(`https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(acc)}`)).map(asRecord);
      if (rows.length === 0) return failure(alphafold, "not_found", `No AlphaFold files found for ${acc}`, false);
      const allFiles = rows.flatMap((row) => {
        const entryId = text(row.entryId) ?? `AF-${acc}-F1`;
        const candidates: Array<{ key: string; format: string; dataType: string; url: string | null }> = [
          { key: "cif", format: "mmCIF", dataType: "coordinates", url: text(row.cifUrl) },
          { key: "pdb", format: "PDB", dataType: "coordinates", url: text(row.pdbUrl) },
          { key: "bcif", format: "BinaryCIF", dataType: "coordinates", url: text(row.bcifUrl) },
          { key: "pae_json", format: "JSON", dataType: "predicted_aligned_error", url: text(row.paeDocUrl) },
          { key: "pae_png", format: "PNG", dataType: "predicted_aligned_error_image", url: text(row.paeImageUrl) },
        ];
        return candidates.flatMap((c) => c.url ? [file({
          id: curie("alphafold.file", `${entryId}:${c.key}`),
          name: c.url.split("/").pop() ?? `${entryId}.${c.key}`,
          format: c.format,
          data_type: c.dataType,
          access: "open",
          url: c.url,
          dataset_id: curie("alphafold", entryId),
          source_extra: { entry_id: entryId, uniprot_accession: acc },
        })] : []);
      });
      const sliced = allFiles.slice(offset, offset + limit);
      return makePage("file", sliced, offset, limit, allFiles.length, alphafold, alphafoldUsage);
    }),
  },
  {
    name: "alphafold_get_summary",
    title: "Get AlphaFold UniProt summary",
    description: "Fetch the AlphaFold DB summary (model count, sequence length, coverage, and confidence range) for one UniProt accession.",
    inputSchema: z.object({
      uniprot_accession: z.string().describe("UniProt accession, e.g. P04637."),
    }),
    call: (client, args) => guard(alphafold, async () => {
      const raw = text(args.uniprot_accession);
      if (!raw) return failure(alphafold, "invalid_input", "uniprot_accession is required", false);
      const acc = cleanAccession(raw);
      const body = asRecord(await client.getJson(`https://alphafold.ebi.ac.uk/api/uniprot/summary/${encodeURIComponent(acc)}.json`));
      const entry = asRecord(body.uniprot_entry);
      const structures = asArray(body.structures).map(asRecord);
      const uniprotAcc = text(entry.ac) ?? acc;
      if (!text(entry.ac) && structures.length === 0) {
        return failure(alphafold, "not_found", `No AlphaFold summary found for ${acc}`, false);
      }
      const firstSummary = asRecord(structures[0]?.summary);
      return makePage("protein", [protein({
        id: curie("uniprot", uniprotAcc),
        name: text(entry.id) ?? uniprotAcc,
        xrefs: structures.slice(0, 10).map((s) => curie("alphafold", text(asRecord(s.summary).model_identifier) ?? `AF-${uniprotAcc}`)),
        source_extra: {
          sequence_length: numberOrNull(entry.sequence_length),
          model_count: structures.length,
          first_model_id: text(firstSummary.model_identifier),
          confidence_type: text(firstSummary.confidence_type),
          confidence_avg_local_score: numberOrNull(firstSummary.confidence_avg_local_score),
          coverage: numberOrNull(firstSummary.coverage),
        },
      })], 0, 1, 1, alphafold, alphafoldUsage);
    }),
  },
];
