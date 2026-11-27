import { curie, ensemblGeneCurie, humanChrom, strandOf } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
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
  sourceOf,
  text,
  variant,
} from "phi-omics-runtime";
import { z } from "zod";

const clinvar = sourceOf("clinvar", "NCBI ClinVar", "https://www.ncbi.nlm.nih.gov/clinvar/", "GRCh38");
const clinvarUsage = {
  summary: "ClinVar variant clinical significance, review status, condition traits, and genomic coordinates from NCBI E-utilities.",
  url: "https://www.ncbi.nlm.nih.gov/home/about/policies/",
};
const gnomad = sourceOf("gnomad", "gnomAD", "https://gnomad.broadinstitute.org/", "v4");
const gnomadUsage = {
  summary: "gnomAD population allele frequencies and gene constraint metrics from the public GraphQL API.",
  url: "https://gnomad.broadinstitute.org/terms",
};
const civic = sourceOf("civic", "CIViC", "https://civicdb.org/", "v2");
const civicUsage = {
  summary: "Clinical Interpretation of Variants in Cancer (CIViC) curated variants and evidence items (CC0 public domain).",
  url: "https://docs.civicdb.org/en/latest/about/citing.html",
};

const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

function ncbiUrl(script: string, params: Record<string, string>): string {
  const url = new URL(`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/${script}`);
  url.searchParams.set("tool", "phi-omics-variants");
  const email = process.env.PHI_LITERATURE_EMAIL || process.env.NCBI_EMAIL || "";
  if (email.trim()) url.searchParams.set("email", email.trim());
  if (process.env.NCBI_API_KEY) url.searchParams.set("api_key", process.env.NCBI_API_KEY);
  for (const [k, v] of Object.entries(params)) url.searchParams.set(k, v);
  return url.toString();
}

function uniqueStrings(values: Array<string | null | undefined>): string[] {
  const seen = new Set<string>();
  const out: string[] = [];
  for (const v of values) {
    if (!v || seen.has(v)) continue;
    seen.add(v);
    out.push(v);
  }
  return out;
}

function mapClinvarSummary(uid: string, row: Record<string, unknown>) {
  const variationSet = asRecord(asArray(row.variation_set)[0]);
  const locs = asArray(variationSet.variation_loc).map(asRecord);
  const grch38 = locs.find((loc) => text(loc.assembly_name) === "GRCh38") ?? locs[0];
  const start = numberOrNull(grch38?.start);
  const stop = numberOrNull(grch38?.stop) ?? start;
  const genes = asArray(row.genes).map(asRecord);
  const firstGene = genes[0];
  const geneIdRaw = text(firstGene?.geneid);
  const geneSymbol = text(firstGene?.symbol);
  const accession = text(row.accession);
  const rsids = asArray(variationSet.variation_xrefs)
    .map(asRecord)
    .filter((x) => text(x.db_source) === "dbSNP" && text(x.db_id))
    .map((x) => `dbsnp:rs${text(x.db_id)!.replace(/^rs/i, "")}`);
  const germ = asRecord(row.germline_classification);
  const onco = asRecord(row.oncogenicity_classification);
  const somatic = asRecord(row.clinical_impact_classification);
  const traitSet = asArray(germ.trait_set).map(asRecord);
  const traits = traitSet.map((t) => text(t.trait_name)).filter((s): s is string => s != null);
  const sig = text(germ.description) ?? text(onco.description) ?? text(somatic.description);
  const xrefs = uniqueStrings([
    accession ? curie("clinvar.accession", accession) : null,
    ...rsids,
  ]);
  return variant({
    id: curie("clinvar", uid),
    gene_id: geneIdRaw ? curie("ncbigene", geneIdRaw) : geneSymbol ? curie("hgnc.symbol", geneSymbol) : null,
    locus: grch38 ? locusOf({
      genomeBuild: text(grch38.assembly_name),
      chrom: humanChrom(text(grch38.chr)),
      start,
      end: stop != null && start != null && stop >= start ? stop : start,
    }) : null,
    consequence: sig,
    xrefs,
    source_extra: {
      accession,
      title: text(row.title),
      obj_type: text(row.obj_type),
      gene_symbol: geneSymbol,
      genes: genes.map((g) => text(g.symbol)).filter(Boolean),
      germline_classification: {
        description: text(germ.description),
        review_status: text(germ.review_status),
        last_evaluated: text(germ.last_evaluated),
        traits,
      },
      oncogenicity_classification: text(onco.description) ? {
        description: text(onco.description),
        review_status: text(onco.review_status),
      } : null,
      clinical_impact_classification: text(somatic.description) ? {
        description: text(somatic.description),
        review_status: text(somatic.review_status),
      } : null,
      hgvs: asArray(row.protein_change).map(text).filter(Boolean),
      canonical_spdi: text(variationSet.canonical_spdi),
    },
  });
}

async function fetchClinvarSummaries(client: JsonClient, ids: string[]) {
  if (ids.length === 0) return [];
  const body = asRecord(await client.getJson(ncbiUrl("esummary.fcgi", {
    db: "clinvar",
    id: ids.join(","),
    retmode: "json",
  })));
  const result = asRecord(body.result);
  const uids = asArray(result.uids).map(text).filter((u): u is string => u != null);
  const ordered = uids.length > 0 ? uids : ids;
  return ordered.flatMap((uid) => {
    const row = asRecord(result[uid]);
    return text(row.uid) || text(row.accession) || text(row.title) ? [mapClinvarSummary(uid, row)] : [];
  });
}

async function gnomadQuery(client: JsonClient, query: string, variables: Record<string, unknown>) {
  const res = asRecord(await client.postJson("https://gnomad.broadinstitute.org/api", { query, variables }));
  const errors = asArray(res.errors).map(asRecord);
  if (errors.length > 0 && !res.data) {
    const msg = text(errors[0]?.message) ?? "gnomAD GraphQL query failed";
    if (/not found/i.test(msg)) return { kind: "not_found" as const, message: msg };
    return { kind: "error" as const, message: msg };
  }
  const data = asRecord(res.data);
  if (errors.length > 0) {
    const msg = text(errors[0]?.message) ?? "gnomAD GraphQL query failed";
    const hasValue = Object.values(data).some((v) => v != null);
    if (!hasValue) {
      return /not found/i.test(msg)
        ? { kind: "not_found" as const, message: msg }
        : { kind: "error" as const, message: msg };
    }
  }
  return { kind: "ok" as const, data };
}

async function civicQuery(client: JsonClient, query: string, variables: Record<string, unknown>) {
  const res = asRecord(await client.postJson("https://civicdb.org/api/graphql", { query, variables }));
  const errors = asArray(res.errors).map(asRecord);
  if (errors.length > 0 && !res.data) {
    return { kind: "error" as const, message: text(errors[0]?.message) ?? "CIViC GraphQL query failed" };
  }
  return { kind: "ok" as const, data: asRecord(res.data) };
}

function normalizeGnomadVariantId(raw: string): string {
  return raw.trim().replace(/^chr/i, "").replace(/:/g, "-");
}

function safeEnsemblCurie(raw: string | null): string | null {
  if (!raw) return null;
  try {
    return ensemblGeneCurie(raw).id;
  } catch {
    return curie("ensembl", raw.replace(/\s+/g, "_"));
  }
}

export const variantTools: OmicsTool[] = [
  {
    name: "clinvar_search_variants",
    title: "Search ClinVar variants",
    description: "Search NCBI ClinVar by gene symbol, free-text term, and optional clinical significance (e.g. pathogenic, benign).",
    inputSchema: z.object({
      gene: z.string().optional().describe("Gene symbol, e.g. BRCA1 or TP53."),
      term: z.string().optional().describe("ClinVar search expression, e.g. NM_000546 or c.818G>A or cardiomyopathy."),
      clinical_significance: z.string().optional().describe("Optional clinical significance filter, e.g. pathogenic, likely_pathogenic."),
      ...pageArgs,
    }),
    call: (client, args) => guard(clinvar, async () => {
      const { offset, limit } = paging(args);
      const geneSym = text(args.gene);
      const freeTerm = text(args.term);
      if (!geneSym && !freeTerm) {
        return failure(clinvar, "invalid_input", "Pass gene or term to search ClinVar", false);
      }
      const clauses: string[] = [];
      if (geneSym) clauses.push(`${geneSym}[gene]`);
      if (freeTerm) clauses.push(`(${freeTerm})`);
      const sig = text(args.clinical_significance);
      if (sig) clauses.push(`"${sig.replace(/_/g, " ")}"[clinical_significance]`);
      const fullTerm = clauses.join(" AND ");
      const searchBody = asRecord(await client.getJson(ncbiUrl("esearch.fcgi", {
        db: "clinvar",
        term: fullTerm,
        retstart: String(offset),
        retmax: String(limit),
        retmode: "json",
      })));
      const esearch = asRecord(searchBody.esearchresult);
      const ids = asArray(esearch.idlist).map(text).filter((id): id is string => id != null);
      const items = await fetchClinvarSummaries(client, ids);
      return makePage("variant", items, offset, limit, numberOrNull(esearch.count), clinvar, clinvarUsage, {
        query: fullTerm,
        query_translation: text(esearch.querytranslation),
      });
    }),
  },
  {
    name: "clinvar_get_variant",
    title: "Get one ClinVar variant",
    description: "Fetch detailed ClinVar summary by numeric variation ID (e.g. 12351) or VCV accession (e.g. VCV000012351).",
    inputSchema: z.object({
      variation_id: z.string().describe("ClinVar numeric ID or VCV accession."),
    }),
    call: (client, args) => guard(clinvar, async () => {
      const raw = text(args.variation_id);
      if (!raw) return failure(clinvar, "invalid_input", "variation_id is required", false);
      let uid = raw.match(/^(?:clinvar:)?(\d+)$/i)?.[1] ?? null;
      if (!uid && /^VCV\d+/i.test(raw)) {
        const searchBody = asRecord(await client.getJson(ncbiUrl("esearch.fcgi", {
          db: "clinvar",
          term: `${raw.split(".")[0]}[accession]`,
          retmax: "1",
          retmode: "json",
        })));
        uid = text(asArray(asRecord(searchBody.esearchresult).idlist)[0]);
      }
      if (!uid) return failure(clinvar, "invalid_input", "variation_id must be numeric or a VCV accession", false);
      const items = await fetchClinvarSummaries(client, [uid]);
      if (items.length === 0) return failure(clinvar, "not_found", `ClinVar variant ${raw} was not found`, false);
      return makePage("variant", items, 0, 1, 1, clinvar, clinvarUsage);
    }),
  },
  {
    name: "gnomad_get_variant",
    title: "Get gnomAD variant allele frequency",
    description: "Fetch population allele frequencies (exome, genome, subpopulations) and transcript consequences from gnomAD v4/v3/v2.1.",
    inputSchema: z.object({
      variant_id: z.string().describe("gnomAD variant ID like 17-7674220-C-T (chr prefix optional)."),
      dataset: z.enum(["gnomad_r4", "gnomad_r3", "gnomad_r2_1"]).optional().describe("Default gnomad_r4."),
      reference_genome: z.enum(["GRCh38", "GRCh37"]).optional().describe("Default GRCh38 for r4/r3 or GRCh37 for r2_1."),
    }),
    call: (client, args) => guard(gnomad, async () => {
      const rawId = text(args.variant_id);
      if (!rawId) return failure(gnomad, "invalid_input", "variant_id is required", false);
      const normalized = normalizeGnomadVariantId(rawId);
      if (!/^[0-9XYMT]+-\d+-[ACGTN]+-[ACGTN]+$/i.test(normalized)) {
        return failure(gnomad, "invalid_input", "variant_id must look like 17-7674220-C-T", false);
      }
      const dataset = (text(args.dataset) ?? "gnomad_r4") as "gnomad_r4" | "gnomad_r3" | "gnomad_r2_1";
      const refGenome = text(args.reference_genome) ?? (dataset === "gnomad_r2_1" ? "GRCh37" : "GRCh38");
      const query = `query ($variantId: String!, $dataset: DatasetId!) {
        variant(variantId: $variantId, dataset: $dataset) {
          variant_id
          chrom
          pos
          ref
          alt
          rsids
          flags
          exome { ac an af homozygote_count filters populations { id ac an } }
          genome { ac an af homozygote_count filters populations { id ac an } }
          transcript_consequences { gene_id gene_symbol transcript_id major_consequence hgvsc hgvsp lof }
        }
      }`;
      const res = await gnomadQuery(client, query, { variantId: normalized, dataset });
      if (res.kind === "not_found") return failure(gnomad, "not_found", res.message, false);
      if (res.kind === "error") return failure(gnomad, "upstream_http", res.message, true);
      const v = asRecord(res.data.variant);
      const vid = text(v.variant_id);
      if (!vid) return failure(gnomad, "not_found", `gnomAD variant ${normalized} was not found in ${dataset}`, false);
      const consequences = asArray(v.transcript_consequences).map(asRecord);
      const primary = consequences[0];
      const pos = numberOrNull(v.pos);
      const refLen = text(v.ref)?.length ?? 1;
      const rsids = asArray(v.rsids).map(text).filter((r): r is string => r != null);
      const row = variant({
        id: curie("gnomad", `${dataset}:${vid}`),
        gene_id: safeEnsemblCurie(text(primary?.gene_id)),
        locus: pos != null ? locusOf({
          genomeBuild: refGenome,
          chrom: humanChrom(text(v.chrom)),
          start: pos,
          end: pos + Math.max(0, refLen - 1),
        }) : null,
        consequence: text(primary?.major_consequence),
        xrefs: rsids.map((r) => curie("dbsnp", r)),
        source_extra: {
          dataset,
          ref: text(v.ref),
          alt: text(v.alt),
          flags: asArray(v.flags),
          exome: v.exome ?? null,
          genome: v.genome ?? null,
          transcript_consequences: consequences.slice(0, 10),
        },
      });
      return makePage("variant", [row], 0, 1, 1, gnomad, gnomadUsage);
    }),
  },
  {
    name: "gnomad_get_gene_constraint",
    title: "Get gnomAD gene constraint",
    description: "Fetch gnomAD gene constraint scores (pLI, LOEUF / oe_lof_upper, missense Z, synonymous Z) for one human gene.",
    inputSchema: z.object({
      gene: z.string().describe("HGNC symbol (e.g. TP53) or Ensembl gene ID (e.g. ENSG00000141510)."),
      reference_genome: z.enum(["GRCh38", "GRCh37"]).optional().describe("Default GRCh38."),
    }),
    call: (client, args) => guard(gnomad, async () => {
      const target = text(args.gene);
      if (!target) return failure(gnomad, "invalid_input", "gene is required", false);
      const refGenome = text(args.reference_genome) ?? "GRCh38";
      const isEnsembl = /^(?:ensembl:)?ENSG\d+/i.test(target);
      const cleanedId = target.replace(/^ensembl:/i, "");
      const query = `query ($geneId: String, $geneSymbol: String, $ref: ReferenceGenomeId!) {
        gene(gene_id: $geneId, gene_symbol: $geneSymbol, reference_genome: $ref) {
          gene_id
          symbol
          name
          chrom
          start
          stop
          strand
          canonical_transcript_id
          omim_id
          gnomad_constraint {
            pLI
            oe_lof
            oe_lof_lower
            oe_lof_upper
            oe_mis
            oe_mis_lower
            oe_mis_upper
            mis_z
            oe_syn
            syn_z
            lof_z
            flags
          }
        }
      }`;
      const res = await gnomadQuery(client, query, {
        geneId: isEnsembl ? cleanedId : null,
        geneSymbol: isEnsembl ? null : target,
        ref: refGenome,
      });
      if (res.kind === "not_found") return failure(gnomad, "not_found", res.message, false);
      if (res.kind === "error") return failure(gnomad, "upstream_http", res.message, true);
      const g = asRecord(res.data.gene);
      const gid = text(g.gene_id);
      if (!gid) return failure(gnomad, "not_found", `Gene ${target} was not found in gnomAD`, false);
      const start = numberOrNull(g.start);
      const stop = numberOrNull(g.stop);
      const row = gene({
        id: safeEnsemblCurie(gid) ?? curie("gnomad.gene", gid),
        symbol: text(g.symbol),
        name: text(g.name),
        taxon_id: "NCBITaxon:9606",
        locus: start != null ? locusOf({
          genomeBuild: refGenome,
          chrom: humanChrom(text(g.chrom)),
          start,
          end: stop != null && stop >= start ? stop : start,
          strand: strandOf(g.strand),
        }) : null,
        xrefs: text(g.omim_id) ? [curie("omim", text(g.omim_id)!)] : [],
        source_extra: {
          canonical_transcript_id: text(g.canonical_transcript_id),
          gnomad_constraint: g.gnomad_constraint ?? null,
        },
      });
      return makePage("gene", [row], 0, 1, 1, gnomad, gnomadUsage);
    }),
  },
  {
    name: "civic_search_variants",
    title: "Search CIViC cancer variants",
    description: "Search curated oncology variants in CIViC v2 by gene symbol and/or variant name (e.g. gene=BRAF, variant=V600E).",
    inputSchema: z.object({
      gene: z.string().optional().describe("Feature/gene symbol, e.g. BRAF, EGFR."),
      variant: z.string().optional().describe("Variant name filter, e.g. V600E, L858R."),
      ...pageArgs,
    }),
    call: (client, args) => guard(civic, async () => {
      const { offset, limit } = paging(args);
      const geneSym = text(args.gene);
      const variantName = text(args.variant);
      if (!geneSym && !variantName) {
        return failure(civic, "invalid_input", "Pass gene or variant to search CIViC", false);
      }
      if (offset + limit > 100) {
        return failure(civic, "invalid_input", "offset + limit must be at most 100 for CIViC search", false);
      }
      const query = `query ($first: Int, $featureName: String, $variantName: String) {
        browseVariants(first: $first, featureName: $featureName, variantName: $variantName) {
          totalCount
          nodes {
            id
            name
            featureName
            link
            featureId
            category
            aliases { name }
            variantTypes { id name }
            evidenceItemCount
          }
        }
      }`;
      const res = await civicQuery(client, query, {
        first: offset + limit,
        ...(geneSym ? { featureName: geneSym } : {}),
        ...(variantName ? { variantName } : {}),
      });
      if (res.kind === "error") return failure(civic, "upstream_http", res.message, true);
      const payload = asRecord(res.data.browseVariants);
      const nodes = asArray(payload.nodes).map(asRecord).slice(offset, offset + limit);
      const items = safeMap(nodes, (row) => {
        const vid = text(row.id) ?? "unknown";
        const sym = text(row.featureName);
        const types = asArray(row.variantTypes).map((vt) => text(asRecord(vt).name)).filter((s): s is string => s != null);
        const aliases = asArray(row.aliases).map((a) => text(asRecord(a).name)).filter((s): s is string => s != null);
        const rsids = aliases.filter((a) => /^RS\d+$/i.test(a)).map((r) => curie("dbsnp", r.toLowerCase()));
        return variant({
          id: curie("civic.variant", vid),
          gene_id: sym ? curie("hgnc.symbol", sym) : null,
          consequence: types[0] ?? null,
          xrefs: uniqueStrings(rsids),
          source_extra: {
            name: text(row.name),
            gene_symbol: sym,
            feature_id: numberOrNull(row.featureId),
            category: text(row.category),
            link: text(row.link) ? `https://civicdb.org${text(row.link)}` : null,
            variant_types: types,
            aliases,
            evidence_item_count: numberOrNull(row.evidenceItemCount),
          },
        });
      });
      return makePage("variant", items, offset, limit, numberOrNull(payload.totalCount), civic, civicUsage);
    }),
  },
  {
    name: "civic_search_evidence",
    title: "Search CIViC clinical evidence",
    description: "Search curated CIViC v2 clinical evidence items linking cancer variants/molecular profiles to therapies, prognosis, or diagnosis.",
    inputSchema: z.object({
      molecular_profile: z.string().optional().describe("Molecular profile name fragment, e.g. BRAF V600E."),
      disease: z.string().optional().describe("Disease name fragment, e.g. Melanoma."),
      evidence_type: z.enum(["PREDICTIVE", "PROGNOSTIC", "DIAGNOSTIC", "PREDISPOSING", "ONCOGENIC", "FUNCTIONAL"]).optional(),
      evidence_level: z.enum(["A", "B", "C", "D", "E"]).optional(),
      ...pageArgs,
    }),
    call: (client, args) => guard(civic, async () => {
      const { offset, limit } = paging(args);
      const mp = text(args.molecular_profile);
      const disease = text(args.disease);
      const evidenceType = text(args.evidence_type);
      const evidenceLevel = text(args.evidence_level);
      if (!mp && !disease && !evidenceType && !evidenceLevel) {
        return failure(civic, "invalid_input", "Pass molecular_profile, disease, evidence_type, or evidence_level", false);
      }
      if (offset + limit > 100) {
        return failure(civic, "invalid_input", "offset + limit must be at most 100 for CIViC search", false);
      }
      const query = `query ($first: Int, $molecularProfileName: String, $diseaseName: String, $evidenceType: EvidenceType, $evidenceLevel: EvidenceLevel) {
        evidenceItems(
          first: $first
          molecularProfileName: $molecularProfileName
          diseaseName: $diseaseName
          evidenceType: $evidenceType
          evidenceLevel: $evidenceLevel
        ) {
          totalCount
          nodes {
            id
            name
            status
            evidenceType
            evidenceLevel
            evidenceDirection
            significance
            evidenceRating
            description
            molecularProfile { id name }
            disease { name doid }
            therapies { name }
            source { citationId sourceType }
          }
        }
      }`;
      const res = await civicQuery(client, query, {
        first: offset + limit,
        ...(mp ? { molecularProfileName: mp } : {}),
        ...(disease ? { diseaseName: disease } : {}),
        ...(evidenceType ? { evidenceType } : {}),
        ...(evidenceLevel ? { evidenceLevel } : {}),
      });
      if (res.kind === "error") return failure(civic, "upstream_http", res.message, true);
      const payload = asRecord(res.data.evidenceItems);
      const nodes = asArray(payload.nodes).map(asRecord).slice(offset, offset + limit);
      const items = safeMap(nodes, (row) => {
        const eid = text(row.id) ?? "unknown";
        const mpObj = asRecord(row.molecularProfile);
        const disObj = asRecord(row.disease);
        const srcObj = asRecord(row.source);
        const doid = text(disObj.doid);
        const citId = text(srcObj.citationId);
        const srcType = text(srcObj.sourceType)?.toUpperCase();
        return variant({
          id: curie("civic.eid", eid),
          consequence: [text(row.evidenceType), text(row.evidenceLevel), text(row.significance)].filter(Boolean).join(":") || null,
          xrefs: uniqueStrings([
            doid ? curie("DOID", doid.replace(/^DOID:/i, "")) : null,
            srcType === "PUBMED" && citId ? curie("pmid", citId) : null,
          ]),
          source_extra: {
            name: text(row.name),
            status: text(row.status),
            evidence_type: text(row.evidenceType),
            evidence_level: text(row.evidenceLevel),
            evidence_direction: text(row.evidenceDirection),
            significance: text(row.significance),
            rating: numberOrNull(row.evidenceRating),
            molecular_profile: text(mpObj.name),
            disease: text(disObj.name),
            therapies: asArray(row.therapies).map((t) => text(asRecord(t).name)).filter(Boolean),
            description: text(row.description),
          },
        });
      });
      return makePage("variant", items, offset, limit, numberOrNull(payload.totalCount), civic, civicUsage);
    }),
  },
];
