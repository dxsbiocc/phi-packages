import { curie } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  compound,
  dataset,
  drugResponse,
  failure,
  guard,
  makePage,
  numberOrNull,
  type OmicsTool,
  paging,
  safeMap,
  sourceOf,
  text,
} from "phi-omics-runtime";
import { z } from "zod";

const chembl = sourceOf("chembl", "EMBL-EBI ChEMBL", "https://www.ebi.ac.uk/chembl/", "REST");
const chemblUsage = {
  summary: "ChEMBL bioactive molecules, physicochemical descriptors, IC50/Ki/Kd bioactivities, and curated mechanisms of action (CC BY-SA 3.0).",
  url: "https://chembl.gitbook.io/chembl-interface-documentation/web-services",
};
const openfda = sourceOf("openfda", "openFDA", "https://open.fda.gov/", "drug");
const openfdaUsage = {
  summary: "openFDA prescription drug labels and FAERS spontaneous adverse event reaction counts. FAERS reports reflect submission counts, not validated incidence rates.",
  url: "https://open.fda.gov/terms/",
};

const pageArgs = {
  offset: z.number().int().nonnegative().optional(),
  limit: z.number().int().min(1).max(50).optional(),
};

function clip(value: string | null, max = 800): string | null {
  if (!value) return null;
  const clean = value.replace(/\s+/g, " ").trim();
  if (!clean) return null;
  return clean.length <= max ? clean : `${clean.slice(0, max - 1).trimEnd()}…`;
}

function firstText(value: unknown): string | null {
  if (Array.isArray(value)) return text(value[0]);
  return text(value);
}

function mapChemblMolecule(row: Record<string, unknown>) {
  const chemblId = text(row.molecule_chembl_id) ?? "CHEMBL0";
  const structs = asRecord(row.molecule_structures);
  const props = asRecord(row.molecule_properties);
  const synonyms = asArray(row.molecule_synonyms)
    .map((s) => text(asRecord(s).molecule_synonym) ?? text(s))
    .filter((s): s is string => s != null);
  return compound({
    id: curie("chembl", chemblId),
    name: text(row.pref_name) ?? chemblId,
    chembl_id: chemblId,
    inchikey: text(structs.standard_inchi_key),
    synonyms: [...new Set(synonyms)].slice(0, 12),
    source_extra: {
      molecule_type: text(row.molecule_type),
      max_phase: numberOrNull(row.max_phase),
      first_approval: numberOrNull(row.first_approval),
      oral: row.oral ?? null,
      canonical_smiles: text(structs.canonical_smiles),
      molecular_formula: text(props.full_molformula),
      molecular_weight: numberOrNull(props.full_mwt),
      alogp: numberOrNull(props.alogp),
      psa: numberOrNull(props.psa),
      hba: numberOrNull(props.hba),
      hbd: numberOrNull(props.hbd),
      ro5_violations: numberOrNull(props.num_ro5_violations),
      qed_weighted: numberOrNull(props.qed_weighted),
    },
  });
}

function fdaTerm(raw: string): string {
  return `"${raw.replace(/"/g, " ").trim()}"`;
}

export const chemistryTools: OmicsTool[] = [
  {
    name: "chembl_search_molecules",
    title: "Search ChEMBL molecules",
    description: "Search ChEMBL compounds by generic name, synonym, or ChEMBL keyword.",
    inputSchema: z.object({
      query: z.string().describe("Drug or compound name, e.g. imatinib, vemurafenib."),
      ...pageArgs,
    }),
    call: (client, args) => guard(chembl, async () => {
      const { offset, limit } = paging(args);
      const q = text(args.query);
      if (!q) return failure(chembl, "invalid_input", "query is required", false);
      const url = new URL("https://www.ebi.ac.uk/chembl/api/data/molecule/search.json");
      url.searchParams.set("q", q);
      url.searchParams.set("limit", String(limit));
      url.searchParams.set("offset", String(offset));
      const body = asRecord(await client.getJson(url.toString()));
      const rows = asArray(body.molecules).map(asRecord);
      const meta = asRecord(body.page_meta);
      return makePage("compound", safeMap(rows, mapChemblMolecule), offset, limit, numberOrNull(meta.total_count), chembl, chemblUsage);
    }),
  },
  {
    name: "chembl_get_molecule",
    title: "Get one ChEMBL molecule",
    description: "Fetch detailed physicochemical properties, SMILES, InChIKey, clinical phase, and synonyms for one ChEMBL compound (e.g. CHEMBL941).",
    inputSchema: z.object({
      chembl_id: z.string().describe("ChEMBL molecule ID, e.g. CHEMBL941."),
    }),
    call: (client, args) => guard(chembl, async () => {
      const raw = text(args.chembl_id)?.replace(/^chembl:/i, "").toUpperCase();
      if (!raw || !/^CHEMBL\d+$/.test(raw)) {
        return failure(chembl, "invalid_input", "chembl_id must look like CHEMBL941", false);
      }
      const row = asRecord(await client.getJson(`https://www.ebi.ac.uk/chembl/api/data/molecule/${raw}.json`));
      if (!text(row.molecule_chembl_id)) {
        return failure(chembl, "not_found", `Molecule ${raw} was not found`, false);
      }
      return makePage("compound", [mapChemblMolecule(row)], 0, 1, 1, chembl, chemblUsage);
    }),
  },
  {
    name: "chembl_get_activities",
    title: "Get ChEMBL bioactivities",
    description: "Fetch experimental IC50/Ki/Kd/EC50 bioactivity records by molecule_chembl_id and/or target_chembl_id.",
    inputSchema: z.object({
      molecule_chembl_id: z.string().optional().describe("ChEMBL compound ID, e.g. CHEMBL941."),
      target_chembl_id: z.string().optional().describe("ChEMBL target ID, e.g. CHEMBL1862."),
      standard_type: z.string().optional().describe("Optional measurement type filter, e.g. IC50, Ki, Kd."),
      ...pageArgs,
    }),
    call: (client, args) => guard(chembl, async () => {
      const { offset, limit } = paging(args);
      const molId = text(args.molecule_chembl_id)?.replace(/^chembl:/i, "").toUpperCase();
      const targetId = text(args.target_chembl_id)?.replace(/^chembl:/i, "").toUpperCase();
      if (!molId && !targetId) {
        return failure(chembl, "invalid_input", "Pass molecule_chembl_id or target_chembl_id", false);
      }
      const url = new URL("https://www.ebi.ac.uk/chembl/api/data/activity.json");
      if (molId) url.searchParams.set("molecule_chembl_id", molId);
      if (targetId) url.searchParams.set("target_chembl_id", targetId);
      const stdType = text(args.standard_type);
      if (stdType) url.searchParams.set("standard_type", stdType);
      url.searchParams.set("limit", String(limit));
      url.searchParams.set("offset", String(offset));
      const body = asRecord(await client.getJson(url.toString()));
      const rows = asArray(body.activities).map(asRecord);
      const meta = asRecord(body.page_meta);
      const items = safeMap(rows, (row) => {
        const aid = text(row.activity_id) ?? "unknown";
        const mId = text(row.molecule_chembl_id);
        const assayId = text(row.assay_chembl_id);
        return drugResponse({
          id: curie("chembl.activity", aid),
          compound_id: mId ? curie("chembl", mId) : null,
          dataset_id: assayId ? curie("chembl.assay", assayId) : null,
          metric: text(row.standard_type),
          value: numberOrNull(row.standard_value),
          unit: text(row.standard_units),
          source_extra: {
            pchembl_value: numberOrNull(row.pchembl_value),
            standard_relation: text(row.standard_relation),
            target_chembl_id: text(row.target_chembl_id),
            target_pref_name: text(row.target_pref_name),
            target_organism: text(row.target_organism),
            assay_type: text(row.assay_type),
            document_chembl_id: text(row.document_chembl_id),
          },
        });
      });
      return makePage("drug_response", items, offset, limit, numberOrNull(meta.total_count), chembl, chemblUsage);
    }),
  },
  {
    name: "chembl_get_mechanisms",
    title: "Get ChEMBL drug mechanisms of action",
    description: "Fetch curated mechanism-of-action and primary target records for one ChEMBL molecule (e.g. CHEMBL941).",
    inputSchema: z.object({
      molecule_chembl_id: z.string().describe("ChEMBL compound ID, e.g. CHEMBL941."),
      ...pageArgs,
    }),
    call: (client, args) => guard(chembl, async () => {
      const { offset, limit } = paging(args);
      const molId = text(args.molecule_chembl_id)?.replace(/^chembl:/i, "").toUpperCase();
      if (!molId || !/^CHEMBL\d+$/.test(molId)) {
        return failure(chembl, "invalid_input", "molecule_chembl_id must look like CHEMBL941", false);
      }
      const parentUrl = new URL("https://www.ebi.ac.uk/chembl/api/data/mechanism.json");
      parentUrl.searchParams.set("parent_molecule_chembl_id", molId);
      parentUrl.searchParams.set("limit", String(limit));
      parentUrl.searchParams.set("offset", String(offset));
      let body = asRecord(await client.getJson(parentUrl.toString()));
      let rows = asArray(body.mechanisms).map(asRecord);
      if (rows.length === 0) {
        const directUrl = new URL("https://www.ebi.ac.uk/chembl/api/data/mechanism.json");
        directUrl.searchParams.set("molecule_chembl_id", molId);
        directUrl.searchParams.set("limit", String(limit));
        directUrl.searchParams.set("offset", String(offset));
        body = asRecord(await client.getJson(directUrl.toString()));
        rows = asArray(body.mechanisms).map(asRecord);
      }
      const meta = asRecord(body.page_meta);
      const items = safeMap(rows, (row) => {
        const mecId = text(row.mec_id) ?? "unknown";
        const mId = text(row.molecule_chembl_id) ?? molId;
        const tgtId = text(row.target_chembl_id);
        return drugResponse({
          id: curie("chembl.mechanism", mecId),
          compound_id: curie("chembl", mId),
          metric: text(row.action_type) ?? "mechanism",
          xrefs: tgtId ? [curie("chembl.target", tgtId)] : [],
          source_extra: {
            mechanism_of_action: text(row.mechanism_of_action),
            target_chembl_id: tgtId,
            binding_site_comment: text(row.binding_site_comment),
            direct_interaction: row.direct_interaction ?? null,
            molecular_mechanism: row.molecular_mechanism ?? null,
            max_phase: numberOrNull(row.max_phase),
          },
        });
      });
      return makePage("drug_response", items, offset, limit, numberOrNull(meta.total_count), chembl, chemblUsage);
    }),
  },
  {
    name: "openfda_search_labels",
    title: "Search openFDA drug labels",
    description: "Search FDA SPL drug labels for indications, boxed warnings, contraindications, adverse reactions, and drug interactions.",
    inputSchema: z.object({
      drug: z.string().describe("Generic or brand drug name, e.g. imatinib, vemurafenib."),
      ...pageArgs,
    }),
    call: (client, args) => guard(openfda, async () => {
      const { offset, limit } = paging(args);
      const drug = text(args.drug);
      if (!drug) return failure(openfda, "invalid_input", "drug is required", false);
      const quoted = fdaTerm(drug);
      const searchExpr = `openfda.generic_name:${quoted}+openfda.brand_name:${quoted}`;
      const url = `https://api.fda.gov/drug/label.json?search=${encodeURIComponent(searchExpr).replace(/%2B/g, "+")}&limit=${limit}&skip=${offset}`;
      const body = asRecord(await client.getJson(url));
      const rows = asArray(body.results).map(asRecord);
      const total = numberOrNull(asRecord(asRecord(body.meta).results).total);
      const items = safeMap(rows, (row, idx) => {
        const metaFda = asRecord(row.openfda);
        const setId = text(row.set_id);
        const rawId = text(row.id) ?? setId ?? `${offset + idx}`;
        const genericNames = asArray(metaFda.generic_name).map(text).filter((s): s is string => s != null);
        const brandNames = asArray(metaFda.brand_name).map(text).filter((s): s is string => s != null);
        return dataset({
          id: curie("openfda.label", rawId.replace(/\s+/g, "_")),
          title: [genericNames[0], brandNames[0]].filter(Boolean).join(" / ") || drug,
          description: clip(firstText(row.indications_and_usage), 900),
          access: "open",
          url: setId ? `https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=${encodeURIComponent(setId)}` : null,
          source_extra: {
            set_id: setId,
            effective_time: text(row.effective_time),
            generic_names: genericNames,
            brand_names: brandNames,
            routes: asArray(metaFda.route).map(text).filter(Boolean),
            boxed_warning: clip(firstText(row.boxed_warning), 700),
            warnings: clip(firstText(row.warnings_and_cautions) ?? firstText(row.warnings), 700),
            contraindications: clip(firstText(row.contraindications), 600),
            adverse_reactions: clip(firstText(row.adverse_reactions), 600),
            drug_interactions: clip(firstText(row.drug_interactions), 600),
          },
        });
      });
      return makePage("dataset", items, offset, limit, total, openfda, openfdaUsage);
    }),
  },
  {
    name: "openfda_adverse_events",
    title: "Count openFDA FAERS adverse events",
    description: "Return top MedDRA adverse event terms from FDA FAERS spontaneous reports for a medicinal product.",
    inputSchema: z.object({
      drug: z.string().describe("Drug name, e.g. imatinib, pembrolizumab."),
      serious_only: z.boolean().optional().describe("Default false. Filter to serious FAERS reports when true."),
      limit: z.number().int().min(1).max(50).optional().describe("Default 20, maximum 50 top MedDRA terms."),
    }),
    call: (client, args) => guard(openfda, async () => {
      const { limit } = paging(args);
      const drug = text(args.drug);
      if (!drug) return failure(openfda, "invalid_input", "drug is required", false);
      const quoted = fdaTerm(drug);
      const searchExpr = args.serious_only === true
        ? `patient.drug.medicinalproduct:${quoted}+AND+serious:1`
        : `patient.drug.medicinalproduct:${quoted}`;
      const url = `https://api.fda.gov/drug/event.json?search=${encodeURIComponent(searchExpr).replace(/%2B/g, "+")}&count=patient.reaction.reactionmeddrapt.exact&limit=${limit}`;
      const body = asRecord(await client.getJson(url));
      const rows = asArray(body.results).map(asRecord);
      const drugSlug = encodeURIComponent(drug.toLowerCase().replace(/\s+/g, "-"));
      const items = safeMap(rows, (row) => {
        const term = text(row.term) ?? "UNKNOWN";
        const count = numberOrNull(row.count);
        return dataset({
          id: curie("openfda.ae", `${drugSlug}:${encodeURIComponent(term)}`),
          title: term,
          description: `FAERS spontaneous reaction reports (${count ?? 0}) for ${drug}.`,
          sample_count: count,
          access: "open",
          url: "https://open.fda.gov/apis/drug/event/",
          source_extra: {
            drug,
            serious_only: args.serious_only === true,
            meddra_pt: term,
            report_count: count,
            caveat: "FAERS counts reflect spontaneous reporting volume and cannot establish causality or incidence.",
          },
        });
      });
      return makePage("dataset", items, 0, limit, items.length, openfda, openfdaUsage);
    }),
  },
];
