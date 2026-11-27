import type { Dataset, Source, ToolResult } from "phi-omics-schemas";
import {
  asArray,
  asRecord,
  dataset,
  failure,
  makePage,
  numberOrNull,
  safeMap,
  text,
  type JsonClient,
  type OmicsTool,
  type ServerModule,
} from "phi-omics-runtime";
import { z } from "zod";

export const ctgovSource: Source = {
  id: "clinicaltrials-gov",
  name: "ClinicalTrials.gov (NLM REST API v2)",
  version: "v2",
  homepage: "https://clinicaltrials.gov/",
  retrieved_at: "1970-01-01T00:00:00.000Z",
};

const USAGE_SUMMARY = [
  "Self-hosted ClinicalTrials.gov MCP server operated by Phi Packages over the US National Library of Medicine (NLM) ClinicalTrials.gov REST API v2.",
  "Returns registered clinical study protocols, eligibility criteria, arms/interventions, primary/secondary endpoints, and posted summary results for research and evidence synthesis.",
  "Clinical trial registry entries are not medical advice or treatment recommendations. Cite ClinicalTrials.gov and the NCT identifier when referencing studies.",
].join(" ");

const API_BASE = "https://clinicaltrials.gov/api/v2";

const STATUS_VALUES = [
  "RECRUITING",
  "NOT_YET_RECRUITING",
  "ACTIVE_NOT_RECRUITING",
  "COMPLETED",
  "ENROLLING_BY_INVITATION",
  "SUSPENDED",
  "TERMINATED",
  "WITHDRAWN",
  "UNKNOWN",
] as const;

const PHASE_VALUES = ["EARLY_PHASE1", "PHASE1", "PHASE2", "PHASE3", "PHASE4", "NA"] as const;

function normalizeNctId(raw: string): string | null {
  const cleaned = raw.trim().replace(/^nct:/i, "").toUpperCase();
  return /^NCT\d{6,10}$/.test(cleaned) ? cleaned : null;
}

function clip(value: string | null | undefined, maxLen = 1200): string | null {
  if (!value) return null;
  const trimmed = value.trim();
  if (!trimmed) return null;
  return trimmed.length > maxLen ? `${trimmed.slice(0, maxLen - 1)}…` : trimmed;
}

function studyToDataset(node: unknown, detailMode: "summary" | "full" | "eligibility" | "results"): Dataset {
  const root = asRecord(node);
  const proto = asRecord(root.protocolSection);
  const ident = asRecord(proto.identificationModule);
  const statusMod = asRecord(proto.statusModule);
  const descMod = asRecord(proto.descriptionModule);
  const condMod = asRecord(proto.conditionsModule);
  const designMod = asRecord(proto.designModule);
  const armsMod = asRecord(proto.armsInterventionsModule);
  const outMod = asRecord(proto.outcomesModule);
  const eligMod = asRecord(proto.eligibilityModule);
  const sponsorMod = asRecord(proto.sponsorCollaboratorsModule);
  const locMod = asRecord(proto.contactsLocationsModule);
  const refsMod = asRecord(proto.referencesModule);
  const resultsSec = asRecord(root.resultsSection);

  const nctId = text(ident.nctId) ?? "NCT00000000";
  const briefTitle = text(ident.briefTitle);
  const officialTitle = text(ident.officialTitle);
  const studyType = text(designMod.studyType) ?? "CLINICAL_TRIAL";
  const phases = asArray(designMod.phases).map(text).filter((p): p is string => Boolean(p));
  const enrollment = numberOrNull(asRecord(designMod.enrollmentInfo).count);
  const conditions = asArray(condMod.conditions).map(text).filter((c): c is string => Boolean(c));
  const interventions = asArray(armsMod.interventions).map((item) => {
    const r = asRecord(item);
    return {
      type: text(r.type),
      name: text(r.name),
      description: detailMode === "summary" ? clip(text(r.description), 240) : text(r.description),
    };
  });

  const pmids = asArray(refsMod.references)
    .map((r) => text(asRecord(r).pmid))
    .filter((p): p is string => p != null && /^\d+$/.test(p));
  const xrefs = Array.from(new Set(pmids.map((p) => `pmid:${p}`)));

  const baseExtra: Record<string, unknown> = {
    nct_id: nctId,
    official_title: officialTitle,
    overall_status: text(statusMod.overallStatus),
    has_results: Boolean(root.hasResults),
    phases,
    conditions,
    interventions: detailMode === "summary" ? interventions.slice(0, 8) : interventions,
    lead_sponsor: text(asRecord(sponsorMod.leadSponsor).name),
    start_date: text(asRecord(statusMod.startDateStruct).date),
    primary_completion_date: text(asRecord(statusMod.primaryCompletionDateStruct).date),
    completion_date: text(asRecord(statusMod.completionDateStruct).date),
    last_update_posted: text(asRecord(statusMod.lastUpdatePostDateStruct).date),
  };

  if (detailMode === "full") {
    baseExtra.detailed_description = clip(text(descMod.detailedDescription), 3000);
    baseExtra.design_info = designMod.designInfo ?? null;
    baseExtra.arm_groups = asArray(armsMod.armGroups).slice(0, 15);
    baseExtra.primary_outcomes = asArray(outMod.primaryOutcomes).slice(0, 15);
    baseExtra.secondary_outcomes = asArray(outMod.secondaryOutcomes).slice(0, 15);
    baseExtra.eligibility_summary = {
      sex: text(eligMod.sex),
      minimum_age: text(eligMod.minimumAge),
      maximum_age: text(eligMod.maximumAge),
      healthy_volunteers: eligMod.healthyVolunteers ?? null,
    };
    baseExtra.locations_sample = asArray(locMod.locations)
      .slice(0, 20)
      .map((l) => {
        const lr = asRecord(l);
        return {
          facility: text(lr.facility),
          city: text(lr.city),
          state: text(lr.state),
          country: text(lr.country),
          status: text(lr.status),
        };
      });
    baseExtra.references = asArray(refsMod.references).slice(0, 20);
  } else if (detailMode === "eligibility") {
    baseExtra.sex = text(eligMod.sex);
    baseExtra.minimum_age = text(eligMod.minimumAge);
    baseExtra.maximum_age = text(eligMod.maximumAge);
    baseExtra.std_ages = asArray(eligMod.stdAges);
    baseExtra.healthy_volunteers = eligMod.healthyVolunteers ?? null;
    baseExtra.study_population = text(eligMod.studyPopulation);
    baseExtra.sampling_method = text(eligMod.samplingMethod);
    baseExtra.eligibility_criteria = clip(text(eligMod.eligibilityCriteria), 8000);
  } else if (detailMode === "results") {
    const outMeasures = asArray(asRecord(resultsSec.outcomeMeasuresModule).outcomeMeasures)
      .slice(0, 10)
      .map((m) => {
        const mr = asRecord(m);
        return {
          type: text(mr.type),
          title: text(mr.title),
          description: clip(text(mr.description), 500),
          time_frame: text(mr.timeFrame),
          unit_of_measure: text(mr.unitOfMeasure),
          groups: asArray(mr.groups).slice(0, 10),
          analyses: asArray(mr.analyses).slice(0, 10),
        };
      });
    const aeMod = asRecord(resultsSec.adverseEventsModule);
    baseExtra.outcome_measures = outMeasures;
    baseExtra.adverse_events_summary = {
      frequency_threshold: text(aeMod.frequencyThreshold),
      time_frame: text(aeMod.timeFrame),
      event_groups: asArray(aeMod.eventGroups).slice(0, 10),
      serious_events_count: asArray(aeMod.seriousEvents).length,
      other_events_count: asArray(aeMod.otherEvents).length,
      top_serious_events: asArray(aeMod.seriousEvents).slice(0, 15),
    };
  }

  return dataset({
    id: `nct:${nctId}`,
    xrefs,
    title: briefTitle ?? officialTitle ?? nctId,
    description: clip(text(descMod.briefSummary), detailMode === "summary" ? 500 : 2000),
    assay: phases.length > 0 ? `${studyType} (${phases.join(", ")})` : studyType,
    taxon_id: "NCBITaxon:9606",
    sample_count: enrollment,
    access: "open",
    url: `https://clinicaltrials.gov/study/${nctId}`,
    source_extra: baseExtra,
  });
}

async function fetchSingleStudy(client: JsonClient, rawNctId: string): Promise<{ nctId: string; study: Record<string, unknown> } | ToolResult> {
  const nctId = normalizeNctId(rawNctId);
  if (!nctId) {
    return failure(ctgovSource, "invalid_input", `Expected a valid NCT identifier such as 'NCT04280705' (got '${rawNctId}').`, false);
  }
  try {
    const body = asRecord(await client.getJson(`${API_BASE}/studies/${nctId}`));
    if (!body.protocolSection) {
      return failure(ctgovSource, "not_found", `Clinical trial '${nctId}' was not found on ClinicalTrials.gov.`, false);
    }
    return { nctId, study: body };
  } catch (err) {
    if (err instanceof Error && "status" in err && (err as { status?: number }).status === 404) {
      return failure(ctgovSource, "not_found", `Clinical trial '${nctId}' was not found on ClinicalTrials.gov.`, false);
    }
    const msg = err instanceof Error ? err.message : String(err);
    return failure(ctgovSource, "upstream_http", msg, true);
  }
}

export const tools: OmicsTool[] = [
  {
    name: "usage_terms",
    title: "Usage terms",
    description: "Return usage terms and NLM ClinicalTrials.gov REST API v2 notes.",
    inputSchema: z.object({}),
    call: async () =>
      makePage("dataset", [], 0, 1, 0, ctgovSource, {
        summary: USAGE_SUMMARY,
        url: "https://clinicaltrials.gov/data-api/about-api",
      }),
  },
  {
    name: "ctgov_search_studies",
    title: "Search ClinicalTrials.gov studies",
    description:
      "Search registered global clinical trials on ClinicalTrials.gov by condition/disease, drug/intervention, general keyword/biomarker, sponsor, location, recruitment status, trial phase, study type, and posted-results availability.",
    inputSchema: z.object({
      condition: z.string().optional().describe("Disease or condition query, e.g. 'hepatocellular carcinoma', 'non-small cell lung cancer'."),
      intervention: z.string().optional().describe("Drug, biologic, device, or procedure query, e.g. 'pembrolizumab', 'sorafenib'."),
      query: z.string().optional().describe("General keyword, target gene, biomarker, or NCT query, e.g. 'KRAS G12C', 'ctDNA'."),
      sponsor: z.string().optional().describe("Lead sponsor or collaborator organization name."),
      location: z.string().optional().describe("Country, state, city, or hospital/facility name."),
      status: z.array(z.enum(STATUS_VALUES)).optional().describe("Filter by recruitment status(es), e.g. ['RECRUITING', 'ACTIVE_NOT_RECRUITING']."),
      phase: z.enum(PHASE_VALUES).optional().describe("Filter by clinical trial phase ('EARLY_PHASE1', 'PHASE1', 'PHASE2', 'PHASE3', 'PHASE4', 'NA')."),
      study_type: z.enum(["INTERVENTIONAL", "OBSERVATIONAL", "EXPANDED_ACCESS"]).optional().describe("Filter by study design type."),
      has_results: z.boolean().optional().describe("Filter to trials with (true) or without (false) posted results on ClinicalTrials.gov."),
      page_token: z.string().optional().describe("Next-page cursor token returned in source_extra.next_page_token."),
      limit: z.number().int().min(1).max(50).optional().describe("Number of studies to return (1..50, default 20)."),
    }),
    call: async (client, args) => {
      const parsed = z
        .object({
          condition: z.string().optional(),
          intervention: z.string().optional(),
          query: z.string().optional(),
          sponsor: z.string().optional(),
          location: z.string().optional(),
          status: z.array(z.enum(STATUS_VALUES)).optional(),
          phase: z.enum(PHASE_VALUES).optional(),
          study_type: z.enum(["INTERVENTIONAL", "OBSERVATIONAL", "EXPANDED_ACCESS"]).optional(),
          has_results: z.boolean().optional(),
          page_token: z.string().optional(),
          limit: z.number().int().min(1).max(50).optional(),
        })
        .safeParse(args);
      if (!parsed.success) {
        return failure(ctgovSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
      }
      const p = parsed.data;
      if (
        !p.condition?.trim() &&
        !p.intervention?.trim() &&
        !p.query?.trim() &&
        !p.sponsor?.trim() &&
        !p.location?.trim() &&
        (!p.status || p.status.length === 0)
      ) {
        return failure(
          ctgovSource,
          "invalid_input",
          "Provide at least one search filter (condition, intervention, query, sponsor, location, or status).",
          false,
        );
      }

      const url = new URL(`${API_BASE}/studies`);
      url.searchParams.set("format", "json");
      url.searchParams.set("countTotal", "true");
      url.searchParams.set("pageSize", String(p.limit ?? 20));
      if (p.condition?.trim()) url.searchParams.set("query.cond", p.condition.trim());
      if (p.intervention?.trim()) url.searchParams.set("query.intr", p.intervention.trim());
      if (p.sponsor?.trim()) url.searchParams.set("query.spons", p.sponsor.trim());
      if (p.location?.trim()) url.searchParams.set("query.locn", p.location.trim());
      if (p.status && p.status.length > 0) url.searchParams.set("filter.overallStatus", p.status.join(","));
      if (p.page_token?.trim()) url.searchParams.set("pageToken", p.page_token.trim());

      const termParts: string[] = [];
      if (p.query?.trim()) termParts.push(`(${p.query.trim()})`);
      if (p.phase) termParts.push(`AREA[Phase]${p.phase}`);
      if (p.study_type) termParts.push(`AREA[StudyType]${p.study_type}`);
      if (p.has_results === true) termParts.push("AREA[ResultsFirstPostDate]RANGE[MIN,MAX]");
      else if (p.has_results === false) termParts.push("NOT AREA[ResultsFirstPostDate]RANGE[MIN,MAX]");
      if (termParts.length > 0) url.searchParams.set("query.term", termParts.join(" AND "));

      try {
        const body = asRecord(await client.getJson(url.toString()));
        const studies = asArray(body.studies);
        const total = numberOrNull(body.totalCount);
        const nextToken = text(body.nextPageToken);
        const items = safeMap(studies, (s) => studyToDataset(s, "summary"));
        return makePage(
          "dataset",
          items,
          0,
          p.limit ?? 20,
          total,
          ctgovSource,
          { summary: USAGE_SUMMARY, url: "https://clinicaltrials.gov/" },
          { next_page_token: nextToken },
        );
      } catch (err) {
        const msg = err instanceof Error ? err.message : String(err);
        return failure(ctgovSource, "upstream_http", msg, true);
      }
    },
  },
  {
    name: "ctgov_get_study",
    title: "Get full ClinicalTrials.gov study protocol",
    description:
      "Retrieve the complete study registration record for one NCT identifier (design, arms, interventions, primary and secondary endpoints, enrollment, sponsors, locations, and linked PubMed citations).",
    inputSchema: z.object({
      nct_id: z.string().min(1).describe("ClinicalTrials.gov identifier, e.g. 'NCT04280705' or 'nct:NCT04280705'."),
    }),
    call: async (client, args) => {
      const parsed = z.object({ nct_id: z.string().min(1) }).safeParse(args);
      if (!parsed.success) return failure(ctgovSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
      const fetched = await fetchSingleStudy(client, parsed.data.nct_id);
      if ("entity_type" in fetched) return fetched;
      const item = studyToDataset(fetched.study, "full");
      return makePage("dataset", [item], 0, 1, 1, ctgovSource, {
        summary: USAGE_SUMMARY,
        url: `https://clinicaltrials.gov/study/${fetched.nctId}`,
      });
    },
  },
  {
    name: "ctgov_get_eligibility",
    title: "Get study inclusion and exclusion criteria",
    description:
      "Retrieve full inclusion/exclusion eligibility criteria text, eligible age range, sex, and healthy-volunteer rules for one NCT clinical trial.",
    inputSchema: z.object({
      nct_id: z.string().min(1).describe("ClinicalTrials.gov identifier, e.g. 'NCT04280705'."),
    }),
    call: async (client, args) => {
      const parsed = z.object({ nct_id: z.string().min(1) }).safeParse(args);
      if (!parsed.success) return failure(ctgovSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
      const fetched = await fetchSingleStudy(client, parsed.data.nct_id);
      if ("entity_type" in fetched) return fetched;
      const item = studyToDataset(fetched.study, "eligibility");
      return makePage("dataset", [item], 0, 1, 1, ctgovSource, {
        summary: USAGE_SUMMARY,
        url: `https://clinicaltrials.gov/study/${fetched.nctId}`,
      });
    },
  },
  {
    name: "ctgov_get_results",
    title: "Get posted trial outcomes and adverse events",
    description:
      "Retrieve posted clinical trial outcome measures, statistical analyses, and serious/other adverse-event summaries from ClinicalTrials.gov resultsSection (for completed trials with posted results).",
    inputSchema: z.object({
      nct_id: z.string().min(1).describe("ClinicalTrials.gov identifier, e.g. 'NCT04280705'."),
    }),
    call: async (client, args) => {
      const parsed = z.object({ nct_id: z.string().min(1) }).safeParse(args);
      if (!parsed.success) return failure(ctgovSource, "invalid_input", parsed.error.issues[0]?.message ?? "invalid input", false);
      const fetched = await fetchSingleStudy(client, parsed.data.nct_id);
      if ("entity_type" in fetched) return fetched;
      if (!fetched.study.hasResults && !fetched.study.resultsSection) {
        return failure(
          ctgovSource,
          "not_found",
          `Study '${fetched.nctId}' does not have posted summary results on ClinicalTrials.gov yet.`,
          false,
        );
      }
      const item = studyToDataset(fetched.study, "results");
      return makePage("dataset", [item], 0, 1, 1, ctgovSource, {
        summary: USAGE_SUMMARY,
        url: `https://clinicaltrials.gov/study/${fetched.nctId}?tab=results`,
      });
    },
  },
];

export const serverModule: ServerModule = {
  name: "clinical-trials",
  version: "1.0.0",
  mcpPath: "/mcp",
  tools,
};
