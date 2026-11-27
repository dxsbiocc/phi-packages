import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

const sampleStudy = {
  hasResults: true,
  protocolSection: {
    identificationModule: {
      nctId: "NCT04280705",
      briefTitle: "Adaptive COVID-19 Treatment Trial (ACTT-1)",
      officialTitle: "A Multicenter, Adaptive, Randomized Blinded Controlled Trial",
    },
    statusModule: {
      overallStatus: "COMPLETED",
      startDateStruct: { date: "2020-02-21" },
      primaryCompletionDateStruct: { date: "2020-04-29" },
      completionDateStruct: { date: "2020-05-21" },
      lastUpdatePostDateStruct: { date: "2021-07-16" },
    },
    descriptionModule: {
      briefSummary: "Randomized controlled trial evaluating remdesivir in hospitalized adults.",
    },
    conditionsModule: {
      conditions: ["COVID-19"],
    },
    designModule: {
      studyType: "INTERVENTIONAL",
      phases: ["PHASE3"],
      enrollmentInfo: { count: 1062, type: "ACTUAL" },
    },
    armsInterventionsModule: {
      interventions: [{ type: "DRUG", name: "Remdesivir", description: "IV infusion" }],
    },
    outcomesModule: {
      primaryOutcomes: [{ measure: "Time to recovery", timeFrame: "Day 1 through Day 29" }],
    },
    eligibilityModule: {
      sex: "ALL",
      minimumAge: "18 Years",
      maximumAge: "99 Years",
      healthyVolunteers: false,
      eligibilityCriteria: "Inclusion Criteria:\n- Adult hospitalized patient\nExclusion Criteria:\n- ALT/AST > 5x ULN",
    },
    sponsorCollaboratorsModule: {
      leadSponsor: { name: "National Institute of Allergy and Infectious Diseases (NIAID)" },
    },
    referencesModule: {
      references: [{ pmid: "32445440", citation: "Beigel JH et al. N Engl J Med 2020." }],
    },
  },
  resultsSection: {
    outcomeMeasuresModule: {
      outcomeMeasures: [
        {
          type: "PRIMARY",
          title: "Time to Recovery",
          timeFrame: "29 days",
          unitOfMeasure: "Days",
        },
      ],
    },
    adverseEventsModule: {
      frequencyThreshold: "5",
      seriousEvents: [{ term: "Respiratory failure" }],
      otherEvents: [{ term: "Anemia" }],
    },
  },
};

const mockClient: JsonClient = {
  async getJson(url: string) {
    if (url.includes("/studies/NCT04280705")) {
      return sampleStudy;
    }
    if (url.includes("/studies?")) {
      return {
        studies: [sampleStudy],
        totalCount: 1,
        nextPageToken: null,
      };
    }
    return {};
  },
  async postJson() {
    return {};
  },
};

test("every clinical-trials tool returns a valid ToolResult schema", async () => {
  assert.equal(tools.length, 5);
  const calls: Record<string, Record<string, unknown>> = {
    usage_terms: {},
    ctgov_search_studies: { condition: "COVID-19", intervention: "Remdesivir", phase: "PHASE3" },
    ctgov_get_study: { nct_id: "NCT04280705" },
    ctgov_get_eligibility: { nct_id: "nct:NCT04280705" },
    ctgov_get_results: { nct_id: "NCT04280705" },
  };
  for (const tool of tools) {
    const res = ToolResultSchema.parse(await tool.call(mockClient, calls[tool.name] ?? {}));
    assert.equal(res.entity_type, "page", tool.name);
    assert.equal(res.source.id, "clinicaltrials-gov");
  }

  const studyPage = ToolResultSchema.parse(
    await tools.find((t) => t.name === "ctgov_get_study")!.call(mockClient, { nct_id: "NCT04280705" }),
  );
  if (studyPage.entity_type === "page") {
    const ds = studyPage.items[0];
    assert.equal(ds?.id, "nct:NCT04280705");
    assert.ok(ds?.xrefs.includes("pmid:32445440"));
  }
});

test("invalid NCT id or empty search returns clean invalid_input error", async () => {
  const badId = ToolResultSchema.parse(
    await tools.find((t) => t.name === "ctgov_get_study")!.call(mockClient, { nct_id: "INVALID123" }),
  );
  assert.equal(badId.entity_type, "error");
  if (badId.entity_type === "error") assert.equal(badId.code, "invalid_input");

  const emptySearch = ToolResultSchema.parse(
    await tools.find((t) => t.name === "ctgov_search_studies")!.call(mockClient, {}),
  );
  assert.equal(emptySearch.entity_type, "error");
  if (emptySearch.entity_type === "error") assert.equal(emptySearch.code, "invalid_input");
});
