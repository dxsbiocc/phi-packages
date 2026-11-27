import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import { UpstreamError, type JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

function molecule() {
  return {
    molecule_chembl_id: "CHEMBL941",
    pref_name: "IMATINIB",
    molecule_type: "Small molecule",
    max_phase: 4,
    first_approval: 2001,
    oral: true,
    molecule_structures: {
      canonical_smiles: "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
      standard_inchi_key: "KTUFNOKKBVMGRW-UHFFFAOYSA-N",
    },
    molecule_properties: {
      full_molformula: "C29H31N7O",
      full_mwt: "493.60",
      alogp: "4.59",
      psa: "86.28",
      hba: 7,
      hbd: 2,
      num_ro5_violations: 0,
      qed_weighted: "0.39",
    },
    molecule_synonyms: [{ molecule_synonym: "Gleevec" }, { molecule_synonym: "STI-571" }],
  };
}

function client(): JsonClient {
  return {
    async getJson(url: string) {
      if (url.includes("chembl/api/data/molecule/search.json")) {
        return { page_meta: { total_count: 1 }, molecules: [molecule()] };
      }
      if (url.includes("chembl/api/data/molecule/CHEMBL941.json")) {
        return molecule();
      }
      if (url.includes("chembl/api/data/activity.json")) {
        return {
          page_meta: { total_count: 1 },
          activities: [{
            activity_id: 1001,
            molecule_chembl_id: "CHEMBL941",
            assay_chembl_id: "CHEMBL829152",
            standard_type: "IC50",
            standard_value: "11.0",
            standard_units: "nM",
            standard_relation: "=",
            pchembl_value: "7.96",
            target_chembl_id: "CHEMBL1862",
            target_pref_name: "Tyrosine-protein kinase ABL1",
            target_organism: "Homo sapiens",
            assay_type: "B",
            document_chembl_id: "CHEMBL1135147",
          }],
        };
      }
      if (url.includes("chembl/api/data/mechanism.json")) {
        return {
          page_meta: { total_count: 1 },
          mechanisms: [{
            mec_id: 215,
            molecule_chembl_id: "CHEMBL941",
            target_chembl_id: "CHEMBL1862",
            action_type: "INHIBITOR",
            mechanism_of_action: "Tyrosine-protein kinase ABL inhibitor",
            direct_interaction: true,
            molecular_mechanism: true,
            max_phase: 4,
          }],
        };
      }
      if (url.includes("api.fda.gov/drug/label.json")) {
        return {
          meta: { results: { total: 1 } },
          results: [{
            id: "lbl-1",
            set_id: "set-123",
            effective_time: "20240101",
            openfda: { generic_name: ["IMATINIB MESYLATE"], brand_name: ["GLEEVEC"], route: ["ORAL"] },
            indications_and_usage: ["Chronic myeloid leukemia."],
            warnings_and_cautions: ["Fluid retention and edema."],
          }],
        };
      }
      if (url.includes("api.fda.gov/drug/event.json")) {
        return {
          results: [{ term: "NAUSEA", count: 4210 }],
        };
      }
      throw new Error(`unmocked ${url}`);
    },
    async postJson() { return {}; },
  };
}

const args: Record<string, Record<string, unknown>> = {
  usage_terms: {},
  chembl_search_molecules: { query: "imatinib", limit: 5 },
  chembl_get_molecule: { chembl_id: "CHEMBL941" },
  chembl_get_activities: { molecule_chembl_id: "CHEMBL941", standard_type: "IC50", limit: 5 },
  chembl_get_mechanisms: { molecule_chembl_id: "CHEMBL941", limit: 5 },
  openfda_search_labels: { drug: "imatinib", limit: 5 },
  openfda_adverse_events: { drug: "imatinib", serious_only: true, limit: 5 },
};

test("every drug-chemistry tool returns the shared schema", async () => {
  const http = client();
  assert.equal(tools.length, 7);
  for (const tool of tools) {
    const result = ToolResultSchema.parse(await tool.call(http, args[tool.name] ?? {}));
    assert.equal(result.entity_type, "page", `${tool.name} ${JSON.stringify(result).slice(0, 200)}`);
    if (result.entity_type === "page" && tool.name !== "usage_terms") {
      assert.ok(result.items.length > 0, tool.name);
    }
  }
});

test("openFDA 404 maps cleanly to not_found and invalid ChEMBL id rejects cleanly", async () => {
  const bad = ToolResultSchema.parse(await tools.find((t) => t.name === "chembl_get_molecule")!.call(client(), { chembl_id: "NOTACHMBL" }));
  assert.equal(bad.entity_type, "error");
  if (bad.entity_type === "error") assert.equal(bad.code, "invalid_input");

  const http = client();
  http.getJson = async () => {
    throw new UpstreamError(404, "404 Not Found", false);
  };
  const notFound = ToolResultSchema.parse(await tools.find((t) => t.name === "openfda_search_labels")!.call(http, { drug: "no_such_drug_xyz" }));
  assert.equal(notFound.entity_type, "error");
  if (notFound.entity_type === "error") assert.equal(notFound.code, "not_found");
});
