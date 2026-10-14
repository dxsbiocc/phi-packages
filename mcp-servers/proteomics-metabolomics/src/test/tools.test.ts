import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

function client(): JsonClient {
  return {
    async getJson(url: string) {
      if (url.includes("/projects/PXD001357/files")) return [{ accession: "abc", fileName: "sample.dat", fileSizeBytes: 10, fileCategory: { name: "RAW" }, projectAccessions: ["PXD001357"], publicFileLocations: [{ value: "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/sample.dat" }] }];
      if (url.includes("/projects/PXD001357")) return { accession: "PXD001357", title: "milk", projectDescription: "dental calculus" };
      if (url.includes("/projects?")) return [{ accession: "PXD001357", title: "milk", projectDescription: "dental calculus" }];
      if (url.includes("search_download.php")) return [{ Gene: "TP53", "Gene synonym": "p53", Ensembl: "ENSG00000141510", Uniprot: ["P04637"] }];
      if (url.includes("proteinatlas.org/ENSG")) return { Gene: "TP53", "Gene description": "tumor protein p53", Ensembl: "ENSG00000141510", Uniprot: ["P04637"], Chromosome: "17", Position: "7661779-7687538", "RNA tissue specificity": "Low tissue specificity" };
      if (url.endsWith("/studies")) return { content: ["MTBLS1", "MTBLS2"] };
      if (url.endsWith("/assays")) return { data: { assays: [{ filename: "a_MTBLS1.txt" }] } };
      if (url.endsWith("/files")) return { private: false, study: [{ file: "s_MTBLS1.txt", type: "metadata", directory: false }] };
      if (url.includes("/studies/MTBLS1")) return { mtblsStudy: { studyStatus: "Public", studyHttpUrl: "https://example.test/MTBLS1", isaInvestigation: { studies: [{ title: "Study", description: "desc" }] } } };
      throw new Error(`unmocked ${url}`);
    },
    async postJson() { return {}; },
  };
}

const args: Record<string, Record<string, unknown>> = {
  usage_terms: {},
  pride_search_projects: { limit: 5 },
  pride_get_project: { accession: "PXD001357" },
  pride_list_files: { accession: "PXD001357", limit: 5 },
  hpa_search_genes: { query: "TP53", limit: 5 },
  hpa_get_gene: { gene: "ENSG00000141510" },
  hpa_list_proteins: { gene: "ENSG00000141510", limit: 5 },
  metabolights_list_studies: { limit: 5 },
  metabolights_get_study: { study_id: "MTBLS1" },
  metabolights_list_files: { study_id: "MTBLS1", limit: 5 },
  metabolights_list_assays: { study_id: "MTBLS1", limit: 5 },
};

test("every proteomics tool returns the shared schema", async () => {
  const http = client();
  for (const tool of tools) {
    const result = ToolResultSchema.parse(await tool.call(http, args[tool.name] ?? {}));
    assert.equal(result.entity_type, "page", `${tool.name} ${JSON.stringify(result).slice(0, 200)}`);
  }
});

test("private MetaboLights files omit URLs", async () => {
  const http = client();
  http.getJson = async () => ({ private: true, study: [{ file: "secret.mzML", directory: false }] });
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "metabolights_list_files")!.call(http, { study_id: "MTBLS1", limit: 5 }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type === "page" && result.items[0]?.entity_type === "file") assert.equal(result.items[0].url, null);
});
