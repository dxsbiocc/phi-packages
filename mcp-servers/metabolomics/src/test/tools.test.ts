import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

function client(): JsonClient {
  return {
    async getJson(url: string) {
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
  metabolights_list_studies: { limit: 5 },
  metabolights_get_study: { study_id: "MTBLS1" },
  metabolights_list_files: { study_id: "MTBLS1", limit: 5 },
  metabolights_list_assays: { study_id: "MTBLS1", limit: 5 },
};

test("every metabolomics tool returns the shared schema", async () => {
  const http = client();
  assert.equal(tools.length, 5);
  for (const tool of tools) {
    const result = ToolResultSchema.parse(await tool.call(http, args[tool.name] ?? {}));
    assert.equal(result.entity_type, "page", `${tool.name} ${JSON.stringify(result).slice(0, 200)}`);
  }
});

test("metabolomics usage terms name this server", async () => {
  const result = ToolResultSchema.parse(await tools[0]!.call(client(), {}));
  assert.equal(result.entity_type, "page");
  if (result.entity_type === "page") assert.equal(result.source.id, "metabolomics");
});

test("private MetaboLights files omit URLs", async () => {
  const http = client();
  http.getJson = async () => ({ private: true, study: [{ file: "secret.mzML", directory: false }] });
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "metabolights_list_files")!.call(http, { study_id: "MTBLS1", limit: 5 }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type === "page" && result.items[0]?.entity_type === "file") assert.equal(result.items[0].url, null);
});
