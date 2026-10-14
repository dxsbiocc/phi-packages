import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import { UpstreamError, type JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

const hmdbXml = `<?xml version="1.0" encoding="UTF-8"?>
<metabolite xmlns="http://www.hmdb.ca">
  <accession>HMDB0000001</accession>
  <status>quantified</status>
  <secondary_accessions><accession>HMDB00001</accession></secondary_accessions>
  <name>1-Methylhistidine</name>
  <description>Derived from dietary anserine.</description>
  <chemical_formula>C7H11N3O2</chemical_formula>
  <iupac_name>(2S)-2-amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid</iupac_name>
  <cas_registry_number>332-80-9</cas_registry_number>
  <smiles>CN1C=NC(C[C@H](N)C(O)=O)=C1</smiles>
  <inchi>InChI=1S/C7H11N3O2</inchi>
  <inchikey>BRMWTNUJHUMWMS-LURJTMIESA-N</inchikey>
  <kegg_id>C01152</kegg_id>
  <pubchem_compound_id>92105</pubchem_compound_id>
  <chebi_id>50599</chebi_id>
  <pathways><pathway><name>Histidine metabolism</name></pathway></pathways>
</metabolite>`;

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
    async getText(url: string) {
      if (url.endsWith("/HMDB0000001.xml")) return hmdbXml;
      if (url.includes("unearth/q")) {
        return `<a href="/metabolites/HMDB0000001">1-Methylhistidine</a><a href="/metabolites/HMDB0000122">D-Glucose</a><a href="/metabolites/HMDB0000001">1-Methylhistidine</a>`;
      }
      throw new Error(`unmocked ${url}`);
    },
  };
}

const args: Record<string, Record<string, unknown>> = {
  usage_terms: {},
  metabolights_list_studies: { limit: 5 },
  metabolights_get_study: { study_id: "MTBLS1" },
  metabolights_list_files: { study_id: "MTBLS1", limit: 5 },
  metabolights_list_assays: { study_id: "MTBLS1", limit: 5 },
  hmdb_get_metabolite: { hmdb_id: "HMDB0000001" },
  hmdb_search_metabolites: { query: "histidine", limit: 5 },
};

test("every metabolomics tool returns the shared schema", async () => {
  const http = client();
  assert.equal(tools.length, 7);
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

test("HMDB XML maps onto the shared metabolite schema", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "hmdb_get_metabolite")!.call(client(), { hmdb_id: "hmdb:HMDB1" }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type !== "page" || result.items[0]?.entity_type !== "metabolite") return;
  assert.equal(result.items[0].id, "hmdb:HMDB0000001");
  assert.equal(result.items[0].name, "1-Methylhistidine");
  assert.equal(result.items[0].formula, "C7H11N3O2");
  assert.equal(result.items[0].inchikey, "BRMWTNUJHUMWMS-LURJTMIESA-N");
  assert.ok(result.items[0].xrefs.includes("CHEBI:50599"));
  assert.ok(result.items[0].xrefs.includes("kegg:C01152"));
  assert.deepEqual(result.items[0].source_extra.secondary_accessions, ["HMDB00001"]);
  assert.equal(result.source.id, "hmdb");
});

test("HMDB search returns metabolite ids from the HTML page", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "hmdb_search_metabolites")!.call(client(), { query: "histidine", offset: 1, limit: 1 }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type !== "page") return;
  assert.equal(result.item_type, "metabolite");
  assert.equal(result.items.length, 1);
  assert.equal(result.items[0]?.id, "hmdb:HMDB0000122");
  assert.equal(result.page.total, 2);
});

test("HMDB rejects an accession that is not an HMDB id", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "hmdb_get_metabolite")!.call(client(), { hmdb_id: "glucose" }));
  assert.equal(result.entity_type, "error");
  if (result.entity_type === "error") assert.equal(result.code, "invalid_input");
});

test("an HMDB Cloudflare challenge is retryable", async () => {
  const http = client();
  http.getText = async () => {
    throw new UpstreamError(403, "403 <title>Just a moment...</title>", false);
  };
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "hmdb_get_metabolite")!.call(http, { hmdb_id: "HMDB0000001" }));
  assert.equal(result.entity_type, "error");
  if (result.entity_type === "error") assert.equal(result.retryable, true);
});

test("private MetaboLights files omit URLs", async () => {
  const http = client();
  http.getJson = async () => ({ private: true, study: [{ file: "secret.mzML", directory: false }] });
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "metabolights_list_files")!.call(http, { study_id: "MTBLS1", limit: 5 }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type === "page" && result.items[0]?.entity_type === "file") assert.equal(result.items[0].url, null);
});
