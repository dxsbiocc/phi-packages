import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

const experiment = { accession: "ENCSR818VNL", assay_title: "Control ChIP-seq", status: "released", description: "chip", biosample_summary: "ovary", biosample_ontology: { term_id: "UBERON:0000992", term_name: "ovary", classification: "tissue" }, award: { rfa: "Roadmap" } };
const encFile = { accession: "ENCFF372VWB", file_format: "tar", file_size: 10, href: "/files/ENCFF372VWB/@@download/ENCFF372VWB.tar.gz", status: "released", output_type: "peaks", restricted: true, dataset: "/experiments/ENCSR818VNL/" };
const biosample = { accession: "ENCBS997DQU", summary: "HEK293", biosample_ontology: { term_id: "EFO:0001184", term_name: "HEK293", classification: "cell line" }, organism: { scientific_name: "Homo sapiens" } };

function client(): JsonClient {
  return {
    async getJson(url: string) {
      if (url.includes("/files/ENCFF")) return encFile;
      if (url.includes("/biosamples/ENCBS")) return biosample;
      if (url.includes("/experiments/ENCSR")) return experiment;
      if (url.includes("type=File")) return { "@graph": [encFile], total: 1 };
      if (url.includes("type=Biosample")) return { "@graph": [biosample], total: 1 };
      return { "@graph": [experiment], total: 1 };
    },
    async postJson() { return {}; },
  };
}

const args: Record<string, Record<string, unknown>> = {
  usage_terms: {},
  encode_search_experiments: { limit: 5 },
  encode_get_experiment: { accession: "ENCSR818VNL" },
  encode_search_roadmap: { limit: 5 },
  encode_search_files: { limit: 5 },
  encode_get_file: { accession: "ENCFF372VWB" },
  encode_search_biosamples: { limit: 5 },
  encode_get_biosample: { accession: "ENCBS997DQU" },
};

test("every epigenomics tool returns the shared schema", async () => {
  const http = client();
  for (const tool of tools) {
    const result = ToolResultSchema.parse(await tool.call(http, args[tool.name] ?? {}));
    assert.equal(result.entity_type, "page", tool.name);
  }
});

test("restricted ENCODE files omit the download URL", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "encode_get_file")!.call(client(), args.encode_get_file));
  assert.equal(result.entity_type, "page");
  if (result.entity_type === "page" && result.items[0]?.entity_type === "file") {
    assert.equal(result.items[0].url, null);
    assert.equal(result.items[0].access, "controlled");
  }
});
