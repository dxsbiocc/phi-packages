import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

function uniprotEntry() {
  return {
    primaryAccession: "P04637",
    uniProtkbId: "P53_HUMAN",
    entryType: "UniProtKB reviewed (Swiss-Prot)",
    organism: { scientificName: "Homo sapiens", taxonId: 9606 },
    proteinDescription: { recommendedName: { fullName: { value: "Cellular tumor antigen p53" } } },
    genes: [{ geneName: { value: "TP53" } }],
    sequence: { length: 393, value: "MEEPQSDPSV" },
    comments: [
      { commentType: "FUNCTION", texts: [{ value: "Acts as a tumor suppressor in many tumor types." }] },
      { commentType: "DISEASE", disease: { diseaseId: "Li-Fraumeni syndrome", acronym: "LFS1", description: "Cancer predisposition." } },
    ],
    features: [
      { type: "Domain", description: "p53 DNA-binding", location: { start: { value: 102 }, end: { value: 292 } } },
      { type: "Binding site", description: "DNA", location: { start: { value: 248 }, end: { value: 248 } } },
    ],
    uniProtKBCrossReferences: [
      { database: "PDB", id: "1TUP" },
      { database: "Ensembl", id: "ENST00000269305", properties: [{ key: "GeneId", value: "ENSG00000141510.18" }] },
    ],
  };
}

function alphafoldEntry() {
  return [{
    entryId: "AF-P04637-F1",
    gene: "TP53",
    taxId: 9606,
    uniprotAccession: "P04637",
    uniprotDescription: "Cellular tumor antigen p53",
    modelEntityId: "AF-P04637-F1",
    globalMetricValue: 81.4,
    fractionPlddtVeryHigh: 0.61,
    fractionPlddtConfident: 0.14,
    fractionPlddtLow: 0.11,
    fractionPlddtVeryLow: 0.14,
    uniprotStart: 1,
    uniprotEnd: 393,
    latestVersion: 4,
    cifUrl: "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v4.cif",
    pdbUrl: "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v4.pdb",
    bcifUrl: "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v4.bcif",
    paeDocUrl: "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-predicted_aligned_error_v4.json",
    paeImageUrl: "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-predicted_aligned_error_v4.png",
  }];
}

function client(): JsonClient {
  return {
    async getJson(url: string) {
      if (url.includes("rest.uniprot.org/uniprotkb/search")) return { results: [uniprotEntry()] };
      if (url.includes("rest.uniprot.org/uniprotkb/P04637.json")) return uniprotEntry();
      if (url.includes("alphafold.ebi.ac.uk/api/prediction/P04637")) return alphafoldEntry();
      if (url.includes("alphafold.ebi.ac.uk/api/uniprot/summary/P04637.json")) {
        return {
          uniprot_entry: { ac: "P04637", id: "P53_HUMAN", sequence_length: 393 },
          structures: [{ summary: { model_identifier: "AF-P04637-F1", confidence_type: "pLDDT", confidence_avg_local_score: 81.4, coverage: 1 } }],
        };
      }
      throw new Error(`unmocked ${url}`);
    },
    async postJson() { return {}; },
  };
}

const args: Record<string, Record<string, unknown>> = {
  usage_terms: {},
  uniprot_search_proteins: { query: "TP53", reviewed_only: true, organism_id: 9606, limit: 5 },
  uniprot_get_protein: { accession: "uniprot:P04637", include_sequence: true },
  alphafold_get_prediction: { uniprot_accession: "P04637" },
  alphafold_list_files: { uniprot_accession: "P04637", limit: 5 },
  alphafold_get_summary: { uniprot_accession: "P04637" },
};

test("every protein-structure tool returns the shared schema", async () => {
  const http = client();
  assert.equal(tools.length, 6);
  for (const tool of tools) {
    const result = ToolResultSchema.parse(await tool.call(http, args[tool.name] ?? {}));
    assert.equal(result.entity_type, "page", `${tool.name} ${JSON.stringify(result).slice(0, 200)}`);
    if (result.entity_type === "page" && tool.name !== "usage_terms") {
      assert.ok(result.items.length > 0, tool.name);
    }
  }
});

test("uniprot_get_protein includes sequence only when requested and alphafold_list_files emits open files", async () => {
  const http = client();
  const withoutSeq = ToolResultSchema.parse(await tools.find((t) => t.name === "uniprot_get_protein")!.call(http, { accession: "P04637" }));
  assert.equal(withoutSeq.entity_type, "page");
  if (withoutSeq.entity_type === "page") {
    assert.equal(withoutSeq.items[0]?.id, "uniprot:P04637");
    assert.equal(withoutSeq.items[0]?.source_extra.sequence, undefined);
  }
  const filesRes = ToolResultSchema.parse(await tools.find((t) => t.name === "alphafold_list_files")!.call(http, { uniprot_accession: "P04637", limit: 5 }));
  assert.equal(filesRes.entity_type, "page");
  if (filesRes.entity_type === "page") {
    assert.equal(filesRes.item_type, "file");
    assert.equal(filesRes.items.length, 5);
    assert.equal(filesRes.items[0]?.entity_type === "file" && filesRes.items[0].access, "open");
  }
});
