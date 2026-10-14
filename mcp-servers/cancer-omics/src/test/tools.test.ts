import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

function client(): JsonClient {
  return {
    async getJson(url: string) {
      if (url.includes("/cases")) return { data: { hits: [{ case_id: "6e3b6b72-142d-4b8d-a462-28a205796e41", submitter_id: "TCGA-44-3918", disease_type: "Adenomas and Adenocarcinomas", primary_site: "Bronchus and lung", demographic: { age_at_index: 60 } }], pagination: { total: 1 } } };
      if (url.includes("/files")) return { data: { hits: [{ file_id: "af0d3d92-f8b5-4612-805b-d6773595262b", file_name: "example.maf.gz", data_type: "Annotated Somatic Mutation", access: "controlled", file_size: 10, cases: [{ submitter_id: "C3L-01032" }] }], pagination: { total: 1 } } };
      if (url.includes("cbioportal.org") && url.includes("/genes/")) return { entrezGeneId: 7157, hugoGeneSymbol: "TP53", type: "protein-coding" };
      if (url.includes("/genes/")) return { data: { gene_id: "ENSG00000141510", symbol: "TP53", name: "tumor protein p53", biotype: "protein_coding", synonyms: ["p53"], gene_chromosome: "17", gene_start: 7661779, gene_end: 7687538, gene_strand: -1 } };
      if (url.includes("/ssms")) return { data: { hits: [{ ssm_id: "40b72e89-d985-5d31-ad42-5f54460b9bc9", chromosome: "chr12", start_position: 109768424, end_position: 109768424, genomic_dna_change: "chr12:g.109768424G>A", mutation_type: "SNP", ncbi_build: "GRCh38" }], pagination: { total: 1 } } };
      if (url.includes("/studies/brca_tcga/samples")) return [{ sampleId: "TCGA-AR-A1AR-01", patientId: "TCGA-AR-A1AR", sampleType: "Primary" }];
      if (url.includes("/studies/brca_tcga/patients")) return [{ patientId: "TCGA-AR-A1AR" }];
      if (url.includes("/studies/brca_tcga")) return { studyId: "brca_tcga", name: "Breast", description: "TCGA", publicStudy: true, readPermission: true, allSampleCount: 10, referenceGenome: "hg19" };
      if (url.includes("/studies")) return [{ studyId: "brca_tcga", name: "Breast", publicStudy: true, readPermission: true, allSampleCount: 10 }];
      throw new Error(`unmocked ${url}`);
    },
    async postJson() {
      return { data: { getPaginatedUIStudy: { total: 1, uiStudies: [{ pdc_study_id: "PDC000759", submitter_id_name: "APOLLO", disease_type: "Breast", analytical_fraction: "Phosphoproteome", program_name: "APOLLO" }] } } };
    },
  };
}

const args: Record<string, Record<string, unknown>> = {
  usage_terms: {},
  gdc_search_cases: { limit: 5 },
  gdc_get_case: { case_id: "TCGA-44-3918" },
  gdc_search_files: { limit: 5 },
  gdc_get_gene: { gene: "ENSG00000141510" },
  gdc_search_mutations: { limit: 5 },
  cbioportal_list_studies: { limit: 5 },
  cbioportal_get_study: { study_id: "brca_tcga" },
  cbioportal_list_samples: { study_id: "brca_tcga", limit: 5 },
  cbioportal_list_patients: { study_id: "brca_tcga", limit: 5 },
  cbioportal_get_gene: { gene: "TP53" },
  pdc_list_studies: { limit: 5 },
  pdc_get_study: { study_id: "PDC000759" },
};

test("every cancer tool returns the shared schema", async () => {
  assert.equal(new Set(tools.map((tool) => tool.name)).size, tools.length);
  const http = client();
  for (const tool of tools) {
    const result = ToolResultSchema.parse(await tool.call(http, args[tool.name] ?? {}));
    assert.equal(result.entity_type, "page", tool.name);
  }
});

test("controlled GDC files do not include a download URL", async () => {
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "gdc_search_files")!.call(client(), { limit: 5 }));
  assert.equal(result.entity_type, "page");
  if (result.entity_type === "page") {
    const item = result.items[0];
    assert.equal(item.entity_type, "file");
    if (item.entity_type === "file") {
      assert.equal(item.access, "controlled");
      assert.equal(item.url, null);
    }
  }
});
