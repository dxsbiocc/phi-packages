import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

const gene = { gencodeId: "ENSG00000141510.18", geneSymbol: "TP53", description: "tumor protein p53", geneType: "protein coding", chromosome: "chr17", start: 7661779, end: 7687538, strand: "-", genomeBuild: "GRCh38/hg38", gencodeVersion: "v39", entrezGeneId: 7157 };
const collection = {
  collection_id: "af893e86-8e9f-41f1-a474-ef05359b1fb7",
  name: "Retina",
  description: "retina collection",
  collection_url: "https://cellxgene.cziscience.com/collections/af893e86-8e9f-41f1-a474-ef05359b1fb7",
  visibility: "PUBLIC",
  datasets: [{
    dataset_id: "ed419b4e-db9b-40f1-8593-68fdf8dfb076",
    title: "Non-neuronal cells",
    explorer_url: "https://cellxgene.cziscience.com/e/ed419b4e-db9b-40f1-8593-68fdf8dfb076.cxg/",
    cell_count: 10,
    assay: [{ label: "10x 3' v3", ontology_term_id: "EFO:0009922" }],
    tissue: [{ label: "fovea centralis", ontology_term_id: "UBERON:0001786" }],
    cell_type: [{ label: "Mueller cell", ontology_term_id: "CL:0000636" }],
    disease: [{ label: "normal", ontology_term_id: "PATO:0000461" }],
    organism: [{ label: "Homo sapiens", ontology_term_id: "NCBITaxon:9606" }],
  }],
};

function client(): JsonClient {
  return {
    async getJson(url: string) {
      if (url.includes("/reference/gene")) return { data: [gene], paging_info: { totalNumberOfItems: 1 } };
      if (url.includes("/tissueSiteDetail")) return { data: [{ tissueSiteDetailId: "Whole_Blood", tissueSiteDetail: "Whole Blood", ontologyId: "UBERON:0013756", datasetId: "gtex_v10", rnaSeqSampleSummary: { totalCount: 800 } }], paging_info: { totalNumberOfItems: 1 } };
      if (url.includes("/medianGeneExpression")) return { data: [{ median: 22.6, tissueSiteDetailId: "Whole_Blood", ontologyId: "UBERON:0013756", datasetId: "gtex_v10", unit: "TPM" }], paging_info: { totalNumberOfItems: 1 } };
      if (url.includes("/association/egene")) return { data: [{ gencodeId: "ENSG00000227232.5", geneSymbol: "WASH7P", tissueSiteDetailId: "Whole_Blood", log2AllelicFoldChange: 0.7, pValue: 1e-10, qValue: 1e-7, datasetId: "gtex_v10" }], paging_info: { totalNumberOfItems: 1 } };
      if (url.includes("/singleTissueEqtl")) return { data: [{ variantId: "chr1_64764_C_T_b38", snpId: "rs769952832", pos: 64764, chromosome: "chr1", nes: 0.4, pValue: 1e-7, datasetId: "gtex_v10" }], paging_info: { totalNumberOfItems: 1 } };
      if (url.includes("/collections/af893e86")) return collection;
      if (url.endsWith("/collections")) return [collection];
      if (url.includes("canonical_marker_genes")) return [{ symbol: "ASCL2", tissue: "small intestine", publication: "10.1016/j.celrep.2021.108819" }];
      if (url.includes("azul.data.humancellatlas.org")) return { hits: [{ projects: [{ projectId: "74b6d569-3b11-42ef-b6b1-a0454522b4a0", projectTitle: "Brain", projectDescription: "cells", dataUseRestriction: "NRES", estimatedCellCount: 10 }], donorOrganisms: [{ genusSpecies: ["Mus musculus"] }], specimens: [{ organ: ["brain"] }] }], pagination: { total: 1 } };
      throw new Error(`unmocked ${url}`);
    },
    async postJson() {
      return { hits: { total: { relation: "eq", value: 1 }, hits: [{ _source: { hubmap_id: "HBM543.RSRV.265", uuid: "abc", entity_type: "Dataset", data_access_level: "public", title: "dataset", donor: { hubmap_id: "HBM624.KWHJ.622" } } }] } };
    },
  };
}

const args: Record<string, Record<string, unknown>> = {
  usage_terms: {},
  gtex_search_genes: { gene: "TP53", limit: 5 },
  gtex_list_tissues: { limit: 5 },
  gtex_gene_expression: { gene: "TP53", limit: 5 },
  gtex_list_egenes: { tissue: "Whole_Blood", limit: 5 },
  gtex_list_eqtls: { gene: "ENSG00000227232.5", tissue: "Whole_Blood", limit: 5 },
  cellxgene_list_collections: { limit: 5 },
  cellxgene_get_collection: { collection_id: "af893e86-8e9f-41f1-a474-ef05359b1fb7" },
  cellxgene_list_datasets: { collection_id: "af893e86-8e9f-41f1-a474-ef05359b1fb7", limit: 5 },
  cellxgene_get_dataset: { collection_id: "af893e86-8e9f-41f1-a474-ef05359b1fb7", dataset_id: "ed419b4e-db9b-40f1-8593-68fdf8dfb076" },
  cellguide_cell_type: { cell_type_id: "CL:0000034" },
  hca_search_projects: { limit: 5 },
  hca_get_project: { project_id: "74b6d569-3b11-42ef-b6b1-a0454522b4a0" },
  hubmap_search_datasets: { limit: 5 },
  hubmap_search_samples: { limit: 5 },
  hubmap_search_donors: { limit: 5 },
};

test("every transcriptomics tool returns the shared schema", async () => {
  assert.equal(new Set(tools.map((tool) => tool.name)).size, tools.length);
  const http = client();
  for (const tool of tools) {
    const result = ToolResultSchema.parse(await tool.call(http, args[tool.name] ?? {}));
    assert.equal(result.entity_type, "page", tool.name);
  }
});

test("controlled-looking HCA restrictions omit the data URL", async () => {
  const http = client();
  const original = http.getJson;
  http.getJson = async (url) => {
    const body = await original(url) as { hits: Array<{ projects: Array<Record<string, unknown>> }> };
    body.hits[0].projects[0].dataUseRestriction = "GRU";
    return body;
  };
  const result = ToolResultSchema.parse(await tools.find((tool) => tool.name === "hca_get_project")!.call(http, args.hca_get_project));
  assert.equal(result.entity_type, "page");
  if (result.entity_type === "page") {
    const item = result.items[0];
    assert.equal(item.entity_type, "dataset");
    if (item.entity_type === "dataset") assert.equal(item.url, null);
  }
});
