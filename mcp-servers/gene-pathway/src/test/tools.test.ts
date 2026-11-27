import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

function client(): JsonClient {
  return {
    async getJson(url: string) {
      if (url.includes("rest.ensembl.org/lookup/")) {
        return {
          id: "ENSG00000141510",
          display_name: "TP53",
          description: "tumor protein p53 [Source:HGNC Symbol;Acc:HGNC:11998]",
          species: "homo_sapiens",
          biotype: "protein_coding",
          assembly_name: "GRCh38",
          seq_region_name: "17",
          start: 7661779,
          end: 7687538,
          strand: -1,
          version: 18,
          canonical_transcript: "ENST00000269305.9",
          Transcript: [{ id: "ENST00000269305", biotype: "protein_coding", is_canonical: 1, start: 7661779, end: 7687538 }],
        };
      }
      if (url.includes("rest.ensembl.org/homology/")) {
        return {
          data: [{
            homologies: [{
              type: "ortholog_one2one",
              taxonomy_level: "Euarchontoglires",
              target: { id: "ENSMUSG00000059552", species: "mus_musculus", protein_id: "ENSMUSP00000104298", perc_id: 77.2 },
            }],
          }],
        };
      }
      if (url.includes("rest.ensembl.org/vep/")) {
        return [{
          id: "ENST00000269305.9:c.818G>A",
          input: "ENST00000269305.9:c.818G>A",
          assembly_name: "GRCh38",
          seq_region_name: "17",
          start: 7673802,
          end: 7673802,
          strand: 1,
          allele_string: "C/T",
          variant_class: "SNV",
          most_severe_consequence: "missense_variant",
          colocated_variants: [{ id: "rs28934576" }],
          transcript_consequences: [{
            gene_id: "ENSG00000141510",
            gene_symbol: "TP53",
            transcript_id: "ENST00000269305",
            biotype: "protein_coding",
            impact: "MODERATE",
            canonical: 1,
            consequence_terms: ["missense_variant"],
            sift_prediction: "deleterious",
            polyphen_prediction: "probably_damaging",
          }],
        }];
      }
      if (url.includes("reactome.org/ContentService/search/query")) {
        return {
          results: [{
            entries: [{
              stId: "R-HSA-3700989",
              name: "Transcriptional Regulation by <span class=\"highlighting\">TP53</span>",
              summation: "<span class=\"highlighting\">TP53</span> pathway.",
              species: ["Homo sapiens"],
              compartmentNames: ["nucleoplasm"],
              isDisease: false,
            }],
          }],
        };
      }
      if (url.includes("reactome.org/ContentService/data/query/enhanced/R-HSA-3700989")) {
        return {
          stId: "R-HSA-3700989",
          displayName: "Transcriptional Regulation by TP53",
          schemaClass: "Pathway",
          species: [{ displayName: "Homo sapiens" }],
          summation: [{ text: "Regulation of p53 target genes." }],
          literatureReference: [{ pubMedIdentifier: 24722188 }],
          hasEvent: [{ stId: "R-HSA-5633007", displayName: "Regulation of TP53 Activity", schemaClass: "Pathway" }],
          goBiologicalProcess: { displayName: "DNA damage response" },
          hasDiagram: true,
        };
      }
      if (url.includes("reactome.org/ContentService/data/mapping/UniProt/P04637/pathways")) {
        return [{
          stId: "R-HSA-3700989",
          displayName: "Transcriptional Regulation by TP53",
          speciesName: "Homo sapiens",
          hasDiagram: true,
          isInDisease: false,
        }];
      }
      throw new Error(`unmocked ${url}`);
    },
    async postJson() { return {}; },
  };
}

const args: Record<string, Record<string, unknown>> = {
  usage_terms: {},
  ensembl_get_gene: { symbol: "TP53" },
  ensembl_get_homologs: { gene_id: "ENSG00000141510", target_species: "mus_musculus", limit: 5 },
  ensembl_vep_hgvs: { variant: "ENST00000269305.9:c.818G>A" },
  reactome_search_pathways: { query: "TP53", limit: 5 },
  reactome_get_pathway: { pathway_id: "R-HSA-3700989" },
  reactome_protein_pathways: { uniprot_accession: "P04637", limit: 5 },
};

test("every gene-pathway tool returns the shared schema", async () => {
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

test("Reactome search strips highlight tags and Ensembl gene strips version", async () => {
  const http = client();
  const geneRes = ToolResultSchema.parse(await tools.find((t) => t.name === "ensembl_get_gene")!.call(http, { symbol: "TP53" }));
  assert.equal(geneRes.entity_type, "page");
  if (geneRes.entity_type === "page") assert.equal(geneRes.items[0]?.id, "ensembl:ENSG00000141510");

  const pathRes = ToolResultSchema.parse(await tools.find((t) => t.name === "reactome_search_pathways")!.call(http, { query: "TP53" }));
  assert.equal(pathRes.entity_type, "page");
  if (pathRes.entity_type === "page" && pathRes.items[0]?.entity_type === "dataset") {
    assert.equal(pathRes.items[0].title, "Transcriptional Regulation by TP53");
  }
});
