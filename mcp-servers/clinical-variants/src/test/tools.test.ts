import assert from "node:assert/strict";
import test from "node:test";
import { ToolResultSchema } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { tools } from "../server.js";

function client(): JsonClient {
  return {
    async getJson(url: string) {
      if (url.includes("esearch.fcgi") && url.includes("clinvar")) {
        return { esearchresult: { count: "1", idlist: ["12351"], querytranslation: "TP53[gene]" } };
      }
      if (url.includes("esummary.fcgi") && url.includes("clinvar")) {
        return {
          result: {
            uids: ["12351"],
            "12351": {
              uid: "12351",
              accession: "VCV000012351",
              title: "NM_000546.6(TP53):c.742C>T (p.Arg248Trp)",
              obj_type: "single nucleotide variant",
              genes: [{ symbol: "TP53", geneid: "7157" }],
              germline_classification: {
                description: "Pathogenic",
                review_status: "reviewed by expert panel",
                last_evaluated: "2024-01-10",
                trait_set: [{ trait_name: "Li-Fraumeni syndrome" }],
              },
              variation_set: [{
                canonical_spdi: "NC_000017.11:7674219:G:A",
                variation_loc: [{ assembly_name: "GRCh38", chr: "17", start: 7674220, stop: 7674220 }],
                variation_xrefs: [{ db_source: "dbSNP", db_id: "121912651" }],
              }],
            },
          },
        };
      }
      throw new Error(`unmocked GET ${url}`);
    },
    async postJson(url: string, body: unknown) {
      const payload = body as { query?: string; variables?: Record<string, unknown> };
      if (url.includes("gnomad.broadinstitute.org/api")) {
        if (payload.query?.includes("variant(")) {
          return {
            data: {
              variant: {
                variant_id: "17-7674220-C-T",
                chrom: "17",
                pos: 7674220,
                ref: "C",
                alt: "T",
                rsids: ["rs121912651"],
                flags: [],
                exome: { ac: 1, an: 1400000, af: 7.1e-7, homozygote_count: 0, filters: [], populations: [] },
                genome: null,
                transcript_consequences: [{
                  gene_id: "ENSG00000141510",
                  gene_symbol: "TP53",
                  transcript_id: "ENST00000269305",
                  major_consequence: "missense_variant",
                  hgvsc: "c.742C>T",
                  hgvsp: "p.Arg248Trp",
                  lof: null,
                }],
              },
            },
          };
        }
        return {
          data: {
            gene: {
              gene_id: "ENSG00000141510",
              symbol: "TP53",
              name: "tumor protein p53",
              chrom: "17",
              start: 7661779,
              stop: 7687538,
              strand: "-",
              canonical_transcript_id: "ENST00000269305",
              omim_id: "191170",
              gnomad_constraint: { pLI: 0.99, oe_lof: 0.12, oe_lof_upper: 0.25 },
            },
          },
        };
      }
      if (url.includes("civicdb.org/api/graphql")) {
        if (payload.query?.includes("browseVariants(")) {
          return {
            data: {
              browseVariants: {
                totalCount: 1,
                nodes: [{
                  id: 12,
                  name: "V600E",
                  featureName: "BRAF",
                  link: "/variants/12",
                  featureId: 5,
                  category: "GENE",
                  aliases: [{ name: "RS113488022" }, { name: "V640E" }],
                  variantTypes: [{ id: 47, name: "Missense Variant" }],
                  evidenceItemCount: 256,
                }],
              },
            },
          };
        }
        return {
          data: {
            evidenceItems: {
              totalCount: 1,
              nodes: [{
                id: 88,
                name: "EID88",
                status: "ACCEPTED",
                evidenceType: "PREDICTIVE",
                evidenceLevel: "A",
                evidenceDirection: "SUPPORTS",
                evidenceRating: 5,
                significance: "SENSITIVITYRESPONSE",
                description: "BRAF V600E melanoma sensitivity to vemurafenib.",
                molecularProfile: { id: 12, name: "BRAF V600E" },
                disease: { name: "Melanoma", doid: "1909" },
                therapies: [{ name: "Vemurafenib" }],
                source: { citationId: "21639808", sourceType: "PUBMED" },
              }],
            },
          },
        };
      }
      throw new Error(`unmocked POST ${url}`);
    },
  };
}

const args: Record<string, Record<string, unknown>> = {
  usage_terms: {},
  clinvar_search_variants: { gene: "TP53", limit: 5 },
  clinvar_get_variant: { variation_id: "VCV000012351" },
  gnomad_get_variant: { variant_id: "chr17-7674220-C-T" },
  gnomad_get_gene_constraint: { gene: "TP53" },
  civic_search_variants: { gene: "BRAF", variant: "V600E", limit: 5 },
  civic_search_evidence: { molecular_profile: "BRAF V600E", limit: 5 },
};

test("every clinical-variants tool returns the shared schema", async () => {
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

test("invalid inputs and missing gnomAD variants return error objects", async () => {
  const http = client();
  const emptyClinvar = ToolResultSchema.parse(await tools.find((t) => t.name === "clinvar_search_variants")!.call(http, {}));
  assert.equal(emptyClinvar.entity_type, "error");
  if (emptyClinvar.entity_type === "error") assert.equal(emptyClinvar.code, "invalid_input");

  const badGnomad = ToolResultSchema.parse(await tools.find((t) => t.name === "gnomad_get_variant")!.call(http, { variant_id: "not-a-variant" }));
  assert.equal(badGnomad.entity_type, "error");
  if (badGnomad.entity_type === "error") assert.equal(badGnomad.code, "invalid_input");

  http.postJson = async () => ({ errors: [{ message: "Variant not found" }], data: { variant: null } });
  const missingGnomad = ToolResultSchema.parse(await tools.find((t) => t.name === "gnomad_get_variant")!.call(http, { variant_id: "1-1-A-G" }));
  assert.equal(missingGnomad.entity_type, "error");
  if (missingGnomad.entity_type === "error") assert.equal(missingGnomad.code, "not_found");
});
