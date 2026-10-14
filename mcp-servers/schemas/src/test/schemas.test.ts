import assert from "node:assert/strict";
import test from "node:test";
import {
  FileSchema,
  PageSchema,
  ToolResultSchema,
  curie,
  ensemblGeneCurie,
  humanChrom,
  ontologyCurie,
  type EntityType,
} from "../index.js";

const source = {
  id: "gtex",
  name: "GTEx",
  version: "v10",
  homepage: "https://gtexportal.org/",
  retrieved_at: "2026-10-14T00:00:00.000Z",
};
const usage = { summary: "GTEx terms apply.", url: "https://gtexportal.org/" };

function page(itemType: EntityType, items: unknown[]) {
  return {
    entity_type: "page",
    item_type: itemType,
    items,
    page: { offset: 0, limit: 20, returned: items.length, total: items.length, next_offset: null },
    source,
    usage,
    source_extra: {},
  };
}

test("CURIE helpers normalize ontology, Ensembl, and chromosomes", () => {
  assert.equal(ontologyCurie("UBERON:0002190"), "UBERON:0002190");
  assert.equal(ontologyCurie("CL_0000636"), "CL:0000636");
  assert.equal(ensemblGeneCurie("ENSG00000141510.18").id, "ensembl:ENSG00000141510");
  assert.equal(ensemblGeneCurie("ensembl:ENSG00000141510").version, null);
  assert.equal(humanChrom("17"), "chr17");
  assert.equal(humanChrom("chrM"), "chrM");
  assert.equal(humanChrom("MT"), "chrM");
  assert.equal(curie("gtex.tissue", "Whole_Blood"), "gtex.tissue:Whole_Blood");
});

test("every entity type fits the shared page envelope", () => {
  const cases: Array<[EntityType, Record<string, unknown>]> = [
    ["gene", { entity_type: "gene", id: "ensembl:ENSG00000141510", xrefs: ["ncbigene:7157"], symbol: "TP53", name: "tumor protein p53", taxon_id: "NCBITaxon:9606", biotype: "protein_coding", synonyms: ["p53"], locus: { genome_build: "GRCh38", chrom: "chr17", start: 7661779, end: 7687538, strand: "-" }, source_extra: {} }],
    ["sample", { entity_type: "sample", id: "hubmap:HBM794.WZTX.678", xrefs: [], name: null, donor_id: "hubmap:HBM624.KWHJ.622", tissue_id: null, cell_type_id: null, disease_ids: [], assay: null, taxon_id: "NCBITaxon:9606", source_extra: {} }],
    ["donor", { entity_type: "donor", id: "gdc.case:6e3b6b72-142d-4b8d-a462-28a205796e41", xrefs: [], sex: null, age_years: 60, taxon_id: "NCBITaxon:9606", disease_ids: [], source_extra: {} }],
    ["dataset", { entity_type: "dataset", id: "pride:PXD001357", xrefs: [], title: "milk", description: null, assay: null, taxon_id: null, tissue_ids: [], disease_ids: [], cell_type_ids: [], sample_count: null, access: "open", url: "https://www.ebi.ac.uk/pride/archive/projects/PXD001357", source_extra: {} }],
    ["file", { entity_type: "file", id: "gdc.file:af0d3d92-f8b5-4612-805b-d6773595262b", xrefs: [], name: "example.maf.gz", format: null, data_type: "Annotated Somatic Mutation", size_bytes: 35992, access: "controlled", url: null, dataset_id: null, source_extra: {} }],
    ["cell_type", { entity_type: "cell_type", id: "CL:0000636", xrefs: [], name: "Mueller cell", description: null, ontology_id: "CL:0000636", parent_ids: [], source_extra: {} }],
    ["tissue", { entity_type: "tissue", id: "gtex.tissue:Adipose_Subcutaneous", xrefs: [], name: "Adipose - Subcutaneous", ontology_id: "UBERON:0002190", source_extra: {} }],
    ["protein", { entity_type: "protein", id: "uniprot:P04637", xrefs: [], name: null, gene_id: "ensembl:ENSG00000141510", gene_symbol: "TP53", taxon_id: "NCBITaxon:9606", source_extra: {} }],
    ["metabolite", { entity_type: "metabolite", id: "CHEBI:15377", xrefs: [], name: "water", formula: "H2O", inchikey: null, source_extra: {} }],
    ["expression", { entity_type: "expression", id: "gtex.expression:ENSG00000141510:Adipose_Subcutaneous:median", xrefs: [], gene_id: "ensembl:ENSG00000141510", tissue_id: "gtex.tissue:Adipose_Subcutaneous", cell_type_id: null, sample_id: null, value: 22.6549, unit: "TPM", statistic: "median", source_extra: {} }],
    ["variant", { entity_type: "variant", id: "gtex.variant:chr1_64764_C_T_b38", xrefs: ["dbsnp:rs769952832"], gene_id: "ensembl:ENSG00000227232", locus: { genome_build: "GRCh38", chrom: "chr1", start: 64764, end: 64764, strand: null }, consequence: null, source_extra: {} }],
  ];
  for (const [itemType, item] of cases) {
    const parsed = ToolResultSchema.parse(page(itemType, [item]));
    assert.equal(parsed.entity_type, "page");
  }
});

test("pages reject mixed entity types and controlled file URLs", () => {
  const gene = { entity_type: "gene", id: "ensembl:ENSG00000141510", xrefs: [], symbol: "TP53", name: null, taxon_id: null, biotype: null, synonyms: [], locus: null, source_extra: {} };
  assert.equal(PageSchema.safeParse(page("tissue", [gene])).success, false);
  const leaked = FileSchema.safeParse({
    entity_type: "file",
    id: "gdc.file:af0d3d92-f8b5-4612-805b-d6773595262b",
    xrefs: [],
    name: "secret.bam",
    format: "BAM",
    data_type: null,
    size_bytes: 10,
    access: "controlled",
    url: "https://api.gdc.cancer.gov/data/af0d3d92-f8b5-4612-805b-d6773595262b",
    dataset_id: null,
    source_extra: {},
  });
  assert.equal(leaked.success, false);
});

test("errors use the shared error object", () => {
  const parsed = ToolResultSchema.parse({
    entity_type: "error",
    code: "not_found",
    message: "missing",
    source,
    retryable: false,
    source_extra: {},
  });
  assert.equal(parsed.entity_type, "error");
});
