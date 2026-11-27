import assert from "node:assert/strict";
import { mkdirSync, symlinkSync, unlinkSync, writeFileSync } from "node:fs";
import { mkdtemp, rm } from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import test from "node:test";
import { DuckDBInstance } from "@duckdb/node-api";
import { ToolResultSchema, type Dataset, type OmicsFile, type Page, type Sample } from "phi-omics-schemas";
import type { JsonClient } from "phi-omics-runtime";
import { serverModule, tools } from "../server.js";

const client: JsonClient = {
  async getJson() {
    return {};
  },
  async postJson() {
    return {};
  },
};

async function callTool(name: string, args: Record<string, unknown> = {}) {
  const tool = tools.find((t) => t.name === name);
  assert.ok(tool, `missing tool ${name}`);
  const parsedArgs = tool.inputSchema.parse(args) as Record<string, unknown>;
  const result = await tool.call(client, parsedArgs);
  return ToolResultSchema.parse(result);
}

async function buildSnapshot(releaseDir: string, extraDataset = false): Promise<void> {
  mkdirSync(releaseDir, { recursive: true });
  const conn = await (await DuckDBInstance.create(":memory:")).connect();
  const copy = async (sql: string, outPath: string) => {
    mkdirSync(path.dirname(outPath), { recursive: true });
    await conn.run(`COPY (${sql}) TO '${outPath}' (FORMAT PARQUET)`);
  };

  const dsValues = [
    `('geo:GSE100001','Single-cell liver HCC atlas','single cell liver hcc atlas','Homo sapiens','Liver','["HCC","liver cancer"]','scRNA-seq','["GPL24676"]',24,'public','CC0','verified_live','h1','2026-11-27T10:00:00+0800','2026-11-27T10:00:00+0800','ok','2026-11-27T10:00:00+0800')`,
    `('gdc_tcga:TCGA-LIHC','TCGA Liver Hepatocellular Carcinoma','tcga liver hepatocellular carcinoma','Homo sapiens','Liver','["HCC"]','RNA-seq','["Illumina"]',420,'controlled','dbGaP','verified_live','h2','2026-11-27T10:00:00+0800','2026-11-27T10:00:00+0800','ok','2026-11-27T10:00:00+0800')`,
  ];
  if (extraDataset) {
    dsValues.push(
      `('pride:PXD000999','Spatial proteomics of cholangiocarcinoma','spatial proteomics of cholangiocarcinoma','Homo sapiens','Bile duct','["CCA"]','proteomics','["Orbitrap"]',36,'public','CC-BY-4.0','metadata_only','h3','2026-11-27T12:00:00+0800','2026-11-27T12:00:00+0800','ok','2026-11-27T12:00:00+0800')`,
    );
  }

  await copy(
    `SELECT * FROM (VALUES ${dsValues.join(",\n")})
     AS t(canonical_id, title, title_norm, organism, tissue, disease_tags, assay, platforms, n_samples, access, license, confidence, content_hash, first_seen_at, updated_at, probe_status, last_checked)`,
    path.join(releaseDir, "datasets.parquet"),
  );

  await copy(
    `SELECT * FROM (VALUES
      (1, 'geo:GSE100001', 'geo', 'GSE100001', 'series', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100001'),
      (2, 'geo:GSE100001', 'ncbi', 'PRJNA900001', 'bioproject', 'https://www.ncbi.nlm.nih.gov/bioproject/PRJNA900001'),
      (3, 'gdc_tcga:TCGA-LIHC', 'gdc_tcga', 'TCGA-LIHC', 'project', 'https://portal.gdc.cancer.gov/projects/TCGA-LIHC')
    ) AS t(id, dataset_id, source, accession, role, external_url)`,
    path.join(releaseDir, "dataset_members.parquet"),
  );

  await copy(
    `SELECT * FROM (VALUES
      ('38000001', '10.1038/s41586-026-00001-1', 'Single-cell atlas of HCC', 'Li X, Wang Y', 'Nature', 2026, 'Comprehensive scRNA-seq study of HCC.', 'https://pubmed.ncbi.nlm.nih.gov/38000001/', '2026-11-27T10:00:00+0800')
    ) AS t(pmid, doi, title, authors, journal, year, abstract, url, fetched_at)`,
    path.join(releaseDir, "publications.parquet"),
  );

  await copy(
    `SELECT * FROM (VALUES
      ('geo:GSE100001', '38000001', 'primary', 'geo_link')
    ) AS t(dataset_id, pmid, relation, source)`,
    path.join(releaseDir, "dataset_publications.parquet"),
  );

  await copy(
    `SELECT * FROM (VALUES
      ('geo:GSM1001', 'geo:GSE100001', 'HCC_Tumor_1', 'Homo sapiens', 'Liver', 'HCC', 'tumor', 'male', '58', '{"stage":"II"}', '["SRR1001"]', 'sh1', '2026-11-27T10:00:00+0800'),
      ('GSM1002', 'geo:GSE100001', 'HCC_Normal_1', 'Homo sapiens', 'Liver', 'HCC', 'adjacent_normal', 'male', '58', '{"stage":"II"}', '[]', 'sh2', '2026-11-27T10:00:00+0800')
    ) AS t(sample_id, dataset_id, title, organism, tissue, disease, condition, sex, age, attributes, run_ids, content_hash, updated_at)`,
    path.join(releaseDir, "samples", "source=geo", "data.parquet"),
  );

  await copy(
    `SELECT * FROM (VALUES
      ('geo:GSE100001:counts', 'geo:GSE100001', NULL::VARCHAR, 'processed', 'counts', 'Raw count matrix', 'GSE100001_counts.tsv.gz', 'tsv', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100001/suppl/GSE100001_counts.tsv.gz', 1048576, 'public', 'CC0', 'geo', NULL::VARCHAR, 'ok', '2026-11-27T10:00:00+0800', 'ah1'),
      ('gdc_tcga:TCGA-LIHC:bam1', 'gdc_tcga:TCGA-LIHC', NULL::VARCHAR, 'raw', 'bam', 'Controlled BAM', 'patient1.bam', 'bam', 'https://api.gdc.cancer.gov/data/secret-uuid', 99999999, 'controlled', 'dbGaP', 'gdc_tcga', NULL::VARCHAR, 'ok', '2026-11-27T10:00:00+0800', 'ah2')
    ) AS t(artifact_id, dataset_id, sample_id, kind, role, title, filename, format, uri, size_bytes, access, license, source, derived_from, probe_status, last_checked, content_hash)`,
    path.join(releaseDir, "dataset_artifacts", "source=geo", "data.parquet"),
  );

  writeFileSync(
    path.join(releaseDir, "_manifest.json"),
    JSON.stringify({
      release: path.basename(releaseDir),
      datasets_count: dsValues.length,
      publications_count: 1,
    }),
  );
}

test("omics-sniffer tools query Parquet snapshots and follow live symlink swaps", async () => {
  const tmp = await mkdtemp(path.join(os.tmpdir(), "phi-omics-sniffer-"));
  const prevEnv = process.env.PHI_OMICS_SNIFFER_DATA_DIR;
  try {
    const servingRoot = path.join(tmp, "serving");
    const rel1 = path.join(servingRoot, "releases", "v_20261127_100000");
    const rel2 = path.join(servingRoot, "releases", "v_20261127_160000");
    await buildSnapshot(rel1, false);
    await buildSnapshot(rel2, true);

    const currentLink = path.join(servingRoot, "current");
    symlinkSync(path.join("releases", "v_20261127_100000"), currentLink);
    process.env.PHI_OMICS_SNIFFER_DATA_DIR = servingRoot;

    assert.equal(await serverModule.ready?.(), null);

    // 0. Usage terms
    const terms = (await callTool("usage_terms")) as Page;
    assert.equal(terms.entity_type, "page");
    assert.match(terms.usage.summary, /v_20261127_100000/);

    // 1. Catalog stats on release 1
    const stats1 = (await callTool("get_catalog_stats")) as Page;
    assert.equal(stats1.page.total, 2);
    assert.equal(stats1.source_extra.snapshot_release, "v_20261127_100000");

    // 2. Search datasets with disease + min_samples filter
    const searchRes = (await callTool("search_omics_datasets", { disease: "HCC", min_samples: 20 })) as Page;
    assert.equal(searchRes.page.total, 2);
    const firstDs = searchRes.items[0] as Dataset;
    assert.equal(firstDs.id, "gdc.tcga:TCGA-LIHC");
    assert.equal(firstDs.access, "controlled");
    assert.equal(firstDs.url, null);

    // 3. Get dataset by bare accession ('GSE100001') and verify xrefs & publications
    const dsPage = (await callTool("get_dataset", { canonical_id: "GSE100001" })) as Page;
    assert.equal(dsPage.items.length, 1);
    const geoDs = dsPage.items[0] as Dataset;
    assert.equal(geoDs.id, "geo:GSE100001");
    assert.ok(geoDs.xrefs.includes("ncbi:PRJNA900001"));
    assert.ok(geoDs.xrefs.includes("pmid:38000001"));

    // 4. Get dataset by normalized CURIE ('gdc.tcga:TCGA-LIHC') round-trips to gdc_tcga:TCGA-LIHC
    const tcgaPage = (await callTool("get_dataset", { canonical_id: "gdc.tcga:TCGA-LIHC" })) as Page;
    assert.equal(tcgaPage.items.length, 1);

    // 5. List samples (including legacy unprefixed GSM1002 normalized to valid CURIE)
    const samplesPage = (await callTool("list_samples", { dataset_id: "GSE100001", condition: "tumor" })) as Page;
    assert.equal(samplesPage.page.total, 1);
    const tumorSample = samplesPage.items[0] as Sample;
    assert.equal(tumorSample.id, "geo:GSM1001");
    assert.ok(tumorSample.xrefs.includes("sra:SRR1001"));

    // 6. List artifacts (open vs controlled URL redaction)
    const openFiles = (await callTool("list_artifacts", { dataset_id: "geo:GSE100001" })) as Page;
    assert.equal(openFiles.items.length, 1);
    assert.match((openFiles.items[0] as OmicsFile).url ?? "", /^https:\/\/ftp\.ncbi\.nlm\.nih\.gov\//);

    const ctrlFiles = (await callTool("list_artifacts", { dataset_id: "gdc.tcga:TCGA-LIHC", source: "geo" })) as Page;
    assert.equal(ctrlFiles.items.length, 1);
    assert.equal((ctrlFiles.items[0] as OmicsFile).access, "controlled");
    assert.equal((ctrlFiles.items[0] as OmicsFile).url, null);

    // 7. Live atomic symlink swap -> release 2 becomes visible immediately without restart
    unlinkSync(currentLink);
    symlinkSync(path.join("releases", "v_20261127_160000"), currentLink);

    const stats2 = (await callTool("get_catalog_stats")) as Page;
    assert.equal(stats2.page.total, 3);
    assert.equal(stats2.source_extra.snapshot_release, "v_20261127_160000");
  } finally {
    if (prevEnv === undefined) delete process.env.PHI_OMICS_SNIFFER_DATA_DIR;
    else process.env.PHI_OMICS_SNIFFER_DATA_DIR = prevEnv;
    await rm(tmp, { recursive: true, force: true });
  }
});

test("missing omics-sniffer Parquet directory returns unavailable error on query tools", async () => {
  const prevEnv = process.env.PHI_OMICS_SNIFFER_DATA_DIR;
  const prevFallback = process.env.OMICS_SNIFFER_SERVING_DIR;
  delete process.env.PHI_OMICS_SNIFFER_DATA_DIR;
  delete process.env.OMICS_SNIFFER_SERVING_DIR;
  try {
    assert.match((await serverModule.ready?.()) ?? "", /PHI_OMICS_SNIFFER_DATA_DIR/);
    const terms = await callTool("usage_terms");
    assert.equal(terms.entity_type, "page");
    const err = await callTool("search_omics_datasets", { query: "HCC" });
    assert.equal(err.entity_type, "error");
    if (err.entity_type === "error") assert.equal(err.code, "unavailable");
  } finally {
    if (prevEnv !== undefined) process.env.PHI_OMICS_SNIFFER_DATA_DIR = prevEnv;
    if (prevFallback !== undefined) process.env.OMICS_SNIFFER_SERVING_DIR = prevFallback;
  }
});
