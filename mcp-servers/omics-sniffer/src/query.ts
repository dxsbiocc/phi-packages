import { existsSync } from "node:fs";
import path from "node:path";
import type { ToolResult } from "phi-omics-schemas";
import { dataset, failure, file, InputError, makePage, safeMap, sample, usagePage } from "phi-omics-runtime";
import {
  bounds,
  dataSource,
  intOrNull,
  mapTaxon,
  normalizeAccess,
  parseJsonSafe,
  partitionFiles,
  queryRows,
  readJsonIfExists,
  resolveServingPaths,
  str,
  toCurie,
  truncateText,
  unavailable,
  type DuckDBValue,
} from "./common.js";

export const USAGE_SUMMARY = [
  "Omics Sniffer indexes public multi-omics dataset, sample, file-artifact, and publication metadata across 30+ repositories (GEO, SRA, ArrayExpress, PRIDE, MassIVE, MetaboLights, Metabolomics Workbench, CELLxGENE, ENCODE, GTEx, GDC/TCGA, cBioPortal, CPTAC, PDC, DepMap, HuBMAP, HPA, ChIP-Atlas, GWAS Catalog, FinnGen, and others).",
  "Queries run read-only over ZSTD-compressed, source-partitioned Parquet snapshots updated incrementally on the server.",
  "Metadata and open repository URLs are returned; controlled-access artifacts expose metadata only without download URLs. Follow each upstream repository's license and citation policy.",
].join(" ");

function defaultDatasetUrl(canonicalId: string): string | null {
  const idx = canonicalId.indexOf(":");
  if (idx <= 0) return null;
  const src = canonicalId.slice(0, idx).toLowerCase();
  const acc = canonicalId.slice(idx + 1);
  if (src === "geo" && acc.startsWith("GSE")) return `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${acc}`;
  if ((src === "ae" || src === "arrayexpress") && acc.startsWith("E-")) return `https://www.ebi.ac.uk/biostudies/arrayexpress/studies/${acc}`;
  if (src === "pride" && acc.startsWith("PXD")) return `https://www.ebi.ac.uk/pride/archive/projects/${acc}`;
  if (src === "metabolights" && acc.startsWith("MTBLS")) return `https://www.ebi.ac.uk/metabolights/${acc}`;
  return null;
}

async function resolveCanonicalDatasetId(snapshotDir: string, inputId: string): Promise<string | null> {
  const trimmed = inputId.trim();
  if (!trimmed) return null;
  const dsPq = path.join(snapshotDir, "datasets.parquet");
  const memPq = path.join(snapshotDir, "dataset_members.parquet");

  // 1. Exact canonical_id match or dot-to-underscore prefix match (e.g. gdc.tcga:X -> gdc_tcga:X)
  const candidates = [trimmed];
  const colon = trimmed.indexOf(":");
  if (colon > 0) {
    const underscoreVariant = `${trimmed.slice(0, colon).replace(/\./g, "_")}:${trimmed.slice(colon + 1)}`;
    if (!candidates.includes(underscoreVariant)) candidates.push(underscoreVariant);
  }
  for (const cand of candidates) {
    const exact = await queryRows(
      "SELECT canonical_id FROM read_parquet(?) WHERE lower(canonical_id) = lower(?) LIMIT 1",
      [dsPq, cand],
    );
    const hit = str(exact[0]?.canonical_id);
    if (hit) return hit;
  }

  // 2. Member accession lookup (e.g. bare GSE164760, E-GEOD-100001, PXD000123, PRJNA...)
  if (existsSync(memPq)) {
    const bareAcc = colon > 0 ? trimmed.slice(colon + 1) : trimmed;
    const byAcc = await queryRows(
      "SELECT dataset_id FROM read_parquet(?) WHERE lower(accession) = lower(?) ORDER BY id ASC LIMIT 1",
      [memPq, bareAcc],
    );
    const hit = str(byAcc[0]?.dataset_id);
    if (hit) return hit;
  }
  return null;
}

export async function usageTerms(): Promise<ToolResult> {
  const paths = resolveServingPaths();
  let summary = USAGE_SUMMARY;
  if (!paths) {
    summary += " PHI_OMICS_SNIFFER_DATA_DIR is not set to an existing Parquet serving directory, so query tools are currently unavailable.";
  } else {
    const manifest = readJsonIfExists(path.join(paths.snapshotDir, "_manifest.json"));
    const rel = str(manifest?.release) ?? path.basename(paths.snapshotDir);
    summary += ` Active Parquet snapshot release: ${rel}.`;
  }
  return usagePage(dataSource, summary, "https://github.com/dxsbiocc/phi-packages");
}

export async function searchOmicsDatasets(args: Record<string, unknown>): Promise<ToolResult> {
  const paths = resolveServingPaths();
  if (!paths) return unavailable("Set PHI_OMICS_SNIFFER_DATA_DIR to the omics-sniffer Parquet serving directory.");

  const { offset, limit } = bounds(args, 50, 20);
  const dsPq = path.join(paths.snapshotDir, "datasets.parquet");
  const memPq = path.join(paths.snapshotDir, "dataset_members.parquet");

  const where: string[] = ["1=1"];
  const filterParams: DuckDBValue[] = [];

  const query = str(args.query);
  const organismArg = str(args.organism);
  const assayArg = str(args.assay);
  const diseaseArg = str(args.disease);
  const tissueArg = str(args.tissue);
  const sourceArg = str(args.source);
  const minSamplesArg = intOrNull(args.min_samples);
  const accessArg = str(args.access);

  if (query) {
    const kw = `%${query.toLowerCase()}%`;
    where.push("(lower(coalesce(d.title, '')) LIKE ? OR lower(d.canonical_id) LIKE ?)");
    filterParams.push(kw, kw);
  }
  if (organismArg) {
    where.push("lower(coalesce(d.organism, '')) LIKE ?");
    filterParams.push(`%${organismArg.toLowerCase()}%`);
  }
  if (assayArg) {
    where.push("lower(coalesce(d.assay, '')) LIKE ?");
    filterParams.push(`%${assayArg.toLowerCase()}%`);
  }
  if (diseaseArg) {
    const kw = `%${diseaseArg.toLowerCase()}%`;
    where.push("(lower(coalesce(d.disease_tags, '')) LIKE ? OR lower(coalesce(d.title, '')) LIKE ?)");
    filterParams.push(kw, kw);
  }
  if (tissueArg) {
    where.push("lower(coalesce(d.tissue, '')) LIKE ?");
    filterParams.push(`%${tissueArg.toLowerCase()}%`);
  }
  if (minSamplesArg != null && minSamplesArg > 0) {
    where.push("coalesce(d.n_samples, 0) >= ?");
    filterParams.push(minSamplesArg);
  }
  if (accessArg) {
    const reqAccess = accessArg.toLowerCase();
    if (reqAccess === "open" || reqAccess === "public") {
      where.push("lower(coalesce(d.access, '')) IN ('public', 'open')");
    } else {
      where.push("lower(coalesce(d.access, '')) = ?");
      filterParams.push(reqAccess);
    }
  }

  let sourceJoin = "";
  const joinParams: DuckDBValue[] = [];
  if (sourceArg && existsSync(memPq)) {
    sourceJoin = "JOIN (SELECT DISTINCT dataset_id FROM read_parquet(?) WHERE lower(source) = lower(?)) m ON m.dataset_id = d.canonical_id";
    joinParams.push(memPq, sourceArg);
  }

  const whereSql = where.join(" AND ");
  const countRows = await queryRows(
    `SELECT count(*)::INTEGER AS n FROM read_parquet(?) d ${sourceJoin} WHERE ${whereSql}`,
    [dsPq, ...joinParams, ...filterParams],
  );
  const total = intOrNull(countRows[0]?.n) ?? 0;

  const found = await queryRows(
    `SELECT d.canonical_id, d.title, d.organism, d.tissue, d.disease_tags,
            d.assay, d.platforms, d.n_samples, d.access, d.license, d.confidence, d.updated_at
     FROM read_parquet(?) d
     ${sourceJoin}
     WHERE ${whereSql}
     ORDER BY coalesce(d.n_samples, 0) DESC, d.canonical_id ASC
     LIMIT ? OFFSET ?`,
    [dsPq, ...joinParams, ...filterParams, limit, offset],
  );

  const items = safeMap(found, (row) => {
    const cid = str(row.canonical_id) ?? "omics:unknown";
    const id = toCurie(cid, "omics.dataset");
    const organism = str(row.organism);
    const access = normalizeAccess(row.access);
    return dataset({
      id,
      title: str(row.title),
      assay: str(row.assay),
      taxon_id: mapTaxon(organism),
      sample_count: intOrNull(row.n_samples),
      access,
      url: defaultDatasetUrl(cid),
      source_extra: {
        canonical_id: cid,
        organism,
        tissue: str(row.tissue),
        disease_tags: parseJsonSafe(str(row.disease_tags)),
        platforms: parseJsonSafe(str(row.platforms)),
        license: str(row.license),
        confidence: str(row.confidence),
        updated_at: str(row.updated_at),
      },
    });
  });

  return makePage("dataset", items, offset, limit, total, dataSource, {
    summary: USAGE_SUMMARY,
    url: "https://github.com/dxsbiocc/phi-packages",
  });
}

export async function getDataset(args: Record<string, unknown>): Promise<ToolResult> {
  const paths = resolveServingPaths();
  if (!paths) return unavailable("Set PHI_OMICS_SNIFFER_DATA_DIR to the omics-sniffer Parquet serving directory.");

  const rawQuery = str(args.canonical_id);
  if (!rawQuery) throw new InputError("canonical_id is required");

  const resolvedCid = await resolveCanonicalDatasetId(paths.snapshotDir, rawQuery);
  if (!resolvedCid) {
    return failure(dataSource, "not_found", `No dataset matched '${rawQuery}' in omics-sniffer.`, false);
  }

  const dsPq = path.join(paths.snapshotDir, "datasets.parquet");
  const memPq = path.join(paths.snapshotDir, "dataset_members.parquet");
  const dpPq = path.join(paths.snapshotDir, "dataset_publications.parquet");
  const pubPq = path.join(paths.snapshotDir, "publications.parquet");

  const dsRows = await queryRows("SELECT * FROM read_parquet(?) WHERE canonical_id = ? LIMIT 1", [dsPq, resolvedCid]);
  const row = dsRows[0];
  if (!row) return failure(dataSource, "not_found", `Dataset '${resolvedCid}' not found.`, false);

  const memberRows = existsSync(memPq)
    ? await queryRows(
        "SELECT source, accession, role, external_url FROM read_parquet(?) WHERE dataset_id = ? ORDER BY id ASC",
        [memPq, resolvedCid],
      )
    : [];

  const pubRows = existsSync(dpPq) && existsSync(pubPq)
    ? await queryRows(
        `SELECT p.pmid, p.doi, p.title, p.authors, p.journal, p.year, p.abstract, p.url,
                dp.relation, dp.source AS link_source
         FROM read_parquet(?) dp
         JOIN read_parquet(?) p ON p.pmid = dp.pmid
         WHERE dp.dataset_id = ?`,
        [dpPq, pubPq, resolvedCid],
      )
    : [];

  const id = toCurie(resolvedCid, "omics.dataset");
  const xrefs: string[] = [];
  let primaryUrl: string | null = defaultDatasetUrl(resolvedCid);

  const members = memberRows.map((m) => {
    const mSource = str(m.source) ?? "omics";
    const mAcc = str(m.accession) ?? "unknown";
    const mCurie = toCurie(`${mSource}:${mAcc}`, "omics");
    if (mCurie !== id && !xrefs.includes(mCurie)) xrefs.push(mCurie);
    const extUrl = str(m.external_url);
    if (!primaryUrl && extUrl) primaryUrl = extUrl;
    return {
      source: mSource,
      accession: mAcc,
      role: str(m.role),
      external_url: extUrl,
    };
  });

  const publications = pubRows.map((p) => {
    const pmid = str(p.pmid);
    if (pmid && /^\d+$/.test(pmid)) {
      const pCurie = `pmid:${pmid}`;
      if (!xrefs.includes(pCurie)) xrefs.push(pCurie);
    }
    return {
      pmid,
      doi: str(p.doi),
      title: str(p.title),
      authors: str(p.authors),
      journal: str(p.journal),
      year: intOrNull(p.year),
      abstract: truncateText(str(p.abstract), 1000),
      url: str(p.url),
      relation: str(p.relation),
      link_source: str(p.link_source),
    };
  });

  const organism = str(row.organism);
  const access = normalizeAccess(row.access);
  const item = dataset({
    id,
    xrefs,
    title: str(row.title),
    assay: str(row.assay),
    taxon_id: mapTaxon(organism),
    sample_count: intOrNull(row.n_samples),
    access,
    url: primaryUrl,
    source_extra: {
      canonical_id: resolvedCid,
      organism,
      tissue: str(row.tissue),
      disease_tags: parseJsonSafe(str(row.disease_tags)),
      platforms: parseJsonSafe(str(row.platforms)),
      license: str(row.license),
      confidence: str(row.confidence),
      probe_status: str(row.probe_status),
      first_seen_at: str(row.first_seen_at),
      updated_at: str(row.updated_at),
      members,
      publications,
    },
  });

  return makePage("dataset", [item], 0, 1, 1, dataSource, {
    summary: USAGE_SUMMARY,
    url: primaryUrl ?? "https://github.com/dxsbiocc/phi-packages",
  });
}

export async function listSamples(args: Record<string, unknown>): Promise<ToolResult> {
  const paths = resolveServingPaths();
  if (!paths) return unavailable("Set PHI_OMICS_SNIFFER_DATA_DIR to the omics-sniffer Parquet serving directory.");

  const rawDs = str(args.dataset_id);
  if (!rawDs) throw new InputError("dataset_id is required");
  const { offset, limit } = bounds(args, 50, 20);
  const tissueArg = str(args.tissue);
  const diseaseArg = str(args.disease);
  const conditionArg = str(args.condition);
  const sourceArg = str(args.source);

  const cid = (await resolveCanonicalDatasetId(paths.snapshotDir, rawDs)) ?? rawDs;
  const inferredPrefix = cid.includes(":") ? cid.slice(0, cid.indexOf(":")) : null;
  const globPath = partitionFiles(paths.snapshotDir, "samples", sourceArg || inferredPrefix);
  if (!globPath) {
    return makePage("sample", [], offset, limit, 0, dataSource, { summary: USAGE_SUMMARY, url: null }, { dataset_id: cid });
  }

  const where = ["dataset_id = ?"];
  const params: DuckDBValue[] = [cid];
  if (tissueArg) {
    where.push("lower(coalesce(tissue, '')) LIKE ?");
    params.push(`%${tissueArg.toLowerCase()}%`);
  }
  if (diseaseArg) {
    where.push("lower(coalesce(disease, '')) LIKE ?");
    params.push(`%${diseaseArg.toLowerCase()}%`);
  }
  if (conditionArg) {
    where.push("lower(coalesce(condition, '')) LIKE ?");
    params.push(`%${conditionArg.toLowerCase()}%`);
  }

  const whereSql = where.join(" AND ");
  const totalRows = await queryRows(
    `SELECT count(*)::INTEGER AS n FROM read_parquet(?, hive_partitioning=true) WHERE ${whereSql}`,
    [globPath, ...params],
  );
  const total = intOrNull(totalRows[0]?.n) ?? 0;

  const found = await queryRows(
    `SELECT sample_id, dataset_id, title, organism, tissue, disease, condition, sex, age, attributes, run_ids, updated_at
     FROM read_parquet(?, hive_partitioning=true)
     WHERE ${whereSql}
     ORDER BY sample_id ASC
     LIMIT ? OFFSET ?`,
    [globPath, ...params, limit, offset],
  );

  const dsCurie = toCurie(cid, "omics.dataset");
  const items = safeMap(found, (row) => {
    const sid = str(row.sample_id) ?? "unknown";
    const organism = str(row.organism);
    const runsParsed = parseJsonSafe(str(row.run_ids));
    const runXrefs = Array.isArray(runsParsed)
      ? runsParsed
          .map((r) => str(r))
          .filter((r): r is string => Boolean(r))
          .map((r) => toCurie(`sra:${r}`, "sra"))
      : [];
    return sample({
      id: toCurie(sid, inferredPrefix ?? "geo"),
      xrefs: [dsCurie, ...runXrefs],
      name: str(row.title) ?? sid,
      taxon_id: mapTaxon(organism),
      source_extra: {
        sample_id: sid,
        dataset_id: cid,
        organism,
        tissue: str(row.tissue),
        disease: str(row.disease),
        condition: str(row.condition),
        sex: str(row.sex),
        age: str(row.age),
        attributes: parseJsonSafe(str(row.attributes)),
        run_ids: runsParsed,
        updated_at: str(row.updated_at),
      },
    });
  });

  return makePage("sample", items, offset, limit, total, dataSource, {
    summary: USAGE_SUMMARY,
    url: defaultDatasetUrl(cid),
  }, { dataset_id: cid });
}

export async function listArtifacts(args: Record<string, unknown>): Promise<ToolResult> {
  const paths = resolveServingPaths();
  if (!paths) return unavailable("Set PHI_OMICS_SNIFFER_DATA_DIR to the omics-sniffer Parquet serving directory.");

  const rawDs = str(args.dataset_id);
  if (!rawDs) throw new InputError("dataset_id is required");
  const { offset, limit } = bounds(args, 50, 20);
  const kindArg = str(args.kind);
  const roleArg = str(args.role);
  const formatArg = str(args.format);
  const sourceArg = str(args.source);

  const cid = (await resolveCanonicalDatasetId(paths.snapshotDir, rawDs)) ?? rawDs;
  const inferredPrefix = cid.includes(":") ? cid.slice(0, cid.indexOf(":")) : null;
  const globPath = partitionFiles(paths.snapshotDir, "dataset_artifacts", sourceArg || inferredPrefix);
  if (!globPath) {
    return makePage("file", [], offset, limit, 0, dataSource, { summary: USAGE_SUMMARY, url: null }, { dataset_id: cid });
  }

  const where = ["dataset_id = ?"];
  const params: DuckDBValue[] = [cid];
  if (kindArg) {
    where.push("lower(coalesce(kind, '')) = lower(?)");
    params.push(kindArg);
  }
  if (roleArg) {
    where.push("lower(coalesce(role, '')) LIKE ?");
    params.push(`%${roleArg.toLowerCase()}%`);
  }
  if (formatArg) {
    where.push("lower(coalesce(format, '')) LIKE ?");
    params.push(`%${formatArg.toLowerCase()}%`);
  }

  const whereSql = where.join(" AND ");
  const totalRows = await queryRows(
    `SELECT count(*)::INTEGER AS n FROM read_parquet(?, hive_partitioning=true) WHERE ${whereSql}`,
    [globPath, ...params],
  );
  const total = intOrNull(totalRows[0]?.n) ?? 0;

  const found = await queryRows(
    `SELECT artifact_id, dataset_id, sample_id, kind, role, title, filename, format, uri, size_bytes, access, license, source
     FROM read_parquet(?, hive_partitioning=true)
     WHERE ${whereSql}
     ORDER BY kind ASC NULLS LAST, filename ASC NULLS LAST, artifact_id ASC
     LIMIT ? OFFSET ?`,
    [globPath, ...params, limit, offset],
  );

  const dsCurie = toCurie(cid, "omics.dataset");
  const items = safeMap(found, (row) => {
    const aid = str(row.artifact_id) ?? "unknown";
    const access = normalizeAccess(row.access);
    const sampleId = str(row.sample_id);
    const xrefs = sampleId ? [toCurie(sampleId, inferredPrefix ?? "geo")] : [];
    return file({
      id: toCurie(aid, "omics.artifact"),
      xrefs,
      name: str(row.filename) ?? str(row.title) ?? aid,
      format: str(row.format),
      data_type: str(row.role) ?? str(row.kind),
      size_bytes: intOrNull(row.size_bytes),
      access,
      url: access === "controlled" ? null : str(row.uri),
      dataset_id: dsCurie,
      source_extra: {
        artifact_id: aid,
        dataset_id: cid,
        sample_id: sampleId,
        kind: str(row.kind),
        role: str(row.role),
        title: str(row.title),
        license: str(row.license),
        source: str(row.source),
      },
    });
  });

  return makePage("file", items, offset, limit, total, dataSource, {
    summary: USAGE_SUMMARY,
    url: defaultDatasetUrl(cid),
  }, { dataset_id: cid });
}

export async function getCatalogStats(): Promise<ToolResult> {
  const paths = resolveServingPaths();
  if (!paths) {
    return usagePage(
      dataSource,
      `${USAGE_SUMMARY} PHI_OMICS_SNIFFER_DATA_DIR is not set to an existing Parquet serving directory, so query tools are currently unavailable.`,
      "https://github.com/dxsbiocc/phi-packages",
    );
  }

  const manifest = readJsonIfExists(path.join(paths.snapshotDir, "_manifest.json"));
  const syncState = readJsonIfExists(path.join(paths.configuredRoot, ".sync_state.json"));
  const dsPq = path.join(paths.snapshotDir, "datasets.parquet");
  const memPq = path.join(paths.snapshotDir, "dataset_members.parquet");

  const dsCountRows = await queryRows("SELECT count(*)::INTEGER AS n FROM read_parquet(?)", [dsPq]);
  const datasetsCount = intOrNull(dsCountRows[0]?.n) ?? 0;

  const bySourceRows = existsSync(memPq)
    ? await queryRows(
        `SELECT source, count(DISTINCT dataset_id)::INTEGER AS datasets
         FROM read_parquet(?)
         GROUP BY source
         ORDER BY datasets DESC, source ASC
         LIMIT 50`,
        [memPq],
      )
    : [];

  const bySource = bySourceRows.map((r) => ({
    source: str(r.source),
    datasets: intOrNull(r.datasets) ?? 0,
  }));

  return makePage(
    "dataset",
    [],
    0,
    1,
    datasetsCount,
    dataSource,
    {
      summary: `${USAGE_SUMMARY} Active snapshot: ${str(manifest?.release) ?? path.basename(paths.snapshotDir)} (${datasetsCount} datasets).`,
      url: "https://github.com/dxsbiocc/phi-packages",
    },
    {
      snapshot_release: str(manifest?.release) ?? path.basename(paths.snapshotDir),
      datasets_count: datasetsCount,
      serving_manifest: manifest,
      sync_state: syncState,
      datasets_by_source: bySource,
    },
  );
}
