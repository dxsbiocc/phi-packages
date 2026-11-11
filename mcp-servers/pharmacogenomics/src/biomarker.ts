import { curie, type ToolResult } from "phi-omics-schemas";
import { dataset as datasetEntity, failure, InputError, makePage, safeMap } from "phi-omics-runtime";
import {
  bounds,
  dataDir,
  dataSource,
  geneKey,
  num,
  rows,
  SILENT_CONSEQUENCES,
  tablePath,
  text,
  tissueFilter,
  unavailable,
  unknownTissue,
  type DuckDBValue,
} from "./common.js";
import { matchSingleCompounds } from "./query.js";
import { metricInfo, moreSensitive } from "./metrics.js";
import { mannWhitney, median, spearman } from "./stats.js";

const FEATURES = ["mutation", "expression", "copy_number", "dependency"] as const;
type Feature = (typeof FEATURES)[number];

const FEATURE_TABLE: Record<Feature, string> = {
  mutation: "mutation",
  expression: "expression",
  copy_number: "cnv",
  dependency: "gene_effect",
};

/** Smallest number of paired models for a rank correlation to be reported. */
const MIN_PAIRED_MODELS = 10;
const DEFAULT_MIN_GROUP = 3;
const SILENT_SQL = SILENT_CONSEQUENCES.map((name) => `'${name}'`).join(", ");

const CAVEATS = "Descriptive screening statistics. Replicate and repeated screens are averaged per model. p-values are uncorrected for the many compound, metric, and dataset combinations, and cell lines are not independent draws from patients. Tissue confounds many associations, so repeat with tissue set before reading one as a biomarker.";

interface ResponseGroup {
  dataset: string;
  metric: string;
  dose: number | null;
  byModel: Map<string, number>;
}

function groupResponses(found: Record<string, unknown>[]): ResponseGroup[] {
  const groups = new Map<string, ResponseGroup>();
  for (const row of found) {
    const model = text(row.model_cvcl);
    const value = num(row.value);
    if (!model || value == null) continue;
    const datasetName = text(row.dataset) ?? "unknown";
    const metric = text(row.metric) ?? "value";
    const dose = num(row.dose);
    const key = `${datasetName}\u0000${metric}\u0000${dose ?? ""}`;
    const group = groups.get(key) ?? { dataset: datasetName, metric, dose, byModel: new Map<string, number>() };
    group.byModel.set(model, value);
    groups.set(key, group);
  }
  return [...groups.values()];
}

interface Association {
  dataset: string;
  metric: string;
  dose: number | null;
  statistics: Record<string, unknown>;
  p: number;
}

interface Skipped {
  dataset: string;
  metric: string;
  dose: number | null;
  reason: string;
}

function mutationAssociation(
  group: ResponseGroup,
  mutated: ReadonlySet<string>,
  sequenced: ReadonlySet<string>,
  minGroup: number,
): Association | Skipped {
  const mutant: number[] = [];
  const wildType: number[] = [];
  for (const [model, value] of group.byModel) {
    if (!sequenced.has(model)) continue;
    (mutated.has(model) ? mutant : wildType).push(value);
  }
  const base = { dataset: group.dataset, metric: group.metric, dose: group.dose };
  if (mutant.length < minGroup || wildType.length < minGroup) {
    return { ...base, reason: `needs at least ${minGroup} models per group, found ${mutant.length} mutant and ${wildType.length} wild type` };
  }
  const test = mannWhitney(mutant, wildType)!;
  const medianMutant = median(mutant)!;
  const medianWild = median(wildType)!;
  const info = metricInfo(group.dataset, group.metric);
  const sensitive = medianMutant === medianWild ? null : moreSensitive(info.direction, medianMutant > medianWild);
  const sensitiveGroup = sensitive === "first" ? "mutant" : sensitive === "second" ? "wild_type" : null;
  return {
    ...base,
    p: test.p,
    statistics: {
      test: "mann_whitney_u",
      n_mutant: mutant.length,
      n_wild_type: wildType.length,
      median_mutant: medianMutant,
      median_wild_type: medianWild,
      median_difference: medianMutant - medianWild,
      rank_biserial: test.effect,
      p_value: test.p,
      direction: info.direction,
      more_sensitive_group: sensitiveGroup,
      reading: "A positive median_difference and rank_biserial mean larger response values in mutant models.",
      note: info.note,
    },
  };
}

function numericAssociation(group: ResponseGroup, feature: ReadonlyMap<string, number>): Association | Skipped {
  const x: number[] = [];
  const y: number[] = [];
  for (const [model, value] of group.byModel) {
    const level = feature.get(model);
    if (level == null) continue;
    x.push(level);
    y.push(value);
  }
  const base = { dataset: group.dataset, metric: group.metric, dose: group.dose };
  if (x.length < MIN_PAIRED_MODELS) return { ...base, reason: `needs at least ${MIN_PAIRED_MODELS} models with both measurements, found ${x.length}` };
  const result = spearman(x, y);
  if (!result) return { ...base, reason: "one of the two measurements is constant across the paired models" };
  const info = metricInfo(group.dataset, group.metric);
  // A rising response value means "first = higher feature" is the larger-value side.
  const sensitive = result.rho === 0 ? null : moreSensitive(info.direction, result.rho > 0);
  const higherFeatureMoreSensitive = sensitive === null ? null : sensitive === "first";
  return {
    ...base,
    p: result.p,
    statistics: {
      test: "spearman",
      n_models: result.n,
      rho: result.rho,
      p_value: result.p,
      direction: info.direction,
      higher_feature_more_sensitive: higherFeatureMoreSensitive,
      reading: "A positive rho means the response value rises with the feature.",
      note: info.note,
    },
  };
}

function isSkipped(value: Association | Skipped): value is Skipped {
  return "reason" in value;
}

async function responseRows(
  dir: string,
  singleIds: string[],
  args: { dataset: string | null; metric: string | null },
  tissueSql: string,
  tissueParams: DuckDBValue[],
): Promise<Record<string, unknown>[]> {
  const responses = tablePath(dir, "drug_response", true)!;
  const described = await rows(dir, "DESCRIBE SELECT * FROM read_parquet(?, hive_partitioning=true)", [responses]);
  const columns = new Set(described.map((row) => text(row.column_name)).filter((item): item is string => item != null));
  const singleOnly = columns.has("combination_id") ? "AND combination_id IS NULL" : "";
  const doseColumn = columns.has("dose") ? "dose" : "NULL::DOUBLE";
  return rows(
    dir,
    `SELECT dataset, metric, ${doseColumn} AS dose, model_cvcl, avg(value) AS value
     FROM read_parquet(?, hive_partitioning=true)
     WHERE list_contains(string_split(?, ','), source_compound_id) ${singleOnly}
       AND coalesce(qc, '') = '' AND model_cvcl IS NOT NULL AND value IS NOT NULL
       AND (? IS NULL OR dataset = ?) AND (? IS NULL OR metric = ?) ${tissueSql}
     GROUP BY 1, 2, 3, 4`,
    [responses, singleIds.join(","), args.dataset, args.dataset, args.metric, args.metric, ...tissueParams],
  );
}

async function featureLevels(
  dir: string,
  table: string,
  key: { id: string | null; symbol: string | null },
  featureDataset: string | null,
  tissueSql: string,
  tissueParams: DuckDBValue[],
): Promise<Map<string, number>> {
  const file = tablePath(dir, table, true)!;
  const found = await rows(
    dir,
    `SELECT model_cvcl, avg(value) AS value FROM read_parquet(?, hive_partitioning=true)
     WHERE ((? IS NOT NULL AND ncbi_gene_id = ?) OR (? IS NOT NULL AND lower(symbol) = lower(?)))
       AND model_cvcl IS NOT NULL AND value IS NOT NULL AND (? IS NULL OR dataset = ?) ${tissueSql}
     GROUP BY 1`,
    [file, key.id, key.id, key.symbol, key.symbol, featureDataset, featureDataset, ...tissueParams],
  );
  const levels = new Map<string, number>();
  for (const row of found) {
    const model = text(row.model_cvcl);
    const value = num(row.value);
    if (model && value != null) levels.set(model, value);
  }
  return levels;
}

async function mutationCalls(
  dir: string,
  key: { id: string | null; symbol: string | null },
  featureDataset: string | null,
  tissueSql: string,
  tissueParams: DuckDBValue[],
): Promise<{ mutated: Set<string>; sequenced: Set<string> }> {
  const file = tablePath(dir, "mutation", true)!;
  const sequenced = await rows(
    dir,
    `SELECT DISTINCT model_cvcl FROM read_parquet(?, hive_partitioning=true)
     WHERE model_cvcl IS NOT NULL AND (? IS NULL OR dataset = ?) ${tissueSql}`,
    [file, featureDataset, featureDataset, ...tissueParams],
  );
  const mutated = await rows(
    dir,
    `SELECT DISTINCT model_cvcl FROM read_parquet(?, hive_partitioning=true)
     WHERE ((? IS NOT NULL AND ncbi_gene_id = ?) OR (? IS NOT NULL AND lower(symbol) = lower(?)))
       AND model_cvcl IS NOT NULL AND lower(coalesce(consequence, '')) NOT IN (${SILENT_SQL})
       AND (? IS NULL OR dataset = ?) ${tissueSql}`,
    [file, key.id, key.id, key.symbol, key.symbol, featureDataset, featureDataset, ...tissueParams],
  );
  const toSet = (list: Record<string, unknown>[]) => new Set(list.map((row) => text(row.model_cvcl)).filter((item): item is string => item != null));
  return { mutated: toSet(mutated), sequenced: toSet(sequenced) };
}

export async function biomarkerAssociation(args: Record<string, unknown>): Promise<ToolResult> {
  const dir = dataDir();
  if (!dir) return unavailable("Set PHI_PGX_DATA_DIR to the processed pharmacogenomics parquet directory.");
  const query = text(args.compound);
  if (!query) throw new InputError("compound is required");
  const geneArg = text(args.gene);
  if (!geneArg) throw new InputError("gene is required");
  const feature = text(args.feature) as Feature | null;
  if (!feature || !FEATURES.includes(feature)) throw new InputError(`feature must be one of ${FEATURES.join(", ")}`);
  const minGroup = args.min_group == null ? DEFAULT_MIN_GROUP : Number(args.min_group);
  if (!Number.isInteger(minGroup) || minGroup < 2 || minGroup > 100) throw new InputError("min_group must be an integer from 2 to 100");
  const { offset, limit } = bounds(args);

  const compounds = tablePath(dir, "compound_xref", false);
  if (!compounds || !tablePath(dir, "drug_response", true)) return unavailable("compound_xref or drug_response parquet is not in PHI_PGX_DATA_DIR.");
  if (!tablePath(dir, FEATURE_TABLE[feature], true)) return unavailable(`${FEATURE_TABLE[feature]} parquet is not in PHI_PGX_DATA_DIR.`);

  const tissue = text(args.tissue);
  const tissueClause = await tissueFilter(dir, tissue);
  if (tissue && tissueClause.unknown) return unknownTissue(tissue);

  const matched = await matchSingleCompounds(dir, compounds, query);
  const singleIds = [...new Set(matched.map((row) => text(row.source_compound_id)).filter((item): item is string => item != null))];
  if (!singleIds.length) return failure(dataSource, "not_found", `No single-agent compound matched ${query}`, false);
  const cids = [...new Set(matched.map((row) => text(row.pubchem_cid)).filter((item): item is string => item != null))];

  const key = geneKey(geneArg);
  const featureDataset = text(args.feature_dataset);
  const responses = groupResponses(await responseRows(
    dir,
    singleIds,
    { dataset: text(args.dataset), metric: text(args.metric) },
    tissueClause.sql,
    tissueClause.params,
  ));

  const associations: Association[] = [];
  const skipped: Skipped[] = [];
  if (feature === "mutation") {
    const calls = await mutationCalls(dir, key, featureDataset, tissueClause.sql, tissueClause.params);
    for (const group of responses) {
      const result = mutationAssociation(group, calls.mutated, calls.sequenced, minGroup);
      if (isSkipped(result)) skipped.push(result);
      else associations.push(result);
    }
  } else {
    const levels = await featureLevels(dir, FEATURE_TABLE[feature], key, featureDataset, tissueClause.sql, tissueClause.params);
    for (const group of responses) {
      const result = numericAssociation(group, levels);
      if (isSkipped(result)) skipped.push(result);
      else associations.push(result);
    }
  }
  associations.sort((left, right) => left.p - right.p || left.dataset.localeCompare(right.dataset) || left.metric.localeCompare(right.metric));

  const page = associations.slice(offset, offset + limit);
  const items = safeMap(page, (association) => datasetEntity({
    id: curie("pgx.association", [association.dataset, association.metric, association.dose == null ? "na" : String(association.dose), feature, geneArg].join("/").replace(/\s+/g, "_")),
    title: `${association.dataset} ${association.metric}${association.dose == null ? "" : ` at dose ${association.dose}`}: ${geneArg} ${feature}`,
    assay: association.metric,
    sample_count: Number(association.statistics.n_models ?? (Number(association.statistics.n_mutant) + Number(association.statistics.n_wild_type))),
    source_extra: {
      drug_dataset: association.dataset,
      metric: association.metric,
      dose: association.dose,
      feature,
      gene: geneArg,
      tissue,
      ...association.statistics,
    },
  }));
  return makePage("dataset", items, offset, limit, associations.length, dataSource, {
    summary: `Association of single-agent response with ${feature === "mutation" ? "non-silent mutation status of" : `${feature.replace("_", " ")} of`} ${geneArg} across cell-line models, one row per response dataset and metric, most significant first. ${CAVEATS}`,
    url: "https://www.cancerrxgene.org/",
  }, {
    pubchem_cids: cids,
    ambiguous_compound: cids.length > 1,
    skipped_groups: skipped.slice(0, 50),
    skipped_group_count: skipped.length,
  });
}
