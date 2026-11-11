import { curie, type ToolResult } from "phi-omics-schemas";
import { dataset as datasetEntity, makePage, safeMap, sliceList } from "phi-omics-runtime";
import { bounds, dataDir, dataSource, num, rows, tablePath, text, unavailable } from "./common.js";
import { metricInfo } from "./metrics.js";

/** The datasets and metrics present in drug_response, with what each metric means. */
export async function listResponseMetrics(args: Record<string, unknown>): Promise<ToolResult> {
  const dir = dataDir();
  if (!dir) return unavailable("Set PHI_PGX_DATA_DIR to the processed pharmacogenomics parquet directory.");
  const responses = tablePath(dir, "drug_response", true);
  if (!responses) return unavailable("drug_response parquet is not in PHI_PGX_DATA_DIR.");
  const dataset = text(args.dataset);
  const { offset, limit } = bounds(args);
  const described = await rows(dir, "DESCRIBE SELECT * FROM read_parquet(?, hive_partitioning=true)", [responses]);
  const columns = new Set(described.map((row) => text(row.column_name)).filter((item): item is string => item != null));
  const classExpr = columns.has("response_class") ? "coalesce(response_class, 'single')" : "'single'";
  const found = await rows(
    dir,
    `SELECT dataset, metric, ${classExpr} AS response_class, any_value(unit) AS unit,
            count(*)::BIGINT AS n, count(DISTINCT model_cvcl)::BIGINT AS models, count(DISTINCT source_compound_id)::BIGINT AS compounds
     FROM read_parquet(?, hive_partitioning=true)
     WHERE (? IS NULL OR dataset = ?)
     GROUP BY 1, 2, 3
     ORDER BY 1, 2, 3`,
    [responses, dataset, dataset],
  );
  const { rows: pageRows, total } = sliceList(found, offset, limit);
  const items = safeMap(pageRows, (row) => {
    const datasetName = text(row.dataset) ?? "unknown";
    const metric = text(row.metric) ?? "value";
    const responseClass = text(row.response_class) ?? "single";
    const info = metricInfo(datasetName, metric);
    return datasetEntity({
      id: curie("pgx.metric", `${datasetName}/${metric}/${responseClass}`),
      title: `${datasetName} ${metric}`,
      assay: metric,
      sample_count: num(row.models),
      source_extra: {
        dataset: datasetName,
        metric,
        response_class: responseClass,
        unit: text(row.unit),
        rows: num(row.n),
        models: num(row.models),
        compounds: num(row.compounds),
        direction: info.direction,
        note: info.note,
      },
    });
  });
  return makePage("dataset", items, offset, limit, total, dataSource, {
    summary: "Every dataset and metric in the drug-response table. direction says whether a lower or higher value is a stronger response, and is unknown where the source does not document it. The same metric name can mean different things in different datasets, so values are only comparable within one dataset and metric.",
    url: "https://www.cancerrxgene.org/",
  });
}
