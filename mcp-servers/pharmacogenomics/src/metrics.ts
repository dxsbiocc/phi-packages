/**
 * What a drug-response metric means, by dataset and metric together.
 *
 * The same metric name does not mean the same thing across datasets. Sanger GDSC and PRISM
 * publish an AUC that falls as a cell line gets more sensitive, while an AUC recomputed by
 * PharmacoGx for an ORCESTRA PharmacoSet may be the area of the response (1 - viability), which
 * rises with sensitivity. Direction is therefore only asserted where the source documents it,
 * and is "unknown" everywhere else rather than guessed.
 */

export type Direction = "lower_is_stronger" | "higher_is_stronger" | "unknown";

export interface MetricInfo {
  direction: Direction;
  note: string | null;
}

/** Potency metrics are concentrations, so a lower value is a stronger response. */
const POTENCY = new Set(["LN_IC50", "IC50", "EC50", "AC50", "GI50", "TGI", "LC50"]);
/**
 * Datasets whose own AUC is the area under the viability curve. ctrpv2 is the CTRPv2.1 PharmacoSet,
 * where auc_recomputed + aac_recomputed is 100 on all 395,263 profiles (percent scale).
 */
const VIABILITY_AUC_DATASETS = new Set(["gdsc1", "gdsc2", "prism_secondary_19q4", "ctrpv2"]);

const AUC_UNKNOWN_NOTE = "AUC direction is not asserted for this dataset: a PharmacoSet AUC can be a viability area (lower is stronger) or a response area (higher is stronger). Check the dataset before comparing it with another one.";

export function metricInfo(dataset: string | null, metric: string | null): MetricInfo {
  const name = (metric ?? "").trim().toUpperCase();
  const source = (dataset ?? "").trim().toLowerCase();
  if (POTENCY.has(name)) return { direction: "lower_is_stronger", note: null };
  if (name === "LFC" || name === "VIABILITY") return { direction: "lower_is_stronger", note: null };
  if (name === "AAC") return { direction: "higher_is_stronger", note: "AAC is the area above the viability curve. Its scale depends on the PharmacoSet build: 0 to 1 in the 2015 CTRPv2 and FIMM PSets, 0 to 100 (percent) in CTRPv2.1." };
  if (name === "AUC") {
    return VIABILITY_AUC_DATASETS.has(source)
      ? { direction: "lower_is_stronger", note: "AUC is the area under the viability curve." }
      : { direction: "unknown", note: AUC_UNKNOWN_NOTE };
  }
  if (name.startsWith("SYNERGY_DELTA_")) return { direction: "unknown", note: "A combination screen delta, not a single-agent response." };
  return { direction: "unknown", note: null };
}

/** Which of two groups is more sensitive, given which has the larger value. */
export function moreSensitive(direction: Direction, firstIsLarger: boolean): "first" | "second" | null {
  if (direction === "unknown") return null;
  const largerIsStronger = direction === "higher_is_stronger";
  return firstIsLarger === largerIsStronger ? "first" : "second";
}
