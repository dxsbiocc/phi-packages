"""What a drug-response metric means, by dataset and metric together.

This mirrors ``mcp-servers/pharmacogenomics/src/metrics.ts``. The same metric name does not
mean the same thing across datasets, so a direction is only asserted where the source
documents it and is ``unknown`` everywhere else. ``tests/test_standardize.py`` compares the
two implementations so they cannot drift apart silently.
"""

from __future__ import annotations

LOWER_IS_STRONGER = "lower_is_stronger"
HIGHER_IS_STRONGER = "higher_is_stronger"
UNKNOWN = "unknown"

# Potency metrics are concentrations, so a lower value is a stronger response.
POTENCY = frozenset({"LN_IC50", "IC50", "EC50", "AC50", "GI50", "TGI", "LC50"})
# Datasets whose own AUC is the area under the viability curve. ctrpv2 is the CTRPv2.1
# PharmacoSet, where auc_recomputed + aac_recomputed is 100 on all 395,263 profiles.
VIABILITY_AUC_DATASETS = frozenset({"gdsc1", "gdsc2", "prism_secondary_19q4", "ctrpv2"})


def direction(dataset: str | None, metric: str | None) -> str:
    """Whether a lower or a higher value of ``metric`` in ``dataset`` is a stronger response."""
    name = (metric or "").strip().upper()
    source = (dataset or "").strip().lower()
    if name in POTENCY or name in {"LFC", "VIABILITY"}:
        return LOWER_IS_STRONGER
    if name == "AAC":
        return HIGHER_IS_STRONGER
    if name == "AUC" and source in VIABILITY_AUC_DATASETS:
        return LOWER_IS_STRONGER
    return UNKNOWN
