# Profiling workflow references

Decision loop: [analysis-playbook.md](analysis-playbook.md) — read it before
chaining tools. Per-tool depth:

- [step-notes/inspect.md](step-notes/inspect.md) — file-kind detection
- [step-notes/merge.md](step-notes/merge.md) — count vs relative estimators, rank filtering
- [step-notes/rarefaction.md](step-notes/rarefaction.md) — depth adequacy, plateau reading
- [step-notes/diversity.md](step-notes/diversity.md) — index assumptions, PERMANOVA caveats
- [step-notes/differential.md](step-notes/differential.md) — compositionality, thresholds
- [step-notes/function.md](step-notes/function.md) — HUMAnN stratification, regrouping
- [step-notes/plot.md](step-notes/plot.md) — reading composition figures, clustering

Upstream tools (fastp, kraken2, bracken, metaphlan, humann3, …) are wrappers,
not plugin scripts: discover them with `wrapper_search` and read their
`wrapper.yaml` before planning a run.

For maintainers: [maintainers/source-map.md](maintainers/source-map.md).