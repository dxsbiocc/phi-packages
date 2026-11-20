# MAG workflow references

Decision loop: [analysis-playbook.md](analysis-playbook.md). Per-tool depth:

- [step-notes/assembly-stats.md](step-notes/assembly-stats.md) — reading QUAST for metagenomes
- [step-notes/bin-qc.md](step-notes/bin-qc.md) — CheckM vs CheckM2, MIMAG tiers
- [step-notes/bin-filter.md](step-notes/bin-filter.md) — threshold choices and pass rates
- [step-notes/mag-table.md](step-notes/mag-table.md) — joins, GTDB gates, dRep clusters

Upstream tools (megahit, spades, metabat2, maxbin2, concoct, dastool,
checkm2, gtdbtk, drep, coverm, quast, and the whole nf-core/mag workflow)
are wrappers: discover with `wrapper_search`, read `wrapper.yaml` before
planning a run.

For maintainers: [maintainers/source-map.md](maintainers/source-map.md).