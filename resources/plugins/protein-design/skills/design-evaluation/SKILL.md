---
name: design-evaluation
description: >-
  Evaluate and diagnose protein design campaign outputs: compute funnel pass
  rates, metric quantiles (pLDDT, ipTM, min_pae_interaction, Kabsch CA RMSD,
  interface contacts, steric clashes, Boltz-2 ligand affinity), and per-parent
  backbone diversity with `pd_funnel_report`; interpret thresholds and decide
  next-round iterative refinement parameters using the metrics reference guide.
  Use after running `local/modules/design-filter`,
  `local/subworkflows/design-postprocess`, or
  `local/workflows/protein-binder-design`.
phi:
  environment: phi:python@1
  attachTo:
    - ProteinDesign
  scripts:
    - name: funnel_report
      description: >-
        Summarize a protein design campaign from DESIGN_FILTER outputs
        (`merged.tsv`, optional `filtered.tsv` and `summary.json`): reports
        total/passing/selected counts, pass rate, numeric quantiles (min,
        median, p90, max) for key metrics, unique parent backbone counts, and
        top ranked candidates.
      run: [python3, ./scripts/funnel_report.py]
      args:
        type: object
        additionalProperties: false
        required: [merged_tsv]
        properties:
          merged_tsv:
            type: string
            format: input-path
            description: Path to merged.tsv produced by local/modules/design-filter.
          filtered_tsv:
            type: string
            format: input-path
            description: Optional path to filtered.tsv produced by local/modules/design-filter.
          summary_json:
            type: string
            format: input-path
            description: Optional path to summary.json produced by local/modules/design-filter.
          top_k:
            type: integer
            description: Number of top candidates to include in the JSON preview (default 10).
      approval: read
      output: ./schemas/funnel-report-result.json
      timeoutSeconds: 60
---

# Design Evaluation & Funnel Diagnostics

Use this skill after running `local/modules/design-filter`,
`local/subworkflows/design-postprocess`, or `local/workflows/protein-binder-design`
to diagnose campaign quality and decide whether to scale, iterate, or package
final candidates.

## Reference guide

- [references/metrics-and-thresholds.md](references/metrics-and-thresholds.md) — Interpretation and starting thresholds for `complex_plddt`, `iptm`, `min_pae_interaction`, `ca_rmsd`, `interface_contacts`, `clash_count`, and Boltz-2 ligand affinity (`affinity_pred_value`, `affinity_probability_binary`).

## Standard evaluation flow

1. Call `pd_funnel_report` on `merged.tsv` (plus `filtered.tsv` and `summary.json` when available).
2. Read the funnel counts (`total_designs`, `passed_filter`, `selected_designs`, `unique_parents_selected`) and per-metric quantiles.
3. Diagnose bottlenecks before spending more GPU compute:
   - Low `passed_filter` with high `ca_rmsd` → tighten MPNN temperature or sample more sequences per backbone.
   - Good `ca_rmsd` but low `iptm` / high `min_pae_interaction` → revisit binder length, hotspots, or target interface definition.
   - Passing designs clustered on only 1–2 parent backbones → increase `n_chunks` / `n_batches` and enforce `max_per_group: 1` or `2`.
