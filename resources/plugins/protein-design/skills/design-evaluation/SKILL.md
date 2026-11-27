---
name: design-evaluation
description: >-
  Evaluate and diagnose protein design campaign outputs: compute funnel pass
  rates, metric quantiles (pLDDT, ipTM, min_pae_interaction, Kabsch CA RMSD,
  interface contacts, steric clashes, Boltz-2 ligand affinity), and per-parent
  backbone diversity with `pd_funnel_report`; run pre-wetlab biophysical and
  chemical-liability triage (`MW`, `pI` vs buffer pH, `eps280` / `A280 0.1%`,
  GRAVY, Guruprasad instability index, aliphatic index, 9-mer hydrophobic
  patch, unpaired Cys, deamidation/isomerization/N-glycosylation motifs) with
  `pd_developability_audit`; and interpret thresholds using the metrics
  reference guide. Use after running `local/modules/design-filter`,
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
    - name: developability_audit
      description: >-
        Run a pure-stdlib pre-wetlab biophysical and chemical-liability audit
        on candidate protein sequences from a FASTA file or a DESIGN_FILTER
        table (`filtered.tsv` / `merged.tsv`), writing `developability_audit.tsv`
        and its `.phi-artifact.json` sidecar into `output_dir`.
      run: [python3, ./scripts/developability_audit.py]
      args:
        type: object
        additionalProperties: false
        required: [sequences_input, output_dir]
        properties:
          sequences_input:
            type: string
            format: input-path
            description: Protein FASTA (.fa/.fasta) or TSV/CSV table containing designed sequences.
          output_dir:
            type: string
            format: project-path
            description: Directory for developability_audit.tsv and its .phi-artifact.json sidecar.
          seq_col:
            type: string
            description: Sequence column name when sequences_input is a TSV/CSV (auto-detected if omitted).
          id_col:
            type: string
            description: Candidate ID column name when sequences_input is a TSV/CSV (auto-detected if omitted).
          buffer_ph:
            type: number
            description: Assay/purification buffer pH for net charge and pI-proximity warning (default 7.4).
          max_hydrophobic_patch:
            type: number
            description: Maximum acceptable 9-aa sliding-window Kyte-Doolittle average before flagging (default 2.2).
      approval: write
      output: ./schemas/developability-audit-result.json
      timeoutSeconds: 120
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
4. Before ordering DNA or gene synthesis, run `pd_developability_audit` on `filtered.tsv` (or a candidate FASTA) to catch `pI` buffer clashes, UV280-blind sequences, unpaired cysteines, hydrophobic patches, and chemical liability motifs.
