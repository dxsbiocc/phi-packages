# Hypothesis & Provenance Ledger Schema (`research_ledger.json` & `PROVENANCE.md`)

Load this reference when initializing, updating, or auditing `research_ledger.json` and `PROVENANCE.md`.

## 1. `research_ledger.json` Canonical Template

```json
{
  "schema_version": "1.0.0",
  "investigation_id": "inv-2026-egfr-tki-persistence",
  "title": "Axl-driven persister survival under osimertinib in EGFR-mutant LUAD",
  "updated_at": "2026-11-27T03:00:00Z",
  "stopping_criteria": {
    "converge_rule": "At least 2 orthogonal modalities (omics + perturbation/wet-lab) support H1 with FDR <= 0.05, H0 batch/cell-cycle artifacts falsified, 0 orphan claims, and 0 unresolved major contradictions.",
    "falsify_rule": "Perturbation rescue/knockdown fails to shift IC50 >= 2-fold (or p > 0.05 across >= 3 biological donors/lines) under matched drug exposure.",
    "current_decision": "CONTINUE",
    "decision_rationale": "Single-cell and qPCR support H1, while proteomics shows lagged AXL induction; pending functional rescue readout."
  },
  "hypotheses": [
    {
      "id": "H0",
      "statement": "Apparent AXL upregulation in osimertinib-treated pools is a cell-cycle arrest / composition artifact rather than cell-intrinsic persister activation.",
      "prior_rationale": "Drug treatment arrests G1/S cells and shifts epithelial vs mesenchymal proportions in bulk lysates.",
      "falsification_threshold": "AXL remains upregulated (log2FC >= 1.0, FDR <= 0.05) within matched G0/G1 single-cell clusters AND in purified epithelial clones.",
      "status": "refuted",
      "evidence_updates": [
        {
          "step": 1,
          "lineage_ref": "LIN-001",
          "direction": "against",
          "delta_log_bf": -2.1,
          "note": "scRNA-seq pseudobulk restricted to G0/G1 tumor cells still shows AXL induction (log2FC = 1.84, FDR = 0.002)."
        }
      ]
    },
    {
      "id": "H1",
      "statement": "Osimertinib induces cell-intrinsic AXL transcriptional activation that sustains persister survival and raises osimertinib IC50 >= 2-fold.",
      "prior_rationale": "Literature reports AXL bypass signaling in EGFR-TKI tolerance (Zhang 2012, PMID:22751098).",
      "falsification_threshold": "Supported if >= 2 modalities show AXL induction (FC >= 2.0, FDR <= 0.05) AND AXL inhibition/KO restores osimertinib sensitivity (IC50 fold-shift >= 2.0, p < 0.05).",
      "status": "inconclusive",
      "evidence_updates": [
        {
          "step": 1,
          "lineage_ref": "LIN-001",
          "direction": "for",
          "delta_log_bf": 1.6,
          "note": "scRNA-seq persister cluster shows AXL upregulation."
        },
        {
          "step": 2,
          "lineage_ref": "LIN-002",
          "direction": "for",
          "delta_log_bf": 1.4,
          "note": "RT-qPCR (wa_qpcr_quant) confirms 3.4x AXL mRNA rise across n=4 independent passages."
        }
      ]
    }
  ],
  "lineage": [
    {
      "claim_id": "LIN-001",
      "manuscript_claim_id": "Result-2.1-Fig1C",
      "input_path": "data/raw/osimertinib_scrna.h5ad",
      "input_sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
      "tool_or_wrapper": "sc_de_markers",
      "tool_version": "single-cell@1.2.0",
      "exact_parameters": {
        "groupby": "treatment_state",
        "condition_test": "persister_day9",
        "condition_ref": "dmso_g0g1",
        "method": "pseudobulk_deseq2",
        "donor_col": "donor_id"
      },
      "artifact_path": "results/sc/persister_vs_dmso_de.csv",
      "sidecar_path": "results/sc/persister_vs_dmso_de.csv.phi-artifact.json",
      "effect_size_ci_p": "log2FC = 1.84, SE = 0.31, FDR = 0.0021",
      "verified_at": "2026-11-26T18:10:00Z"
    },
    {
      "claim_id": "LIN-002",
      "manuscript_claim_id": "Result-2.2-Fig2A",
      "input_path": "data/wetlab/axl_qpcr_day9.csv",
      "input_sha256": "9f86d081884c7d659a2feaa0c55ad015a3bf4f1b2b0b822cd15d6c15b0f00a08",
      "tool_or_wrapper": "wa_qpcr_quant",
      "tool_version": "wetlab-assay@1.0.0",
      "exact_parameters": {
        "method": "pfaffl",
        "ref_genes": ["GAPDH", "ACTB"],
        "biological_unit": "passage_id",
        "alpha": 0.05
      },
      "artifact_path": "results/wetlab/axl_qpcr_summary.json",
      "sidecar_path": "results/wetlab/axl_qpcr_summary.json.phi-artifact.json",
      "effect_size_ci_p": "FC = 3.42, 95% CI [2.38, 4.91], Welch p = 0.0048, n = 4 passages",
      "verified_at": "2026-11-26T21:45:00Z"
    }
  ],
  "contradictions": [
    {
      "contradiction_id": "CON-001",
      "claim_axis": "AXL abundance at 24 h post-osimertinib",
      "source_a": "LIN-002 (wa_qpcr_quant, 24 h mRNA FC = 2.9, p = 0.009)",
      "source_b": "LIN-003 (proteomics DIA, 24 h AXL protein log2FC = 0.18, FDR = 0.42)",
      "discordance_type": "cross-modality (mRNA vs total protein)",
      "adjudication_axes": {
        "sample_size_power": "qPCR n=4 passages; DIA n=3 replicates (underpowered for <1.5x protein change)",
        "batch_confounding": "Matched passage blocks in both assays",
        "assay_modality": "Transcriptional surge precedes steady-state receptor accumulation; surface/phospho-AXL not resolved in whole-cell DIA",
        "biological_unit": "Independent passage (consistent)",
        "bias_risk": "Low in both; timing mismatch"
      },
      "verdict": "reconciled_by_context",
      "resolution_note": "72 h and Day-9 timepoints + phospho-AXL ELISA (wa_immunoassay_curve) scheduled to verify delayed protein/activity accumulation."
    }
  ]
}
```

## 2. `PROVENANCE.md` Human-Audit Layout

Render a concise Markdown mirror alongside `research_ledger.json` whenever the ledger is updated:

1. **Header & Stopping-Gate Banner**: Investigation ID, timestamp, current decision (`CONVERGE` / `CONTINUE` / `STOP`), and blocking items.
2. **Hypothesis Status Board**:

| ID | Mechanism / Claim | Falsification Threshold | Status | Net Evidence (`delta_log_bf`) |
| --- | --- | --- | --- | --- |
| `H0` | Cell-cycle/composition artifact | Remains DE (`log2FC >= 1`, `FDR <= 0.05`) in G0/G1-matched clones | `refuted` | `-2.1` |
| `H1` | Cell-intrinsic AXL persister resistance | `>= 2` modalities + AXL-i rescue `>= 2x` IC50 shift (`p < 0.05`) | `inconclusive` | `+3.0` (awaiting rescue) |

3. **Data-to-Claim Lineage Table**:

| Claim ID | Manuscript Target | Input SHA-256 (prefix) | Tool & Key Params | Output Artifact | Verified Effect (`CI`, `P/FDR`) |
| --- | --- | --- | --- | --- | --- |
| `LIN-001` | `Result-2.1-Fig1C` | `e3b0c442` | `sc_de_markers` (`pseudobulk_deseq2`, `donor_id`) | `results/sc/persister_vs_dmso_de.csv` | `log2FC = 1.84`, `SE = 0.31`, `FDR = 0.0021` |
| `LIN-002` | `Result-2.2-Fig2A` | `9f86d081` | `wa_qpcr_quant` (`pfaffl`, `GAPDH/ACTB`, `passage_id`) | `results/wetlab/axl_qpcr_summary.json` | `FC = 3.42`, `95% CI [2.38, 4.91]`, `p = 0.0048` |

4. **Cross-Literature & Cross-Modality Contradiction Matrix**:

| ID | Axis | Source A vs Source B | Adjudication (Power / Modality / Bias) | Verdict | Next Action |
| --- | --- | --- | --- | --- | --- |
| `CON-001` | 24 h AXL level | mRNA `FC=2.9` vs DIA protein `log2FC=0.18` | Kinetic translation lag + whole-lysate DIA dilution (`n=3`) | `reconciled_by_context` | Test 72 h `wa_immunoassay_curve` + `pgx_hts_qc` rescue |

## 3. Audit Checklist Before Handing Off to `manuscript-writing`

- [ ] Every numeric value in Results/Abstract maps to a `LIN-*` entry in `research_ledger.json`.
- [ ] All `input_sha256` values match current files on disk (zero `[STALE LINEAGE]` flags).
- [ ] `H0` null/artifact explanations have documented falsification entries, not implicit dismissal.
- [ ] Every known conflicting study (`PMID`/`DOI`) or discordant modality appears in `contradictions[]` with an explicit `verdict`.
- [ ] `stopping_criteria.current_decision == "CONVERGE"` before writing definitive causal conclusions.
