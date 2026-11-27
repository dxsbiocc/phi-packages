---
name: assay-workflow
description: >-
  Quantitative wet-lab assay and biochemical kinetics workflow: MIQE-compliant
  RT-qPCR relative quantification (Pfaffl, geNorm M and Andersen reference-gene
  stability, NTC/SD QC, biological-replicate statistics); weighted 4PL/5PL
  immunoassay calibration (AICc selection, blank LOD, FDA/EMA LLOQ/ULOQ,
  range-censored unknown interpolation); steady-state Michaelis-Menten and Hill
  enzyme kinetics with progress-curve v0 extraction, AICc reversible inhibitor
  mechanism selection, and Cheng-Prusoff Ki; and global 1:1 Langmuir SPR/BLI
  binding kinetics (kon, koff, KD, t1/2, Rmax) with steady-state KD_eq
  cross-validation. Use for qPCR, RT-qPCR, Pfaffl, geNorm, ELISA, 4PL, 5PL,
  LLOQ, Michaelis-Menten, enzyme kinetics, SPR, BLI, Langmuir, 定量PCR, 酶动力学,
  表面等离子共振.
phi:
  environment: plugin:wa
  attachTo:
    - WetlabAssay
  scripts:
    - name: qpcr_quant
      description: >-
        MIQE-compliant RT-qPCR relative quantification (2^-ddCq and Pfaffl
        efficiency-corrected), geNorm M and Andersen reference-gene stability
        ranking with geometric-mean normalization, NTC and technical replicate SD
        QC guards, and biological-replicate statistical testing with BH FDR.
      run: [python, ./scripts/wa_qpcr_quant.py]
      args:
        type: object
        additionalProperties: false
        required: [cq_table, reference_genes, control_group, output_dir]
        properties:
          cq_table:
            type: string
            format: input-path
            description: Long-format qPCR Cq TSV/CSV with sample, group, gene, and cq columns (technical replicate rows supported; Undetermined/NaN allowed for NTCs).
          reference_genes:
            type: array
            items:
              type: string
            description: Reference (housekeeping) gene names used for geometric-mean normalization and stability ranking.
          control_group:
            type: string
            description: Calibrator / control condition name for delta-delta Cq and group comparisons.
          output_dir:
            type: string
            format: project-path
            description: Directory for sample quantities, group comparisons, reference stability tables, and summary figure.
          standard_curve:
            type: string
            format: input-path
            description: Optional TSV/CSV with gene, dilution_or_copy, and cq columns to estimate gene-specific slope, R^2, and amplification efficiency E.
          efficiencies:
            type: array
            items:
              type: string
            description: Optional explicit per-gene amplification efficiencies as GENE=E pairs (e.g. GAPDH=1.98; default 2.0 unless standard_curve provided).
          ntc_group:
            type: string
            description: Group or sample label identifying No-Template Control wells (default NTC).
      approval: write
      output: ./schemas/qpcr-quant-result.json
      timeoutSeconds: 600
    - name: immunoassay_curve
      description: >-
        Weighted 4PL/5PL immunoassay standard-curve calibration (ELISA, Luminex,
        MSD, BCA/Bradford) with AICc model selection, blank LOD, FDA/EMA
        LLOQ/ULOQ validation (80-120% recovery, <=20% CV), and range-censored
        unknown sample interpolation.
      run: [python, ./scripts/wa_immunoassay_curve.py]
      args:
        type: object
        additionalProperties: false
        required: [standards, output_dir]
        properties:
          standards:
            type: string
            format: input-path
            description: Standard curve TSV/CSV with concentration (0 allowed for blanks) and signal columns.
          output_dir:
            type: string
            format: project-path
            description: Directory for standard recovery table, interpolated unknowns table, and calibration figure.
          unknowns:
            type: string
            format: input-path
            description: Optional unknown samples TSV/CSV with sample, signal, and optional dilution_factor (default 1.0) and group columns.
          model:
            type: string
            enum: [auto, 4PL, 5PL]
            description: Calibration curve family (default auto selects between 4PL and 5PL via AICc).
          weighting:
            type: string
            enum: [1_y2, 1_y, none]
            description: Nonlinear least-squares residual weighting scheme (default 1_y2 for heteroscedastic immunoassays).
      approval: write
      output: ./schemas/immunoassay-curve-result.json
      timeoutSeconds: 600
    - name: enzyme_kinetics
      description: >-
        Steady-state Michaelis-Menten and Hill enzyme kinetics, optional
        progress-curve initial-velocity (v0) extraction, AICc-ranked global
        reversible inhibitor mechanism discrimination (Competitive,
        Uncompetitive, Noncompetitive, Mixed), kcat/Km, and Cheng-Prusoff Ki.
      run: [python, ./scripts/wa_enzyme_kinetics.py]
      args:
        type: object
        additionalProperties: false
        required: [kinetics_data, output_dir]
        properties:
          kinetics_data:
            type: string
            format: input-path
            description: TSV/CSV with substrate_conc, velocity, and optional inhibitor_conc (or time, signal, substrate_conc when from_progress_curves is true).
          output_dir:
            type: string
            format: project-path
            description: Directory for kinetics model comparison, velocity fit table, and 2-panel Michaelis-Menten/Lineweaver-Burk figure.
          from_progress_curves:
            type: boolean
            description: Set true when kinetics_data contains time and signal progress curves requiring automated initial-velocity (v0) extraction (default false).
          enzyme_conc:
            type: number
            description: Optional total enzyme concentration [E]t in the same molar units as substrate_conc to compute kcat and kcat/Km.
          ic50:
            type: number
            description: Optional measured IC50 value (requires cheng_prusoff_substrate) to compute Cheng-Prusoff Ki.
          cheng_prusoff_substrate:
            type: number
            description: Substrate concentration [S] used in the IC50 assay for Cheng-Prusoff Ki calculation.
      approval: write
      output: ./schemas/enzyme-kinetics-result.json
      timeoutSeconds: 600
    - name: binding_kinetics
      description: >-
        Global 1:1 Langmuir SPR/BLI binding kinetics (kon, koff, KD, t1/2, Rmax)
        across multi-concentration sensorgrams with reduced chi^2, mass-transport
        and slow-off-rate diagnostics, and steady-state equilibrium KD_eq
        isotherm cross-validation.
      run: [python, ./scripts/wa_binding_kinetics.py]
      args:
        type: object
        additionalProperties: false
        required: [sensorgrams, output_dir]
        properties:
          sensorgrams:
            type: string
            format: input-path
            description: Reference-subtracted SPR/BLI sensorgram TSV/CSV with time (s), response, concentration (conc_m/conc_nm/conc), and optional phase (assoc/dissoc).
          output_dir:
            type: string
            format: project-path
            description: Directory for kinetic summary, fitted sensorgram traces, and 2-panel sensorgram/isotherm figure.
          dissociation_start:
            type: number
            description: Time in seconds when buffer wash / dissociation begins (required if phase column is absent).
          conc_unit:
            type: string
            enum: [nM, uM, M]
            description: Concentration unit of the input concentration column when not explicitly named conc_m (default nM).
          local_rmax:
            type: boolean
            description: Fit per-concentration local Rmax instead of the default globally shared Rmax (default false).
      approval: write
      output: ./schemas/binding-kinetics-result.json
      timeoutSeconds: 600
---

# Wet-Lab Assay Quantification & Biochemical Kinetics Workflow

Read [references/analysis-playbook.md](references/analysis-playbook.md) for the mathematical formulations, MIQE/FDA/EMA QC thresholds, and diagnostic guardrails.

## Tool Overview

```
RT-qPCR raw Cq table (+ std curve)   ──wa_qpcr_quant────────▶ geNorm/Andersen stability, Pfaffl FC, biol-rep stats
ELISA / BCA standards (+ unknowns)   ──wa_immunoassay_curve─▶ Weighted 4PL/5PL (AICc), LOD/LLOQ/ULOQ, censored unknowns
Enzyme velocities / progress curves  ──wa_enzyme_kinetics───▶ v0 extraction, MM/Hill, AICc inhibitor discrimination, Ki
SPR / BLI multi-conc sensorgrams     ──wa_binding_kinetics──▶ Global 1:1 Langmuir (kon, koff, KD), steady-state KD_eq
```

## Hard Red Lines

1. **No technical pseudoreplication in qPCR**: statistical tests run strictly on biological-replicate $\Delta C_q$ values after technical replicate aggregation (`sd_cq > 0.5` flagged) and NTC contamination checks.
2. **No silent immunoassay extrapolation**: unknown concentrations outside `[LLOQ, ULOQ]` (`80%–120%` standard recovery, `<= 20%` CV) are explicitly labeled `below_LOD`, `below_LLOQ`, or `above_ULOQ`.
3. **No unweighted Lineweaver-Burk parameter fitting**: enzyme kinetics parameters are fitted exclusively via nonlinear least squares; Lineweaver-Burk plots serve solely as diagnostic geometry overlays.
4. **Always validate SPR/BLI kinetic $K_D$**: cross-check against steady-state $K_{D,\text{eq}}$ ($>3\times$ discrepancy flagged) and inspect mass-transport ($k_{\text{on}} > 10^7\text{ M}^{-1}\text{s}^{-1}$) and slow-off-rate ($k_{\text{off}} < 10^{-5}\text{ s}^{-1}$) diagnostics.
