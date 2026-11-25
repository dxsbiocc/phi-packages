# pharmacogenomics · 0.1.0

Pharmacogenomics, high-throughput drug screening, and combination therapy synergy plugin (`toolPrefix: pgx`, agent `Pharmacogenomics`) in a locked offline Python environment (`plugin:pgx`: Python 3.13 + `pandas/numpy/scipy/scikit-learn/matplotlib-base`).

## Relationship with `pharmacogenomics` MCP connector

- **Upstream database lookup (`resources/connectors/mcp-servers/pharmacogenomics`)**: queries public GDSC, PRISM/DepMap, PharmGKB, and CPIC reference records.
- **Local quantitative modeling & cohort integration (this plugin, 3 `pgx_*` tools)**:
  - `pgx_dose_response`: fit 4-parameter log-logistic Hill curves ($E_0, E_{\infty}, \text{EC}_{50}, h$) across compound–cell-line screens, compute relative $\text{EC}_{50}$, absolute $\text{IC}_{50}$, normalized trapezoidal `AUC` & `AAC`, $R^2$, and flag `beyond_tested_range` whenever $\text{IC}_{50} > \max(\text{dose})$.
  - `pgx_synergy`: evaluate 2D drug-combination viability matrices under **HSA**, **Bliss Independence**, **Loewe Additivity**, and **ZIP (Zero Interaction Potency)** models (`> 10` synergistic, `-10..10` additive, `< -10` antagonistic), exporting full synergy matrices and 2-panel landscape heatmaps.
  - `pgx_biomarker`: screen somatic mutations/CNVs (MAF or binary matrix via Welch $t$ / Mann-Whitney + Cohen's $d$) or continuous expression/pathway scores (Spearman $\rho$) against drug sensitivity (`IC50`/`LN_IC50`/`AUC`/`AAC`) with explicit sensitivity polarity enforcement and BH FDR control.
