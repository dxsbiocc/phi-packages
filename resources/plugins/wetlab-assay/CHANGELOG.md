# Changelog — wetlab-assay

## 0.1.0 (2026-11-26)

Initial release of the wet-lab assay quantification and biochemical kinetics plugin (`toolPrefix: wa`, agent `WetlabAssay`, private locked environment `plugin:wa`).

- `wa_qpcr_quant`: MIQE-compliant RT-qPCR relative quantification (Pfaffl & $2^{-\Delta\Delta C_q}$), Vandesompele geNorm $M_j$ and Andersen $\rho_j$ reference stability ranking, geometric-mean multi-reference normalization, NTC/SD QC guards, and biological-replicate ANOVA/Welch testing with BH FDR.
- `wa_immunoassay_curve`: Weighted 4PL/5PL standard-curve calibration with $\text{AIC}_c$ selection, blank LOD, bioanalytical LLOQ/ULOQ bounds (`80%–120%` recovery, `<= 20%` CV), and range-censored unknown interpolation.
- `wa_enzyme_kinetics`: Progress-curve $v_0$ extraction, nonlinear Michaelis-Menten & Hill fitting, $\text{AIC}_c$ global reversible inhibitor discrimination (`Competitive`, `Uncompetitive`, `Noncompetitive`, `Mixed`), and Cheng-Prusoff $K_i$.
- `wa_binding_kinetics`: Global 1:1 Langmuir SPR/BLI sensorgram fitting ($k_{\text{on}}, k_{\text{off}}, K_D, t_{1/2}, R_{\max}$) with steady-state $K_{D,\text{eq}}$ cross-validation and mass-transport / slow-dissociation QC warnings.
