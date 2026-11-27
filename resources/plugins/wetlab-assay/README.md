# wetlab-assay · 0.1.0

Quantitative wet-lab assay analysis and biochemical kinetics plugin (`toolPrefix: wa`, agent `WetlabAssay`) running in a locked MIT/BSD `conda-forge` environment (`plugin:wa`: Python 3.13 + `pandas`/`numpy`/`scipy`/`scikit-learn`/`matplotlib-base`).

## Tools (`wa_*`)

- `wa_qpcr_quant`: MIQE-compliant RT-qPCR relative quantification ($2^{-\Delta\Delta C_q}$ and Pfaffl efficiency-corrected), geNorm ($M_j$) & Andersen ($\rho_j$) reference-gene stability ranking with geometric-mean normalization, NTC contamination & technical replicate SD (`> 0.5`) QC, and biological-replicate statistical testing (Welch $t$-test or one-way ANOVA + Tukey/Dunnett + BH FDR).
- `wa_immunoassay_curve`: ELISA / Luminex / MSD / Bradford-BCA calibration using weighted (`1/y^2`, `1/y`, or unweighted) 4PL and 5PL nonlinear least squares with $\text{AIC}_c$ model selection, blank-based LOD ($\bar{y}_{\text{blank}} + 3.3\,\sigma_{\text{blank}}$), FDA/EMA LLOQ/ULOQ validation (`80%–120%` recovery, `<= 20%` CV), and strict out-of-range censoring (`below_LOD`, `below_LLOQ`, `above_ULOQ`) for diluted unknowns.
- `wa_enzyme_kinetics`: Steady-state Michaelis-Menten and Hill kinetics, automated initial-velocity ($v_0$) extraction from progress curves ($R^2 \ge 0.98$, $<15\%$ substrate depletion warning), $\text{AIC}_c$-ranked reversible inhibitor mechanism discrimination (`Competitive`, `Uncompetitive`, `Noncompetitive`, `Mixed`), $k_{\text{cat}}$ / $k_{\text{cat}}/K_m$, and Cheng-Prusoff $K_i$.
- `wa_binding_kinetics`: Global 1:1 Langmuir SPR/BLI association/dissociation kinetics ($k_{\text{on}}, k_{\text{off}}, K_D, t_{1/2}, R_{\max}$, reduced $\chi^2$, mass-transport & slow-off-rate diagnostics) paired with steady-state equilibrium isotherm $K_{D,\text{eq}}$ cross-validation.
