# Wet-Lab Assay & Biochemical Kinetics Playbook

## 1. MIQE RT-qPCR Quantification (`wa_qpcr_quant`)

- **Technical Replicate QC**: Aggregate technical wells per `(sample, gene)` as mean $C_q$ and sample SD (`sd_cq`). Any `(sample, gene)` with `sd_cq > 0.5` triggers a MIQE precision warning.
- **NTC Contamination QC**: Inspect `ntc_group` (default `NTC`) wells per gene. Contamination is flagged whenever `NTC_Cq < 35.0` or `NTC_Cq - min(sample_Cq) < 5.0` (while `Undetermined` / `NaN` NTC wells pass).
- **Amplification Efficiency**: From a dilution/copy-number standard curve ($C_q$ vs $\log_{10}(\text{quantity})$), fit linear regression $C_q = m \log_{10} Q + b$:
  $$E = 10^{-1/m}$$
  MIQE acceptable range is $90\%\text{–}110\%$ efficiency, i.e., $E \in [1.90, 2.10]$ (or $R^2 \ge 0.98$). When neither `standard_curve` nor `efficiencies` is supplied for a gene, $E = 2.0$ is assumed.
- **Reference Gene Stability**:
  - **Vandesompele et al. (2002) geNorm $M_j$**: For each pair of reference genes $(j, k)$, compute log-transformed relative quantity ratios $A_{jk, i} = \log_2(Q_{ji} / Q_{ki}) = -\Delta C_{q,ji}\log_2 E_j + \Delta C_{q,ki}\log_2 E_k$ across samples $i$; $V_{jk} = \operatorname{SD}_i(A_{jk, i})$, and $M_j = \frac{1}{R - 1}\sum_{k \ne j} V_{jk}$. Lower $M_j$ indicates higher stability ($M < 0.5$ optimal for homogeneous samples, $M < 1.0$ acceptable for heterogeneous tissues).
  - **Andersen et al. (2004) Variance Stability $\rho_j$**: Decompose log2-quantities $y_{jgi} = \log_2(Q_{jgi})$ across groups $g$ and samples $i$ into group deviations and sample-level residuals relative to the multi-reference consensus; compute intra-group variance $\hat{\sigma}_{jg}^2$ and inter-group variance $\hat{\tau}_j^2$, yielding stability index $\rho_j$ combining inter-group bias and intra-group dispersion (lower = more stable).
- **Pfaffl Multi-Reference Normalization & Statistical Testing**:
  For each biological replicate $s$ and target gene $t$, define calibrator mean $\bar{C}_{q,\text{ctrl}}$ in `control_group`. Compute geometric-mean reference factor:
  $$\text{NF}_s = \left(\prod_{r=1}^{R} E_r^{\bar{C}_{q,r,\text{ctrl}} - C_{q,r,s}}\right)^{1/R}$$
  and efficiency-corrected relative expression:
  $$\text{Ratio}_{t,s} = \frac{E_t^{\bar{C}_{q,t,\text{ctrl}} - C_{q,t,s}}}{\text{NF}_s}$$
  When all $E = 2.0$, $\text{Ratio}_{t,s}$ equals $2^{-\Delta\Delta C_q}$. Statistical comparisons are performed strictly on biological replicates on the $\log_2(\text{Ratio}_{t,s})$ (equivalent $-\Delta\Delta C_q$) scale: Welch's $t$-test for 2 groups, or one-way ANOVA followed by `scipy.stats.tukey_hsd` / Dunnett-style vs `control_group` for $\ge 3$ groups, with Benjamini-Hochberg FDR correction.

## 2. Immunoassay Calibration (`wa_immunoassay_curve`)

- **Blank & LOD**: Wells with `concentration == 0` yield $\bar{y}_{\text{blank}}$ and $\sigma_{\text{blank}}$; $\text{LOD}_{\text{signal}} = \bar{y}_{\text{blank}} + 3.3\,\sigma_{\text{blank}}$ (or $+3.3\,\sigma_{\text{res}}$ fallback if blanks are absent or zero-variance), inverted onto the fitted curve.
- **Weighted 4PL & 5PL Models**:
  - **4PL**: $y(x) = D + \frac{A - D}{1 + (x/C)^B}$
  - **5PL**: $y(x) = D + \frac{A - D}{(1 + (x/C)^B)^S}$
  Fitted on positive standards via multi-start bounded `scipy.optimize.least_squares` with weights $w_i \in \{1/\hat{y}_i^2,\; 1/\hat{y}_i,\; 1\}$. When `model="auto"`, select between 4PL and 5PL using small-sample corrected $\text{AIC}_c = n \ln(\text{WRSS}/n) + 2k + \frac{2k(k+1)}{n - k - 1}$.
- **Bioanalytical Validation**: Back-calculate every standard well (`%Recovery = 100 * back_calc / nominal`) and level `%CV`. `LLOQ` and `ULOQ` are the lowest and highest standard levels meeting `80% <= %Recovery <= 120%` and `%CV <= 20%`. Unknowns are diluted-adjusted (`conc = back_calc * dilution_factor`) and strictly tagged `in_quantifiable_range`, `below_LOD`, `below_LLOQ`, or `above_ULOQ`.

## 3. Enzyme Kinetics & Inhibition (`wa_enzyme_kinetics`)

- **Progress Curves (`from_progress_curves`)**: For each curve, identify the early linear phase maximizing length while maintaining $R^2 \ge 0.98$, extracting slope $v_0$; flag curves exceeding $15\%$ substrate depletion.
- **Uninhibited Kinetics**: Fit **Michaelis-Menten** $v = \frac{V_{\max}[S]}{K_m + [S]}$ and **Hill** $v = \frac{V_{\max}[S]^{n_H}}{K_{0.5}^{n_H} + [S]^{n_H}}$ via nonlinear least squares, reporting SEs, 95% CIs, $R^2$, $\text{AIC}_c$, and optionally $k_{\text{cat}} = V_{\max}/[E]_t$ and $k_{\text{cat}}/K_m$.
- **Global Reversible Inhibition Discrimination**: Fit all 4 models simultaneously across `[S]` and `[I]`:
  - **Competitive**: $v = \frac{V_{\max}[S]}{K_m(1 + [I]/K_i) + [S]}$
  - **Uncompetitive**: $v = \frac{V_{\max}[S]}{K_m + [S](1 + [I]/K_i')}$
  - **Pure Noncompetitive**: $v = \frac{V_{\max}[S]}{(K_m + [S])(1 + [I]/K_i)}$
  - **Mixed**: $v = \frac{V_{\max}[S]}{K_m(1 + [I]/K_i) + [S](1 + [I]/K_i')}$
  Rank models by $\text{AIC}_c$ and Akaike weights $w_m = \frac{e^{-\Delta_m/2}}{\sum_j e^{-\Delta_j/2}}$. Optionally compute Cheng-Prusoff $K_i = \frac{\text{IC}_{50}}{1 + [S]/K_m}$.

## 4. Global 1:1 Langmuir SPR / BLI Kinetics (`wa_binding_kinetics`)

- **1:1 Langmuir Model** (association $t \in [t_0, t_d]$, dissociation $t > t_d$):
  $$R(t) = \frac{k_{\text{on}} C R_{\max}}{k_{\text{on}} C + k_{\text{off}}}\left(1 - e^{-(k_{\text{on}} C + k_{\text{off}})(t - t_0)}\right), \quad t \le t_d$$
  $$R(t) = R(t_d)\,e^{-k_{\text{off}}(t - t_d)}, \quad t > t_d$$
  Fitted globally across all concentrations via `least_squares` with Jacobian-derived 95% CIs, kinetic $K_D = k_{\text{off}}/k_{\text{on}}$, $t_{1/2} = \ln(2)/k_{\text{off}}$, and reduced $\chi^2 = \text{RSS}/\text{df}$.
- **Steady-State Equilibrium Cross-Validation**: Fit pre-dissociation plateau responses $R_{\text{eq}}(C) = \frac{C\,R_{\max,\text{eq}}}{K_{D,\text{eq}} + C}$ when $\ge 4$ concentrations are available; warn if $K_D$ vs $K_{D,\text{eq}}$ differs by $>3\times$, $k_{\text{on}} > 10^7\text{ M}^{-1}\text{s}^{-1}$ (mass transport), or $k_{\text{off}} < 10^{-5}\text{ s}^{-1}$ (insufficient dissociation).
