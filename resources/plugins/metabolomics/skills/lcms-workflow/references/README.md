# LC-MS workflow references

Decision loop: [analysis-playbook.md](analysis-playbook.md) — read it before
chaining tools. Per-tool depth:

- [step-notes/inspect.md](step-notes/inspect.md) — vendor formats, RT units, QC detection
- [step-notes/qc.md](step-notes/qc.md) — RSD/detection filtering, QC design
- [step-notes/preprocess.md](step-notes/preprocess.md) — imputation, normalization, transform
- [step-notes/multivariate.md](step-notes/multivariate.md) — PCA/PLS-DA, VIP, permutation Q2
- [step-notes/univariate.md](step-notes/univariate.md) — test choice, FDR + FC, power
- [step-notes/annotate.md](step-notes/annotate.md) — isotopes/adducts/losses, MSI levels
- [step-notes/ms2-match.md](step-notes/ms2-match.md) — MSP/MGF parsing, cosine, library provenance
- [step-notes/pathway.md](step-notes/pathway.md) — hypergeometric vs mummichog-style permutation
- [step-notes/plot.md](step-notes/plot.md) — heatmap reading, clustering choices

Upstream tools (msconvert, xcms, MZmine, MS-DIAL) have **no wrapper** in this
repository; the playbook §1 explains what to ask the user for instead.

For maintainers: [maintainers/source-map.md](maintainers/source-map.md).