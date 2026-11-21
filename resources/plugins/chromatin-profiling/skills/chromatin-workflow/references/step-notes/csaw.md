Purpose: window-based differential binding with csaw — the preferred route for
broad marks (H3K27me3, H3K36me3, H3K9me3) whose signal spreads over kilobases
rather than sharp peaks. Principle: csaw counts reads in sliding windows
(windowCounts, width/spacing), filters background windows (filterWindows),
normalizes for composition bias from large bins (normFactors/TMM), fits a
quasi-likelihood glm (edgeR estimateDisp/glmQLFit/glmQLFTest), then merges adjacent
windows into regions and combines window p-values (mergeWindows/combineTests,
Simes). Peak-based DiffBind can miss diffuse changes this catches. Interpretation:
n_windows -> n_regions -> n_significant at FDR; the region table is FDR-sorted.
Pitfalls: csaw has no osx-arm64 bioconda build — on darwin-arm64 the step returns
an error naming the gap (use DiffBind, or linux-x64/darwin-x64, or a project env);
never silently fall back and call it csaw. Needs a sample sheet with a Condition
column and indexed BAM; with >2 conditions coef tests the last level vs the first
(a note says so) — set explicit contrasts for others. fragment_size extends
single-end reads to the estimated fragment length.