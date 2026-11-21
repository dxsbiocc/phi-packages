Purpose: differential binding with DiffBind + edgeR — the alternative to DESeq2
for TF/narrow-mark ChIP and CUT&Tag. Principle: edgeR's exact-test/QL framework is
often more sensitive with few replicates, while DESeq2 is more conservative; both
run on the same DiffBind consensus count matrix, so switching method does not
change the peaks being tested. Interpretation: dba.contrast sets the comparison
(explicit contrast or DBA_CONDITION categories with minMembers);
dba.analyze(DBA_EDGER) fits; dba.report(th=fdr, fold=lfc, bUsePval=FALSE) returns
FDR-sorted peaks; the MA plot uses the edgeR method. The choice between DESeq2 and
edgeR is the analyst's — state it. Pitfalls: same red lines as DESeq2 — replicates
are the experimental unit (min_members>=2), raw counts with internal normalization
only, and a null result is not "no regulation". With few replicates edgeR may call
more peaks; that is sensitivity, not a license to over-claim.