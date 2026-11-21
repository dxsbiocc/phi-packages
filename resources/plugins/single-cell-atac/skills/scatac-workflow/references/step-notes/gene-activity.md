Purpose: aggregate accessibility over gene bodies and promoters into a
gene-by-cell matrix that proxies expression. Principle: make_gene_matrix sums
fragment accessibility across each gene's body plus an upstream promoter window
(gene_upstream bp), using a gene annotation; the result is a separate gene_adata
object. Interpretation: gene activity lets you score marker panels and transfer
scRNA labels — it is the bridge between accessibility and cell identity.
Pitfalls: gene activity is a PROXY, not measured RNA — it is sparse and biased by
gene length and local accessibility, so use it for annotation/marker
visualization, never differential-expression claims; it needs a gene annotation
(GTF/GFF3), without gene_anno the step reports a gap rather than downloading a
~50MB annotation silently; gene_name_key must match the annotation's attribute
(gene_name) or genes come out as IDs. include_gene_body and gene_upstream change
what "activity" means — say what you used.