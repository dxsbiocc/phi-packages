# Gene activity aggregates accessibility over each gene body plus its promoter
# (upstream bp) into a gene-by-cell matrix — a PROXY for expression, not RNA. It lets
# scATAC clusters be annotated with the same marker panels used for scRNA, and lets
# label transfer from a scRNA reference work in a shared gene space. It needs a gene
# annotation and the fragments, so without gene_anno (or for a 10x matrix input with
# no fragments) the cell reports the gap and continues. The result is a separate
# object (gene_adata); gene_name_key keeps gene symbols available for panel scoring.
if not gene_anno:
    analysis_notes.append("gene activity skipped: no gene_anno (pass a GTF/GFF3 to enable it)")
    info("gene activity skipped: gene_anno is null")
elif not os.path.exists(gene_anno):
    analysis_notes.append("gene activity skipped: gene_anno not found: " + gene_anno)
    info("gene_anno not found: " + gene_anno)
else:
    try:
        os.makedirs(tables_dir, exist_ok=True)
        _gene_path = os.path.join(tables_dir, "gene_activity.h5ad")
        gene_adata = snap.pp.make_gene_matrix(
            adata, gene_anno, inplace=False, file=_gene_path,
            id_type="gene", upstream=gene_upstream, include_gene_body=gene_body,
            gene_name_key="gene_name",
        )
        analysis_steps.append("gene-activity")
        print(f"gene activity: {gene_adata.n_obs} cells x {gene_adata.n_vars} genes -> {_gene_path}")
    except Exception as error:
        analysis_notes.append("gene activity skipped: " + str(error))
        info("gene activity skipped: " + str(error))