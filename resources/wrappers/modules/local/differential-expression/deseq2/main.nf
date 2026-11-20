process DESEQ2_DIFFERENTIAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned:ccbe9c69fe5b6749' :
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned@sha256:e500cffd28caf431b76344610bd213a30b536e2d2fa58418255ec00238e3ab12' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target), val(formula), val(comparison)
    tuple val(meta2), path(samplesheet), path(counts)
    tuple val(meta3), path(control_genes_file)
    tuple val(meta4), path(transcript_lengths_file)

    output:
    tuple val(meta), path("*.deseq2.results.tsv")              , emit: results
    tuple val(meta), path("*.deseq2.qc.tsv")                   , emit: research_qc
    tuple val(meta), path("*.deseq2.dispersion.png")           , emit: dispersion_plot_png
    tuple val(meta), path("*.deseq2.dispersion.pdf")           , emit: dispersion_plot_pdf
    tuple val(meta), path("*.dds.rld.rds")                     , emit: rdata
    tuple val(meta), path("*.deseq2.sizefactors.tsv")          , emit: size_factors
    tuple val(meta), path("*.normalised_counts.tsv")           , emit: normalised_counts
    tuple val(meta), path("*.rlog.tsv")                        , optional: true, emit: rlog_counts
    tuple val(meta), path("*.vst.tsv")                         , optional: true, emit: vst_counts
    tuple val(meta), path("*.deseq2.model.txt")                , emit: model
    tuple val(meta), path("*.R_sessionInfo.log")               , emit: session_info
    tuple val(meta), path("*.deseq2.de_genes_filtered.tsv")    , emit: de_genes
    tuple val(meta), path("*.deseq2.de_genes_filtered_up.tsv") , emit: de_genes_up
    tuple val(meta), path("*.deseq2.de_genes_filtered_down.tsv"), emit: de_genes_down
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq2_differential.R'

    stub:
    """
    # Header-only results: no evaluable rows, so the QC sidecar is 'invalid'
    # with evaluated=hits=0 and the MD5 of these exact stub bytes.
    results_file=${meta.id}.deseq2.results.tsv
    printf '${meta.gene_id_col ?: 'gene_id'}\\tbaseMean\\tlog2FoldChange\\tlfcSE\\tstat\\tpvalue\\tpadj\\n' > \$results_file
    for suffix in de_genes_filtered de_genes_filtered_up de_genes_filtered_down; do
        printf '${meta.gene_id_col ?: 'gene_id'}\\tbaseMean\\tlog2FoldChange\\tlfcSE\\tstat\\tpvalue\\tpadj\\n' > ${meta.id}.deseq2.\$suffix.tsv
    done
    touch ${meta.id}.deseq2.dispersion.png
    touch ${meta.id}.deseq2.dispersion.pdf
    touch ${meta.id}.dds.rld.rds
    touch ${meta.id}.deseq2.sizefactors.tsv
    touch ${meta.id}.normalised_counts.tsv
    touch ${meta.id}.rlog.tsv
    touch ${meta.id}.deseq2.model.txt
    touch ${meta.id}.R_sessionInfo.log

    result_md5=\$(Rscript -e "cat(unname(tools::md5sum('\$results_file')))")
    printf 'schemaVersion\\tdisposition\\tcontrast\\tgeneIdColumn\\tpadjThreshold\\tlog2fcThreshold\\tevaluated\\thits\\tresultMd5\\n' > ${meta.id}.deseq2.qc.tsv
    printf '1\\tinvalid\\t%s\\t%s\\t%s\\t%s\\t0\\t0\\t%s\\n' '${meta.contrast ?: 'stub'}' '${meta.gene_id_col ?: 'gene_id'}' '${meta.padj_threshold != null ? meta.padj_threshold : 0.05}' '${meta.log2fc_threshold != null ? meta.log2fc_threshold : 1}' "\$result_md5" >> ${meta.id}.deseq2.qc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
