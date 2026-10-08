// Adapted from nf-core/rnaseq's DESEQ2_QC module (Harshil Patel, Gavin
// Kelly; MIT license) -- rewritten to Nextflow's `template` mechanism to
// match this family's other processes (deseq2, limma, timeseries). The
// optional MultiQC fragments are retained for the full RNA-seq workflow.
process DESEQ2_QC {
    label "process_medium"

    // Shared family image (resources/wrappers/images/differential-expression-r/) --
    // built and tagged locally, not published to any registry yet, so no
    // singularity/apptainer branch here (see that directory's Dockerfile).
    conda "${moduleDir}/../../../../images/differential-expression-r/environment.yml"
    container 'phi/differential-expression-r:1.0.0'

    input:
    path counts
    path pca_header_multiqc
    path clustering_header_multiqc

    output:
    path "*.pdf"              , optional:true, emit: pdf
    path "*.RData"             , optional:true, emit: rdata
    path "*.pca.vals.txt"      , optional:true, emit: pca_txt
    path "*pca.vals_mqc.tsv"   , optional:true, emit: pca_multiqc
    path "*.sample.dists.txt"  , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    path "*.log"               , optional:true, emit: log
    path "size_factors"        , optional:true, emit: size_factors
    path "versions.yml"        , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq2_qc.R'

    stub:
    prefix = task.ext.prefix ?: "deseq2"
    """
    touch ${prefix}.dds.RData
    touch ${prefix}.pca.vals.txt
    touch sample.pca.vals_mqc.tsv
    touch ${prefix}.plots.pdf
    touch ${prefix}.sample.dists.txt
    touch sample.sample.dists_mqc.tsv
    touch R_sessionInfo.log

    mkdir size_factors
    touch size_factors/${prefix}.size_factors.RData
    for i in `head $counts -n 1 | cut -f3-`;
    do
        touch size_factors/\${i}.size_factors.RData
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
