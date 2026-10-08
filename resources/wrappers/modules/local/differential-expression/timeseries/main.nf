process VARIANCEPARTITION_DREAM {
    tag "${meta.id}"
    label 'process_high'

    // Shared family image (resources/wrappers/images/differential-expression-r/) --
    // built and tagged locally, not published to any registry yet, so no
    // singularity/apptainer branch here (see that directory's Dockerfile).
    conda "${moduleDir}/../../../../images/differential-expression-r/environment.yml"
    container 'phi/differential-expression-r:1.0.0'

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target), val(formula), val(comparison)
    tuple val(meta2), path(samplesheet), path(counts)

    output:
    tuple val(meta), path("*.dream.results.tsv")        , emit: results
    tuple val(meta), path("*.dream.model.txt")          , emit: model
    tuple val(meta), path("*.normalised_counts.tsv")    , emit: normalised_counts, optional: true
    path "versions.yml"                                 , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'dream.R'

    stub:
    """
    touch "${meta.id}.dream.results.tsv"
    touch "${meta.id}.dream.model.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
        r-variancepartition: \$(Rscript -e "library(variancePartition); cat(as.character(packageVersion('variancePartition')))")
    END_VERSIONS
    """
}
