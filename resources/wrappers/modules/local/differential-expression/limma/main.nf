process LIMMA_DIFFERENTIAL {
    tag "$meta.id"
    label 'process_single'

    // Shared family image (resources/wrappers/images/differential-expression-r/) --
    // built and tagged locally, not published to any registry yet, so no
    // singularity/apptainer branch here (see that directory's Dockerfile).
    conda "${moduleDir}/../../../../images/differential-expression-r/environment.yml"
    container 'phi/differential-expression-r:1.0.0'

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target), val(formula), val(comparison)
    tuple val(meta2), path(samplesheet), path(intensities)

    output:
    tuple val(meta), path("*.limma.results.tsv")          , emit: results
    tuple val(meta), path("*.limma.mean_difference.png")  , emit: md_plot
    tuple val(meta), path("*.MArrayLM.limma.rds")         , emit: rdata
    tuple val(meta), path("*.limma.model.txt")            , emit: model
    tuple val(meta), path("*.R_sessionInfo.log")          , emit: session_info
    tuple val(meta), path("*.normalised_counts.tsv")      , emit: normalised_counts, optional: true
    path "versions.yml"                                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'limma_de.R'

    stub:
    prefix              = task.ext.prefix   ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    library(limma)
    a <- file("${prefix}.limma.results.tsv", "w")
    close(a)
    a <- file("${prefix}.limma.mean_difference.png", "w")
    close(a)
    a <- file("${prefix}.MArrayLM.limma.rds", "w")
    close(a)
    a <- file("${prefix}.normalised_counts.tsv", "w")
    close(a)
    a <- file("${prefix}.limma.model.txt", "w")
    close(a)
    a <- file("${prefix}.R_sessionInfo.log", "w")
    close(a)
    ## VERSIONS FILE
    r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
    limma.version <- as.character(packageVersion('limma'))
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', r.version),
            paste('    bioconductor-limma:', limma.version)
        ),
        'versions.yml'
    )
    """
}
