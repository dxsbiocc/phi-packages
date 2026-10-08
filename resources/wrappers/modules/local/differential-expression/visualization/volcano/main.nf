// New local process (no upstream module to vendor) -- statistical
// classification adapted from bioSkills's differential-expression/
// de-visualization example, styling (diverging red/blue direction colors,
// restrained theme) adapted from resources/skills/omics-visualization's
// scatter/volcano template, using only packages already in this family's
// shared image (resources/wrappers/images/differential-expression-r/)
// rather than that template's own ggrepel/ggprism dependencies.
process DE_VOLCANO {
    label "process_single"

    // Shared family image (resources/wrappers/images/differential-expression-r/) --
    // built and tagged locally, not published to any registry yet, so no
    // singularity/apptainer branch here (see that directory's Dockerfile).
    conda "${moduleDir}/../../../../../images/differential-expression-r/environment.yml"
    container 'phi/differential-expression-r:1.0.0'

    input:
    path results

    output:
    path "*.volcano.pdf"       , optional:true, emit: pdf
    path "*.volcano.labels.tsv", optional:true, emit: labels
    path "versions.yml"        , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'volcano.R'

    stub:
    prefix = task.ext.prefix ?: "de"
    """
    touch ${prefix}.volcano.pdf
    touch ${prefix}.volcano.labels.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        r-ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
    """
}
