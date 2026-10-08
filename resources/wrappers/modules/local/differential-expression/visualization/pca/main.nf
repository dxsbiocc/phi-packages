// New local process (no upstream module to vendor) -- PCA computation
// adapted from bioSkills's differential-expression/de-visualization
// pca_plot.R example and this family's own ../../qc/templates/deseq2_qc.R
// (plotPCA_vst helper). Styling (categorical colorblind-safe palette,
// restrained theme) adapted from resources/skills/omics-visualization's
// scatter/group template, using only packages already in this family's
// shared image (resources/wrappers/images/differential-expression-r/).
process DE_PCA {
    label "process_single"

    // Shared family image (resources/wrappers/images/differential-expression-r/) --
    // built and tagged locally, not published to any registry yet, so no
    // singularity/apptainer branch here (see that directory's Dockerfile).
    conda "${moduleDir}/../../../../../images/differential-expression-r/environment.yml"
    container 'phi/differential-expression-r:1.0.0'

    input:
    path matrix
    path samplesheet

    output:
    path "*.pca.pdf"     , optional:true, emit: pdf
    path "*.pca.vals.tsv", optional:true, emit: vals
    path "versions.yml"  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'pca.R'

    stub:
    prefix = task.ext.prefix ?: "de"
    """
    touch ${prefix}.pca.pdf
    touch ${prefix}.pca.vals.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        r-ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
    """
}
