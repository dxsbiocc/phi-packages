// New local process (no upstream module to vendor) -- top-DE-gene heatmap
// adapted from bioSkills's differential-expression/de-visualization
// heatmap.R example. Styling (diverging z-score palette) adapted from
// resources/skills/omics-visualization's heatmap/cluster_basic template's
// intent, implemented with pheatmap (already in this family's shared
// image, resources/wrappers/images/differential-expression-r/) rather
// than that template's own ComplexHeatmap dependency, which the image
// does not have.
process DE_HEATMAP {
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
    path "*.heatmap.pdf", optional:true, emit: pdf
    path "versions.yml" , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'heatmap.R'

    stub:
    prefix = task.ext.prefix ?: "de"
    """
    touch ${prefix}.heatmap.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        r-pheatmap: \$(Rscript -e "cat(as.character(packageVersion('pheatmap')))")
    END_VERSIONS
    """
}
