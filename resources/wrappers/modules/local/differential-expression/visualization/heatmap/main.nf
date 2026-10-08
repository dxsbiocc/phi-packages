// New local process (no upstream module to vendor) -- top-DE-gene heatmap
// adapted from bioSkills's differential-expression/de-visualization
// heatmap.R example. Styling (diverging z-score palette) adapted from
// resources/skills/omics-visualization's heatmap/cluster_basic template's
// intent, implemented with pheatmap (already in this family's shared
// dependencies declared in this module's environment.yml) rather
// than that template's own ComplexHeatmap dependency, which the image
// does not have.
process DE_HEATMAP {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned:ccbe9c69fe5b6749' :
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned@sha256:e500cffd28caf431b76344610bd213a30b536e2d2fa58418255ec00238e3ab12' }"

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
