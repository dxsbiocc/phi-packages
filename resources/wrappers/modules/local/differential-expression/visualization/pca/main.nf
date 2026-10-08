// New local process (no upstream module to vendor) -- PCA computation
// adapted from bioSkills's differential-expression/de-visualization
// pca_plot.R example and this family's own ../../qc/templates/deseq2_qc.R
// (plotPCA_vst helper). Styling (categorical colorblind-safe palette,
// restrained theme) adapted from resources/skills/omics-visualization's
// scatter/group template, using only packages already in this family's
// declared module-local dependencies.
process DE_PCA {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned:ccbe9c69fe5b6749' :
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned@sha256:e500cffd28caf431b76344610bd213a30b536e2d2fa58418255ec00238e3ab12' }"

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
