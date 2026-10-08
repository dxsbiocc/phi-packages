// New local process (no upstream module to vendor) -- statistical
// classification adapted from bioSkills's differential-expression/
// de-visualization example, styling (diverging red/blue direction colors,
// restrained theme) adapted from resources/skills/omics-visualization's
// scatter/volcano template, using only packages already in this family's
// declared module-local dependencies
// rather than that template's own ggrepel/ggprism dependencies.
process DE_VOLCANO {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned:ccbe9c69fe5b6749' :
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned@sha256:e500cffd28caf431b76344610bd213a30b536e2d2fa58418255ec00238e3ab12' }"

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
