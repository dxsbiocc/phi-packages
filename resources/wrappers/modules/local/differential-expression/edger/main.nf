// New local process (no upstream module to vendor -- nf-core/modules has
// no standalone edgeR differential-expression module, only limma/
// differential's own limma-voom path for RNA-seq counts). Statistical
// approach (TMM normalization, glmQLFit/glmQLFTest -- the modern
// quasi-likelihood F-test edgeR default) follows bioSkills's
// differential-expression/edger-basics skill.
process EDGER_DIFFERENTIAL {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned:ccbe9c69fe5b6749' :
        'community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned@sha256:e500cffd28caf431b76344610bd213a30b536e2d2fa58418255ec00238e3ab12' }"

    input:
    path counts
    path samplesheet

    output:
    path "*.edger.results.tsv"      , optional:true, emit: results
    path "*.edger.bcv.pdf"          , optional:true, emit: bcv_plot
    path "*.normalised_counts.tsv"  , optional:true, emit: normalised_counts
    path "*.edger.model.txt"        , optional:true, emit: model
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'edger_differential.R'

    stub:
    prefix = task.ext.prefix ?: "de"
    """
    touch ${prefix}.edger.results.tsv
    touch ${prefix}.edger.bcv.pdf
    touch ${prefix}.normalised_counts.tsv
    touch ${prefix}.edger.model.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """
}
