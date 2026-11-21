process LEFSE_PLOT_RES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/lefse:1.1.2--pyhdfd78af_0' :
        'quay.io/biocontainers/lefse:1.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(res)

    output:
    tuple val(meta), path("${prefix}.png"), emit: plot
    path "versions.yml"                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    lefse_plot_res.py \\
        ${res} \\
        ${prefix}.png \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lefse: 1.1.2
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lefse: 1.1.2
    END_VERSIONS
    
    """
}
