process LEFSE_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/lefse:1.1.2--pyhdfd78af_0' :
        'quay.io/biocontainers/lefse:1.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(formatted)

    output:
    tuple val(meta), path("${prefix}.res"), emit: results
    path "versions.yml"                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    lefse_run.py \\
        ${formatted} \\
        ${prefix}.res \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lefse: 1.1.2
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}.res
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lefse: 1.1.2
    END_VERSIONS
    
    """
}
