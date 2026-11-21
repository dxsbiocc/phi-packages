process LEFSE_FORMAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/lefse:1.1.2--pyhdfd78af_0' :
        'quay.io/biocontainers/lefse:1.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(table)

    output:
    tuple val(meta), path("${prefix}.lefse"), emit: formatted
    path "versions.yml"                     , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    lefse_format_input.py \\
        ${table} \\
        ${prefix}.lefse \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lefse: 1.1.2
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}.lefse
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lefse: 1.1.2
    END_VERSIONS
    
    """
}
