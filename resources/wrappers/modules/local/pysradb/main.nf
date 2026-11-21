process PYSRADB {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/pysradb:2.5.1--pyhdfd78af_0' :
        'quay.io/biocontainers/pysradb:2.5.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(command), val(query)

    output:
    tuple val(meta), path("${prefix}*"), emit: results
    path "versions.yml"                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    pysradb \\
        ${command} \\
        ${query} \\
        --saveto ${prefix}.tsv \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysradb: 2.5.1
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysradb: 2.5.1
    END_VERSIONS
    
    """
}
