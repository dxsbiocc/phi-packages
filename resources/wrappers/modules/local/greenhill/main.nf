process GREENHILL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/greenhill:1.1.0--h663a4a6_3' :
        'quay.io/biocontainers/greenhill:1.1.0--h663a4a6_3' }"

    input:
    tuple val(meta), path(hic)
    path(aware)
    path(bubble)
    path(cph)
    path(long)
    val(hic_type)

    output:
    tuple val(meta), path("${prefix}*"), emit: results
    path "versions.yml"                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = hic_type == 'paired' ? 'HIC' : 'hic'
    def aware_opt = aware ? "-c ${aware}" + (bubble ? " -b ${bubble.join(' ')}" : '') : ''
    def cph_opt = cph ? "-cph ${cph}" : ''
    def long_opt = long ? "-p ${long.join(' ')}" : ''
    
    """
    greenhill \\
        ${aware_opt} \\
        ${cph_opt} \\
        ${long_opt} \\
        -t ${task.cpus} \\
        ${mode} ${hic} \\
        -o ${prefix} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        greenhill: 1.1.0
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}.scaffolds.fa
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        greenhill: 1.1.0
    END_VERSIONS
    
    """
}
