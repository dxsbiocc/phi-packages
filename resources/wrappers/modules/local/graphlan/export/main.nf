process GRAPHLAN_EXPORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/export2graphlan:0.22--py_0' :
        'quay.io/biocontainers/export2graphlan:0.22--py_0' }"

    input:
    tuple val(meta), path(taxonomy)

    output:
    tuple val(meta), path("${prefix}.tree.xml") , emit: tree
    tuple val(meta), path("${prefix}.annot.xml"), emit: annotation
    path "versions.yml"                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    export2graphlan.py \\
        -i ${taxonomy} \\
        --tree ${prefix}.tree.xml \\
        --annotation ${prefix}.annot.xml \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        export2graphlan: 0.22
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}.tree.xml
    touch ${meta.id}.annot.xml
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        export2graphlan: 0.22
    END_VERSIONS
    
    """
}
