process GRAPHLAN_ANNOTATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/graphlan:1.1.3--2' :
        'quay.io/biocontainers/graphlan:1.1.3--2' }"

    input:
    tuple val(meta), path(tree), path(annotation)

    output:
    tuple val(meta), path("${prefix}.xml"), emit: xml
    path "versions.yml"                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    graphlan_annotate.py \\
        --annot ${annotation} \\
        ${args} \\
        ${tree} \\
        ${prefix}.xml
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphlan: 1.1.3
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}.xml
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphlan: 1.1.3
    END_VERSIONS
    
    """
}
