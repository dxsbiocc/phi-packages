process ISESCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/isescan:1.7.3--hab16a5f_1' :
        'quay.io/biocontainers/isescan:1.7.3--hab16a5f_1' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("*.tsv")  , emit: tsv
    tuple val(meta), path("*.fasta"), emit: fasta, optional: true
    tuple val(meta), path("*.gff")  , emit: gff, optional: true
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    isescan.py \\
        --seqfile ${contigs} \\
        --output . \\
        --nthread ${task.cpus} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isescan: 1.7.3
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isescan: 1.7.3
    END_VERSIONS
    
    """
}
