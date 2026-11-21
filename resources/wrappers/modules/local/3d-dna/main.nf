process THREEDDNA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/3d-dna:201008--hdfd78af_0' :
        'quay.io/biocontainers/3d-dna:201008--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(mnd)

    output:
    tuple val(meta), path("*.Final.fasta")           , emit: final_fasta, optional: true
    tuple val(meta), path("*.final.assembly.fasta")  , emit: assembly, optional: true
    tuple val(meta), path("*.hic")                   , emit: hic, optional: true
    path "versions.yml"                              , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    3d-dna \\
        ${args} \\
        ${fasta} \\
        ${mnd}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        3d-dna: 201008
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}.Final.fasta
    touch ${meta.id}.final.assembly.fasta
    touch ${meta.id}.hic
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        3d-dna: 201008
    END_VERSIONS
    
    """
}
