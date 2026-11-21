process DRAMV {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/dram:1.5.0--pyhdfd78af_0' :
        'quay.io/biocontainers/dram:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(affi)
    path(config)
    val(distill)

    output:
    tuple val(meta), path("annotate/annotations.tsv"), emit: annotations
    tuple val(meta), path("distill")                 , emit: distill, optional: true
    path "versions.yml"                              , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config_opt = config ? "--config_loc ${config}" : ''
    
    """
    DRAM-v.py annotate \\
        --threads ${task.cpus} \\
        -i ${fasta} \\
        -v ${affi} \\
        -o annotate \\
        ${config_opt} \\
        ${args}

    if [ "${distill}" = "true" ]; then
        DRAM-v.py distill \\
            -i annotate/annotations.tsv \\
            -o distill \\
            ${config_opt}
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: 1.5.0
    END_VERSIONS
    """
    stub:
    """
    mkdir -p annotate distill
    touch annotate/annotations.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: 1.5.0
    END_VERSIONS
    
    """
}
