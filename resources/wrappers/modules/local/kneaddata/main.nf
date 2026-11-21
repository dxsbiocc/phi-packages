process KNEADDATA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/kneaddata:0.12.4--pyhdfd78af_0' :
        'quay.io/biocontainers/kneaddata:0.12.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    tuple val(meta), path("*_kneaddata.fastq"), emit: reads, optional: true
    tuple val(meta), path("*contam.fastq")    , emit: contaminant, optional: true
    tuple val(meta), path("kneaddata.txt")    , emit: counts
    tuple val(meta), path("*.log")            , emit: log
    path "versions.yml"                       , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def read_list = reads instanceof List ? reads : [reads]
    def input_cmd = read_list.size() == 2 ? "--input1 ${read_list[0]} --input2 ${read_list[1]}" : "--unpaired ${read_list[0]}"
    
    """
    kneaddata \\
        ${input_cmd} \\
        --reference-db ${db} \\
        -t ${task.cpus} \\
        --output . \\
        ${args}

    kneaddata_read_count_table \\
        --input . \\
        --output kneaddata.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kneaddata: 0.12.4
    END_VERSIONS
    """
    stub:
    """
    touch sample_kneaddata.fastq
    touch kneaddata.txt
    touch sample.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kneaddata: 0.12.4
    END_VERSIONS
    
    """
}
