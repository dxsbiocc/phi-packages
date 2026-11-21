process VIRSORTER2 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/virsorter:2.2.4--pyhdfd78af_2' :
        'quay.io/biocontainers/virsorter:2.2.4--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(contigs)
    path(db_dir)
    val(mode)
    val(include_groups)
    val(min_length)

    output:
    tuple val(meta), path("results/final-viral-combined.fa"), emit: combined, optional: true
    tuple val(meta), path("results/final-viral-score.tsv")  , emit: score, optional: true
    tuple val(meta), path("results/for-dramv")              , emit: dramv, optional: true
    path "versions.yml"                                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def groups_opt = include_groups ? "--include-groups ${include_groups}" : ''
    def db_opt = db_dir ? "--db-dir ${db_dir}" : ''
    def setup = ''
    if (db_dir) {
        setup = "if [ ! -f ${db_dir}/.vs2_setup.done ]; then virsorter setup -d ${db_dir} -j ${task.cpus} && touch ${db_dir}/.vs2_setup.done; fi; virsorter config --init-source --db-dir=${db_dir}"
    }
    
    """
    ${setup}

    virsorter run \\
        -w results \\
        -i ${contigs} \\
        --min-length ${min_length} \\
        -j ${task.cpus} \\
        ${groups_opt} \\
        ${db_opt} \\
        ${args} \\
        ${mode}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virsorter: 2.2.4
    END_VERSIONS
    """
    stub:
    """
    mkdir -p results/for-dramv
    touch results/final-viral-combined.fa
    touch results/final-viral-score.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virsorter: 2.2.4
    END_VERSIONS
    
    """
}
