process GNINA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://gnina/gnina:latest' :
        'gnina/gnina:latest' }"

    input:
    tuple val(meta), path(receptor), path(ligands), path(autobox_ligand)
    val center_x
    val center_y
    val center_z
    val size_x
    val size_y
    val size_z
    val cnn_scoring
    val exhaustiveness
    val num_modes
    val cpu_threads

    output:
    tuple val(meta), path("*.docked.sdf.gz")     , emit: poses
    tuple val(meta), path("*.score_table.tsv")   , emit: score_table
    tuple val(meta), path("*.gnina_summary.json") , emit: summary
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ncpu   = cpu_threads ?: task.cpus
    def box_arg = autobox_ligand ? "--autobox_ligand ${autobox_ligand}" : "--center_x ${center_x} --center_y ${center_y} --center_z ${center_z} --size_x ${size_x} --size_y ${size_y} --size_z ${size_z}"
    """
    gnina \\
        --receptor ${receptor} \\
        --ligand ${ligands} \\
        ${box_arg} \\
        --cnn_scoring ${cnn_scoring} \\
        --exhaustiveness ${exhaustiveness} \\
        --num_modes ${num_modes} \\
        --cpu ${ncpu} \\
        --out ${prefix}.docked.sdf.gz \\
        ${args}

    python3 ${moduleDir}/bin/parse_gnina_sdf.py \\
        --sdf ${prefix}.docked.sdf.gz \\
        --output_tsv ${prefix}.score_table.tsv \\
        --output_json ${prefix}.gnina_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gnina: \$(gnina --version 2>&1 | head -n 1 | sed 's/gnina[[:space:]]*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "" | gzip > ${prefix}.docked.sdf.gz
    printf "ligand_id\\tpose_rank\\tscore\\tminimized_affinity\\tcnn_score\\tcnn_affinity\\tcnn_vs\\tsource_file\\n" > ${prefix}.score_table.tsv
    echo '{"tool":"gnina","n_ligands":0}' > ${prefix}.gnina_summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gnina: 1.3
    END_VERSIONS
    """
}
