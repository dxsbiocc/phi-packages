process AUTODOCK_VINA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/autodock-vina:1.2.5--h4ac6f70_1' :
        'quay.io/biocontainers/autodock-vina:1.2.5--h4ac6f70_1' }"

    input:
    tuple val(meta), path(receptor_pdbqt), path(ligands_dir)
    val center_x
    val center_y
    val center_z
    val size_x
    val size_y
    val size_z
    val exhaustiveness
    val num_modes
    val cpu_per_ligand

    output:
    tuple val(meta), path("*_docked/*.pdbqt")    , emit: poses
    tuple val(meta), path("*.score_table.tsv")   , emit: score_table
    tuple val(meta), path("*.vina_summary.json") , emit: summary
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ncpu   = cpu_per_ligand ?: Math.min(task.cpus as int, 4)
    """
    mkdir -p ${prefix}_docked
    shopt -s nullglob
    for lig in ${ligands_dir}/*.pdbqt; do
        bname=\$(basename "\$lig" .pdbqt)
        vina \\
            --receptor ${receptor_pdbqt} \\
            --ligand "\$lig" \\
            --center_x ${center_x} \\
            --center_y ${center_y} \\
            --center_z ${center_z} \\
            --size_x ${size_x} \\
            --size_y ${size_y} \\
            --size_z ${size_z} \\
            --exhaustiveness ${exhaustiveness} \\
            --num_modes ${num_modes} \\
            --cpu ${ncpu} \\
            --out "${prefix}_docked/\${bname}_out.pdbqt" \\
            ${args} \\
            > "${prefix}_docked/\${bname}.vina.log" 2>&1
    done

    python3 ${moduleDir}/bin/parse_vina_results.py \\
        --docked_dir "${prefix}_docked" \\
        --output_tsv "${prefix}.score_table.tsv" \\
        --output_json "${prefix}.vina_summary.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autodock-vina: \$(vina --version 2>&1 | sed -n 's/.*Vina //p' | head -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_docked
    touch ${prefix}_docked/stub_ligand_out.pdbqt
    printf "ligand_id\\tpose_rank\\tscore\\taffinity_kcal_mol\\trmsd_lb\\trmsd_ub\\tpose_file\\n" > ${prefix}.score_table.tsv
    echo '{"tool":"autodock_vina","n_ligands":0}' > ${prefix}.vina_summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autodock-vina: 1.2.5
    END_VERSIONS
    """
}
