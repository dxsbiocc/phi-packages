process FPOCKET {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fpocket:4.0--he881be0_2' :
        'quay.io/biocontainers/fpocket:4.0--he881be0_2' }"

    input:
    tuple val(meta), path(receptor_pdb)
    val min_alpha_spheres
    val box_padding

    output:
    tuple val(meta), path("*_out")                  , emit: out_dir
    tuple val(meta), path("*.pockets.tsv")          , emit: pockets_table
    tuple val(meta), path("*.top_pocket_box.json")  , emit: top_pocket_box
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fpocket \\
        -f ${receptor_pdb} \\
        -i ${min_alpha_spheres} \\
        ${args}

    stem=\$(basename ${receptor_pdb} .pdb)
    out_dir="\${stem}_out"
    if [ ! -d "\$out_dir" ]; then
        out_dir=\$(ls -d *_out | head -n 1)
    fi

    python3 ${moduleDir}/bin/parse_fpocket_info.py \\
        --fpocket_out "\$out_dir" \\
        --padding ${box_padding} \\
        --output_tsv ${prefix}.pockets.tsv \\
        --output_json ${prefix}.top_pocket_box.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fpocket: \$(fpocket --version 2>&1 | head -n 1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_out/pockets
    printf "pocket_id\\tscore\\tdruggability_score\\tn_alpha_spheres\\ttotal_sasa\\tvolume\\tcenter_x\\tcenter_y\\tcenter_z\\tsize_x\\tsize_y\\tsize_z\\n" > ${prefix}.pockets.tsv
    echo '{"tool":"fpocket","top_pocket":null}' > ${prefix}.top_pocket_box.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fpocket: 4.0
    END_VERSIONS
    """
}
