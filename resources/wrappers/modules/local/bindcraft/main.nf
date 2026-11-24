process BINDCRAFT {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/martinpacesa/bindcraft:latest' :
        'ghcr.io/martinpacesa/bindcraft:latest' }"

    input:
    tuple val(meta), path(starting_pdb), path(settings_json), path(filters_json), path(advanced_json)
    val(binder_name)
    val(chains)
    val(target_hotspot_residues)
    val(min_length)
    val(max_length)
    val(number_of_final_designs)

    output:
    tuple val(meta), path("bindcraft_structures/*.pdb")     , emit: structures, optional: true
    tuple val(meta), path("${prefix}.bindcraft_designs.tsv"), emit: metrics_tsv
    tuple val(meta), path("${prefix}.bindcraft_designs.json"), emit: metrics_json
    path "versions.yml"                                     , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def pdb_opt = starting_pdb ? "--starting-pdb ${starting_pdb}" : ''
    def set_opt = settings_json ? "--settings-json ${settings_json}" : ''
    def hot_opt = target_hotspot_residues ? "--target-hotspot-residues '${target_hotspot_residues}'" : ''
    def filt_arg = filters_json ? "--filters ${filters_json}" : '--filters "${BINDCRAFT_DIR}/settings_filters/default_filters.json"'
    def adv_arg = advanced_json ? "--advanced ${advanced_json}" : '--advanced "${BINDCRAFT_DIR}/settings_advanced/default_4stage_multimer.json"'
    """
    BINDCRAFT_DIR=\${BINDCRAFT_HOME:-/opt/BindCraft}

    python3 "${moduleDir}/bin/prepare_bindcraft_inputs.py" \\
        ${set_opt} \\
        ${pdb_opt} \\
        --design-path bindcraft_raw \\
        --binder-name '${binder_name}' \\
        --chains '${chains}' \\
        ${hot_opt} \\
        --min-length ${min_length} \\
        --max-length ${max_length} \\
        --num-designs ${number_of_final_designs} \\
        --out-json staged_bindcraft_settings.json

    python3 -u "\${BINDCRAFT_DIR}/bindcraft.py" \\
        --settings staged_bindcraft_settings.json \\
        ${filt_arg} \\
        ${adv_arg} \\
        ${args}

    python3 "${moduleDir}/bin/compile_bindcraft_outputs.py" \\
        --design-path bindcraft_raw \\
        --step-name bindcraft \\
        --out-structures-dir bindcraft_structures \\
        --out-tsv ${prefix}.bindcraft_designs.tsv \\
        --out-json ${prefix}.bindcraft_designs.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bindcraft: "1.x"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p bindcraft_structures
    touch bindcraft_structures/${prefix}_stub.pdb
    printf "design_id\\tparent_design_id\\tstep\\tsequence\\tstructure_path\\tconfidence_score\\tcomplex_plddt\\tiptm\\n" > ${prefix}.bindcraft_designs.tsv
    printf "[]\\n" > ${prefix}.bindcraft_designs.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bindcraft: "1.x"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
