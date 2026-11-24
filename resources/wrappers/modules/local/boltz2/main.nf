process BOLTZ2 {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/jwohlwend/boltz:latest' :
        'ghcr.io/jwohlwend/boltz:latest' }"

    input:
    tuple val(meta), path(inputs), path(sequences_tsv)
    path(cache_dir)
    val(target_sequence)
    val(ligand_ccd)
    val(ligand_smiles)
    val(predict_affinity)
    val(use_msa_server)
    val(use_potentials)
    val(recycling_steps)
    val(diffusion_samples)
    val(output_format)

    output:
    tuple val(meta), path("boltz2_structures/*.{cif,pdb}"), emit: structures
    tuple val(meta), path("${prefix}.boltz2_metrics.tsv") , emit: metrics_tsv
    tuple val(meta), path("${prefix}.boltz2_metrics.json"), emit: metrics_json
    path "versions.yml"                                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def inputs_opt = inputs ? "--inputs ${inputs}" : ''
    def seq_tsv_opt = sequences_tsv ? "--sequences-tsv ${sequences_tsv}" : ''
    def target_opt = target_sequence ? "--target-sequence '${target_sequence}'" : ''
    def ccd_opt = ligand_ccd ? "--ligand-ccd '${ligand_ccd}'" : ''
    def smiles_opt = ligand_smiles ? "--ligand-smiles '${ligand_smiles}'" : ''
    def aff_opt = predict_affinity ? '--predict-affinity' : ''
    def msa_yaml_mode = use_msa_server ? 'server' : 'empty'
    def cache_opt = cache_dir ? "--cache \$(readlink -f ${cache_dir})" : ''
    def msa_flag = use_msa_server ? '--use_msa_server' : ''
    def pot_flag = use_potentials ? '--use_potentials' : ''
    """
    python3 "${moduleDir}/bin/prepare_boltz2_inputs.py" \\
        ${inputs_opt} \\
        ${seq_tsv_opt} \\
        ${target_opt} \\
        ${ccd_opt} \\
        ${smiles_opt} \\
        ${aff_opt} \\
        --msa-mode ${msa_yaml_mode} \\
        --out-dir boltz_inputs \\
        --out-manifest boltz_manifest.json

    boltz predict boltz_inputs \\
        --out_dir boltz_raw \\
        ${cache_opt} \\
        --recycling_steps ${recycling_steps} \\
        --diffusion_samples ${diffusion_samples} \\
        --output_format ${output_format} \\
        ${msa_flag} \\
        ${pot_flag} \\
        --override \\
        ${args}

    python3 "${moduleDir}/bin/compile_boltz2_outputs.py" \\
        --out-dir boltz_raw \\
        --manifest boltz_manifest.json \\
        --step-name boltz2 \\
        --out-structures-dir boltz2_structures \\
        --out-tsv ${prefix}.boltz2_metrics.tsv \\
        --out-json ${prefix}.boltz2_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        boltz: \$(boltz --version 2>/dev/null || echo "2.x")
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p boltz2_structures
    touch boltz2_structures/${prefix}_stub.${output_format == 'pdb' ? 'pdb' : 'cif'}
    printf "design_id\\tparent_design_id\\tstep\\tstructure_path\\tconfidence_score\\tcomplex_plddt\\tiptm\\n" > ${prefix}.boltz2_metrics.tsv
    printf "[]\\n" > ${prefix}.boltz2_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        boltz: "2.x"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
