process LIGANDMPNN {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/rosettacommons/ligandmpnn:latest' :
        'ghcr.io/rosettacommons/ligandmpnn:latest' }"

    input:
    tuple val(meta), path(structures)
    path(weights_dir)
    val(model_type)
    val(num_seq_per_pdb)
    val(batch_size)
    val(temperature)
    val(seed)
    val(output_mode)
    val(fixed_residues)
    val(chains_to_design)

    output:
    tuple val(meta), path("${prefix}.sequences.tsv")   , emit: sequences_tsv
    tuple val(meta), path("${prefix}.sequences.fasta") , emit: sequences_fasta
    tuple val(meta), path("${prefix}.mpnn_metrics.json"), emit: metrics_json
    tuple val(meta), path("mpnn_raw/packed/*.pdb")     , emit: packed_structures, optional: true
    path "versions.yml"                                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def n_batches = Math.max(1, (int) Math.ceil((num_seq_per_pdb as double) / (batch_size as double)))
    def pack_opt = output_mode == 'structures' ? '--pack_side_chains 1 --number_of_packs_per_design 1' : '--pack_side_chains 0'
    def fixed_opt = fixed_residues ? "--fixed_residues '${fixed_residues}'" : ''
    def chains_opt = chains_to_design ? "--chains_to_design '${chains_to_design}'" : ''
    def ckpt_env = weights_dir ? "export LIGANDMPNN_WEIGHTS_DIR=\$(readlink -f ${weights_dir})" : ''
    """
    ${ckpt_env}

    python3 "${moduleDir}/bin/prepare_ligandmpnn_inputs.py" \\
        --structures ${structures} \\
        --out-multi-json pdb_path_multi.json \\
        --out-manifest parent_manifest.json

    MPNN_RUN=\${LIGANDMPNN_RUN_PY:-/opt/LigandMPNN/run.py}
    python3 "\${MPNN_RUN}" \\
        --model_type "${model_type}" \\
        --seed ${seed} \\
        --pdb_path_multi pdb_path_multi.json \\
        --out_folder mpnn_raw \\
        --batch_size ${batch_size} \\
        --number_of_batches ${n_batches} \\
        --temperature ${temperature} \\
        ${pack_opt} \\
        ${fixed_opt} \\
        ${chains_opt} \\
        ${args}

    python3 "${moduleDir}/bin/compile_ligandmpnn_outputs.py" \\
        --out-dir mpnn_raw \\
        --manifest parent_manifest.json \\
        --model-type "${model_type}" \\
        --step-name ligandmpnn \\
        --out-tsv ${prefix}.sequences.tsv \\
        --out-fasta ${prefix}.sequences.fasta \\
        --out-json ${prefix}.mpnn_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ligandmpnn: "1.0.0"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "design_id\\tparent_design_id\\tstep\\tmodel_type\\tsequence\\n" > ${prefix}.sequences.tsv
    printf ">stub_mpnn_1\\nACDEFGHIKLMNPQRSTVWY\\n" > ${prefix}.sequences.fasta
    printf "[]\\n" > ${prefix}.mpnn_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ligandmpnn: "1.0.0"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
