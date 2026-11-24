process RFD3 {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/rosettacommons/foundry:latest' :
        'ghcr.io/rosettacommons/foundry:latest' }"

    input:
    tuple val(meta), path(inputs_json), path(input_structure)
    path(weights_dir)
    val(contig)
    val(ligand)
    val(n_batches)
    val(diffusion_batch_size)
    val(num_timesteps)
    val(low_memory_mode)
    val(chunk_id)

    output:
    tuple val(meta), path("rfd3_raw/**/*.{cif,cif.gz,pdb}") , emit: structures
    tuple val(meta), path("${prefix}.rfd3_designs.tsv")     , emit: metrics_tsv
    tuple val(meta), path("${prefix}.rfd3_designs.json")    , emit: metrics_json
    path "versions.yml"                                     , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def json_opt = inputs_json ? "--inputs-json ${inputs_json}" : ''
    def struct_opt = input_structure ? "--input-structure ${input_structure}" : ''
    def contig_opt = contig ? "--contig '${contig}'" : ''
    def ligand_opt = ligand ? "--ligand '${ligand}'" : ''
    def ckpt_env = weights_dir ? "export FOUNDRY_CHECKPOINT_DIRS=\$(readlink -f ${weights_dir})" : ''
    def low_mem = low_memory_mode ? 'true' : 'false'
    """
    ${ckpt_env}

    python3 "${moduleDir}/bin/prepare_rfd3_inputs.py" \\
        ${json_opt} \\
        ${struct_opt} \\
        ${contig_opt} \\
        ${ligand_opt} \\
        --design-name '${prefix}' \\
        --chunk-id ${chunk_id} \\
        --out-json staged_rfd3_inputs.json

    rfd3 design \\
        out_dir=rfd3_raw \\
        inputs=staged_rfd3_inputs.json \\
        n_batches=${n_batches} \\
        diffusion_batch_size=${diffusion_batch_size} \\
        inference_sampler.num_timesteps=${num_timesteps} \\
        low_memory_mode=${low_mem} \\
        skip_existing=false \\
        ${args}

    python3 "${moduleDir}/bin/compile_rfd3_outputs.py" \\
        --out-dir rfd3_raw \\
        --step-name rfd3 \\
        --out-tsv ${prefix}.rfd3_designs.tsv \\
        --out-json ${prefix}.rfd3_designs.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rfd3: \$(rfd3 --version 2>/dev/null || echo "foundry")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p rfd3_raw
    touch rfd3_raw/${prefix}_0_model_0.cif
    printf "design_id\\tparent_design_id\\tstep\\tstructure_path\\n" > ${prefix}.rfd3_designs.tsv
    printf "[]\\n" > ${prefix}.rfd3_designs.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rfd3: foundry
    END_VERSIONS
    """
}
