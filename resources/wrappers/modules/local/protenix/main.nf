process PROTENIX {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/bytedance/protenix:latest' :
        'ghcr.io/bytedance/protenix:latest' }"

    input:
    tuple val(meta), path(inputs), path(sequences_tsv)
    path(cache_dir)
    val(model_name)
    val(target_sequence)
    val(ligand_ccd)
    val(ligand_smiles)
    val(seeds)
    val(recycling_steps)
    val(diffusion_steps)
    val(diffusion_samples)
    val(use_msa)
    val(dtype)

    output:
    tuple val(meta), path("protenix_structures/*.cif")      , emit: structures
    tuple val(meta), path("${prefix}.protenix_metrics.tsv") , emit: metrics_tsv
    tuple val(meta), path("${prefix}.protenix_metrics.json"), emit: metrics_json
    path "versions.yml"                                     , emit: versions, topic: versions

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
    def cache_env = cache_dir ? "export PROTENIX_ROOT_DIR=\$(readlink -f ${cache_dir})" : ''
    def msa_bool = use_msa ? 'true' : 'false'
    """
    ${cache_env}

    python3 "${moduleDir}/bin/prepare_protenix_inputs.py" \\
        ${inputs_opt} \\
        ${seq_tsv_opt} \\
        ${target_opt} \\
        ${ccd_opt} \\
        ${smiles_opt} \\
        --out-json staged_protenix.json \\
        --out-manifest protenix_manifest.json

    protenix pred \\
        --input staged_protenix.json \\
        --out_dir protenix_raw \\
        --model_name "${model_name}" \\
        --seeds "${seeds}" \\
        --cycle ${recycling_steps} \\
        --step ${diffusion_steps} \\
        --sample ${diffusion_samples} \\
        --dtype "${dtype}" \\
        --use_msa ${msa_bool} \\
        ${args}

    python3 "${moduleDir}/bin/compile_protenix_outputs.py" \\
        --out-dir protenix_raw \\
        --manifest protenix_manifest.json \\
        --step-name protenix \\
        --out-structures-dir protenix_structures \\
        --out-tsv ${prefix}.protenix_metrics.tsv \\
        --out-json ${prefix}.protenix_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        protenix: \$(protenix --version 2>/dev/null || echo "${model_name}")
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p protenix_structures
    touch protenix_structures/${prefix}_stub.cif
    printf "design_id\\tparent_design_id\\tstep\\tstructure_path\\tconfidence_score\\tranking_score\\tptm\\tiptm\\tplddt\\tcomplex_plddt\\n" > ${prefix}.protenix_metrics.tsv
    printf "[]\\n" > ${prefix}.protenix_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        protenix: "${model_name}"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
