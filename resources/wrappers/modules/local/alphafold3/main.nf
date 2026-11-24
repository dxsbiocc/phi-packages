process ALPHAFOLD3 {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/google-deepmind/alphafold3:latest' :
        'ghcr.io/google-deepmind/alphafold3:latest' }"

    input:
    tuple val(meta), path(inputs), path(sequences_tsv)
    path(model_dir)
    path(db_dir)
    val(target_sequence)
    val(ligand_ccd)
    val(ligand_smiles)
    val(model_seeds)
    val(run_data_pipeline)
    val(num_diffusion_samples)

    output:
    tuple val(meta), path("af3_structures/*.cif")       , emit: structures
    tuple val(meta), path("${prefix}.af3_metrics.tsv")  , emit: metrics_tsv
    tuple val(meta), path("${prefix}.af3_metrics.json") , emit: metrics_json
    path "versions.yml"                                 , emit: versions, topic: versions

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
    def prep_data_flag = run_data_pipeline ? '--run-data-pipeline' : ''
    def db_opt = db_dir ? "--db_dir=\$(readlink -f ${db_dir})" : ''
    def af3_data_flag = run_data_pipeline ? '--run_data_pipeline' : '--norun_data_pipeline'
    """
    python3 "${moduleDir}/bin/prepare_af3_inputs.py" \\
        ${inputs_opt} \\
        ${seq_tsv_opt} \\
        ${target_opt} \\
        ${ccd_opt} \\
        ${smiles_opt} \\
        --model-seeds '${model_seeds}' \\
        ${prep_data_flag} \\
        --out-dir af3_inputs \\
        --out-manifest af3_manifest.json

    AF3_RUN=\${ALPHAFOLD3_RUN_PY:-/app/alphafold/run_alphafold.py}
    python3 "\${AF3_RUN}" \\
        --input_dir=af3_inputs \\
        --output_dir=af3_raw \\
        --model_dir=\$(readlink -f ${model_dir}) \\
        ${db_opt} \\
        ${af3_data_flag} \\
        --run_inference \\
        --num_diffusion_samples=${num_diffusion_samples} \\
        --force_output_dir \\
        ${args}

    python3 "${moduleDir}/bin/compile_af3_outputs.py" \\
        --out-dir af3_raw \\
        --manifest af3_manifest.json \\
        --step-name alphafold3 \\
        --out-structures-dir af3_structures \\
        --out-tsv ${prefix}.af3_metrics.tsv \\
        --out-json ${prefix}.af3_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alphafold3: "3.0"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p af3_structures
    touch af3_structures/${prefix}_stub.cif
    printf "design_id\\tparent_design_id\\tstep\\tstructure_path\\tconfidence_score\\tranking_score\\tptm\\tiptm\\tplddt\\tcomplex_plddt\\tmin_pae_interaction\\n" > ${prefix}.af3_metrics.tsv
    printf "[]\\n" > ${prefix}.af3_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alphafold3: "3.0"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
