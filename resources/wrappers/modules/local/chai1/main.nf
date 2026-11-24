process CHAI1 {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/chaidiscovery/chai-lab:latest' :
        'ghcr.io/chaidiscovery/chai-lab:latest' }"

    input:
    tuple val(meta), path(inputs), path(sequences_tsv)
    path(cache_dir)
    val(target_sequence)
    val(ligand_smiles)
    val(use_msa_server)
    val(use_templates_server)
    val(num_trunk_recycles)
    val(num_diffn_samples)
    val(seed)

    output:
    tuple val(meta), path("chai1_structures/*.cif")      , emit: structures
    tuple val(meta), path("${prefix}.chai1_metrics.tsv") , emit: metrics_tsv
    tuple val(meta), path("${prefix}.chai1_metrics.json"), emit: metrics_json
    path "versions.yml"                                  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def inputs_opt = inputs ? "--inputs ${inputs}" : ''
    def seq_tsv_opt = sequences_tsv ? "--sequences-tsv ${sequences_tsv}" : ''
    def target_opt = target_sequence ? "--target-sequence '${target_sequence}'" : ''
    def smiles_opt = ligand_smiles ? "--ligand-smiles '${ligand_smiles}'" : ''
    def cache_env = cache_dir ? "export CHAI_DOWNLOADS_DIR=\$(readlink -f ${cache_dir})" : ''
    def msa_flag = use_msa_server ? '--use-msa-server' : ''
    def tmpl_flag = use_templates_server ? '--use-templates-server' : ''
    """
    ${cache_env}

    python3 "${moduleDir}/bin/prepare_chai1_inputs.py" \\
        ${inputs_opt} \\
        ${seq_tsv_opt} \\
        ${target_opt} \\
        ${smiles_opt} \\
        --out-dir chai1_inputs \\
        --out-manifest chai1_manifest.json

    mkdir -p chai1_raw
    for fa in chai1_inputs/*.fasta; do
        job_name=\$(basename "\${fa}" .fasta)
        job_out="chai1_raw/\${job_name}"
        rm -rf "\${job_out}"
        mkdir -p "\${job_out}"
        chai-lab fold \\
            ${msa_flag} \\
            ${tmpl_flag} \\
            --num-trunk-recycles ${num_trunk_recycles} \\
            --num-diffn-samples ${num_diffn_samples} \\
            --seed ${seed} \\
            ${args} \\
            "\${fa}" \\
            "\${job_out}"
    done

    python3 "${moduleDir}/bin/compile_chai1_outputs.py" \\
        --out-dir chai1_raw \\
        --manifest chai1_manifest.json \\
        --step-name chai1 \\
        --out-structures-dir chai1_structures \\
        --out-tsv ${prefix}.chai1_metrics.tsv \\
        --out-json ${prefix}.chai1_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chai-lab: \$(python3 -c "import importlib.metadata; print(importlib.metadata.version('chai_lab'))" 2>/dev/null || echo "0.x")
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p chai1_structures
    touch chai1_structures/${prefix}_stub.cif
    printf "design_id\\tparent_design_id\\tstep\\tstructure_path\\tconfidence_score\\taggregate_score\\tptm\\tiptm\\n" > ${prefix}.chai1_metrics.tsv
    printf "[]\\n" > ${prefix}.chai1_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chai-lab: "0.x"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
