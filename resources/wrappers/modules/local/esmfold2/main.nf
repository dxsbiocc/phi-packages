process ESMFOLD2 {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/rosettacommons/esmfold2:latest' :
        'ghcr.io/rosettacommons/esmfold2:latest' }"

    input:
    tuple val(meta), path(fasta), path(sequences_tsv)
    path(model_dir)
    val(num_recycles)
    val(chunk_size)

    output:
    tuple val(meta), path("esmfold2_structures/*.{pdb,cif}"), emit: structures
    tuple val(meta), path("${prefix}.esmfold2_metrics.tsv") , emit: metrics_tsv
    tuple val(meta), path("${prefix}.esmfold2_metrics.json"), emit: metrics_json
    path "versions.yml"                                     , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def fa_opt = fasta ? "--fasta ${fasta}" : ''
    def tsv_opt = sequences_tsv ? "--sequences-tsv ${sequences_tsv}" : ''
    def model_env = model_dir ? "export ESMFOLD2_MODEL_DIR=\$(readlink -f ${model_dir})" : ''
    def chunk_opt = chunk_size && chunk_size > 0 ? "--chunk-size ${chunk_size}" : ''
    """
    ${model_env}

    python3 "${moduleDir}/bin/prepare_esmfold2_inputs.py" \\
        ${fa_opt} \\
        ${tsv_opt} \\
        --out-fasta staged_inputs.fasta \\
        --out-manifest esmfold2_manifest.json

    ESMFOLD_CLI=\${ESMFOLD2_CLI:-esm-fold}
    \${ESMFOLD_CLI} \\
        -i staged_inputs.fasta \\
        -o esmfold2_raw \\
        --num-recycles ${num_recycles} \\
        ${chunk_opt} \\
        ${args}

    python3 "${moduleDir}/bin/compile_esmfold2_outputs.py" \\
        --out-dir esmfold2_raw \\
        --manifest esmfold2_manifest.json \\
        --step-name esmfold2 \\
        --out-structures-dir esmfold2_structures \\
        --out-tsv ${prefix}.esmfold2_metrics.tsv \\
        --out-json ${prefix}.esmfold2_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esmfold2: "2.0"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p esmfold2_structures
    touch esmfold2_structures/${prefix}_stub.pdb
    printf "design_id\\tparent_design_id\\tstep\\tstructure_path\\tconfidence_score\\tplddt\\tcomplex_plddt\\tptm\\n" > ${prefix}.esmfold2_metrics.tsv
    printf "[]\\n" > ${prefix}.esmfold2_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esmfold2: "2.0"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
