process PROTEIN_METRICS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://python:3.12-slim' :
        'python:3.12-slim' }"

    input:
    tuple val(meta), path(structures)
    path(reference)
    val(binder_chains)
    val(target_chains)
    val(contact_cutoff)
    val(hbond_cutoff)
    val(clash_cutoff)

    output:
    tuple val(meta), path("${prefix}.protein_metrics.tsv") , emit: metrics_tsv
    tuple val(meta), path("${prefix}.protein_metrics.json"), emit: metrics_json
    path "versions.yml"                                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ref_opt = reference ? "--reference ${reference}" : ''
    def binder_opt = binder_chains ? "--binder-chains '${binder_chains}'" : ''
    def target_opt = target_chains ? "--target-chains '${target_chains}'" : ''
    """
    python3 "${moduleDir}/bin/compute_protein_metrics.py" \\
        --structures ${structures} \\
        ${ref_opt} \\
        ${binder_opt} \\
        ${target_opt} \\
        --contact-cutoff ${contact_cutoff} \\
        --hbond-cutoff ${hbond_cutoff} \\
        --clash-cutoff ${clash_cutoff} \\
        --out-tsv ${prefix}.protein_metrics.tsv \\
        --out-json ${prefix}.protein_metrics.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        protein_metrics: 1.0.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "design_id\\tstructure_path\\tca_rmsd\\tinterface_contacts\\tclash_count\\n" > ${prefix}.protein_metrics.tsv
    printf "[]\\n" > ${prefix}.protein_metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        protein_metrics: 1.0.0
    END_VERSIONS
    """
}
