process DESIGN_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://python:3.12-slim' :
        'python:3.12-slim' }"

    input:
    tuple val(meta), path(tables)
    val(id_column)
    val(join_mode)
    val(expr)
    val(rank_by)
    val(rank_direction)
    val(group_by)
    val(max_per_group)
    val(top_n)

    output:
    tuple val(meta), path("${prefix}.merged.tsv")  , emit: merged_tsv
    tuple val(meta), path("${prefix}.filtered.tsv"), emit: filtered_tsv
    tuple val(meta), path("${prefix}.summary.json"), emit: summary_json
    path "versions.yml"                            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def expr_opt = expr ? "--expr \"${expr.replace('"', '\\"')}\"" : ''
    def rank_opt = rank_by ? "--rank-by '${rank_by}' --rank-direction '${rank_direction}'" : ''
    def group_opt = group_by ? "--group-by '${group_by}' --max-per-group ${max_per_group}" : ''
    def top_opt = top_n && top_n > 0 ? "--top-n ${top_n}" : ''
    """
    python3 "${moduleDir}/bin/merge_and_filter_designs.py" \\
        --tables ${tables} \\
        --id-column '${id_column}' \\
        --join-mode '${join_mode}' \\
        ${expr_opt} \\
        ${rank_opt} \\
        ${group_opt} \\
        ${top_opt} \\
        --out-merged ${prefix}.merged.tsv \\
        --out-filtered ${prefix}.filtered.tsv \\
        --out-summary ${prefix}.summary.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        design_filter: 1.0.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "${id_column}\\tpassed_filter\\n" > ${prefix}.merged.tsv
    printf "${id_column}\\tpassed_filter\\n" > ${prefix}.filtered.tsv
    printf "{\\"total_designs\\": 0, \\"passed_filter\\": 0, \\"selected_designs\\": 0}\\n" > ${prefix}.summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        design_filter: 1.0.0
    END_VERSIONS
    """
}
