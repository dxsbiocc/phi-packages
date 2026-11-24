process FOLDSEEK_EASYCLUSTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/foldseek:9.427df8a--pl5321hb365157_1' :
        'quay.io/biocontainers/foldseek:9.427df8a--pl5321hb365157_1' }"

    input:
    tuple val(meta), path(structures)
    val(coverage)
    val(tmscore_threshold)

    output:
    tuple val(meta), path("${prefix}.foldseek_cluster.tsv")     , emit: metrics_tsv
    tuple val(meta), path("${prefix}.foldseek_cluster.json")    , emit: metrics_json
    tuple val(meta), path("${prefix}.foldseek_cluster_raw.tsv") , emit: raw_tsv
    path "versions.yml"                                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p query_staged tmp_foldseek
    for f in ${structures}; do
        if [ -d "\$f" ]; then
            cp -r "\$f"/* query_staged/ 2>/dev/null || true
        else
            cp "\$f" query_staged/
        fi
    done

    foldseek easy-cluster \\
        query_staged \\
        cluster_out \\
        tmp_foldseek \\
        --threads ${task.cpus} \\
        -c ${coverage} \\
        --tmscore-threshold ${tmscore_threshold} \\
        ${args}
    cp cluster_out_clu.tsv ${prefix}.foldseek_cluster_raw.tsv

    python3 "${moduleDir}/bin/compile_foldseek_cluster.py" \\
        --input-tsv ${prefix}.foldseek_cluster_raw.tsv \\
        --step-name foldseek_easycluster \\
        --out-tsv ${prefix}.foldseek_cluster.tsv \\
        --out-json ${prefix}.foldseek_cluster.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek version 2>/dev/null || echo "9.x")
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.foldseek_cluster_raw.tsv
    printf "design_id\\tstep\\tcluster_id\\tcluster_representative\\tis_representative\\tcluster_size\\n" > ${prefix}.foldseek_cluster.tsv
    printf "[]\\n" > ${prefix}.foldseek_cluster.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: "9.x"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
