process FOLDSEEK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/foldseek:9.427df8a--pl5321hb365157_1' :
        'quay.io/biocontainers/foldseek:9.427df8a--pl5321hb365157_1' }"

    input:
    tuple val(meta), path(structures)
    path(target_db)
    val(mode)
    val(evalue)
    val(coverage)
    val(tmscore_threshold)
    val(novelty_tm_cutoff)

    output:
    tuple val(meta), path("${prefix}.foldseek.tsv")     , emit: metrics_tsv
    tuple val(meta), path("${prefix}.foldseek.json")    , emit: metrics_json
    tuple val(meta), path("${prefix}.foldseek_raw.tsv") , emit: raw_tsv
    path "versions.yml"                                 , emit: versions, topic: versions

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

    if [ "${mode}" = "search" ]; then
        foldseek easy-search \\
            query_staged \\
            ${target_db} \\
            ${prefix}.foldseek_raw.tsv \\
            tmp_foldseek \\
            --threads ${task.cpus} \\
            -e ${evalue} \\
            --format-output query,target,fident,alnlen,evalue,bits,alntmscore,qtmscore,ttmscore,lddt \\
            ${args}
    elif [ "${mode}" = "cluster" ]; then
        foldseek easy-cluster \\
            query_staged \\
            cluster_out \\
            tmp_foldseek \\
            --threads ${task.cpus} \\
            -c ${coverage} \\
            --tmscore-threshold ${tmscore_threshold} \\
            ${args}
        cp cluster_out_clu.tsv ${prefix}.foldseek_raw.tsv
    else
        echo "Unsupported Foldseek mode: ${mode}" >&2
        exit 1
    fi

    python3 "${moduleDir}/bin/compile_foldseek_outputs.py" \\
        --mode "${mode}" \\
        --input-tsv ${prefix}.foldseek_raw.tsv \\
        --query-dir query_staged \\
        --novelty-tm-cutoff ${novelty_tm_cutoff} \\
        --step-name foldseek \\
        --out-tsv ${prefix}.foldseek.tsv \\
        --out-json ${prefix}.foldseek.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek version 2>/dev/null || echo "9.x")
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.foldseek_raw.tsv
    printf "design_id\\tstep\\tmode\\n" > ${prefix}.foldseek.tsv
    printf "[]\\n" > ${prefix}.foldseek.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: "9.x"
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
