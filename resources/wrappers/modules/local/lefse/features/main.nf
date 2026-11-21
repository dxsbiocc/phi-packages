process LEFSE_PLOT_FEATURES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/lefse:1.1.2--pyhdfd78af_0' :
        'quay.io/biocontainers/lefse:1.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(dataset), path(res)

    output:
    tuple val(meta), path("feature_plots/*.png"), emit: plots
    path "versions.yml"                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    mkdir -p feature_plots

    lefse_plot_features.py \\
        ${dataset} \\
        ${res} \\
        feature_plots \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lefse: 1.1.2
    END_VERSIONS
    """
    stub:
    """
    mkdir -p feature_plots
    touch feature_plots/feature.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lefse: 1.1.2
    END_VERSIONS
    
    """
}
