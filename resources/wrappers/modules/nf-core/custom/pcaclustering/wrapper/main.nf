#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_PCACLUSTERING } from '../main.nf'


params.features   = null
params.algorithm  = 'kmeans'
params.n_clusters = 3
params.dbscan_eps = '0.5'
params.dbscan_min_samples = 5
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], file(params.features, checkIfExists: true)])

    CUSTOM_PCACLUSTERING(in_ch, params.algorithm, params.n_clusters, params.dbscan_eps, params.dbscan_min_samples)
}
