#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { SNPCLUSTERING } from '../main.nf'

params.vcf = null
params.vcf_index = null
params.region = null
params.n_pcs = 10
params.use_approx = false
params.algorithm = "kmeans"
params.n_clusters = 3
params.dbscan_eps = 0.5
params.dbscan_min_samples = 5
params.sample_id = "sample"
params.outdir = "results"

workflow {
    vcf_ch = Channel.value([[id: params.sample_id], file(params.vcf, checkIfExists: true), file(params.vcf_index, checkIfExists: true)])

    SNPCLUSTERING(vcf_ch, [], [], params.region ?: null, params.n_pcs, params.use_approx, params.algorithm, params.n_clusters, params.dbscan_eps, params.dbscan_min_samples)
}
