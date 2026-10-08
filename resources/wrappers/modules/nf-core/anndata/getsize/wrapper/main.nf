#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { ANNDATA_GETSIZE } from '../main.nf'

params.anndata   = null
params.size_type = 'cells'
params.outdir    = null

workflow {
    anndata_ch = Channel
        .fromPath(params.anndata, checkIfExists: true)
        .map { h5ad -> [[id: h5ad.simpleName], h5ad] }

    ANNDATA_GETSIZE(anndata_ch, params.size_type)
}
