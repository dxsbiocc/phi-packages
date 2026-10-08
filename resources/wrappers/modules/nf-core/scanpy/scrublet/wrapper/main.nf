#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SCANPY_SCRUBLET } from '../main.nf'

params.anndata   = null
params.batch_col = 'sample'
params.outdir    = null

workflow {
    anndata_ch = Channel
        .fromPath(params.anndata, checkIfExists: true)
        .map { h5ad -> [[id: h5ad.simpleName], h5ad] }

    SCANPY_SCRUBLET(anndata_ch, params.batch_col)
}
