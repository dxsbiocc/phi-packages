#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CELLBENDER_REMOVEBACKGROUND } from '../main.nf'

params.anndata = null
params.epochs  = 5
params.outdir  = null

workflow {
    anndata_ch = Channel
        .fromPath(params.anndata, checkIfExists: true)
        .map { h5ad -> [[id: h5ad.simpleName], h5ad] }

    CELLBENDER_REMOVEBACKGROUND(anndata_ch)
}
