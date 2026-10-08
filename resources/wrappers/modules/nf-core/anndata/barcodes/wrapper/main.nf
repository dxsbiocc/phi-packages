#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { ANNDATA_BARCODES } from '../main.nf'

params.anndata  = null
params.barcodes = null
params.outdir   = null

workflow {
    // meta.id must not equal the input h5ad's own basename: the module
    // writes its output as "${meta.id}.h5ad" next to the input, and Nextflow
    // excludes anything matching an input filename from the output glob, so
    // a matching name makes the run fail with "Missing output file(s)".
    anndata_ch = Channel
        .fromPath(params.anndata, checkIfExists: true)
        .map { h5ad -> [[id: "${h5ad.simpleName}_subset"], h5ad] }

    barcodes_ch = Channel.fromPath(params.barcodes, checkIfExists: true)

    ANNDATA_BARCODES(anndata_ch.combine(barcodes_ch))
}
