#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SCANPY_HASHSOLO } from '../main.nf'

params.anndata              = null
params.cell_hashing_columns = '0,1,2,3,4,5,6,7,8,9'
params.outdir                = null

workflow {
    columns = params.cell_hashing_columns.split(',') as List

    data_ch = Channel
        .fromPath(params.anndata, checkIfExists: true)
        .map { h5ad -> [[id: h5ad.simpleName], h5ad, columns] }

    SCANPY_HASHSOLO(data_ch)
}
