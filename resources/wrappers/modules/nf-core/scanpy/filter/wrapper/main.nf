#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SCANPY_FILTER } from '../main.nf'

params.anndata             = null
params.min_genes           = 20
params.min_cells           = 20
params.min_counts_gene     = 50
params.min_counts_cell     = 50
params.max_mito_percentage = 50
params.symbol_col          = 'index'
params.outdir              = null

workflow {
    anndata_ch = Channel
        .fromPath(params.anndata, checkIfExists: true)
        .map { h5ad -> [[id: h5ad.simpleName], h5ad] }

    SCANPY_FILTER(
        anndata_ch,
        params.min_genes,
        params.min_cells,
        params.min_counts_gene,
        params.min_counts_cell,
        params.max_mito_percentage,
        params.symbol_col
    )
}
