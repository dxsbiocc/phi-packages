#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored limma differential-expression
// module at ../main.nf (nf-core/modules `limma/differential`, as used by
// nf-core/differentialabundance for both microarray/intensity data and,
// with use_voom enabled, RNA-seq counts via limma-voom). See
// docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { LIMMA_DIFFERENTIAL } from '../main.nf'

params.counts             = null
params.samplesheet        = null
params.contrast_variable  = null
params.reference_level    = null
params.target_level       = null
params.formula            = null
params.comparison         = null
params.sample_id_col      = 'sample'
params.gene_id_col        = 'gene_id'
params.use_voom           = true
params.prefix             = 'de'
params.outdir             = null

workflow {
    contrast_ch = Channel.of(
        tuple(
            [id: params.prefix, variable: params.contrast_variable],
            params.contrast_variable,
            params.reference_level,
            params.target_level,
            params.formula,
            params.comparison
        )
    )

    matrix_ch = Channel.of(
        tuple(
            [id: params.prefix],
            file(params.samplesheet, checkIfExists: true),
            file(params.counts, checkIfExists: true)
        )
    )

    LIMMA_DIFFERENTIAL(contrast_ch, matrix_ch)
}
