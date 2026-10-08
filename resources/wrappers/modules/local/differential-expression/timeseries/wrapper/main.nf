#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored dream module at ../main.nf
// (nf-core/modules `variancepartition/dream`, as used by
// nf-core/differentialabundance for repeated-measures / mixed-model designs
// — e.g. time-course data with samples nested within subject, via a
// formula with a random-effect term such as '~ time + (1 | subject)'). See
// docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { VARIANCEPARTITION_DREAM } from '../main.nf'

params.counts             = null
params.samplesheet        = null
params.contrast_variable  = null
params.reference_level    = null
params.target_level       = null
params.formula            = null
params.comparison         = null
params.sample_id_col      = 'sample'
params.apply_voom         = true
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

    VARIANCEPARTITION_DREAM(contrast_ch, matrix_ch)
}
