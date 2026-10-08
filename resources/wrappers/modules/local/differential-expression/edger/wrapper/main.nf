#!/usr/bin/env nextflow
// Thin agent-facing adapter over the EDGER_DIFFERENTIAL process at
// ../main.nf. See docs/design/phi-wrapper-agent-composition-design.md
// section 1.
nextflow.enable.dsl = 2

include { EDGER_DIFFERENTIAL } from '../main.nf'

params.counts            = null
params.samplesheet       = null
params.contrast_variable = 'null'
params.reference_level   = 'null'
params.target_level      = 'null'
params.formula           = 'null'
params.comparison        = 'null'
params.sample_id_col     = 'sample'
params.gene_id_col       = 'gene_id'
params.prefix             = 'de'
params.outdir             = null

workflow {
    counts_ch      = Channel.fromPath(params.counts, checkIfExists: true)
    samplesheet_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)

    EDGER_DIFFERENTIAL(counts_ch, samplesheet_ch)
}
