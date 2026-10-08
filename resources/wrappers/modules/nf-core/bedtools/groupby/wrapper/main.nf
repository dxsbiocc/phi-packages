#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEDTOOLS_GROUPBY } from '../main.nf'


params.bed        = null
params.summary_col = 5
params.outdir     = null

workflow {
    bed_ch = Channel.value([[id: 'result'], file(params.bed, checkIfExists: true)])

    BEDTOOLS_GROUPBY(bed_ch, params.summary_col)
}
