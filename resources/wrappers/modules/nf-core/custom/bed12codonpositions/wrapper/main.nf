#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_BED12CODONPOSITIONS } from '../main.nf'


params.bed12      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], file(params.bed12, checkIfExists: true)])

    CUSTOM_BED12CODONPOSITIONS(in_ch)
}
