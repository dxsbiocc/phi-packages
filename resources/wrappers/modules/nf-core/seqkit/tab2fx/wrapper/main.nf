#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SEQKIT_TAB2FX } from '../main.nf'


params.text       = null
params.out_ext    = 'fa.zst'
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], file(params.text, checkIfExists: true)])

    SEQKIT_TAB2FX(in_ch, params.out_ext)
}
