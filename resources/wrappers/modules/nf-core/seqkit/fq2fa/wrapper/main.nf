#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SEQKIT_FQ2FA } from '../main.nf'


params.fastq      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result', single_end: false], file(params.fastq, checkIfExists: true)])

    SEQKIT_FQ2FA(in_ch)
}
