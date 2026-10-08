#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SEQKIT_GREP } from '../main.nf'

params.sequence = null
params.pattern  = null
params.outdir   = null

workflow {
    sequence_ch = Channel
        .fromPath(params.sequence, checkIfExists: true)
        .map { f -> [[id: "${f.simpleName}_grep", single_end: false], f] }
    pattern_ch  = Channel.value([])

    SEQKIT_GREP(sequence_ch, pattern_ch, 'fa')
}
