#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { DEEPTOOLS_PLOTFINGERPRINT } from '../main.nf'

params.bam    = null
params.bai    = null
params.outdir = null
// Referenced directly by the vendored module's script block (not passed
// as a process input) — must be defined here or it evaluates to null.
params.fragment_size = 0

workflow {
    input_ch = Channel.value([[id: 'test', single_end: false], [file(params.bam, checkIfExists: true)], [file(params.bai, checkIfExists: true)]])

    DEEPTOOLS_PLOTFINGERPRINT(input_ch)
}
