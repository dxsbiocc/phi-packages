#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_MULTIQCCUSTOMBIOTYPE } from '../main.nf'

params.count  = null
params.header = null
params.outdir = null

workflow {
    count_ch = Channel
        .fromPath(params.count, checkIfExists: true)
        .map { count -> [[id: count.simpleName], count] }

    header_ch = Channel.value([[:], file(params.header, checkIfExists: true)])

    CUSTOM_MULTIQCCUSTOMBIOTYPE(count_ch, header_ch)
}
