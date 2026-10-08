#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GAWK } from '../main.nf'


params.input      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], [file(params.input, checkIfExists: true)]])

    GAWK(in_ch, [], false)
}
