#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_GENETICMAPCONVERT } from '../main.nf'


params.map_file   = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result', chr: 'chr21'], file(params.map_file, checkIfExists: true)])

    CUSTOM_GENETICMAPCONVERT(in_ch)
}
