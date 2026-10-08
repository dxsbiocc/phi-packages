#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { MINIMAC4_COMPRESSREF } from '../main.nf'


params.ref        = null
params.ref_csi    = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'input'], file(params.ref, checkIfExists: true), file(params.ref_csi, checkIfExists: true)])

    MINIMAC4_COMPRESSREF(in_ch)
}
