#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { MUSE_SUMP } from '../main.nf'


params.call_txt   = null
params.dbsnp      = null
params.dbsnp_tbi  = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.call_txt, checkIfExists: true), file(params.dbsnp, checkIfExists: true), file(params.dbsnp_tbi, checkIfExists: true)])

    MUSE_SUMP(in_ch)
}
