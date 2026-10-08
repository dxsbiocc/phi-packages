#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { MUSE_CALL } from '../main.nf'


params.tumor      = null
params.tumor_bai  = null
params.normal     = null
params.normal_bai = null
params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.tumor, checkIfExists: true), file(params.tumor_bai, checkIfExists: true),
                           file(params.normal, checkIfExists: true), file(params.normal_bai, checkIfExists: true), file(params.fasta, checkIfExists: true)])

    MUSE_CALL(in_ch)
}
