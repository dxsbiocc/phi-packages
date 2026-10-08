#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { ASCAT } from '../main.nf'


params.normal     = null
params.normal_bai = null
params.tumor      = null
params.tumor_bai  = null
params.alleles    = null
params.loci       = null
params.gc         = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', single_end: false], file(params.normal, checkIfExists: true), file(params.normal_bai, checkIfExists: true),
                           file(params.tumor, checkIfExists: true), file(params.tumor_bai, checkIfExists: true)])

    ASCAT(in_ch, [file(params.alleles, checkIfExists: true)], [file(params.loci, checkIfExists: true)], [], [], [file(params.gc, checkIfExists: true)], [])
}
