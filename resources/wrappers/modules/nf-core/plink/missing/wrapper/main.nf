#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PLINK_MISSING } from '../main.nf'


params.bed        = null
params.bim        = null
params.fam        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.bed, checkIfExists: true), file(params.bim, checkIfExists: true), file(params.fam, checkIfExists: true)])

    PLINK_MISSING(in_ch)
}
