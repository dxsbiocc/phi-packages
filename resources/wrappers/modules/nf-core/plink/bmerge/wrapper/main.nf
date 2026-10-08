#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PLINK_BMERGE } from '../main.nf'


params.bed        = null
params.bim        = null
params.fam        = null
params.bed2       = null
params.bim2       = null
params.fam2       = null
params.outdir     = null

workflow {
    in0 = Channel.value([[id: 'sample'], file(params.bed, checkIfExists: true), file(params.bim, checkIfExists: true), file(params.fam, checkIfExists: true)])
    in1 = Channel.value([[id: 'sample'], file(params.bed2, checkIfExists: true), file(params.bim2, checkIfExists: true), file(params.fam2, checkIfExists: true)])

    PLINK_BMERGE(in0, in1)
}
