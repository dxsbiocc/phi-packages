#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PLINK2_FILTER } from '../main.nf'


params.pgen       = null
params.pvar       = null
params.psam       = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.pgen, checkIfExists: true), file(params.pvar, checkIfExists: true), file(params.psam, checkIfExists: true)])
    PLINK2_FILTER(in_ch)

}
