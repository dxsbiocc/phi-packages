#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { HMMCOPY_MAPCOUNTER } from '../main.nf'
include { HMMCOPY_GENERATEMAP } from '../../generatemap/main.nf'

params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.fasta, checkIfExists: true)])
    HMMCOPY_GENERATEMAP(in_ch)
    HMMCOPY_MAPCOUNTER(HMMCOPY_GENERATEMAP.out.bigwig)
}
