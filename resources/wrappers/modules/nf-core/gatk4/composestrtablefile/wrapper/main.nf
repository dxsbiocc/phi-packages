#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_COMPOSESTRTABLEFILE } from '../main.nf'


params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))

    GATK4_COMPOSESTRTABLEFILE(fasta_ch, fai_ch, dict_ch)
}
