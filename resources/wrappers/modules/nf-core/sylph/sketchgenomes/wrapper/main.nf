#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SYLPH_SKETCHGENOMES } from '../main.nf'


params.fasta1     = null
params.fasta2     = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], [file(params.fasta1, checkIfExists: true), file(params.fasta2, checkIfExists: true)]])

    SYLPH_SKETCHGENOMES(in_ch)
}
