#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BWAMEME_MEM } from '../main.nf'
include { BWAMEME_INDEX } from '../../index/main.nf'

params.reads1     = null
params.reads2     = null
params.fasta      = null
params.outdir     = null

workflow {
    reads_ch = Channel.value([[id: 'sample', single_end: false], [file(params.reads1, checkIfExists: true), file(params.reads2, checkIfExists: true)]])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])

    BWAMEME_INDEX(fasta_ch)
    BWAMEME_MEM(reads_ch, BWAMEME_INDEX.out.index, fasta_ch, false)
}
