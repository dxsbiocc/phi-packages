#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BWA_SAMSE } from '../main.nf'
include { BWA_INDEX } from '../../index/main.nf'
include { BWA_ALN } from '../../aln/main.nf'

params.reads      = null
params.fasta      = null
params.outdir     = null

workflow {
    reads_ch = Channel.value([[id: 'sample', single_end: true], [file(params.reads, checkIfExists: true)]])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])

    BWA_INDEX(fasta_ch)
    BWA_ALN(reads_ch, BWA_INDEX.out.index)
    BWA_SAMSE(reads_ch.join(BWA_ALN.out.sai), BWA_INDEX.out.index)
}
