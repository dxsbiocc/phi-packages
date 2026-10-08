#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { DRAGMAP_ALIGN } from '../main.nf'
include { DRAGMAP_HASHTABLE } from '../../hashtable/main.nf'

params.reads1     = null
params.reads2     = null
params.fasta      = null
params.outdir     = null

workflow {
    reads_ch = Channel.value([[id: 'sample', single_end: false], [file(params.reads1, checkIfExists: true), file(params.reads2, checkIfExists: true)]])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])

    DRAGMAP_HASHTABLE(fasta_ch)
    DRAGMAP_ALIGN(reads_ch, DRAGMAP_HASHTABLE.out.hashmap, fasta_ch, false)
}
