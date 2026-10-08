#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PARAPHASE } from '../main.nf'
include { HTSLIB_BGZIPTABIX } from '../../htslib/bgziptabix/main.nf'

params.bam        = null
params.bai        = null
params.fasta      = null
params.gene       = 'PRODH'
params.outdir     = null

workflow {
    HTSLIB_BGZIPTABIX(Channel.value([[id: 'test_ref'], file(params.fasta, checkIfExists: true), [], []]), 'decompress', false, 'fa')
    in_ch = Channel.value([[id: 'sample', single_end: true], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])

    PARAPHASE(in_ch, HTSLIB_BGZIPTABIX.out.output, [[:], []])
}
