#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CNVNATOR_CNVNATOR } from '../main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_RD } from '../../cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_HIST } from '../../cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_STAT } from '../../cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_PARTITION } from '../../cnvnator/main.nf'

params.bam        = null
params.bai        = null
params.outdir     = null

workflow {
    bam_ch   = Channel.value([[id: 'sample', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    empty2   = Channel.value([[:], []])
    empty3   = Channel.value([[:], [], []])

    CNVNATOR_RD(bam_ch, empty2, empty2, empty2, 'rd')
    CNVNATOR_HIST(empty3, CNVNATOR_RD.out.root, empty2, empty2, 'his')
    CNVNATOR_STAT(empty3, CNVNATOR_HIST.out.root, empty2, empty2, 'stat')
    CNVNATOR_PARTITION(empty3, CNVNATOR_STAT.out.root, empty2, empty2, 'partition')
    CNVNATOR_CNVNATOR(empty3, CNVNATOR_PARTITION.out.root, empty2, empty2, 'call')
}
