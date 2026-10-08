#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SVABA } from '../main.nf'
include { BWA_INDEX } from '../../bwa/index/main.nf'

params.bam        = null
params.bai        = null
params.fasta      = null
params.fai        = null
params.dbsnp      = null
params.dbsnp_tbi  = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), [], []])
    fasta_ch = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'fasta'], file(params.fai, checkIfExists: true)])
    dbsnp_ch = Channel.value([[id: 'dbsnp'], file(params.dbsnp, checkIfExists: true)])
    dbsnp_tbi_ch = Channel.value([[id: 'dbsnp'], file(params.dbsnp_tbi, checkIfExists: true)])
    BWA_INDEX(fasta_ch)

    SVABA(bam_ch, fasta_ch, fai_ch, BWA_INDEX.out.index, dbsnp_ch, dbsnp_tbi_ch, Channel.value([[], []]))
}
