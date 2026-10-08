#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_PRINTSVEVIDENCE } from '../main.nf'
include { GATK4_COLLECTSVEVIDENCE } from '../../collectsvevidence/main.nf'

params.bam        = null
params.bai        = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))
    GATK4_COLLECTSVEVIDENCE(Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), [], []]),
                            fasta_ch, fai_ch, dict_ch)
    ev_ch = GATK4_COLLECTSVEVIDENCE.out.paired_end_evidence.combine(GATK4_COLLECTSVEVIDENCE.out.paired_end_evidence_index, by: 0)
                .map { meta, file, index -> [[id: 'printed'], file, index] }.groupTuple()
    GATK4_PRINTSVEVIDENCE(ev_ch, [], fasta_ch, fai_ch, dict_ch)
}
