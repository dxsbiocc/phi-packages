#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_GERMLINECNVCALLER } from '../main.nf'
include { GATK4_COLLECTREADCOUNTS } from '../../collectreadcounts/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../../determinegermlinecontigploidy/main.nf'
include { GATK4_BEDTOINTERVALLIST } from '../../bedtointervallist/main.nf'

params.bam        = null
params.bai        = null
params.bam2       = null
params.bai2       = null
params.intervals  = null
params.ploidy_priors = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])
    bed_ch = Channel.value(file(params.intervals, checkIfExists: true))
    samples = Channel.of(
        [[id: 'test', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)],
        [[id: 'test2', single_end: false], file(params.bam2, checkIfExists: true), file(params.bai2, checkIfExists: true)])
    GATK4_COLLECTREADCOUNTS(samples.combine(bed_ch), fasta_ch, fai_ch, dict_ch)
    counts_ch = GATK4_COLLECTREADCOUNTS.out.tsv.map { meta, tsv -> [[id: 'cohort'], tsv] }.groupTuple()
    GATK4_DETERMINEGERMLINECONTIGPLOIDY(counts_ch.combine(bed_ch).map { meta, counts, bed -> [meta, counts, bed, []] },
                                        [[], []], file(params.ploidy_priors, checkIfExists: true))
    GATK4_BEDTOINTERVALLIST(bed_ch.map { bed -> [[id: 'cohort'], bed] }, dict_ch.map { meta, dict -> [meta, dict] })
    GATK4_GERMLINECNVCALLER(counts_ch.combine(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls)
                                     .combine(GATK4_BEDTOINTERVALLIST.out.interval_list)
                                     .map { meta, counts, meta2, calls, meta3, regions -> [meta, counts, regions, calls, []] })
}
