#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BISMARK_SUMMARY } from '../main.nf'
include { BISMARK_GENOMEPREPARATION } from '../../genomepreparation/main.nf'
include { BISMARK_ALIGN } from '../../align/main.nf'
include { BISMARK_DEDUPLICATE } from '../../deduplicate/main.nf'
include { BISMARK_METHYLATIONEXTRACTOR } from '../../methylationextractor/main.nf'

params.reads      = null
params.fasta      = null
params.outdir     = null

workflow {
    reads_ch = Channel.value([[id: 'sample', single_end: true], file(params.reads, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    BISMARK_GENOMEPREPARATION(fasta_ch)
    BISMARK_ALIGN(reads_ch, fasta_ch, BISMARK_GENOMEPREPARATION.out.index)
    BISMARK_DEDUPLICATE(BISMARK_ALIGN.out.bam)
    BISMARK_METHYLATIONEXTRACTOR(BISMARK_ALIGN.out.bam, BISMARK_GENOMEPREPARATION.out.index)
    BISMARK_SUMMARY(BISMARK_ALIGN.out.bam.collect { meta, bam -> bam }, BISMARK_ALIGN.out.report.collect { meta, report -> report },
                    BISMARK_DEDUPLICATE.out.report.collect { meta, report -> report }, BISMARK_METHYLATIONEXTRACTOR.out.report.collect { meta, report -> report },
                    BISMARK_METHYLATIONEXTRACTOR.out.mbias.collect { meta, txt -> txt })
}
