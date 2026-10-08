#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BISMARK_COVERAGE2CYTOSINE } from '../main.nf'
include { BISMARK_GENOMEPREPARATION } from '../../genomepreparation/main.nf'
include { BISMARK_METHYLATIONEXTRACTOR } from '../../methylationextractor/main.nf'

params.bam        = null
params.fasta      = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'sample', single_end: false], file(params.bam, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    BISMARK_GENOMEPREPARATION(fasta_ch)
    BISMARK_METHYLATIONEXTRACTOR(bam_ch, BISMARK_GENOMEPREPARATION.out.index)
    BISMARK_COVERAGE2CYTOSINE(BISMARK_METHYLATIONEXTRACTOR.out.coverage, fasta_ch, BISMARK_GENOMEPREPARATION.out.index)
}
