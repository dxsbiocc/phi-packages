#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { TRGT_GENOTYPE } from '../main.nf'
include { GUNZIP } from '../../../gunzip/main.nf'
include { SAMTOOLS_FAIDX } from '../../../samtools/faidx/main.nf'

params.bam        = null
params.bai        = null
params.karyotype  = 'XX'
params.fasta_gz   = null
params.repeats    = 'chr22\\t18890357\\t18890451\\tID=TEST;MOTIFS=AT;STRUC=(AT)n'
params.outdir     = null

workflow {
    GUNZIP(Channel.value([[id: 'chr22'], file(params.fasta_gz, checkIfExists: true)]))
    SAMTOOLS_FAIDX(GUNZIP.out.gunzip.combine(Channel.of([[]])), false)
    bam_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), params.karyotype])
    repeats_ch = Channel.of(params.repeats).collectFile(name: 'repeats.bed', newLine: false).map { file -> [[id: 'chr22'], file] }

    TRGT_GENOTYPE(bam_ch, GUNZIP.out.gunzip, SAMTOOLS_FAIDX.out.fai, repeats_ch)
}
