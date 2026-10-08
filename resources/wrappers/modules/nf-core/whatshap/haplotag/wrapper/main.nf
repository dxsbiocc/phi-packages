#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { WHATSHAP_HAPLOTAG } from '../main.nf'
include { SAMTOOLS_FAIDX } from '../../../samtools/faidx/main.nf'
include { WHATSHAP_PHASE } from '../../phase/main.nf'

params.vcf        = null
params.tbi        = null
params.bam        = null
params.bai        = null
params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true),
                          file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), []])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true), []])
    SAMTOOLS_FAIDX(fasta_ch, false)
    ref_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)]).join(SAMTOOLS_FAIDX.out.fai)

    WHATSHAP_PHASE(in_ch, ref_ch)
    tag_in = Channel.value([[id: 'sample']]).join(WHATSHAP_PHASE.out.vcf).join(WHATSHAP_PHASE.out.tbi)
        .map { meta, vcf, tbi -> [meta, vcf, tbi, file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)] }
    fasta_only_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_only_ch   = Channel.value([[id: 'genome']]).join(SAMTOOLS_FAIDX.out.fai)

    WHATSHAP_HAPLOTAG(tag_in, fasta_only_ch, fai_only_ch, false)
}
