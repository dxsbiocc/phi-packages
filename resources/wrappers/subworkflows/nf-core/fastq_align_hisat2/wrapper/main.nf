#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored fastq_align_hisat2 subworkflow
// at ../main.nf. A pre-built HISAT2 index (`index`) is used as-is; without
// one, hisat2/build builds a splice-aware index from `fasta` + `gtf` first.
// Splice sites always come from `gtf`.
nextflow.enable.dsl = 2

include { HISAT2_EXTRACTSPLICESITES } from '../../../../modules/nf-core/hisat2/extractsplicesites/main.nf'
include { HISAT2_BUILD              } from '../../../../modules/nf-core/hisat2/build/main.nf'
include { FASTQ_ALIGN_HISAT2        } from '../main.nf'

params.reads_1 = null
params.reads_2 = null
params.index   = null
params.fasta   = null
params.gtf     = null
params.outdir  = null

workflow {
    if (!params.index && !params.fasta) {
        error "fastq-align-hisat2 needs either `index` (a pre-built HISAT2 index directory) or `fasta` to build one."
    }

    ch_gtf = Channel.fromPath(params.gtf, checkIfExists: true)
    HISAT2_EXTRACTSPLICESITES(ch_gtf.map { gtf -> [[id: gtf.baseName], gtf] })

    if (params.index) {
        ch_index = Channel.value([[id: 'hisat2_index'], file(params.index, checkIfExists: true)])
    } else {
        ch_build_input = Channel.fromPath(params.fasta, checkIfExists: true)
            .combine(ch_gtf)
            .combine(HISAT2_EXTRACTSPLICESITES.out.txt.map { _meta, txt -> txt })
            .map { fasta, gtf, splicesites -> [[id: fasta.baseName], fasta, gtf, splicesites] }
        HISAT2_BUILD(ch_build_input, '1.GB')
        ch_index = HISAT2_BUILD.out.index
    }

    def single_end = !params.reads_2
    def read1      = file(params.reads_1, checkIfExists: true)
    def reads      = [read1]
    if (!single_end) reads << file(params.reads_2, checkIfExists: true)
    def sample_id  = single_end ? read1.simpleName : read1.simpleName.replaceAll(/[._-]R?1$/, '')
    reads_ch = channel.of([[id: sample_id, single_end: single_end], reads])

    ch_fasta_fai = channel.value([[:], [], []])

    FASTQ_ALIGN_HISAT2(
        reads_ch,
        ch_index,
        HISAT2_EXTRACTSPLICESITES.out.txt,
        ch_fasta_fai,
        false
    )
}
