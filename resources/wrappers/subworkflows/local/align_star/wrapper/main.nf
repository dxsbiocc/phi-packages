#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored align_star subworkflow at
// ../main.nf. A pre-built STAR index (`index`) is used as-is; without one,
// star/genomegenerate builds it from `fasta` + `gtf` first.
nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../../../../modules/nf-core/star/genomegenerate/main.nf'
include { ALIGN_STAR } from '../main.nf'

params.reads_1              = null
params.reads_2              = null
params.index                = null
params.fasta                = null
params.gtf                  = null
params.star_ignore_sjdbgtf  = false
params.skip_markduplicates  = false
params.outdir               = null

workflow {
    if (!params.index && !params.fasta) {
        error "align-star needs either `index` (a pre-built STAR index directory) or `fasta` to build one."
    }

    gtf_ch = Channel.value([[id: 'genome'], [file(params.gtf, checkIfExists: true)]])

    if (params.index) {
        index_ch = Channel.value([[id: 'genome'], file(params.index, checkIfExists: true)])
    } else {
        STAR_GENOMEGENERATE(Channel.value([[id: 'genome'], [file(params.fasta, checkIfExists: true)]]), gtf_ch)
        index_ch = STAR_GENOMEGENERATE.out.index
    }

    def read1 = file(params.reads_1, checkIfExists: true)
    reads_ch = Channel.of([
        [id: read1.simpleName.replaceAll(/[._-]R?1$/, ''), single_end: false],
        [read1, file(params.reads_2, checkIfExists: true)]
    ])

    // samtools stats only uses the FASTA as an optional reference.
    fasta_fai_ch = params.fasta
        ? Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true), []])
        : Channel.value([[id: 'genome'], [], []])

    ALIGN_STAR(
        reads_ch,
        index_ch,
        gtf_ch,
        params.star_ignore_sjdbgtf,
        fasta_fai_ch,
        params.skip_markduplicates
    )
}
