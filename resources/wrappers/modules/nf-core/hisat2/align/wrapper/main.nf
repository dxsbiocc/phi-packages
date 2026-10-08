#!/usr/bin/env nextflow
// Thin agent-facing adapter over ../main.nf (HISAT2_ALIGN). A pre-built HISAT2
// index (`index`) is used as-is; without one, ../../build/main.nf builds a
// splice-aware index from `fasta` + `gtf` first. Splice sites always come from
// `gtf` (../../extractsplicesites/main.nf) — cheap, and HISAT2_ALIGN needs them
// either way.
nextflow.enable.dsl = 2

include { HISAT2_EXTRACTSPLICESITES } from '../../extractsplicesites/main.nf'
include { HISAT2_BUILD              } from '../../build/main.nf'
include { HISAT2_ALIGN              } from '../main.nf'

params.reads_1 = null
params.reads_2 = null
params.index   = null
params.fasta   = null
params.gtf     = null
params.outdir  = null

workflow {
    if (!params.index && !params.fasta) {
        error "hisat2-align needs either `index` (a pre-built HISAT2 index directory) or `fasta` to build one."
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
        ch_index = HISAT2_BUILD.out.index.first()
    }

    def read1 = file(params.reads_1, checkIfExists: true)
    reads_ch = channel.of([
        [id: read1.simpleName.replaceAll(/[._-]R?1$/, ''), single_end: false],
        [read1, file(params.reads_2, checkIfExists: true)]
    ])

    HISAT2_ALIGN(reads_ch, ch_index, HISAT2_EXTRACTSPLICESITES.out.txt.first(), false)
}
