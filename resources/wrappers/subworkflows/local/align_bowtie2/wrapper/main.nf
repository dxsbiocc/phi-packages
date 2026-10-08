#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored align_bowtie2 subworkflow at
// ../main.nf. A pre-built Bowtie2 index (`index`) is used as-is; without one,
// bowtie2/build builds it from `fasta` first.
nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../../modules/nf-core/bowtie2/build/main.nf'
include { ALIGN_BOWTIE2 } from '../main.nf'

params.reads_1        = null
params.reads_2        = null
params.index          = null
params.fasta          = null
params.save_unaligned = false
params.outdir         = null

workflow {
    if (!params.index && !params.fasta) {
        error "align-bowtie2 needs either `index` (a pre-built Bowtie2 index directory) or `fasta` to build one."
    }

    if (params.index) {
        index_ch = Channel.value(file(params.index, checkIfExists: true))
    } else {
        fasta_ch = Channel
            .fromPath(params.fasta, checkIfExists: true)
            .map { fasta -> [[id: fasta.baseName], fasta] }
        BOWTIE2_BUILD(fasta_ch)
        index_ch = BOWTIE2_BUILD.out.index.map { meta, index -> index }
    }

    def read1 = file(params.reads_1, checkIfExists: true)
    reads_ch = Channel.of([
        [id: read1.simpleName.replaceAll(/[._-]R?1$/, ''), single_end: false],
        [read1, file(params.reads_2, checkIfExists: true)]
    ])

    fasta_fai_ch = Channel.value([[:], [], []])

    ALIGN_BOWTIE2(reads_ch, index_ch, fasta_fai_ch)
}
