#!/usr/bin/env nextflow
// Thin agent-facing adapter over ../main.nf (STAR_ALIGN). A pre-built STAR
// index (`index`) is used as-is; without one, ../genomegenerate/main.nf builds
// it from `fasta` + `gtf` first — the only way to align against a genome that
// has no index yet. This still composes only vendored modules directly,
// matching the composition rule in
// docs/design/phi-wrapper-agent-composition-design.md section 6.
nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../../genomegenerate/main.nf'
include { STAR_ALIGN } from '../main.nf'

params.index  = null
params.fasta  = null
params.gtf    = null
params.reads  = null
params.outdir = null

// fromFilePairs names each sample by the file name up to the first `{` or `[` of the
// glob, so a glob that picks samples with braces (`S{01,02}_*_R{1,2}.fq.gz`) lumps
// them all into one "sample". Stop with the reason instead of running on merged files.
def checkReadGroup(sample, reads) {
    def files = reads instanceof List ? reads : [reads]
    if (files.size() > 2) {
        error "The reads glob put ${files.size()} files into one sample '${sample}': ${files*.name.join(', ')}. " +
            "A sample is named by the file name up to the first { or [ in the glob, so braces or brackets before the last * merge samples. " +
            "Use a glob like /data/*_R{1,2}.fastq.gz, or run each group of samples separately."
    }
    [sample, reads]
}

workflow {
    if (!params.index && !params.fasta) {
        error "star-align needs either `index` (a pre-built STAR index directory) or `fasta` to build one."
    }

    gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true).map { f -> [[id: f.simpleName], f] }

    if (params.index) {
        index_ch = Channel.value([[id: 'star_index'], file(params.index, checkIfExists: true)])
    } else {
        fasta_ch = Channel.fromPath(params.fasta, checkIfExists: true).map { f -> [[id: f.simpleName], f] }
        STAR_GENOMEGENERATE(fasta_ch, gtf_ch)
        index_ch = STAR_GENOMEGENERATE.out.index.first()
    }

    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { sample, reads -> checkReadGroup(sample, reads) }
        .map { sample, reads -> [[id: sample, single_end: reads instanceof List ? reads.size() == 1 : true], reads] }

    STAR_ALIGN(reads_ch, index_ch, gtf_ch.first(), false)
}
