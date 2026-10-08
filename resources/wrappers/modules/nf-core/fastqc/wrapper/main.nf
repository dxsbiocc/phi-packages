#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — turns
// wrapper.yaml's `reads`/`outdir` params into the module's channel shape
// and lets the module's own publishDir (declared in this file's
// nextflow.config) land the outputs. Does not modify the module itself.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { FASTQC } from '../main.nf'

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
    Channel
        .fromFilePairs(params.reads, size: -1)
        .map { sample, reads -> checkReadGroup(sample, reads) }
        .map { sample, reads -> [[id: sample], reads] }
        .set { reads_ch }

    FASTQC(reads_ch)
}
