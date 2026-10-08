#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Composes untar directly (vendored module, not another wrapper) so callers
// supply the db as a plain tar.gz instead of pre-extracting it themselves.
nextflow.enable.dsl = 2

include { UNTAR } from '../../../untar/main.nf'
include { KRAKEN2_KRAKEN2 } from '../main.nf'

params.reads                 = null
params.db                    = null
params.save_output_fastqs    = false
params.save_reads_assignment = false
params.outdir                = null

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
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { sample, reads -> checkReadGroup(sample, reads) }
        .map { sample, reads -> [[id: sample, single_end: reads instanceof List ? reads.size() == 1 : true], reads] }

    db_archive_ch = Channel
        .fromPath(params.db, checkIfExists: true)
        .map { db -> [[id: 'kraken2_db'], db] }

    UNTAR(db_archive_ch)

    KRAKEN2_KRAKEN2(reads_ch, UNTAR.out.untar.map { meta, dir -> dir }, params.save_output_fastqs, params.save_reads_assignment)
}
