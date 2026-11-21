#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_TRIM_FASTP_FASTQC } from '../main.nf'

params.reads = null
params.adapter_fasta = null
params.save_trimmed_fail = false
params.discard_trimmed_pass = false
params.save_merged = false
params.skip_fastp = false
params.skip_fastqc = false
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads, params.adapter_fasta ? file(params.adapter_fasta, checkIfExists: true) : []] }

    FASTQ_TRIM_FASTP_FASTQC(reads_ch, params.save_trimmed_fail, params.discard_trimmed_pass, params.save_merged, params.skip_fastp, params.skip_fastqc)
}
