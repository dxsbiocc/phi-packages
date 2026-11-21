#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_QC_STATS } from '../main.nf'

params.reads = null
params.skip_fastqc = false
params.skip_seqfu_check = false
params.skip_seqfu_stats = false
params.skip_seqkit_stats = false
params.skip_seqtk_comp = false
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }

    FASTQ_QC_STATS(reads_ch, params.skip_fastqc, params.skip_seqfu_check, params.skip_seqfu_stats, params.skip_seqkit_stats, params.skip_seqtk_comp)
}
