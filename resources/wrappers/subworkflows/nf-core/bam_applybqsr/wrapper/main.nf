#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { BAM_APPLYBQSR } from '../main.nf'

params.bam = null
params.bai = null
params.recal = null
params.fasta = null
params.fai = null
params.dict = null
params.intervals = null
params.output_suffix = "bam"
params.sample_id = "sample"
params.outdir = "results"

workflow {
    input_ch = Channel.value([[id: params.sample_id], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), file(params.recal, checkIfExists: true)])
    references_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true), file(params.dict, checkIfExists: true)])
    intervals_ch = params.intervals ? Channel.value([file(params.intervals, checkIfExists: true), 1]) : Channel.value([[], 0])

    BAM_APPLYBQSR(input_ch, references_ch, intervals_ch, params.output_suffix)
}
