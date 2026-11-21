#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { BAM_CREATE_SOM_PON_GATK } from '../main.nf'

params.bam = null
params.bai = null
params.fasta = null
params.fai = null
params.dict = null
params.pon_norm = "pon"
params.intervals_gendb = null
params.intervals = null
params.sample_id = "normal"
params.outdir = "results"

workflow {
    input_ch = Channel.value([[id: params.sample_id], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'reference'], file(params.fai, checkIfExists: true), []])
    dict_ch = Channel.value([[id: 'reference'], file(params.dict, checkIfExists: true)])
    gendb_ch = params.intervals_gendb ? Channel.value(file(params.intervals_gendb, checkIfExists: true)) : Channel.empty()
    intervals_num_ch = params.intervals ? Channel.value([file(params.intervals, checkIfExists: true), 1]) : Channel.value([[], 0])

    BAM_CREATE_SOM_PON_GATK(input_ch, fasta_ch, fai_ch, dict_ch, params.pon_norm, gendb_ch, intervals_num_ch)
}
