#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { BAM_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING_STRELKA } from '../main.nf'

params.normal_cram = null
params.normal_crai = null
params.tumor_cram = null
params.tumor_crai = null
params.fasta = null
params.fai = null
params.dict = null
params.manta_vcf = null
params.manta_tbi = null
params.intervals = null
params.intervals_tbi = null
params.sample_id = "tumor_normal"
params.outdir = "results"

workflow {
    cram_ch = Channel.value([[id: params.sample_id], file(params.normal_cram, checkIfExists: true), file(params.normal_crai, checkIfExists: true), file(params.tumor_cram, checkIfExists: true), file(params.tumor_crai, checkIfExists: true), params.manta_vcf ? file(params.manta_vcf, checkIfExists: true) : [], params.manta_tbi ? file(params.manta_tbi, checkIfExists: true) : []])
    dict_ch = params.dict ? Channel.value([[id: 'reference'], file(params.dict, checkIfExists: true)]) : Channel.empty()
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch = Channel.value(file(params.fai, checkIfExists: true))
    intervals_ch = params.intervals ? Channel.value([file(params.intervals, checkIfExists: true), file(params.intervals_tbi, checkIfExists: true), 1]) : Channel.value([[], [], 0])

    BAM_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING_STRELKA(cram_ch, dict_ch, fasta_ch, fai_ch, intervals_ch)
}
