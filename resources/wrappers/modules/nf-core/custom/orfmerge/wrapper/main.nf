#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_ORFMERGE } from '../main.nf'
include { CUSTOM_ORFNORMALISE } from '../../orfnormalise/main.nf'

params.ribotish   = null
params.ribocode   = null
params.gtf        = null
params.outdir     = null

workflow {
    gtf_ch = Channel.value([[id: 'reference'], file(params.gtf, checkIfExists: true)])
    calls = Channel.of(
        [[id: 'sample1', caller: 'ribotish'], file(params.ribotish, checkIfExists: true), 'ribotish'],
        [[id: 'sample1', caller: 'ribocode'], file(params.ribocode, checkIfExists: true), 'ribocode'],
        [[id: 'sample2', caller: 'ribotish'], file(params.ribotish, checkIfExists: true), 'ribotish'])

    CUSTOM_ORFNORMALISE(calls, gtf_ch)
    ch_beds = CUSTOM_ORFNORMALISE.out.bed12.map { meta, bed -> bed }.collect().map { beds -> ['cohort', beds] }
    ch_tsvs = CUSTOM_ORFNORMALISE.out.tsv.map { meta, tsv -> tsv }.collect().map { tsvs -> ['cohort', tsvs] }
    CUSTOM_ORFMERGE(ch_beds.combine(ch_tsvs, by: 0).map { k, beds, tsvs -> [[id: 'cohort'], beds, tsvs] })
}
