#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { VCF_EXTRACT_RELATE_SOMALIER } from '../main.nf'

params.vcfs = null
params.fasta = null
params.fai = null
params.somalier_sites = null
params.ped = null
params.sample_groups = null
params.cohort = "cohort"
params.common_id = "cohort"
params.outdir = "results"

workflow {
    vcfs_ch = Channel
        .fromPath(params.vcfs, checkIfExists: true)
        .map { f -> [[id: f.simpleName, cohort: params.cohort], f, [], 1] }
    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'reference'], file(params.fai, checkIfExists: true)])
    sites_ch = Channel.value(file(params.somalier_sites, checkIfExists: true))
    peds_ch = Channel.value([[id: 'ped'], file(params.ped, checkIfExists: true)])
    groups_ch = params.sample_groups ? Channel.value(file(params.sample_groups, checkIfExists: true)) : []

    VCF_EXTRACT_RELATE_SOMALIER(vcfs_ch, fasta_ch, fai_ch, sites_ch, peds_ch, groups_ch, params.common_id ?: null)
}
