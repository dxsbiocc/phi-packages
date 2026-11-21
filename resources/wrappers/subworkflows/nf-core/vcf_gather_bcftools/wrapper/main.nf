#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { VCF_GATHER_BCFTOOLS } from '../main.nf'

params.vcfs = null
params.vcf_indexes = null
params.common_meta = "id"
params.sort = false
params.outdir = "results"

workflow {
    vcfs_ch = Channel
        .fromPath(params.vcfs, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f, [], 1] }

    VCF_GATHER_BCFTOOLS(vcfs_ch, params.common_meta.tokenize(','), params.sort)
}
