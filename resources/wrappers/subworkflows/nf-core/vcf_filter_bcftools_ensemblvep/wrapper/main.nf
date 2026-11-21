#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { VCF_FILTER_BCFTOOLS_ENSEMBLVEP } from '../main.nf'

params.vcf = null
params.filter_vep_feature_file = null
params.filter_with_bcftools = true
params.filter_with_filter_vep = false
params.outdir = "results"

workflow {
    vcf_ch = Channel
        .fromPath(params.vcf, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f] }
    feature_ch = params.filter_vep_feature_file ? Channel.fromPath(params.filter_vep_feature_file, checkIfExists: true).map { f -> [[id: f.simpleName], f] } : Channel.empty()

    VCF_FILTER_BCFTOOLS_ENSEMBLVEP(vcf_ch, feature_ch, params.filter_with_bcftools, params.filter_with_filter_vep)
}
