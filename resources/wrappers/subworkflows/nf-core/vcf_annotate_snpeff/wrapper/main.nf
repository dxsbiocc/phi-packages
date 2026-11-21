#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { VCF_ANNOTATE_SNPEFF } from '../main.nf'

params.vcf = null
params.snpeff_db = "GRCh38.99"
params.snpeff_cache = null
params.outdir = "results"

workflow {
    vcf_ch = Channel
        .fromPath(params.vcf, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f] }
    cache_ch = params.snpeff_cache ? Channel.value(file(params.snpeff_cache, checkIfExists: true)) : Channel.empty()

    VCF_ANNOTATE_SNPEFF(vcf_ch, params.snpeff_db, cache_ch)
}
