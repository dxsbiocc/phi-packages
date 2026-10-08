#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// No pre-built cache is passed — snpEff downloads the requested database
// from its own server on first use, matching what happens whenever a
// caller doesn't already have a cache (the module's own test instead
// pulls a pre-built cache from a private S3 bucket this environment has
// no credentials for).
nextflow.enable.dsl = 2

include { SNPEFF_SNPEFF } from '../main.nf'

params.vcf    = null
params.db     = 'NC_045512.2'
params.outdir = null

workflow {
    vcf_ch = Channel
        .fromPath(params.vcf, checkIfExists: true)
        .map { vcf -> [[id: vcf.simpleName], vcf] }

    cache_ch = Channel.value([[:], []])

    SNPEFF_SNPEFF(vcf_ch, params.db, cache_ch)
}
