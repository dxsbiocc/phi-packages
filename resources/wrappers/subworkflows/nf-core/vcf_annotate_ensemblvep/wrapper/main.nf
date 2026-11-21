#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { VCF_ANNOTATE_ENSEMBLVEP } from '../main.nf'

params.vcf = null
params.custom_files = null
params.fasta = null
params.genome = "GRCh38"
params.species = "homo_sapiens"
params.cache_version = "110"
params.cache = null
params.extra_files = null
params.gtf = null
params.gtf_tbi = null
params.outdir = "results"

workflow {
    vcf_ch = Channel
        .fromPath(params.vcf, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f, params.custom_files ? file(params.custom_files, checkIfExists: true) : []] }
    fasta_ch = params.fasta ? Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)]) : Channel.empty()
    cache_ch = params.cache ? Channel.value([[id: 'cache'], file(params.cache, checkIfExists: true)]) : Channel.empty()
    extra_ch = params.extra_files ? Channel.value(file(params.extra_files, checkIfExists: true)) : Channel.empty()
    gtf_ch = params.gtf ? Channel.value([file(params.gtf, checkIfExists: true), file(params.gtf_tbi, checkIfExists: true)]) : Channel.empty()

    VCF_ANNOTATE_ENSEMBLVEP(vcf_ch, fasta_ch, params.genome, params.species, params.cache_version, cache_ch, extra_ch, gtf_ch)
}
