#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { VCF_ANNOTATE_ENSEMBLVEP_SNPEFF } from '../main.nf'

params.vcf = null
params.vcf_tbi = null
params.custom_files = null
params.fasta = null
params.vep_genome = "GRCh38"
params.vep_species = "homo_sapiens"
params.vep_cache_version = "110"
params.vep_cache = null
params.vep_extra_files = null
params.snpeff_db = "GRCh38.99"
params.snpeff_cache = null
params.tools_to_use = "ensemblvep,snpeff"
params.sites_per_chunk = 10000
params.vep_gtf = null
params.vep_gtf_tbi = null
params.outdir = "results"

workflow {
    vcf_ch = Channel
        .fromPath(params.vcf, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f, file(params.vcf_tbi, checkIfExists: true), params.custom_files ? file(params.custom_files, checkIfExists: true) : []] }
    fasta_ch = params.fasta ? Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)]) : Channel.empty()
    vep_cache_ch = params.vep_cache ? Channel.value(file(params.vep_cache, checkIfExists: true)) : Channel.empty()
    vep_extra_ch = params.vep_extra_files ? Channel.value(file(params.vep_extra_files, checkIfExists: true)) : Channel.empty()
    snpeff_cache_ch = params.snpeff_cache ? Channel.value(file(params.snpeff_cache, checkIfExists: true)) : Channel.empty()
    vep_gtf_ch = params.vep_gtf ? Channel.value([file(params.vep_gtf, checkIfExists: true), file(params.vep_gtf_tbi, checkIfExists: true)]) : Channel.empty()

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF(vcf_ch, fasta_ch, params.vep_genome, params.vep_species, params.vep_cache_version, vep_cache_ch, vep_extra_ch, params.snpeff_db, snpeff_cache_ch, params.tools_to_use.tokenize(','), params.sites_per_chunk, vep_gtf_ch)
}
