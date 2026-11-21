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
    if (params.vep_gtf && !params.vep_gtf_tbi) {
        error "vcf-annotate-ensemblvep-snpeff: 'vep_gtf_tbi' is required with 'vep_gtf'."
    }
    if (params.tools_to_use.tokenize(',')*.trim().contains('ensemblvep') && !params.vep_cache && !params.vep_gtf) {
        error "vcf-annotate-ensemblvep-snpeff: VEP needs 'vep_cache' or 'vep_gtf' (+ 'vep_gtf_tbi'); set 'tools_to_use' to 'snpeff' only if VEP is not wanted."
    }

    vcf_ch = Channel
        .fromPath(params.vcf, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f, file(params.vcf_tbi, checkIfExists: true), params.custom_files ? file(params.custom_files, checkIfExists: true) : []] }
    fasta_ch = Channel.value(params.fasta ? [[id: 'reference'], file(params.fasta, checkIfExists: true)] : [[], []])
    vep_cache_ch = Channel.value(params.vep_cache ? file(params.vep_cache, checkIfExists: true) : [])
    vep_extra_ch = Channel.value(params.vep_extra_files ? [file(params.vep_extra_files, checkIfExists: true)] : [])
    snpeff_cache_ch = Channel.value(params.snpeff_cache ? file(params.snpeff_cache, checkIfExists: true) : [])
    vep_gtf_ch = Channel.value(params.vep_gtf ? [file(params.vep_gtf, checkIfExists: true), file(params.vep_gtf_tbi, checkIfExists: true)] : [[], []])

    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF(vcf_ch, fasta_ch, params.vep_genome, params.vep_species, params.vep_cache_version, vep_cache_ch, vep_extra_ch, params.snpeff_db, snpeff_cache_ch, params.tools_to_use.tokenize(','), params.sites_per_chunk, vep_gtf_ch)
}
