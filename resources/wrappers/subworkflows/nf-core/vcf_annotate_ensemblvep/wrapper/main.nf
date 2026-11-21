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
    if (!params.cache && !params.gtf) {
        error "vcf-annotate-ensemblvep: provide 'cache' (an Ensembl VEP cache) or 'gtf' (+ 'gtf_tbi') for offline annotation."
    }
    if (params.gtf && !params.gtf_tbi) {
        error "vcf-annotate-ensemblvep: 'gtf_tbi' is required with 'gtf'."
    }

    vcf_ch = Channel
        .fromPath(params.vcf, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f, params.custom_files ? file(params.custom_files, checkIfExists: true) : []] }
    fasta_ch = Channel.value(params.fasta ? [[id: 'reference'], file(params.fasta, checkIfExists: true)] : [[], []])
    cache_ch = Channel.value(params.cache ? [[id: 'cache'], file(params.cache, checkIfExists: true)] : [[], []])
    extra_ch = Channel.value(params.extra_files ? [file(params.extra_files, checkIfExists: true)] : [])
    gtf_ch = Channel.value(params.gtf ? [file(params.gtf, checkIfExists: true), file(params.gtf_tbi, checkIfExists: true)] : [[], []])

    VCF_ANNOTATE_ENSEMBLVEP(vcf_ch, fasta_ch, params.genome, params.species, params.cache_version, cache_ch, extra_ch, gtf_ch)
}
