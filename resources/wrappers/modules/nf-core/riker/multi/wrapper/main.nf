#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { RIKER_MULTI } from '../main.nf'

params.bam = null
params.bai = null
params.error_vcf = null
params.error_vcf_idx = null
params.error_intervals = null
params.gcbias_exclude_intervals = null
params.hybcap_baits = null
params.hybcap_targets = null
params.rna_gene_model = null
params.rna_ribosomal_intervals = null
params.wgs_intervals = null
params.fasta = null
params.fai = null
params.args = null
params.outdir = null

workflow {
    def opt = { p -> p ? file(p, checkIfExists: true) : [] }
    def dirlist = { p -> p ? file(p, checkIfExists: true).listFiles().findAll { it.isFile() }.toList() : [] }

    RIKER_MULTI(
        Channel.value([
            [id: 'riker', single_end: false],
            file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true),
            opt(params.error_vcf), opt(params.error_vcf_idx), opt(params.error_intervals),
            opt(params.gcbias_exclude_intervals),
            dirlist(params.hybcap_baits), dirlist(params.hybcap_targets),
            opt(params.rna_gene_model), opt(params.rna_ribosomal_intervals), opt(params.wgs_intervals)
        ]),
        Channel.value([[id: 'reference'], opt(params.fasta), opt(params.fai)])
    )
}
