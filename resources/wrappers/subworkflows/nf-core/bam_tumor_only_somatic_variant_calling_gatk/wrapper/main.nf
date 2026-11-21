#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK } from '../main.nf'

params.bam = null
params.bai = null
params.fasta = null
params.fai = null
params.dict = null
params.germline_resource = null
params.germline_resource_tbi = null
params.alleles = null
params.alleles_tbi = null
params.panel_of_normals = null
params.panel_of_normals_tbi = null
params.interval_file = null
params.sample_id = "tumor"
params.outdir = "results"

workflow {
    input_ch = Channel.value([[id: params.sample_id], [file(params.bam, checkIfExists: true)], [file(params.bai, checkIfExists: true)], []])
    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'reference'], file(params.fai, checkIfExists: true), []])
    dict_ch = Channel.value([[id: 'reference'], file(params.dict, checkIfExists: true)])
    alleles_ch = params.alleles ? Channel.value(file(params.alleles, checkIfExists: true)) : Channel.value([])
    alleles_tbi_ch = params.alleles_tbi ? Channel.value(file(params.alleles_tbi, checkIfExists: true)) : Channel.value([])
    germline_ch = Channel.value(file(params.germline_resource, checkIfExists: true))
    germline_tbi_ch = Channel.value(file(params.germline_resource_tbi, checkIfExists: true))
    pon_ch = params.panel_of_normals ? Channel.value(file(params.panel_of_normals, checkIfExists: true)) : Channel.value([])
    pon_tbi_ch = params.panel_of_normals_tbi ? Channel.value(file(params.panel_of_normals_tbi, checkIfExists: true)) : Channel.value([])
    interval_ch = params.interval_file ? Channel.value(file(params.interval_file, checkIfExists: true)) : Channel.value([])

    BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK(input_ch, fasta_ch, fai_ch, dict_ch, alleles_ch, alleles_tbi_ch, germline_ch, germline_tbi_ch, pon_ch, pon_tbi_ch, interval_ch)
}
