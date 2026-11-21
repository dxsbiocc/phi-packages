#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { BAM_VARIANT_CALLING_MPILEUP_BCFTOOLS } from '../main.nf'

params.bam = null
params.bai = null
params.posfile = null
params.fasta = null
params.fai = null
params.panel_id = "panel"
params.chr = "chr22"
params.sample_merge_key = "id"
params.sample_merge_value = "all_samples"
params.region_gather_keys = "panel_id,id"
params.sort_region_gather = false
params.annotate = false
params.outdir = "results"

workflow {
    bam_ch = Channel.value([[id: params.sample_merge_value], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    posfile_ch = Channel.value([[panel_id: params.panel_id, chr: params.chr], file(params.posfile, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])

    BAM_VARIANT_CALLING_MPILEUP_BCFTOOLS(bam_ch, posfile_ch, fasta_ch, [params.sample_merge_key], [params.sample_merge_value], params.region_gather_keys.tokenize(','), params.sort_region_gather, params.annotate)
}
