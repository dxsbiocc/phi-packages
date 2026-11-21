#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { BAM_TELOMERE_ESTIMATION } from '../main.nf'

params.bam = null
params.bai = null
params.data_type = "short"
params.bed = null
params.control_bam = null
params.control_bai = null
params.fasta = null
params.fai = null
params.cytoband = null
params.run_telomerehunter = false
params.length_estimator = null
params.sample_id = "sample"
params.outdir = "results"

workflow {
    reads_ch = Channel.value([[id: params.sample_id], params.data_type, file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), params.bed ? file(params.bed, checkIfExists: true) : []])
    control_ch = params.control_bam ? Channel.value([[id: params.sample_id], file(params.control_bam, checkIfExists: true), file(params.control_bai, checkIfExists: true)]) : Channel.empty()
    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])
    cytoband_ch = Channel.value(params.cytoband ? [file(params.cytoband, checkIfExists: true)] : [[]])

    BAM_TELOMERE_ESTIMATION(reads_ch, control_ch, fasta_ch, cytoband_ch, params.run_telomerehunter, params.length_estimator ?: null)
}
