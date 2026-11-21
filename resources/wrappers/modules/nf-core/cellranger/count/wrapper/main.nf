#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { CELLRANGER_COUNT } from '../main.nf'

params.fastqs = null
params.reference = null
params.args = null
params.outdir = null

workflow {
    reads_ch = Channel.fromPath("${params.fastqs}/*.fastq.gz").collect()
    ref_ch   = Channel.value([[id: 'reference'], file(params.reference, checkIfExists: true)])

    CELLRANGER_COUNT(reads_ch.map { r -> [[id: 'cellranger', single_end: false], r] }, ref_ch)
}
