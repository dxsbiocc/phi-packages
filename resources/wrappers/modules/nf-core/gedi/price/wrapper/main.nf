#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { GEDI_PRICE } from '../main.nf'

params.bams = null
params.index = null
params.args = null
params.outdir = null

workflow {
    bams_ch = Channel.fromPath("${params.bams}/*.bam").collect()
    bais_ch = Channel.fromPath("${params.bams}/*.bam.bai").collect()

    GEDI_PRICE(
        bams_ch.combine(bais_ch).map { b, bi -> [[id: 'price', single_end: true], b, bi] },
        Channel.value([[id: 'index'], file(params.index, checkIfExists: true)])
    )
}
