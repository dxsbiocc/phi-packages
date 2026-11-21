#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { COMEBIN_RUNCOMEBIN } from '../main.nf'

params.assembly = null
params.bams = null
params.args = null
params.outdir = null

workflow {
    bams_ch = Channel.fromPath("${params.bams}/*.bam").collect()

    COMEBIN_RUNCOMEBIN(bams_ch.map { b ->
        [[id: 'comebin'], file(params.assembly, checkIfExists: true), b]
    })
}
