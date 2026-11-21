#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { ASHLAR } from '../main.nf'

params.images = null
params.opt_dfp = null
params.opt_ffp = null
params.args = null
params.outdir = null

workflow {
    images_ch = Channel.fromPath("${params.images}/*").collect()
    def dfp = params.opt_dfp ? file(params.opt_dfp, checkIfExists: true) : []
    def ffp = params.opt_ffp ? file(params.opt_ffp, checkIfExists: true) : []

    ASHLAR(images_ch.map { imgs -> [[id: 'ashlar'], imgs] }, dfp, ffp)
}
