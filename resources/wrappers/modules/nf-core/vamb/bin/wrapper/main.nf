#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { VAMB_BIN } from '../main.nf'

params.assembly = null
params.abundance_tsv = null
params.bams = null
params.taxonomy = null
params.args = null
params.outdir = null

workflow {
    def bams = params.bams
        ? file(params.bams, checkIfExists: true).listFiles().findAll { it.isFile() && it.name.endsWith('.bam') }.toList()
        : []
    def opt = { p -> p ? file(p, checkIfExists: true) : [] }

    VAMB_BIN(Channel.value([
        [id: 'vamb'],
        file(params.assembly, checkIfExists: true),
        opt(params.abundance_tsv),
        bams,
        opt(params.taxonomy)
    ]))
}
