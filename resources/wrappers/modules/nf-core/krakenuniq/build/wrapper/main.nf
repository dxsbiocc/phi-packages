#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { KRAKENUNIQ_BUILD } from '../main.nf'

params.library = null
params.taxonomy = null
params.seqid2taxid = null
params.keep_intermediate = null
params.args = null
params.outdir = null

workflow {
    lib_ch = Channel.fromPath("${params.library}/*").collect()
    tax_ch = Channel.fromPath("${params.taxonomy}/*").collect()

    KRAKENUNIQ_BUILD(
        lib_ch.combine(tax_ch).map { lib, tax ->
            [[id: 'krakenuniq'], lib, tax, file(params.seqid2taxid, checkIfExists: true)]
        },
        params.keep_intermediate
    )
}
