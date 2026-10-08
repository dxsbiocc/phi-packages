#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { ARRIBA_VISUALISATION } from '../main.nf'


params.fusions    = null
params.gtf        = null
params.protein_domains = null
params.cytobands  = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', single_end: false], [], [], file(params.fusions, checkIfExists: true)])
    gtf_ch = Channel.value([[id: 'gtf'], file(params.gtf, checkIfExists: true)])
    domains_ch = Channel.value([[id: 'protein_domains'], file(params.protein_domains, checkIfExists: true)])
    cytobands_ch = Channel.value([[id: 'cytobands'], file(params.cytobands, checkIfExists: true)])

    ARRIBA_VISUALISATION(in_ch, gtf_ch, domains_ch, cytobands_ch)
}
