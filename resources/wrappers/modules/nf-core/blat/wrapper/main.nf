#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BLAT } from '../main.nf'

params.query   = null
params.subject = null
params.outdir  = null

workflow {
    query_ch   = Channel.value([[id: 'query'], file(params.query, checkIfExists: true)])
    subject_ch = Channel.value([[id: 'subject'], file(params.subject, checkIfExists: true)])

    BLAT(query_ch, subject_ch)
}
