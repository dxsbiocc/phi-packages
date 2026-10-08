#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_GTFFILTER } from '../main.nf'

params.gtf    = null
params.fasta  = null
params.outdir = null

workflow {
    // meta.id must not equal the input GTF's own basename: the module writes
    // its output as "${meta.id}.gtf" next to the input, and a matching name
    // would truncate the input file while it is still being read.
    gtf_ch = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .map { gtf -> [[id: "${gtf.simpleName}_filtered"], gtf] }

    fasta_ch = params.fasta
        ? Channel.fromPath(params.fasta, checkIfExists: true).map { fasta -> [[id: fasta.simpleName], fasta] }
        : Channel.value([[:], []])

    CUSTOM_GTFFILTER(gtf_ch, fasta_ch)
}
