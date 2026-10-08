#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../main.nf'

params.fasta  = null
params.gtf    = null
params.outdir = null

workflow {
    ch_fasta = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> [[id: fasta.baseName], fasta] }
    ch_gtf = Channel
        .fromPath(params.gtf, checkIfExists: true)
        .map { gtf -> [[id: gtf.baseName], gtf] }

    STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
}
