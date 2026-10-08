#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1.
//
// Composes HISAT2_EXTRACTSPLICESITES + HISAT2_BUILD (matching nf-core's own
// module test setup for HISAT2_BUILD) so a splice-aware index is the
// default, not a bare-genome one — HISAT2_BUILD's own `gtf`/`splicesites`
// inputs are optional, but dropping them silently produces a lower-quality
// index for RNA-seq use, which is the point of using HISAT2 at all.
// `hisat2_memory_input` is set to a low '1.GB' threshold, not left empty:
// the module's own logic is `avail_mem >= hisat2_build_memory` decides
// splice-aware — a *higher* threshold makes splice-aware *less* likely to
// trigger, not more (confirmed by a real run: passing '' fell through to
// Integer.MAX_VALUE and silently built a non-splice-aware index instead).
// A low threshold means any real task.memory clears it easily.
nextflow.enable.dsl = 2

include { HISAT2_EXTRACTSPLICESITES } from '../../extractsplicesites/main.nf'
include { HISAT2_BUILD              } from '../main.nf'

params.fasta  = null
params.gtf    = null
params.outdir = null

workflow {
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
    ch_gtf   = Channel.fromPath(params.gtf, checkIfExists: true)

    HISAT2_EXTRACTSPLICESITES(ch_gtf.map { gtf -> [[id: gtf.baseName], gtf] })

    ch_build_input = ch_fasta
        .combine(ch_gtf)
        .combine(HISAT2_EXTRACTSPLICESITES.out.txt.map { _meta, txt -> txt })
        .map { fasta, gtf, splicesites -> [[id: fasta.baseName], fasta, gtf, splicesites] }

    HISAT2_BUILD(ch_build_input, '1.GB')
}
