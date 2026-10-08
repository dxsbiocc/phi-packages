#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1.
//
// Composes STRINGTIE_STRINGTIE + STRINGTIE_MERGE internally: merge's whole
// point is combining multiple samples' assembled transcript GTFs into one,
// so a runnable default needs at least one real assembled GTF to merge —
// there's nothing to point this wrapper's default params.json at otherwise.
// A single-sample merge is a real (if minimal) invocation, not a stand-in.
nextflow.enable.dsl = 2

include { STRINGTIE_STRINGTIE } from '../../stringtie/main.nf'
include { STRINGTIE_MERGE     } from '../main.nf'

params.bam    = null
params.gtf    = null
params.outdir = null

workflow {
    ch_gtf = file(params.gtf, checkIfExists: true)

    Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.baseName], bam, []] }
        .set { bam_ch }

    STRINGTIE_STRINGTIE(bam_ch, [], ch_gtf)

    ch_assembled = STRINGTIE_STRINGTIE.out.transcript_gtf
        .map { _meta, gtf -> gtf }
        .collect()
        .map { gtfs -> [[id: 'merged'], gtfs] }

    STRINGTIE_MERGE(ch_assembled, [[:], ch_gtf])
}
