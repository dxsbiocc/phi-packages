#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PARAGRAPH_MULTIGRMPY } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.cram       = null
params.crai       = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    manifest = Channel.of('id\tpath\tdepth\tread length\ntest\ttest.paired_end.sorted.cram\t0.73\t150').collectFile(name: 'manifest.txt', newLine: true)
    in_ch = Channel.of([[id: 'sample'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true), file(params.cram, checkIfExists: true), file(params.crai, checkIfExists: true)]).combine(manifest)
    fasta_ch = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'fasta_fai'], file(params.fai, checkIfExists: true)])

    PARAGRAPH_MULTIGRMPY(in_ch, fasta_ch, fai_ch)
}
