#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_NORM } from '../main.nf'

params.vcf    = null
params.fasta  = null
params.outdir = null

workflow {
    // meta.id must not equal the input VCF's own basename: the module
    // writes its output as "${meta.id}.vcf.gz" next to the input, and
    // Nextflow excludes anything matching an input filename from the
    // output glob, so a matching name makes the run fail with
    // "Missing output file(s)".
    vcf_ch = Channel
        .fromPath(params.vcf, checkIfExists: true)
        .map { vcf -> [[id: "${vcf.simpleName}_norm"], vcf, []] }

    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])

    BCFTOOLS_NORM(vcf_ch, fasta_ch)
}
