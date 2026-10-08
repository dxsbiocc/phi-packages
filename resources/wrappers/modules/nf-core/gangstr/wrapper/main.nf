#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GANGSTR } from '../main.nf'


params.cram       = null
params.crai       = null
params.regions    = null
params.fasta      = null
params.fai        = null
params.outdir     = null

workflow {
    bed_ch = Channel.of(file(params.regions, checkIfExists: true)).map { bed ->
        def content = bed.text.replace('\n', '').tokenize('\t')
        content.addAll(['5', 'CGCGC'])
        [content.join('\t')]
    }.collectFile(newLine: true) { content -> ['genome.bed', content[0]] }
    in_ch = Channel.value([[id: 'sample'], file(params.cram, checkIfExists: true), file(params.crai, checkIfExists: true)]).combine(bed_ch)

    GANGSTR(in_ch, file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true))
}
