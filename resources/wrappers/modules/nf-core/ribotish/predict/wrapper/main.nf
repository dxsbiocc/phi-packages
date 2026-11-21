#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { RIBOTISH_PREDICT } from '../main.nf'

params.bam_ribo = null
params.bai_ribo = null
params.bam_ti = null
params.bai_ti = null
params.fasta = null
params.gtf = null
params.candidate_orfs = null
params.para_ribo = null
params.para_ti = null
params.reference_gtf = null
params.args = null
params.outdir = null

workflow {
    def opt2 = { a, b -> (a && b) ? [file(a, checkIfExists: true), file(b, checkIfExists: true)] : [[], []] }
    def opt1 = { p -> p ? [file(p, checkIfExists: true)] : [[]] }

    RIBOTISH_PREDICT(
        Channel.value([[id: 'ribo', single_end: true, strandedness: 'forward'],
                       file(params.bam_ribo, checkIfExists: true), file(params.bai_ribo, checkIfExists: true)]),
        Channel.value([[id: 'ti']] + opt2(params.bam_ti, params.bai_ti)),
        Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.gtf, checkIfExists: true)]),
        Channel.value([[id: 'orfs']] + opt1(params.candidate_orfs)),
        Channel.value([[id: 'para_ribo']] + opt1(params.para_ribo)),
        Channel.value([[id: 'para_ti']] + opt1(params.para_ti)),
        Channel.value([[id: 'secondary']] + opt1(params.reference_gtf))
    )
}
