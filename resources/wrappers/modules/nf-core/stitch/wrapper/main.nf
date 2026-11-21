#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { STITCH } from '../main.nf'

params.crams = null
params.crais = null
params.posfile = null
params.input = null
params.genetic_map = null
params.rdata = null
params.samplename = null
params.chromosome_name = null
params.start = null
params.end = null
params.buffer = null
params.k = null
params.n_gen = null
params.fasta = null
params.fai = null
params.seed = null
params.args = null
params.outdir = null

workflow {
    ch_reads = Channel.fromPath("${params.crams}/*").collect()
        .map { fs -> fs.findAll { it.name.endsWith('.cram') || it.name.endsWith('.bam') } }
    ch_idx = Channel.fromPath("${params.crais}/*").collect()
    ch_list = ch_reads.map { fs -> fs*.name.join('\n') }.collectFile(name: 'cramlist.txt', newLine: true)
    def opt = { p -> p ? file(p, checkIfExists: true) : [] }

    ch_main = ch_reads.combine(ch_idx).combine(ch_list).map { reads, idx, list ->
        [[id: 'stitch'], reads, idx, list,
         opt(params.samplename),
         file(params.posfile, checkIfExists: true),
         opt(params.input),
         opt(params.genetic_map),
         opt(params.rdata),
         params.chromosome_name, params.start, params.end, params.buffer, params.k, params.n_gen]
    }

    STITCH(ch_main,
           Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)]),
           params.seed)
}
