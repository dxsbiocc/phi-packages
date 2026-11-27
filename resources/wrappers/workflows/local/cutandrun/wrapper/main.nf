nextflow.enable.dsl = 2

include { CUTANDRUN } from '../main.nf'

workflow {
    ch_rows = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)

    ch_reads = ch_rows
        .filter { row -> row.fastq_1 }
        .map { row ->
            def single_end = (row.single_end?.toString()?.toLowerCase() == 'true') || !row.fastq_2
            def meta = [id: row.sample, single_end: single_end]
            def reads = single_end
                ? file(row.fastq_1, checkIfExists: true)
                : [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
            tuple(meta, reads)
        }

    ch_seacr = ch_rows
        .filter { row -> row.bedgraph }
        .map { row ->
            def meta = [id: row.sample]
            def ctrl = row.ctrl_bedgraph ? file(row.ctrl_bedgraph, checkIfExists: true) : []
            tuple(meta, file(row.bedgraph, checkIfExists: true), ctrl)
        }

    ch_macs2 = ch_rows
        .filter { row -> row.bam }
        .map { row ->
            def single_end = (row.single_end?.toString()?.toLowerCase() == 'true')
            def meta = [id: row.sample, single_end: single_end]
            def ctrl = row.ctrl_bam ? file(row.ctrl_bam, checkIfExists: true) : []
            tuple(meta, file(row.bam, checkIfExists: true), ctrl)
        }

    ch_index = params.bowtie2_index
        ? Channel.of(tuple([id: 'index'], file(params.bowtie2_index, checkIfExists: true)))
        : Channel.empty()

    ch_fasta = params.fasta
        ? Channel.of(tuple([id: 'fasta'], file(params.fasta, checkIfExists: true)))
        : Channel.of(tuple([id: 'none'], []))

    def do_align = (params.run_align == true) && (params.bowtie2_index != null)

    CUTANDRUN(
        ch_reads,
        ch_index,
        ch_fasta,
        ch_seacr,
        ch_macs2,
        params.caller ?: 'seacr',
        params.seacr_threshold ?: 0.01,
        params.macs2_gsize ?: 'hs',
        do_align
    )
}
