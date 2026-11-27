//
// CUT&RUN / CUT&Tag Alignment & Peak Calling Workflow
//
// Composes vendored nf-core modules:
//   - BOWTIE2_ALIGN (paired/single-end FASTQ alignment -> coordinate-sorted BAM)
//   - SEACR_CALLPEAK (CUT&RUN/CUT&Tag bedGraph peak calling vs IgG control or top threshold)
//   - MACS2_CALLPEAK (optional MACS2 narrow/broad peak calling from BAM)
//
// Called peaks feed directly into chromatin-profiling (chip_inspect / chip_consensus / chip_super_enhancer).
//

include { BOWTIE2_ALIGN  } from '../../../modules/nf-core/bowtie2/align/main.nf'
include { SEACR_CALLPEAK } from '../../../modules/nf-core/seacr/callpeak/main.nf'
include { MACS2_CALLPEAK } from '../../../modules/nf-core/macs2/callpeak/main.nf'

workflow CUTANDRUN {
    take:
    ch_reads         // channel: [ val(meta), path(reads) ] (optional)
    ch_bowtie2_index // channel: [ val(meta), path(index_dir) ] (optional)
    ch_fasta         // channel: [ val(meta), path(fasta) ] (optional)
    ch_seacr_input   // channel: [ val(meta), path(bedgraph), path(ctrl_bedgraph_or_empty) ]
    ch_macs2_input   // channel: [ val(meta), path(ip_bam), path(ctrl_bam_or_empty) ]
    caller           // val: string ('seacr' | 'macs2' | 'both')
    seacr_threshold  // val: number (used when ctrl_bedgraph is empty)
    macs2_gsize      // val: string (e.g. 'hs', 'mm', '2.7e9')
    run_align        // val: boolean

    main:
    ch_bam         = Channel.empty()
    ch_seacr_peaks = Channel.empty()
    ch_macs2_peaks = Channel.empty()

    if (run_align) {
        BOWTIE2_ALIGN(
            ch_reads,
            ch_bowtie2_index.first(),
            ch_fasta.first(),
            false,
            true
        )
        ch_bam = BOWTIE2_ALIGN.out.bam
    }

    if (caller in ['seacr', 'both']) {
        SEACR_CALLPEAK(
            ch_seacr_input,
            seacr_threshold
        )
        ch_seacr_peaks = SEACR_CALLPEAK.out.bed
    }

    if (caller in ['macs2', 'both']) {
        MACS2_CALLPEAK(
            ch_macs2_input,
            macs2_gsize
        )
        ch_macs2_peaks = MACS2_CALLPEAK.out.peak
    }

    emit:
    bam         = ch_bam
    seacr_peaks = ch_seacr_peaks
    macs2_peaks = ch_macs2_peaks
}
