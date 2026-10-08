//
// Alignment with STAR
//
// Sentieon (commercial) and Parabricks (GPU-only) acceleration branches
// removed — see docs/design/phi-wrapper-agent-composition-design.md
// discussion: edge hardware/licensing paths, never exercised or verified,
// not worth the vendored module weight. Plain STAR_ALIGN is now the only
// path.
include { STAR_ALIGN                                } from '../../../modules/nf-core/star/align'
include { BAM_SORT_STATS_SAMTOOLS                   } from '../../nf-core/bam_sort_stats_samtools'


//
// Function that parses and returns the alignment rate from the STAR log output
//
def getStarPercentMapped(_params, align_log) {
    def percent_aligned = 0
    def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
    align_log.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            percent_aligned = matcher[0][1].toFloat()
        }
    }

    return percent_aligned
}

workflow ALIGN_STAR {
    take:
    reads                // channel: [ val(meta), [ reads ] ]
    index                // channel: [ val(meta), [ index ] ]
    gtf                  // channel: [ val(meta), [ gtf ] ]
    star_ignore_sjdbgtf  // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    fasta_fai            // channel: [ val(meta), path(fasta), path(fai) ]
    skip_markduplicates  // boolean: whether to skip marking duplicates

    main:

    //
    // Map reads with STAR
    //
    STAR_ALIGN(reads, index, gtf, star_ignore_sjdbgtf)

    ch_orig_bam = STAR_ALIGN.out.bam
    ch_log_final = STAR_ALIGN.out.log_final
    ch_log_out = STAR_ALIGN.out.log_out
    ch_log_progress = STAR_ALIGN.out.log_progress
    ch_bam_sorted = STAR_ALIGN.out.bam_sorted
    ch_bam_transcript = STAR_ALIGN.out.bam_transcript
    ch_fastq = STAR_ALIGN.out.fastq
    ch_tab = STAR_ALIGN.out.tab
    ch_percent_mapped = ch_log_final.map { meta, log -> [ meta, getStarPercentMapped(params, log) ] }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS(ch_orig_bam, fasta_fai)

    emit:
    orig_bam = ch_orig_bam                          // channel: [ val(meta), bam            ]
    log_final = ch_log_final                        // channel: [ val(meta), log_final      ]
    log_out = ch_log_out                            // channel: [ val(meta), log_out        ]
    log_progress = ch_log_progress                  // channel: [ val(meta), log_progress   ]
    bam_sorted = ch_bam_sorted                      // channel: [ val(meta), bam_sorted     ]
    bam_transcript = ch_bam_transcript              // channel: [ val(meta), bam_transcript ]
    fastq = ch_fastq                                // channel: [ val(meta), fastq          ]
    tab = ch_tab                                    // channel: [ val(meta), tab            ]
    bam = BAM_SORT_STATS_SAMTOOLS.out.bam           // channel: [ val(meta), [ bam ] ]
    index = BAM_SORT_STATS_SAMTOOLS.out.index       // channel: [ val(meta), [ index ] ]
    stats = BAM_SORT_STATS_SAMTOOLS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    percent_mapped = ch_percent_mapped              // channel: [ val(meta), percent_mapped ]
}
