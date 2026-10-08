include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main'


workflow BAM_STRINGTIE_MERGE {
    take:
    bam_sorted // channel: [ meta, bam ]
    ch_chrgtf // channel: [ meta, gtf ]

    main:

    // modules/nf-core/stringtie/stringtie (post-3.26.0) added long-read BAM
    // support: input is now (meta,srbam,lrbam) + a `mode` flag list. No
    // lrbam, no mode — reproduces the old version's plain short-read-only
    // assembly behaviour exactly.
    STRINGTIE_STRINGTIE(
        bam_sorted.map { meta, bam -> [meta, bam, []] },
        [],
        ch_chrgtf.map { _meta, gtf -> [gtf] },
    )

    STRINGTIE_STRINGTIE.out.transcript_gtf
        .map { _meta, gtf -> gtf }
        .toSortedList { a, b -> a.name <=> b.name }
        .map { gtfs -> [ [id: 'stringtie_merge'], gtfs ] }
        .set { collected_gtfs }

    STRINGTIE_MERGE(
        collected_gtfs,
        ch_chrgtf
    )

    emit:
    stringtie_gtf = STRINGTIE_MERGE.out.merged_gtf // channel: [ meta, gtf ]
}
