// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules


//
// Run HUMAnN:functional profiling from metagenomic or metatranscriptomic FASTQ data
//

include {   HUMANN_DOWNLOADCHOCOPHLANDB     } from 'modules/local/humann/downloadchocophlandb/main.nf'
include {   HUMANN_DOWNLOADUNIREFDB         } from 'modules/local/humann/downloadunirefdb/main.nf'
include {   HUMANN_DOWNLOADUTILITYMAPPINGDB } from 'modules/local/humann/downloadutilitymappingdb/main.nf'
include {   CAT_FASTQ                       } from 'modules/nf-core/cat/fastq/main.nf'
include {   HUMANN_HUMANN                   } from 'modules/local/humann/humann/main.nf'
include {   HUMANN_JOINTABLES               } from 'modules/local/humann/join_tables/main.nf'
include {   HUMANN_RENORMTABLE              } from 'modules/local/humann/renorm_table/main.nf'
include {   HUMANN_REGROUPTABLE             } from 'modules/local/humann/regroup_table/main.nf'
include {   HUMANN_RENAMETABLE              } from 'modules/local/humann/rename_table/main.nf'

workflow HUMANNSUBWORKFLOW {

    take:
    // TODO nf-core: edit input (take) channels
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SAMTOOLS_SORT ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

