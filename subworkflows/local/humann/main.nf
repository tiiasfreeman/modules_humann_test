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

workflow HUMANN {

    take:
    ch_reads
    ch_metaphlan_profile

    main:

    ch_versions = Channel.empty()

    HUMANN_DOWNLOADCHOCOPHLANDB ( )
    ch_versions = ch_versions.mix(HUMANN_DOWNLOADCHOCOPHLANDB.out.versions.first())

    HUMANN_DOWNLOADUNIREFDB ( )
    ch_versions = ch_versions.mix(HUMANN_DOWNLOADUNIREFDB.out.versions.first())

    HUMANN_HUMANN (ch_reads, ch_metaphlan_profile, HUMANN_DOWNLOADCHOCOPHLANDB.out.chocophlan_db, HUMANN_DOWNLOADUNIREFDB.out.uniref_db)
    ch_versions = ch_versions.mix(HUMANN_HUMANN.out.versions.first())

    versions = ch_versions

}
