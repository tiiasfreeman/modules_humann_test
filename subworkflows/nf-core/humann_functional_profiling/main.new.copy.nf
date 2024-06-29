//
// Run HUMAnN:functional profiling from metagenomic or metatranscriptomic FASTQ data
//

include {   HUMANN_DOWNLOADCHOCOPHLANDB                             } from '/workspace/modules_humann_test/modules/nf-core/humann/downloadchocophlandb/main.nf'
include {   HUMANN_DOWNLOADUNIREFDB                                 } from '/workspace/modules_humann_test/modules/nf-core/humann/downloadunirefdb/main.nf'
include {   HUMANN_DOWNLOADUTILITYMAPPINGDB                         } from '/workspace/modules_humann_test/modules/nf-core/humann/downloadutilitymappingdb/main.nf'

include {   HUMANN_HUMANN                                           } from '/workspace/modules_humann_test/modules/nf-core/humann/humann/main.nf'

include {   HUMANN_RENORMTABLE  as HUMANN_RENORM_PATHABUNDANCE      } from '/workspace/modules_humann_test/modules/nf-core/humann/renorm_table/main.nf'

//include {   HUMANN_RENORMTABLE  as HUMANN_RENORM_GENEFAMILIES       } from '/workspace/modules_humann_test/modules/nf-core/humann/renorm_table/main.nf'
//include {   HUMANN_REGROUPTABLE                                     } from '/workspace/modules_humann_test/modules/nf-core/humann/regroup_table/main.nf'
//include {   HUMANN_RENAMETABLE                                      } from '/workspace/modules_humann_test/modules/nf-core/humann/rename_table/main.nf'

include {   HUMANN_JOINTABLES   as HUMANN_JOIN_PATHABUNDANCE        } from '/workspace/modules_humann_test/modules/nf-core/humann/join_tables/main.nf'
//include {   HUMANN_JOINTABLES   as HUMANN_JOIN_GENEFAMILIES         } from '/workspace/modules_humann_test/modules/nf-core/humann/join_tables/main.nf'
//include {   HUMANN_JOINTABLES   as HUMANN_JOIN_PATHCOVERAGE         } from '/workspace/modules_humann_test/modules/nf-core/humann/join_tables/main.nf'


workflow HUMANN_FUNCTIONAL_PROFILING {

    take:

    val_chocophlan_db_version           //   value: version of ChocoPhlan database to download,         options; ["DEMO", "full", "ec_filtered"]
    val_uniref_db_version               //   value: version of UniRef database to download,             options; ["DEMO_diamond", "uniref90_diamond", "uniref90_ec_filtered_diamond"]
    val_utility_mapping_db_version      //   value: version of utility mapping database to download,    options; ["full"]
    ch_reads                            // channel: [ val(meta), [fastq list]]                          options; single-end fastq ["single_end:true"], paired-end fastq ["single_end:false"]
    ch_metaphlan_profile                // channel: [ val(meta), [profile]]
    val_renorm_units                    //   value: renormalization scheme units,                       options; copies per million ["cpm"], relative abundance ["relab"], reads per kilobase ["rpk"]
    val_group_type                      //   value: regrouping data by function options,                options; ["uniref90_rxn"]
    val_name_type                       //   value: renaming table features,                            options; ["kegg-orthology", "ec", "metacyc-rxn", "metacyc-pwy", "pfam", "eggnog", "go" or "infogo1000"]
    val_multi_sample_join               //   value: merge multi-sample final output tables,             options; ["Y","N"]

    main:

    ch_versions = Channel.empty()

//
// Downloading ChocoPhlan, UniRef and utility mapping databases
//

    HUMANN_DOWNLOADCHOCOPHLANDB (val_chocophlan_db_version)                                 //change to "full", delete variable
    ch_versions = ch_versions.mix(HUMANN_DOWNLOADCHOCOPHLANDB.out.versions.first())

    HUMANN_DOWNLOADUNIREFDB (val_uniref_db_version)                                         //variable should remain, depends on mapping types
    ch_versions = ch_versions.mix(HUMANN_DOWNLOADUNIREFDB.out.versions.first())

    HUMANN_DOWNLOADUTILITYMAPPINGDB (val_utility_mapping_db_version)                        //change to "full", delete variable
    ch_versions = ch_versions.mix(HUMANN_DOWNLOADUTILITYMAPPINGDB.out.versions.first())

//
// Running HUMAnN module
//

    ch_reads.each { tuple_sequence ->
        HUMANN_HUMANN (
            tuple_sequence,
            ch_metaphlan_profile,
            HUMANN_DOWNLOADCHOCOPHLANDB.out.chocophlan_db,
            HUMANN_DOWNLOADUNIREFDB.out.uniref_db) }
        ch_versions = ch_versions.mix(HUMANN_HUMANN.out.versions.first())

    def ch_humann_pathabundance = HUMANN_HUMANN.out.pathabundance
    def ch_humann_genefamilies = HUMANN_HUMANN.out.genefamilies
    def ch_humann_pathcoverage = HUMANN_HUMANN.out.pathcoverage


//
// Customizing HUMAnN output files: Pathabundance
//

    if (val_renorm_units == "rpk") {
        println "Pathabundance file already in reads per kilobase (rpk), skipping renormalization."

    } else {
        HUMANN_RENORM_PATHABUNDANCE (
            "pathabundance",
            ch_humann_pathabundance,
            val_renorm_units)
        ch_versions = ch_versions.mix(HUMANN_RENORM_PATHABUNDANCE.out.versions.first())
    }

    def ch_humann_renorm_pathabundance = HUMANN_RENORM_PATHABUNDANCE.out.humann_output_renorm.collect()
        .flatMap { collected_tuples -> collected_tuples.collect { it[1] } }

//
// Customizing HUMAnN output files: Genefamilies
//


//
// For use in multi-sample analysis, join multiple output tables.
//

        HUMANN_JOIN_PATHABUNDANCE (
            ch_humann_renorm_pathabundance,
            "pathabundance")
        ch_versions = ch_versions.mix(HUMANN_JOIN_PATHABUNDANCE.out.versions.first())

    emit:

    pathabundance_merged    = HUMANN_JOIN_PATHABUNDANCE.out.humann_output_merged               // channel: [ val(meta), [humann_output_merged]]; merged pathabundace data from multi-sample fastq data
    humann_pathabundance    = ch_humann_pathabundance
    humann_genefamilies     = ch_humann_genefamilies
    humann_pathcoverage     = ch_humann_pathcoverage

    versions                = ch_versions                                                   // channel: [ versions.yml ]
    humann_log              = HUMANN_HUMANN.out.log
}
