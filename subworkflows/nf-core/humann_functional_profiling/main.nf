//
// Run HUMAnN:functional profiling from metagenomic or metatranscriptomic FASTQ data
//

include {   HUMANN_DOWNLOADCHOCOPHLANDB                             } from '/workspace/modules_humann_test/modules/local/humann/downloadchocophlandb/main.nf'
include {   HUMANN_DOWNLOADUNIREFDB                                 } from '/workspace/modules_humann_test/modules/local/humann/downloadunirefdb/main.nf'
include {   HUMANN_DOWNLOADUTILITYMAPPINGDB                         } from '/workspace/modules_humann_test/modules/local/humann/downloadutilitymappingdb/main.nf'

include {   CAT_FASTQ                                               } from '/workspace/modules_humann_test/modules/nf-core/cat/fastq/main.nf'

include {   HUMANN_HUMANN                                           } from '/workspace/modules_humann_test/modules/local/humann/humann/main.nf'

include {   HUMANN_RENORMTABLE  as HUMANN_RENORM_PATHABUNDANCE      } from '/workspace/modules_humann_test/modules/local/humann/renorm_table/main.nf'

include {   HUMANN_RENORMTABLE  as HUMANN_RENORM_GENEFAMILIES       } from '/workspace/modules_humann_test/modules/local/humann/renorm_table/main.nf'
include {   HUMANN_REGROUPTABLE                                     } from '/workspace/modules_humann_test/modules/local/humann/regroup_table/main.nf'
include {   HUMANN_RENAMETABLE                                      } from '/workspace/modules_humann_test/modules/local/humann/rename_table/main.nf'
// See notes for including HUMANN_JOINTABLES at bottom

workflow HUMANN_FUNCTIONAL_PROFILING {

    take:

    val_chocophlan_db_version           //   value: version of ChocoPhlan database to download,         options; ["DEMO", "full", "ec_filtered"]
    val_uniref_db_version               //   value: version of UniRef database to download,             options; ["uniref90_diamond", "uniref90_ec_filtered_diamond"]
    val_utility_mapping_db_version      //   value: version of utility mapping database to download,    options; ["full"]
    ch_fastq                            // channel: [ val(meta), [fastq list]]
    ch_metaphlan_profile                // channel: [ val(meta), [profile]]
    val_renorm_units                    //   value: renormalization scheme units,                       options; copies per million ["cpm"], relative abundance ["relab"], reads per kilobase ["rpk"]
        // put in no renormalization if val_renorm_units = "rpk" --> not actual option for HUMANN_RENORMTABLE units
    val_group_type                      //   value:                                                     options; ["uniref90_rxn"]
        // "regroup" currently under construction, only one option (if extra time, other files already added to db, just need coded in)
    val_name_type                       //   value: renaming table features,                            options; ["kegg-orthology", "ec", "metacyc-rxn", "metacyc-pwy", "pfam", "eggnog", "go" or "infogo1000"]

    main:

    ch_versions = Channel.empty()

//
// Downloading ChocoPhlan, UniRef and utility mapping databases
//

    HUMANN_DOWNLOADCHOCOPHLANDB (val_chocophlan_db_version)
    ch_versions = ch_versions.mix(HUMANN_DOWNLOADCHOCOPHLANDB.out.versions.first())

    HUMANN_DOWNLOADUNIREFDB (val_uniref_db_version)
    ch_versions = ch_versions.mix(HUMANN_DOWNLOADUNIREFDB.out.versions.first())

    HUMANN_DOWNLOADUTILITYMAPPINGDB (val_utility_mapping_db_version)
    ch_versions = ch_versions.mix(HUMANN_DOWNLOADUTILITYMAPPINGDB.out.versions.first())

//
// Formatting FASTQ data as single .fastq file for input into HUMAnN
//

    CAT_FASTQ (ch_fastq)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

//
// Running HUMAnN module
//

    HUMANN_HUMANN (
        CAT_FASTQ.out.reads,
        ch_metaphlan_profile,
        HUMANN_DOWNLOADCHOCOPHLANDB.out.chocophlan_db,
        HUMANN_DOWNLOADUNIREFDB.out.uniref_db)
    ch_versions = ch_versions.mix(HUMANN_HUMANN.out.versions.first())

//
// Customizing HUMAnN output files: Pathabundance
//

    if (val_renorm_units == "rpk") {

        println "Pathabundance file already in reads per kilobase (rpk), skipping renormalization."

    } else {

    HUMANN_RENORM_PATHABUNDANCE (
        "pathabundance",
        HUMANN_HUMANN.out.pathabundance,
        val_renorm_units)
    ch_versions = ch_versions.mix(HUMANN_RENORM_PATHABUNDANCE.out.versions.first())

    }

//
// Customizing HUMAnN output files: Genefamilies
//

    if (val_renorm_units == "rpk") {

        println "Genefamilies file already in reads per kilobase (rpk), skipping renormalization."

    HUMANN_REGROUPTABLE (
        HUMANN_HUMANN.out.genefamilies,
        val_renorm_units,
        val_group_type,
        HUMANN_DOWNLOADUTILITYMAPPINGDB.out.utility_mapping_db)
    ch_versions = ch_versions.mix(HUMANN_REGROUPTABLE.out.versions.first())

    } else {

    HUMANN_RENORM_GENEFAMILIES (
        "genefamilies",
        HUMANN_HUMANN.out.genefamilies,
        val_renorm_units)
    ch_versions = ch_versions.mix(HUMANN_RENORM_GENEFAMILIES.out.versions.first())

    HUMANN_REGROUPTABLE (
        HUMANN_RENORM_GENEFAMILIES.out.humann_output_renorm,
        val_renorm_units,
        val_group_type,
        HUMANN_DOWNLOADUTILITYMAPPINGDB.out.utility_mapping_db)
    ch_versions = ch_versions.mix(HUMANN_REGROUPTABLE.out.versions.first())
    }

    HUMANN_RENAMETABLE (
        HUMANN_REGROUPTABLE.out.humann_output_regroup,
        val_renorm_units,
        val_group_type,
        val_name_type,
        HUMANN_DOWNLOADUTILITYMAPPINGDB.out.utility_mapping_db)
    ch_versions = ch_versions.mix(HUMANN_RENAMETABLE.out.versions.first())


    // Conditional Outputs
    if (val_renorm_units == "rpk") {
        emit: pathabundance_merged = HUMANN_HUMANN.out.pathabundance                        // channel: [ val(meta), [pathabundance]]; pathabundance data without renorm
    } else {
        emit: pathabundance_merged = HUMANN_RENORM_PATHABUNDANCE.out.humann_output_renorm   // channel: [ val(meta), [humann_output_renorm]]; renorm pathabundance data
    }

    emit:

    genefamilies_merged     = HUMANN_RENAMETABLE.out.humann_output_rename   // channel: [ val(meta), [humann_output_rename]]; genefamilies data after renorm (if specified), regroup, rename
    pathcoverage_merged     = HUMANN_HUMANN.out.pathcoverage                // channel: [ val(meta), [pathcoverage]]; pathcoverage data without adjustment

    versions                = ch_versions                                   // channel: [ versions.yml ]
    humann_log              = HUMANN_HUMANN.out.log                         // channel: [ val(meta), [log] ]
}



// Notes for later addition of HUMANN_JOINTABLES module for multi-sample processing (make conditional)

    //include {   HUMANN_JOINTABLES   as HUMANN_JOIN_PATHABUNDANCE        } from '/workspace/modules_humann_test/modules/local/humann/join_tables/main.nf'
    //include {   HUMANN_JOINTABLES   as HUMANN_JOIN_GENEFAMILIES         } from '/workspace/modules_humann_test/modules/local/humann/join_tables/main.nf'
    //include {   HUMANN_JOINTABLES   as HUMANN_JOIN_PATHCOVERAGE         } from '/workspace/modules_humann_test/modules/local/humann/join_tables/main.nf'

    //take: val_merge_fastqs                    //   value: [Y/N], for multi-fastq analysis; [Y] to combine fastq files before running HUMAnN, [N] to run HUMAnN analysis separately then join output tables

//
// Customizing HUMAnN output files: Pathcoverage
//
