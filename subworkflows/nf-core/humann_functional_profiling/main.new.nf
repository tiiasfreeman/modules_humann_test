//
// Run HUMAnN:functional profiling from metagenomic or metatranscriptomic FASTQ data
//

include {   HUMANN_DOWNLOADCHOCOPHLANDB                             } from '/workspace/modules_humann_test/modules/nf-core/humann/downloadchocophlandb/main.nf'
include {   HUMANN_DOWNLOADUNIREFDB                                 } from '/workspace/modules_humann_test/modules/nf-core/humann/downloadunirefdb/main.nf'
include {   HUMANN_DOWNLOADUTILITYMAPPINGDB                         } from '/workspace/modules_humann_test/modules/nf-core/humann/downloadutilitymappingdb/main.nf'

include {   CAT_FASTQ                                               } from '/workspace/modules_humann_test/modules/nf-core/cat/fastq/main.nf'

include {   HUMANN_HUMANN                                           } from '/workspace/modules_humann_test/modules/nf-core/humann/humann/main.nf'

include {   HUMANN_RENORMTABLE  as HUMANN_RENORM_PATHABUNDANCE      } from '/workspace/modules_humann_test/modules/nf-core/humann/renorm_table/main.nf'

include {   HUMANN_RENORMTABLE  as HUMANN_RENORM_GENEFAMILIES       } from '/workspace/modules_humann_test/modules/nf-core/humann/renorm_table/main.nf'
include {   HUMANN_REGROUPTABLE                                     } from '/workspace/modules_humann_test/modules/nf-core/humann/regroup_table/main.nf'
include {   HUMANN_RENAMETABLE                                      } from '/workspace/modules_humann_test/modules/nf-core/humann/rename_table/main.nf'

include {   HUMANN_JOINTABLES   as HUMANN_JOIN_PATHABUNDANCE        } from '/workspace/modules_humann_test/modules/nf-core/humann/join_tables/main.nf'
include {   HUMANN_JOINTABLES   as HUMANN_JOIN_GENEFAMILIES         } from '/workspace/modules_humann_test/modules/nf-core/humann/join_tables/main.nf'
include {   HUMANN_JOINTABLES   as HUMANN_JOIN_PATHCOVERAGE         } from '/workspace/modules_humann_test/modules/nf-core/humann/join_tables/main.nf'


workflow HUMANN_FUNCTIONAL_PROFILING {

    take:

    val_chocophlan_db_version           //   value: version of ChocoPhlan database to download,         options; ["DEMO", "full", "ec_filtered"]
    val_uniref_db_version               //   value: version of UniRef database to download,             options; ["DEMO_diamond", "uniref90_diamond", "uniref90_ec_filtered_diamond"]
    val_utility_mapping_db_version      //   value: version of utility mapping database to download,    options; ["full"]
    ch_fastq                            // channel: [ val(meta), [fastq list]]                          options; single-end fastq ["single_end:true"], paired-end fastq ["single_end:false"]
    ch_metaphlan_profile                // channel: [ val(meta), [profile]]
    val_renorm_units                    //   value: renormalization scheme units,                       options; copies per million ["cpm"], relative abundance ["relab"], reads per kilobase ["rpk"]
    val_group_type                      //   value: regrouping data by function options,                options; ["uniref90_rxn"]
        // (if extra time, other files already added to db, just need coded in)
    val_name_type                       //   value: renaming table features,                            options; ["kegg-orthology", "ec", "metacyc-rxn", "metacyc-pwy", "pfam", "eggnog", "go" or "infogo1000"]
    val_multi_sample_join               //   value: merge multi-sample final output tables,             options; ["Y","N"]


    main:

    ch_versions = Channel.empty()

    if (val_renorm_units == "rpk" && val_multi_sample_join == "Y") {
            fail "Error: Cannot merge multiple sample HUMAnN output files when val_renorm_units = "rpk". Please re-enter inputs where val_renorm_units = "cpm"."
        } else if (val_renorm_units == "relab" && val_multi_sample_join == "Y") {
            fail "Error: Cannot merge multiple sample HUMAnN output files when val_renorm_units = "relab". Please re-enter inputs where val_renorm_units = "cpm"."
        }

//
// Downloading ChocoPhlan, UniRef and utility mapping databases
//

    HUMANN_DOWNLOADCHOCOPHLANDB (val_chocophlan_db_version)                             //change to "full", delete variable
    ch_versions = ch_versions.mix(HUMANN_DOWNLOADCHOCOPHLANDB.out.versions.first())

    HUMANN_DOWNLOADUNIREFDB (val_uniref_db_version)                                     //variable should remain, depends on mapping types
    ch_versions = ch_versions.mix(HUMANN_DOWNLOADUNIREFDB.out.versions.first())

    HUMANN_DOWNLOADUTILITYMAPPINGDB ("val_utility_mapping_db_version")                  //change to "full", delete variable
    ch_versions = ch_versions.mix(HUMANN_DOWNLOADUTILITYMAPPINGDB.out.versions.first())

//
// Formatting FASTQ **paired-end data** as .fastq file for input into HUMAnN
//
    // Note: HUMAnN documentation suggests all paired-end files are concatenated into single
    //      FASTQ or FASTA file because HUMAnN doesn't handle single paired-end files well

    def sample_sequences = ch_fastq.flatMap{ metadata, fastq_files -> fastq_files}

    sample_sequences.collectFile(name: 'concatenated_fastq_files.fastq') {path -> path} into {concatenated_fastq_files}

//
// Running HUMAnN module
//

    HUMANN_HUMANN (
        concatenated_fastq_files,
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

        println "Proceeding to regrouping genefamilies table..."

    // Regrouping genefamilies table using functional group type...
    HUMANN_REGROUPTABLE (
        HUMANN_HUMANN.out.genefamilies,
        val_renorm_units,
        val_group_type,
        HUMANN_DOWNLOADUTILITYMAPPINGDB.out.utility_mapping_db)
    ch_versions = ch_versions.mix(HUMANN_REGROUPTABLE.out.versions.first())

    } else {

    // Renormalizing genefamilies table using renorm units...
    HUMANN_RENORM_GENEFAMILIES (
        "genefamilies",
        HUMANN_HUMANN.out.genefamilies,
        val_renorm_units)
    ch_versions = ch_versions.mix(HUMANN_RENORM_GENEFAMILIES.out.versions.first())

    // Regrouping genefamilies table using functional group type...
    HUMANN_REGROUPTABLE (
        HUMANN_RENORM_GENEFAMILIES.out.humann_output_renorm,
        val_renorm_units,
        val_group_type,
        HUMANN_DOWNLOADUTILITYMAPPINGDB.out.utility_mapping_db)
    ch_versions = ch_versions.mix(HUMANN_REGROUPTABLE.out.versions.first())
    }

    // Renaming table features using name type...
    HUMANN_RENAMETABLE (
        HUMANN_REGROUPTABLE.out.humann_output_regroup,
        val_renorm_units,
        val_group_type,
        val_name_type,
        HUMANN_DOWNLOADUTILITYMAPPINGDB.out.utility_mapping_db)
    ch_versions = ch_versions.mix(HUMANN_RENAMETABLE.out.versions.first())
}

//
// For use in multi-sample analysis, join multiple output tables.
//
    if (val_multi_sample_join == "N") {

        println "No output tables were merged, individual files will appear for each sample."

    // Conditional Outputs
        if (val_renorm_units == "rpk") {
            emit: pathabundance_output = HUMANN_HUMANN.out.pathabundance                        // channel: [ val(meta), [pathabundance]]; pathabundance data without renorm
        } else {
            emit: pathabundance_output = HUMANN_RENORM_PATHABUNDANCE.out.humann_output_renorm   // channel: [ val(meta), [humann_output_renorm]]; renorm pathabundance data
        }

    emit:

    genefamilies_output     = HUMANN_RENAMETABLE.out.humann_output_rename   // channel: [ val(meta), [humann_output_rename]]; genefamilies data after renorm (if specified), regroup, rename
    pathcoverage_output     = HUMANN_HUMANN.out.pathcoverage                // channel: [ val(meta), [pathcoverage]]; pathcoverage data without adjustment

    versions                = ch_versions                                   // channel: [ versions.yml ]
    humann_log              = HUMANN_HUMANN.out.log                         // channel: [ val(meta), [log] ]

    } else if (val_multi_sample_join == "Y") {

    // Merging genefamilies output files from multiple samples...
    HUMANN_JOIN_GENEFAMILIES (
        HUMANN_RENAMETABLE.out.humann_output_rename.map{ [ [id:'genefamilies'], it[1] ] }.groupTuple( sort: "deep" ),
        "genefamilies")
    ch_versions = ch_versions.mix(HUMANN_JOIN_GENEFAMILIES.out.versions.first())

    // Merging pathcoverage output files from multiple samples...
    HUMANN_JOIN_PATHCOVERAGE (
        HUMANN_HUMANN.out.pathcoverage.map{ [ [id:'pathcoverage'], it[1] ] }.groupTuple( sort: "deep" ),
        "pathcoverage")
    ch_versions = ch_versions.mix(HUMANN_JOIN_PATHCOVERAGE.out.versions.first())

    // Merging pathabundance output files from multiple samples...
        if (val_renorm_units == "rpk") {
            fail "Error: Cannot merge pathabundance files when val_renorm_units = "rpk"."
        } else {

        HUMANN_JOIN_PATHABUNDANCE (
                HUMANN_RENORM_PATHABUNDANCE.out.humann_output_renorm.map{ [ [id:'pathabundance'], it[1] ] }.groupTuple( sort: "deep" ),
                "pathabundance")
        ch_versions = ch_versions.mix(HUMANN_JOIN_PATHABUNDANCE.out.versions.first())

        println "Output tables with the same name_type were merged, both individual and merged files will appear in HUMAnN output."
        }

    emit:

    genefamilies_merged = HUMANN_JOIN_GENEFAMILIES.out.humann_output_merged                 // channel: [ val(meta), [humann_output_merged]]; merged genefamilies data from multi-sample fastq data
    pathabundance_merged = HUMANN_JOIN_PATHABUNDANCE.out.humann_output_merged               // channel: [ val(meta), [humann_output_merged]]; merged pathabundace data from multi-sample fastq data
    genefamilies_merged = HUMANN_JOIN_PATHCOVERAGE.out.humann_output_merged                 // channel: [ val(meta), [humann_output_merged]]; merged pathcoverage data from multi-sample fastq data

    versions                = ch_versions                                                   // channel: [ versions.yml ]
    humann_log              = HUMANN_HUMANN.out.log                                         // channel: [ val(meta), [log] ]
    }
}

//Other notes to use

// Create a channel from a collection of output files from different processes
//def combinedChannel = Channel.from([PROCESS1.out.file, PROCESS2.out.file, PROCESS3.out.file])

//output: dir 'humann_outputdir'
//mkdir -p humann_outputdir

//input: file inputFile from myChannel.collect()
//output: file "humann_outputdir/${inputFile.name}" into

//humann_output_ch = HUMANN_RENAMETABLE.out.humann_output_rename
    //.join(HUMANN_HUMANN.out.pathcoverage)
    //.join(HUMANN_HUMANN.out.pathabundance)
    //.join(HUMANN_RENORM_PATHABUNDANCE.out.humann_output_renorm)

//for phython shell scripts: foreach fastq.file in input_file { process()}
