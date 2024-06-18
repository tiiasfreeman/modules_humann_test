//
// Run PanPhlAn: Metagenomic profiling to identify gene composition of individual strains
//

include {   PANPHLAN_DOWNLOADPANGENOME              }   from '/workspace/modules_humann_test/modules/nf-core/panphlan/downloadpangenome/main.nf'
include {   PANPHLAN_MAP                            }   from '/workspace/modules_humann_test/modules/nf-core/panphlan/panphlan_map/main.nf'
include {   PANPHLAN_PROFILING                      }   from '/workspace/modules_humann_test/modules/nf-core/panphlan/panphlan_profiling/main.nf'

workflow PANPHLAN_GENOMIC_PROFILING {

take:
    val_species_name                        //   value: "Genus_species";                    name of species for analysis
    ch_sample_sequences                     // channel: val(meta), path(list of files);     meta: ID is final name for run file, fastq sequence(s) for samples
    val_sample_count                        //   value: "single", "multiple";               number of sample sequence files being assessed

    main:

    ch_versions = Channel.empty()

    // Download species-specific pangenome, indexes and annotations
    PANPHLAN_DOWNLOADPANGENOME (val_species_name)
    ch_versions = ch_versions.mix(PANPHLAN_DOWNLOADPANGENOME.out.versions)

        println("Running PanPhlAn for multiple sample files.")

        // Separate individual sequencing samples to run through PANPLAN_MAP
        def sample_sequences = ch_sample_sequences
            .flatMap{ metadata, files -> files.collect {file -> [metadata, file] } }

        // Individually run each sample sequence through PANPHLAN_MAP
        sample_sequences.each { tuple_sequence ->
            PANPHLAN_MAP (
                PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
                PANPHLAN_DOWNLOADPANGENOME.out.indexes,
                tuple_sequence,
                val_species_name) }
            ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())

        //Create channel with all PANPHLAN_MAP maps for PANPHLAN_PROFILING
        def ch_ combined_mapping_dir = PANPHLAN_MAP.out.sample_map
            .map{ [ [id:'all_samples'], it[1] ] }.groupTuple( sort: { it.getName() }  )
            .view()

        //Use maps to create a absence/presence profile matrix with PANPHLAN_PROFILING
        PANPHLAN_PROFILING (
            ch_combined_mapping_dir,
            PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
            val_species_name)
        ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())

        emit:
            pangenome               = PANPHLAN_DOWNLOADPANGENOME.out.pangenome          // channel: $prefix/{species_name}_pangenome.tsv; pangenome used for selected species
            mapping_directory       = ch_combined_mapping_dir                           // channel: val(meta), path("${species_name}_map_dir"); mapping directory built for selected samples
            panphlan_profile        = PANPHLAN_PROFILING.out.profile_matrix             // channel: val(meta), path("*.tsv"); final profile matrix
            versions                = ch_versions                                       // channel: [ versions.yml ]
}
