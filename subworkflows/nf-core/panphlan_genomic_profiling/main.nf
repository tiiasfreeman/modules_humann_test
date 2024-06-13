//
// Run PanPhlAn: Metagenomic profiling to identify gene composition of individual strains
//

include {   PANPHLAN_DOWNLOADPANGENOME              }   from '/workspace/modules_humann_test/modules/nf-core/panphlan/downloadpangenome/main.nf'
include {   PANPHLAN_MAP                            }   from '/workspace/modules_humann_test/modules/nf-core/panphlan/panphlan_map/main.nf'
include {   PANPHLAN_PROFILING                      }   from '/workspace/modules_humann_test/modules/nf-core/panphlan/panphlan_profiling/main.nf'

workflow PANPHLAN_GENOMIC_PROFILING {

take:
    val_species_name                    //   value: "Genus_species";                    name of species for analysis
    ch_sample_sequences                  // channel: val(meta), path(list of files);     meta: ID is final name for run file, fastq sequence(s) for samples
    //val_single_end                      //   value: "[True/False]";                     [True] single-end samples, [False] paired-end samples

    main:

    ch_versions = Channel.empty()

    PANPHLAN_DOWNLOADPANGENOME (val_species_name)
    ch_versions = ch_versions.mix(PANPHLAN_DOWNLOADPANGENOME.out.versions.first())

    def sample_sequences = ch_sample_sequences
        .flatMap{ metadata, files -> files.collect {file -> [metadata, file] } }
        //.flatten()

    sample_sequences.each { tuple_sequence ->
        PANPHLAN_MAP (
            PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
            PANPHLAN_DOWNLOADPANGENOME.out.indexes,
            tuple_sequence,
            val_species_name)
        }
    ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())

    panphlan_merged_maps_txt = PANPHLAN_MAP.out.mapping_dir
        .map{ it[1].toString() }
        .collectFile(name: 'panphlan_merged_maps.txt') { paths -> paths.join('\n')}
    panphlan_merged_maps_txt.view()

    PANPHLAN_PROFILING (
        panphlan_merged_maps_txt,
        PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
        val_species_name)
    ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())

emit:

    pangenome               = PANPHLAN_DOWNLOADPANGENOME.out.pangenome          // channel: $prefix/{species_name}_pangenome.tsv; pangenome used for selected species
    mapping_directory       = PANPHLAN_MAP.out.mapping_dir                      // channel: val(meta), path("${species_name}_map"); mapping directory built for selected species
    panphlan_merged_maps_txt
    panphlan_profile        = PANPHLAN_PROFILING.out.profile_matrix             // channel: val(meta), path("*.tsv"); final profile matrix

    versions                = ch_versions                                       // channel: [ versions.yml ]
}
