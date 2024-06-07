//
// Run PanPhlAn: Metagenomic profiling to identify gene composition of individual strains
//

include {   PANPHLAN_DOWNLOADPANGENOME              }   from '/workspace/modules_humann_test/modules/nf-core/panphlan/downloadpangenome/main.nf'
include {   SEQKIT_PAIR                             }   from '/workspace/modules_humann_test/modules/nf-core/seqkit/pair/main.nf'
include {   PANPHLAN_MAP                            }   from '/workspace/modules_humann_test/modules/nf-core/panphlan/panphlan_map/main.nf'
include {   PANPHLAN_PROFILING                      }   from '/workspace/modules_humann_test/modules/nf-core/panphlan/panphlan_profiling/main.nf'

workflow PANPHLAN_GENOMIC_PROFILING {

take:
    val_species_name                    //   value: "Genus_species";                    name of species for analysis
    ch_sample_sequence                  // channel: val(meta), path(list of files);     meta: ID is final name for run file, fastq sequence(s) for samples
    val_single_end                      //   value: "[True/False]";                     [True] single-end samples, [False] paired-end samples

    main:

    ch_versions = Channel.empty()

    PANPHLAN_DOWNLOADPANGENOME (val_species_name)
    ch_versions = ch_versions.mix(PANPHLAN_DOWNLOADPANGENOME.out.versions.first())

    if (val_single_end == "true") {

        println "All samples are single ended, skipping concatenation."

        PANPHLAN_MAP (
            PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
            PANPHLAN_DOWNLOADPANGENOME.out.indexes,
            ch_sample_sequence,
            val_species_name)
        ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())

    } else if (val_single_end == "false") {

        SEQKIT_PAIR ( ch_sample_sequence )
        ch_versions = ch_versions.mix(SEQKIT_PAIR.out.versions.first())

        SEQKIT_PAIR.out.reads.collect { meta, reads ->
            def paired_ends_merged = file("${val_species_name}.paired_ends_merged.fastq.gz")

            // Concatenate the files using Groovy and shell commands
            def paired_reads = SEQKIT_PAIR.out.reads.collect { it.toString() }.join(' ')
            def concatCommand = "cat ${paired_reads} > ${paired_ends_merged}"

            def process = ["bash", "-c", concatCommand].execute()
            process.waitFor()

            // Return the species name and the concatenated output file
            [val_species_name, paired_ends_merged]
        }.set { concatenated_reads }

        println "Paired end samples were concatenated."

        PANPHLAN_MAP (
            PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
            PANPHLAN_DOWNLOADPANGENOME.out.indexes,
            concatenated_reads,
            val_species_name)
        ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())
    }

    PANPHLAN_PROFILING (
        PANPHLAN_MAP.out.mapping_dir,
        PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
        val_species_name)
    ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())

emit:

    pangenome               = PANPHLAN_DOWNLOADPANGENOME.out.pangenome          // channel: $prefix/{species_name}_pangenome.tsv; pangenome used for selected species
    mapping_directory       = PANPHLAN_MAP.out.mapping_dir                      // channel: val(meta), path("${species_name}_map"); mapping directory built for selected species
    panphlan_profile        = PANPHLAN_PROFILING.out.profile_matrix             // channel: val(meta), path("*.tsv"); final profile matrix

    versions                = ch_versions                                       // channel: [ versions.yml ]
}
