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



test("Bacteroides_fragilis - combined single end - fails") {
        when {
            workflow {
                """
                input[0] = "Bacteroides_fragilis"
                input[1] = ([
                    [ id:'test', single_end:false ], // meta map
                        file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test1_1.fastq.gz", checkIfExists:true),
                        file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test1_2.fastq.gz", checkIfExists:true),
                        file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test2_1.fastq.gz", checkIfExists:true),
                        file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test2_2.fastq.gz", checkIfExists:true)
                ])
                input[2] = "false"
                """
            }
        }

        then {
            assertAll(
                { assert workflow.success},
                { assert snapshot(workflow.out).match()}
            )
        }
    }


    test("Staphylococcus aureus - paired end - works") {
        when {
            workflow {
                """
                input[0] = "Staphylococcus_aureus"
                input[1] = Channel.of([
                    [ id:'test', single_end:false ], // meta map
                        file("/workspace/modules_humann_test/modules/local/data/Staph_aureus_SRR29085994.fastq.gz")
                ])
                input[2] = "false"
                """
            }
        }

        then {
            assertAll(
                { assert workflow.success},
                { assert snapshot(workflow.out).match()}
            )
        }
    }


    test("Eubacterium rectale sample fastq - paired_end") {

        when {
            workflow {
                """
                input[0] = "Eubacterium_rectale"
                input[1] = Channel.of ([
                            [ id:"test1", single_end:false ],
                                    file("/workspace/modules_humann_test/panphlan_tutorial/samples_fastq/CCMD34381688ST-21-0.fastq")
                            ])
                """
            }
        }

        then {
            assertAll(
                { assert workflow.success},
                { assert snapshot(workflow.out).match()}
            )
        }
    }


    PANPHLAN_DOWNLOADPANGENOME (val_species_name)
    ch_versions = ch_versions.mix(PANPHLAN_DOWNLOADPANGENOME.out.versions.first())

    //def sample_sequence = ch_sample_sequences.split{ it }

    //def sample_sequences = ch_sample_sequences
        //.map{ tuple ->
            //def (metadata,files) = tuple
        //files.collect {file -> tuple(metadata, file) } }
        //.flatten()

    def sample_sequence = ch_sample_sequences.split()

    PANPHLAN_MAP (
        PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
        PANPHLAN_DOWNLOADPANGENOME.out.indexes,
        sample_sequence,
        val_species_name)
    ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())

    panphlan_merged_maps_txt = PANPHLAN_MAP.out.mapping_dir.map{ it[1].toString() }.collectFile(name: 'panphlan_merged_maps.txt') { paths -> paths.join('\n')}
    panphlan_merged_maps_txt.view()

    PANPHLAN_PROFILING (
        panphlan_merged_maps_txt,
        PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
        val_species_name)
    ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())


    panphlan_merged_maps_txt = PANPHLAN_MAP.out.mapping_dir
        .map{ it[1].toString() }
        .collectFile(name: 'panphlan_merged_maps.txt') { paths -> paths.join('\n')}
    panphlan_merged_maps_txt.view()


    def merged_mapping_dirs = PANPHLAN_MAP.out.mapping_dir
        .collectFile(name: 'merged_map_dir.tsv') { files ->
            files.collect {it.text}.join('\n')
        }
    merged_mapping_dirs.view()


def merged_mapping_dir = combined_mapping_dir.map { maps ->
        def combined_maps = file("merged_mapping_dir.tsv")
        def header = maps[0].text.readLines().head()
        def content = maps.collect { it.text.readLines().drop(1) }.flatten()

        combined_maps.text = [header, *content].join('\n')
        return combined_maps
    }




// Copy of
if else refined version of the subworkflow

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

    def sample_count = sample_sequences.count()

    if (sample_count < 1) {

        println("Please provide sample sequence to analyze.")

    } else if (sample_count = 1) {

        PANPHLAN_MAP (
            PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
            PANPHLAN_DOWNLOADPANGENOME.out.indexes,
            ch_sample_sequences,
            val_species_name)
        ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())

        def ch_mapping_dir = PANPHLAN_MAP.out.sample_map
            .view()

        PANPHLAN_PROFILING (
            ch_mapping_dir,
            PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
            val_species_name)
        ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())

        emit:
            pangenome               = PANPHLAN_DOWNLOADPANGENOME.out.pangenome          // channel: $prefix/{species_name}_pangenome.tsv; pangenome used for selected species
            mapping_directory       = ch_mapping_dir                                    // channel: val(meta), path("${species_name}_map"); mapping directory built for selected species
            panphlan_profile        = PANPHLAN_PROFILING.out.profile_matrix             // channel: val(meta), path("*.tsv"); final profile matrix
            versions                = ch_versions                                       // channel: [ versions.yml ]

    } else if (sample_count > 1) {

        sample_sequences.each { tuple_sequence ->
            PANPHLAN_MAP (
                PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
                PANPHLAN_DOWNLOADPANGENOME.out.indexes,
                tuple_sequence,
                val_species_name
            ) }
        ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())

        def ch_combined_mapping_dir = PANPHLAN_MAP.out.sample_map
            .map{ [ [id:'all_samples'], it[1] ] }.groupTuple( sort: { it.getName() }  )
            .view()

        PANPHLAN_PROFILING (
            ch_combined_mapping_dir,
            PANPHLAN_DOWNLOADPANGENOME.out.pangenome,
            val_species_name)
        ch_versions = ch_versions.mix(PANPHLAN_MAP.out.versions.first())

        emit:
            pangenome               = PANPHLAN_DOWNLOADPANGENOME.out.pangenome          // channel: $prefix/{species_name}_pangenome.tsv; pangenome used for selected species
            mapping_directory       = ch_combined_mapping_dir                           // channel: val(meta), path("${species_name}_map"); mapping directory built for selected species
            panphlan_profile        = PANPHLAN_PROFILING.out.profile_matrix             // channel: val(meta), path("*.tsv"); final profile matrix
            versions                = ch_versions                                       // channel: [ versions.yml ]
    }
}

test("Bacteroides_fragilis - single end, multiple samples - fails") {
        when {
            workflow {
                """
                input[0] = "Bacteroides_fragilis"
                input[1] = Channel.of([ [ id:'test', single_end:false ], // meta map
                        [file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test1_1.fastq.gz", checkIfExists:true),
                        file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test1_2.fastq.gz", checkIfExists:true),
                        file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test2_1.fastq.gz", checkIfExists:true),
                        file("https://github.com/nf-core/test-datasets/raw/modules/data/genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test2_2.fastq.gz", checkIfExists:true)]
                    ])
                input[2] = "multiple"
                """
            }
        }

        then {
            assertAll(
                { assert workflow.success},
                { assert snapshot(workflow.out).match()}
            )
        }
    }
