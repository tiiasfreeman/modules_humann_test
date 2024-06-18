process PANPHLAN_PROFILING {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e3818ba00c0f15bcda92c714b4fa500c626067fe:cae224e6965f9abd8b08407231aff3bba632f1de-0':
        'biocontainers/mulled-v2-e3818ba00c0f15bcda92c714b4fa500c626067fe:cae224e6965f9abd8b08407231aff3bba632f1de-0' }"

    input:
    tuple val(meta), path(sample_map)
    path(pangenome)
    val(species_name)

    output:
    tuple val(meta), path("*.tsv"), emit: profile_matrix
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${species_name}_map_dir
    mv *_${species_name}_*.tsv ${species_name}_map_dir

    panphlan_profiling.py \\
        --i_dna ${species_name}_map_dir \\
        --o_matrix ${species_name}_profile.tsv \\
        -p ${pangenome} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panphlan: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${species_name}_map_dir
    mv *_${species_name}_*.tsv ${species_name}_map_dir

    touch ${species_name}_profile.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panphlan: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
