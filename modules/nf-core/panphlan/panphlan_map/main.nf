process PANPHLAN_MAP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e3818ba00c0f15bcda92c714b4fa500c626067fe:cae224e6965f9abd8b08407231aff3bba632f1de-0':
        'biocontainers/mulled-v2-e3818ba00c0f15bcda92c714b4fa500c626067fe:cae224e6965f9abd8b08407231aff3bba632f1de-0' }"

    input:
    path(pangenome)
    path(indexes)
    tuple val(meta), path(sample_sequence)
    val(species_name)

    output:
    tuple val(meta), path("${species_name}_map")    , emit: mapping_dir
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    panphlan_map.py \\
        -p $pangenome \\
        --indexes ${indexes}/${species_name} \\
        -i ${sample_sequence} \\
        -o ${meta.id}_${species_name}.tsv \\
        ${args}

    mkdir ${species_name}_map
    mv ${meta.id}_${species_name}.tsv ${species_name}_map

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panphlan: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panphlan: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
