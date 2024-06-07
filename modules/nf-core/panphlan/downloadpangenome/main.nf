process PANPHLAN_DOWNLOADPANGENOME {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e3818ba00c0f15bcda92c714b4fa500c626067fe:cae224e6965f9abd8b08407231aff3bba632f1de-0':
        'biocontainers/mulled-v2-e3818ba00c0f15bcda92c714b4fa500c626067fe:cae224e6965f9abd8b08407231aff3bba632f1de-0' }"

    input:
    val(species_name)

    output:
    path("indexes")                       , emit: indexes
    path("${prefix}/*_pangenome.tsv")     , emit: pangenome
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${species_name}"
    """
    panphlan_download_pangenome.py \\
        -i ${species_name} \\
        -o . \\
        ${args}

    mkdir indexes
    mv ${prefix}/*.bt2 indexes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panphlan: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panphlan: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
