process HUMANN_REGROUPTABLE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0':
        'biocontainers/humann:3.8--pyh7cba7a3_0' }"


    input:
    tuple val(meta), path(regroup_input)
    val renorm_units
    val group_type
    path utility_mapping_db


    output:
    tuple val(meta), path("*_genefamilies-${renorm_units}-${group_type}.tsv")   , emit: humann_output_regroup
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    humann_regroup_table \\
        --input ${regroup_input} \\
        --groups ${group_type} \\
        --precision 3 \\
        --ungrouped Y \\
        --protected Y \\
        --output ${prefix}_genefamilies-${renorm_units}-${group_type}.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1 | sed 's/^.*humann //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_genefamilies-${renorm_units}-${group_type}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1 | sed 's/^.*humann //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
