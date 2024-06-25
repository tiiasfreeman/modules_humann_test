process HUMANN_RENORMTABLE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0':
        'biocontainers/humann:3.8--pyh7cba7a3_0' }"


    input:
    val humann_output_type
        // "genefamilies" OR "pathabundance"
    tuple val(meta), path(humann_output)
    val renorm_units
        // renormalization scheme; copies per million [cpm] OR relative abundance [relab]

    output:
    tuple val(meta), path("*_${humann_output_type}-${renorm_units}.tsv")     , emit: humann_output_renorm
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    humann_renorm_table \\
        --input ${humann_output} \\
        --output ${prefix}_${humann_output_type}-${renorm_units}.tsv \\
        --units ${renorm_units} \\
        --mode community \\
        --special y \\
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
    touch ${prefix}_${humann_output_type}-${renorm_units}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1 | sed 's/^.*humann //; s/Using.*\$//' ))
    END_VERSIONS
    """
}

