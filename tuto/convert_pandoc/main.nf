// Process to convert markdown to various formats using Pandoc
// Demonstrates: dynamic output format based on meta.output_format

process CONVERT_PANDOC {
    container 'pandoc/latex:3.1-ubuntu'
    tag "${meta.id}"
    cpus 4
    memory 4.GB

    input:
    tuple val(meta), path(markdown)

    output:
    tuple val(meta), path("${meta.id}.${meta.output_format}"), emit: document

    script:
    def output_format = meta.output_format ?: 'pdf'
    """
    pandoc ${markdown} \\
        -o ${meta.id}.${output_format} \\
        --from markdown \\
        --to ${output_format}
    """

    stub:
    """
    touch ${meta.id}.${meta.output_format ?: 'pdf'}
    """
}
