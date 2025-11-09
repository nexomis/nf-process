// Process to write greeting markdown files
// Demonstrates: using meta attributes for dynamic content

process WRITE_GREETING {
    container 'ubuntu:22.04'
    tag "${meta.id}"
    cpus 4
    memory 4.GB

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}.md"), emit: markdown

    script:
    def title = meta.title ?: 'Mr/Ms'
    def level = meta.level ?: 'beginner'
    """
    cat > ${meta.id}.md <<EOF
# Hello, ${title} ${meta.id}!

Welcome to the tutorial for **${level}** level.

Best regards,  
The Tutorial Team
EOF
    """

    stub:
    """
    touch ${meta.id}.md
    """
}
