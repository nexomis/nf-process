process KALLISTO_INDEX {
    tag "$meta.id"
    cpus 2
    memory 8.GB

    container "quay.io/biocontainers/kallisto:0.50.1--h6de1650_2"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.idx") , emit: idx

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    """
    kallisto index \\
        $args \\
        -i ${prefix}.idx \\
        ${fasta}
    """

    stub:
    def prefix = meta.id
    """
    touch ${prefix}.idx
    """
}
