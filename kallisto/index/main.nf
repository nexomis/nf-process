process KALLISTO_INDEX {

    label 'cpu_x2'
    label 'mem_8G'

    container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/kallisto:0.50.1--h6de1650_2"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.idx") , emit: idx

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    """
    kallisto index \\
        -i ${prefix}.idx \\
        ${fasta}
    """

    stub:
    """
    touch ${prefix}.idx
    """
}
