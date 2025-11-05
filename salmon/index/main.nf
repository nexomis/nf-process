process SALMON_INDEX {
    container "quay.io/biocontainers/salmon:1.10.3--haf24da9_3"
    tag "$meta.id"
    label 'cpu_x4'
    label 'mem_16G'

    input:
    tuple val(meta), path(fasta, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.label ?: meta.id}.idx", arity: 1, type: "dir"), emit: idx

    script:
    """
    salmon index \\
        -t $fasta \\
        -i ${meta.label ?: meta.id}.idx \\
        -p $task.cpus \\
        ${task.ext.args ?: ''}
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/info.json
    """
}
