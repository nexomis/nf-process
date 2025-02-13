process SALMON_QUANT {
    container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/salmon:1.10.3--haf24da9_3"

    label 'cpu_med'
    label 'mem_8G'

    input:
    tuple val(meta), path(reads, stageAs: 'input_reads/*', arity: 1..2)
    tuple val(meta2), path(index, stageAs: 'input_index/*')

    output:
    tuple val(meta), path("${meta.label ?: meta.id}/quant.sf.gz"), emit: quant
    tuple val(meta), path("${meta.label ?: meta.id}/aux_info/${meta.label ?: meta.id}.meta_info.json"), emit: log
    tuple val(meta), path("${meta.label ?: meta.id}/${meta.label ?: meta.id}.lib_format_counts.json"), emit: lib_count
    tuple val(meta), path("${meta.label ?: meta.id}/libParams/${meta.label ?: meta.id}.flenDist.txt"), emit: flenDist

    script:
    def reads_args = (reads.size() == 1) ? "-1 ${reads[0]} -2 ${reads[1]}" : "-r ${reads[0]}"
    """
    salmon quant \\
        -i $index \\
        ${meta.args_salmon ?: '-l A'} \\
        $reads_args \\
        -o ${meta.label ?: meta.id} \\
        -p $task.cpus \\
        ${task.ext.args ?: ''}
    
    mv ${meta.label ?: meta.id}/aux_info/meta_info.json ${meta.label ?: meta.id}/aux_info/${meta.label ?: meta.id}.meta_info.json
    mv ${meta.label ?: meta.id}/lib_format_counts.json ${meta.label ?: meta.id}/${meta.label ?: meta.id}.lib_format_counts.json
    mv ${meta.label ?: meta.id}/libParams/flenDist.txt ${meta.label ?: meta.id}/libParams/${meta.label ?: meta.id}.flenDist.txt

    gzip ${meta.label ?: meta.id}/quant.sf
    """

    stub:
    """
    mkdir -p ${meta.label ?: meta.id}
    touch ${meta.label ?: meta.id}/quant.sf
    """
}
