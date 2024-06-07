def getKallistoStrandArgs (reads_orientation) {  // create a method only to including this command in the script AND in the stub part, while minimizing code repetition
    def list_reads_orientation = ["unstranded", "fr-stranded"  , "rf-stranded"  ]
    def list_kallisto_strands =  [""          , "--fr-stranded", "--rf-stranded"]
    def hit_index = list_reads_orientation.indexOf(reads_orientation)
    if (hit_index >= 0) {
        return list_kallisto_strands[hit_index]
    } else {
        throw new IllegalArgumentException("Invalid value for reads_orientation: '$reads_orientation'. Expected one of: $list_reads_orientation.")
    }
}

//--fragment-length: 0 if UTR, unless if PE and recquired if SE
//--single : ??

process KALLISTO_QUANT {
    container 'quay.io/biocontainers/kallisto:0.50.1--h6de1650_2'
    tag "$smpl_name"

    label 'cpu_medium'
    label 'mem_8G'

    input:
    tuple val(smpl_name), path(reads, arity: 1..2)
    path index
    val(reads_orientation)

    output:
    tuple val(smpl_name), path(smpl_name)


    script:
    def kallisto_strand = getKallistoStrandArgs(reads_orientation)
    """
    #!/usr/bin/bash

    kallisto quant --threads ${task.cpus} \\
        --index ${index} \\
        --output-dir ${smpl_name} \\
        ${kallisto_strand} \\
        ${task.ext.args ?: ''} \\
        ${reads}
    """

    stub:
    def kallisto_strand = getKallistoStrandArgs(reads_orientation)
    """
    #!/usr/bin/bash

    mkdir ${smpl_name}
    touch ${smpl_name}/abundance.tsv ${smpl_name}/abundance.h5 ${smpl_name}/abundance.tsv ${smpl_name}/run_info.json
    """
}

