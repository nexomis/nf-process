process KALLISTO_QUANT {
    container 'ghcr.io/nexomis/kallisto:v0.50.1'
    tag "$smpl_name"

    label 'cpu_medium'
    label 'mem_8G'

    input:
    tuple val(smpl_name), path(reads, arity: 1..2)
    path index

    output:
    tuple val(smpl_name), path(smpl_name)

    //add $strand (not recquired, in case of unstranded library ?)
    //--fragment-length: 0 if UTR, unless if PE, recquired if SE
    //--single : ??

    script:

    """
    #!/usr/bin/bash

    kallisto quant --threads ${task.cpus} \\
        --index ${index} \\
        --output-dir ${smpl_name} \\
        ${task.ext.args ?: ''} \\
        ${reads}
    """

    stub:
    """
    #!/usr/bin/bash

    mkdir ${smpl_name}
    touch ${smpl_name}/abundance.tsv ${smpl_name}/abundance.h5 ${smpl_name}/abundance.tsv ${smpl_name}/run_info.json
    """
}

