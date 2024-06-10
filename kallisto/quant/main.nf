def formatKallistoArgs (strand, sgl_overhang, reads, mean_FlD, sd_FlD) {  // create a method only to including this command in the script AND in the stub part, while minimizing code repetition
    // strand
    def list_strand = ["unstranded", "fr-stranded"  , "rf-stranded"  ]
    def list_kallisto_strand =  [""          , "--fr-stranded", "--rf-stranded"]
    def hit_index = list_strand.indexOf(strand)
    def strand_args
    if (hit_index >= 0) {
        strand_args = list_kallisto_strand[hit_index]
    } else {
        throw new IllegalArgumentException("Invalid value for strand: '$strand'. Expected one of: $list_strand.")
    }

    // sgl_overhand
    if (!(sgl_overhang instanceof Boolean)) {
        throw new IllegalArgumentException("Invalid value for sgl_overhang: $sgl_overhang. Expected a boolean value (true or false).")
    }
    def sgl_overhang_args = sgl_overhang ? "--single-overhang" : ""

    // FlD
    def FlD_args = ""
    def mean_FlD_args = ""
    def sd_FlD_args = ""
    if (reads.size() == 1) {
        if (!mean_FlD) {
            mean_FlD = 230
        }
        if (!sd_FlD) {
            sd_FlD = 60
        }
        FlD_args = "--fragment-length $mean_FlD --sd $sd_FlD"
    } else if (reads.size() == 2) {
        if (mean_FlD) {
            mean_FlD_args = "--fragment-length $mean_FlD"
        }
        if (sd_FlD) {
            sd_FlD_args = "--sd $sd_FlD"
        }
        FlD_args = "$mean_FlD_args $sd_FlD_args".trim()
    } else {
        throw new IllegalArgumentException("Invalid nb of file for $smpl_name: $reads. Expected 1 or 2 file(s) path(s)")
    }
    
    //all_args = "$strand_args $sgl_overhang_args $FlD_args".trim()
    return "$strand_args $sgl_overhang_args $FlD_args".trim()
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
    tuple val(strand), val(sgl_overhang), val(mean_FlD), val(sd_FlD)

    output:
    tuple val(smpl_name), path(smpl_name)


    script:
    def args = formatKallistoArgs (strand, sgl_overhang, reads, mean_FlD, sd_FlD)
    """
    #!/usr/bin/bash

    kallisto quant --threads ${task.cpus} \\
        --index ${index} \\
        --output-dir ${smpl_name} \\
        ${args} \\
        ${task.ext.args ?: ''} \\
        ${reads}
    """

    stub:
    def args = formatKallistoArgs (strand, sgl_overhang, reads, mean_FlD, sd_FlD)
    """
    #!/usr/bin/bash

    mkdir ${smpl_name}
    touch ${smpl_name}/abundance.tsv ${smpl_name}/abundance.h5 ${smpl_name}/abundance.tsv ${smpl_name}/run_info.json
    """
}

