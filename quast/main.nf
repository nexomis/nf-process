process QUAST {
  container "staphb/quast:5.2.0"
  // container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/quast:5.2.0--py312pl5321hc60241a_4"  // ~700MB !
  //"staphb/quast:5.2.0-slim" 80MB instead of 362MB but seems not compatible with bam ...: 
    //Running Reads analyzer...
    //Compiling BWA (details are in /usr/local/lib/python3.10/dist-packages/quast_libs/bwa/make.log and make.err)
    //[Errno 13] Permission denied: '/usr/local/lib/python3.10/dist-packages/quast_libs/bwa/make.log'

  label 'cpu_med'
  label 'mem_med'
  
  input:
  tuple val(meta) , path(assembly)   // assembled scaffold (can be multiple files)
  tuple val(meta_bam) , path(bam)    // [optionnel] sorted bam: cleanedreads on assembled scaffold (need here 'bai' ?)
  path(ref_fa)   // [optionnel] reference genome  (can be mutliple files ?)
  // include spades ??

  output:
  tuple val(meta), path("${meta.id}/report.html", type: 'file') , emit: html
  tuple val(meta), path("${meta.id}/report.tsv", type: 'file')  , emit: tsv
  tuple val(meta), path("${meta.id}", type: 'dir')              , emit: quast_output  // all output
  tuple val(meta), path("${meta.id}.log", type: 'file')         , emit: log

  script:
  // quast option: --rna-finding ?  --glimmer ?  --features(gff) ? --bam ?
  """
  quast.py \\
    --output-dir ${meta.id} \\
    --threads $task.cpus \\
    ${bam ? "--bam ${bam}" : ""} \\
    ${ref_fa ? "-r ${ref_fa}" : ""} \\
    ${task.ext.args ?: ''} \\
    ${assembly.join(' ')} \\
    2> ${meta.id}.log
  """

  stub:
  """
  #!/usr/bin/bash
  mkdir ${meta.id}
  touch ${meta.id}.log ${meta.id}/report.html ${meta.id}/report.tsv ${meta.id}/report.tex ${meta.id}/icarus.html ${meta.id}/report.pdf 
  """
}
