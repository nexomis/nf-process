process SAM_BAM_SORT_IDX {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/samtools:1.20--h50ea8bc_1"

  label 'cpu_high'
  label 'mem_high'

  input:
  tuple val(meta), path(sam, arity: 1, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.id}.bam") , emit: bam
  tuple val(meta), path("${meta.id}.bam.bai") , emit: bai
  //tuple val(meta), path("${meta.id}.bam", "${meta.id}.bai") , emit: bam_bai
  // output_dir usefull ?

  script:
  sort_bam = task.ext.sort_bam ?: 'true'
  // TODO: add 2 boolean parameters in otpion :
  // for status of sorting and indexing steps (and determine automatically the need to convert step using extension of input file)
  // probably more appropriate to use internal process parameters instead of external task parameters, as it would be useful to manage this automatically in the worflow (useful to manage the possibility of overwriting the option with external parameters??).

  """
  #!/usr/bin/bash

  #sam_bn=\$(echo ${sam} | sed "s/\\.sam\$//")
  if ${sort_bam}; then
    samtools view ${task.ext.args_view ?: ''} \\
      -bS ${sam} | \\
        samtools sort ${task.ext.args_sort ?: ''} \\
        -o ${sam.baseName}.bam \\
        -@ $task.cpus
  else
    samtools view ${task.ext.args_view ?: ''} \\
      -bS ${sam} \\
      -o ${sam.baseName}.bam
  fi

  samtools index ${task.ext.args_index ?: ''} ${sam.baseName}.bam
  """

  stub:

  """
  #!/usr/bin/bash

  touch ${sam.baseName}.bam ${sam.baseName}.bam.bai
  """
}