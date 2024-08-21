process SAM_BAM_SORT_IDX {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/samtools:1.20--h50ea8bc_1"

  label 'cpu_high'
  label 'mem_high'

  input:
  tuple val(meta), path(sam, arity: 1, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.id}.bam") , emit: bam
  tuple val(meta), path("${meta.id}.bai") , emit: bai

  script:
  sort_bam = ${task.ext.sort_bam ?: 'true'}

  """
  #!/usr/bin/bash

  #sam_bn=\$(echo ${sam} | sed "s/\\.sam\$//")
  if ${sort_bam}; then
    samtools view -bS ${sam} | samtools sort -o ${sam.baseName}.bam
  else
    samtools view -bS ${sam} -o ${sam.baseName}.bam
  fi

  samtools index ${sam.baseName}.bam
  """

  stub:

  """
  #!/usr/bin/bash

  touch ${sam.baseName}.bam ${sam.baseName}.bai
  """
}