process SAM_BAM_SORT_IDX {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/samtools:1.20--h50ea8bc_1"

  label 'cpu_high'
  label 'mem_2G_per_cpu'

  input:
  tuple val(meta), path(sam, arity: 1, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}.bam"), path("${meta.label ?: meta.id}.bam.bai")

  script:
  sort_bam = task.ext.sort_bam ?: 'true'
  """
  #!/usr/bin/bash

  samtools view ${task.ext.args_samtools_view ?: ''} \\
    -bS ${sam} | \\
    samtools sort ${task.ext.args_samtools_sort ?: ''} \\
      -o ${meta.label ?: meta.id}.bam \\
      -@ $task.cpus

  samtools index ${task.ext.args_samtools_index ?: ''} ${meta.label ?: meta.id}.bam
  """

  stub:
  """
  #!/usr/bin/bash

  touch ${meta.label ?: meta.id}.bam ${meta.label ?: meta.id}.bam.bai
  """
}