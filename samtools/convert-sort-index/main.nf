// todo: make sort and flagstat optionall ?

process SAM_BAM_SORT_IDX {
  container "quay.io/biocontainers/samtools:1.20--h50ea8bc_1"

  cpus 16
  memory 2.GB * task.cpus

  input:
  tuple val(meta), path(sam, arity: 1, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}.bam"), path("${meta.label ?: meta.id}.bam.bai") , optional:false , emit: bam_bai
  tuple val(meta), path("${meta.label ?: meta.id}.flagstats", type: 'file') , optional:false , emit: flagstat

  script:
  //sort_bam = task.ext.sort_bam ?: 'true'
  //flagstat = task.ext.flagstat ?: 'true'

  """
  #!/usr/bin/bash

  samtools view ${task.ext.args_samtools_view ?: ''} \\
    -bS ${sam} | \\
    samtools sort ${task.ext.args_samtools_sort ?: ''} \\
      -o ${meta.label ?: meta.id}.bam \\
      -@ $task.cpus

  samtools index ${task.ext.args_samtools_index ?: ''} ${meta.label ?: meta.id}.bam

  # interest option '-O/--output-fmt'='json' or 'tsv'
  samtools flagstat \\
    -@ $task.cpus \\
    ${meta.label ?: meta.id}.bam \\
    ${task.ext.args_samtools_flagstat ?: ''} \\
    > ${meta.label ?: meta.id}.flagstats
  """

  stub:
  """
  #!/usr/bin/bash

  touch ${meta.label ?: meta.id}.bam ${meta.label ?: meta.id}.bam.bai ${meta.label ?: meta.id}.flagstats
  """
}