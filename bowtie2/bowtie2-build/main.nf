
process BOWTIE2_BUILD {
  container "quay.io/biocontainers/bowtie2:2.5.4--he20e202_3"
  tag "$meta.id"
  cpus 16
  memory 8.GB // 4GB max of foot print: https://hpc.nih.gov/apps/bowtie2.html

  input:
  tuple val(meta), path(fasta, arity: 1, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}/", type: 'dir')     , emit: idx

  script:

  """
  #!/usr/bin/bash

  mkdir ${meta.label ?: meta.id}/
  bowtie2-build --threads $task.cpus \\
    ${task.ext.args ?: ''} \\
    $fasta \\
    ${meta.label ?: meta.id}/index
  """

  stub:

  """
  #!/usr/bin/bash

  mkdir ${meta.label ?: meta.id}/
  touch ${meta.label ?: meta.id}}/${meta.label ?: meta.id}.rev.1.bt2
  """
}