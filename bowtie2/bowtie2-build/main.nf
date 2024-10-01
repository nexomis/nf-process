
process BOWTIE2_BUILD {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/bowtie2:2.5.4--he20e202_3"

  label 'cpu_high'
  label 'mem_8G'

  input:
  tuple val(meta), path(fasta, arity: 1, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.id}/", type: 'dir')     , emit: idx

  script:

  """
  #!/usr/bin/bash

  mkdir ${meta.id}/
  bowtie2-build --threads $task.cpus \\
    ${task.ext.args ?: ''} \\
    $fasta \\
    ${meta.id}/${meta.id} \\

  """

  stub:

  """
  #!/usr/bin/bash

  mkdir ${meta.id}/
  touch ${meta.id}/${meta.id}.rev.1.bt2
  """
}