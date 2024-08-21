
process BOWTIE2_BUILD {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/bowtie2:2.5.4--he20e202_3"

  label 'cpu_high'
  label 'mem_8G'

  input:
  tuple val(meta), path(fasta, arity: 1, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.id}") , emit: idx
  path "${meta}.log"                  , emit: log

  script:

  """
  #!/usr/bin/bash

  mkdir ${meta.id}/
  bowtie2-build ${task.ext.args ?: ''} --threads $task.cpus $fasta ${meta.id}/${meta.id} 2>${meta.id}.log
  """

  stub:

  """
  #!/usr/bin/bash

  mkdir ${meta.id}/
  touch mkdir ${meta.id}/${meta.id}.rev.1.bt2 ${meta.id}.log
  """
}