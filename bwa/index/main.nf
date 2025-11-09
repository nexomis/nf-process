process BWA_INDEX {
  container "quay.io/biocontainers/bwa:0.7.18--he4a0461_1"
  tag "$meta.id"
  cpus 1
  label 'mem_8GB'

  input:
  tuple val(meta), path(fasta, arity: 1, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.label ?: meta.id}/", type: 'dir'), emit: idx

  script:
  """
  #!/usr/bin/bash

  mkdir -p ${meta.label ?: meta.id}/
  
  bwa index \\
    -p ${meta.label ?: meta.id}/index \\
    ${task.ext.args ?: ''} \\
    $fasta
  """

  stub:
  """
  #!/usr/bin/bash

  mkdir -p ${meta.label ?: meta.id}/
  touch ${meta.label ?: meta.id}/index.bwt
  """
}
