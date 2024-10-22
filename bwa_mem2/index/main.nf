// IN PROGRESS !
process BWA_MEM2_INDEX {
  container "quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_3"  // ????

  label 'cpu_high'
  label 'mem_8GB'

  input:
  tuple val(meta), path(fasta, arity: 1, stageAs: 'input_raw/*')

  output:
  tuple val(meta), path("${meta.id}/", type: 'dir'), emit: idx

  script:
  """
  #!/usr/bin/bash

  mkdir -p ${meta.id}/
  
  bwa-mem2 index -t $task.cpus $fasta ${meta.id}/index
  """

  stub:
  """
  #!/usr/bin/bash

  mkdir -p ${meta.id}/
  touch ${meta.id}/index.bwt.2bit.64
  """
}
