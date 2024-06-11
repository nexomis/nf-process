

process SPRING_DECOMPRESS {
  container 'ghcr.io/nexomis/spring:1.1.1'

  label 'cpu_low'
  label 'mem_12G'

  input:
  tuple val(meta), path(spring_file, arity: 1)

  output:
  tuple val(meta), path("${meta.id}*.fq.gz", arity: 1..2)

  script:
  """
  #!/usr/bin/bash

  spring -d -g -t ${task.cpus} -i ${meta.id}.spring -o ${meta.id}.fq.gz
  if [[ -e ${meta.id}.fq.gz.1 && -e ${meta.id}.fq.gz.2 ]]; then
    # Move (rename) the files
    mv ${meta.id}.fq.gz.1 ${meta.id}_R1.fq.gz
    mv ${meta.id}.fq.gz.2 ${meta.id}_R2.fq.gz
  fi
  """

  stub:
  """
  #!/usr/bin/bash

  touch ${meta.id}.fq.gz

  """
}
