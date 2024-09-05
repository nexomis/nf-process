

process SPRING_DECOMPRESS {

  container 'ghcr.io/nexomis/spring:1.1.1'

  label 'cpu_low'
  label 'mem_12G'

  input:
  tuple val(meta), path(spring_file, arity: 1)

  output:
  tuple val(meta), path("${meta.id}*.fq.gz", arity: 1..2)

  script:
  def R1 = meta.id + "_R1.fq.gz"
  def R2 = meta.id + "_R2.fq.gz"
  meta.remove("read_type")
  """
  #!/usr/bin/bash

  spring -d -g -t ${task.cpus} -i $spring_file -o $R1 $R2
  """

  stub:
  """
  #!/usr/bin/bash

  touch ${meta.id}.fq.gz

  """
}
