process SPRING_COMPRESS {
  container 'ghcr.io/nexomis/spring:1.1.1'

  label 'cpu_low'
  label 'mem_med'

  input:
  tuple val(meta), path(files, arity: 1..2, stageAs: "inputs/*")

  output:
  tuple val(meta), path("*.spring", arity: 1), emit: spring

  script:
  def gz_extensions = ['gz', 'gzip', 'z']
  meta.read_type = "spring"
  """
  #!/usr/bin/bash

  spring -q ${params.quality_mode} -t ${task.cpus} -c -i ${files} -o ${meta.id}.spring \\
    ${task.ext.args ? task.ext.args : ''} \\
    ${(files[0].getExtension() in gz_extensions) ? '-g' : ''}
  """

  stub:
  """
  #!/usr/bin/bash

  touch ${meta.id}.spring ${meta.id}.spring${files.size() > 1 ? '.1' : '' }.fq ${files.size() > 1 ? meta.id + '.spring.2.fq' : '' }

  """
}