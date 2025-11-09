process GUNZIP {
  container 'quay.io/biocontainers/pigz:2.8'
  tag "$meta.id"
  cpus 4
  label 'mem_8Gu'
  
  input:
  tuple val(meta), path(files, arity: 1..2, stageAs: "ungz/*")

  output:
  tuple val(meta), path("ungz/*", arity: 1..2)

  script:
  """
  cd ungz
  pigz -f -p ${task.cpus} --no-time -d *
  """
  
  stub:
  """
  touch ${files[0]}
  ${(files.size() == 2)? 'touch ' + files[1].toString() : ''}
  """

}