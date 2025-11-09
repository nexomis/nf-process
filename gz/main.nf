process GZ {
  container 'quay.io/biocontainers/pigz:2.8'
  tag "$meta.id"
  cpus 4
  memory 8.GB
  
  input:
  tuple val(meta), path(files, arity: 1..2)

  output:
  tuple val(meta), path("*.gz")

  script:
  def cmd = "pigz -f -p $task.cpus --keep --no-time "
  """
  $cmd $files
  """
  
  stub:
  """
  touch ${files[0]}
  ${(files.size() == 2)? 'touch ' + files[1].toString() : ''}
  """

}