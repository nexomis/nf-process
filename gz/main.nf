process GZ {
  container 'quay.io/biocontainers/pigz:2.8'

  label 'cpu_low'
  label 'mem_2G_per_cpu'
  
  input:
  tuple val(meta), path(files, arity: 1..2)

  output:
  tuple val(meta), path("*.gz")

  script:
  def cmd = "pigz -f -p " + task.cpus + " --keep --no-time "
  """
  $cmd $files
  """
  
  stub:
  """
  touch ${files[0]}
  ${(files.size() === 2)? 'touch ' + files[1].toString() : ''}
  """

}